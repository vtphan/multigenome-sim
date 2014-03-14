package main

import (
   "flag"
   "log"
   "fmt"
   "os"
   "bytes"
   "io/ioutil"
   "bufio"
   "strings"
   "strconv"
)

var SEQ []byte

func ReadSequence(file string) {
   f, err := os.Open(file)
   if err != nil {
      panic(err)
   }
   defer f.Close()

   if file[len(file)-6:] == ".fasta" {
      scanner := bufio.NewScanner(f)
      byte_array := make([]byte, 0)
      for scanner.Scan() {
         line := scanner.Bytes()
         if len(line)>0 && line[0] != '>' {
            byte_array = append(byte_array, bytes.Trim(line,"\n\r ")...)
         }
      }
      SEQ = append(byte_array, byte('$'))
   } else {
      byte_array, err := ioutil.ReadFile(file)
      if err != nil {
         panic(err)
      }
      SEQ = append(bytes.Trim(byte_array, "\n\r "), byte('$'))
   }
}


func print_read_substr(read, substr string, pos int, is_error bool) {
   fmt.Println(read, "index", pos, "is_error:", is_error)
   for k:=0; k<len(substr); k++ {
      if k != pos {
         fmt.Printf(" ")
      } else {
         fmt.Printf("%s", string(substr[pos]))
      }
   }
   fmt.Println(" (genome substring)")
}

func check_for_invalid_symbol_in_genome(s []byte) {
   for i:=0; i<len(s); i++ {
      base := s[i]
      if base != 'A' && base != 'C' && base != 'G' && base != 'T'  && base != '$' && base != 'N' {
         fmt.Println("\tGenome contains an invalid nucleotide.", base)
         return
      }
      if base == '$' {
         fmt.Println("\tGenome contains $ at position", i)
      }
   }
}

func check_for_invalid_symbol_in_read(s string) {
   for i:=0; i<len(s); i++ {
      base := s[i]
      if base != 'A' && base != 'C' && base != 'G' && base != 'T' && base != 'N' {
         fmt.Println("\tRead contains an invalid nucleotide:", base, s)
         return
      }
   }
}

func verify(read string, positions []int, error_positions map[int]bool) {
   check_for_invalid_symbol_in_read(read)
   n := len(read)
   for i:=0; i<len(positions); i++ {
      pos := positions[i]
      substr := string(SEQ[pos: pos+n])
      for j:=0; j<len(read); j++ {
         _, is_error := error_positions[j]

         if substr[j]!='N' && ((is_error && read[j] == substr[j]) || (!is_error && read[j] != substr[j])) {
            print_read_substr(read, substr, j, is_error)
         }
      }
   }
}

func main () {
   var genome_file = flag.String("s", "", "Specify a file containing the genome.")
   var read_file = flag.String("r", "", "Specify a file containing the reads.")
   flag.Parse()
   if *genome_file == "" || *read_file == "" {
      log.Fatalln("Must provide both genome file (-s option) and erad file (-r option).")
   }
   ReadSequence(*genome_file)
   fmt.Println("Checking invalid symbols in genome...")
   check_for_invalid_symbol_in_genome(SEQ)

   fmt.Println("Checking invalid reads...")
   f, err := os.Open(*read_file)
   if err != nil {
      log.Fatalln("Unable to read", *read_file)
   }
   scanner := bufio.NewScanner(f)
   count := 0
   for scanner.Scan() {
      // line := strings.Trim(scanner.Text(), "\n\r\t ")
      line := scanner.Text()
      items := strings.Split(line, " ")
      read := items[0]
      occ, _ := strconv.Atoi(items[1])
      positions := make([]int, occ)
      for i:=0; i<occ; i++ {
         positions[i], _ = strconv.Atoi(items[2+i])
      }
      error_count, _ := strconv.Atoi(items[2+occ])
      error_positions := make(map[int]bool)
      for i:=0; i<error_count; i++ {
         p, _ := strconv.Atoi(items[3+occ+i])
         error_positions[p] = true
      }
      verify(read, positions, error_positions)
      count++
   }
   err = scanner.Err()
   if err != nil {
      fmt.Println("May not process all reads to due error:", err)
   }

   fmt.Printf("Finish verifying %d reads (from %s) in genome %s (len=%d)\n", count, *read_file, *genome_file, len(SEQ))
}