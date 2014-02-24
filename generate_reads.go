/*
   Copyright 2013 Vinhthuy Phan
   Usage:  go run generate_reads.go --help

Expected input format: fasta format file.

Each line of the output has the following format:
   R  N p1 … pN  E q1 … qE
where
   + R is the content of the read
   + N is the number of occurrences of this read in the genome
   + p1 to pN are the locations of the read in the genome
   + E is the number of errors
   + q1 to qE is the locations of the errors in the read.  If E is equal to 0, then this list is empty.

Example:

ATTTAGGTACTACTAACTCGTGTGTTGCAGTTTTCGAAAATGAAAAAGTCGCGTTATAGAAAATTCAGAACGTGCCCATACTACTCACCCTTCTATAATT 1 25 2 51 70

   This means the read occurs at 1 genome location: 25, with 2 errors at read locations: 51 and 70.
*/
package main

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"os"
	"sort"
	"flag"
	"log"
	"encoding/gob"
	"bufio"
	"path"
	"math/rand"
	"time"
)

var Debug bool

//-----------------------------------------------------------------------------
// Global variables: sequence (SEQ), suffix array (SA), BWT, FM index (C, OCC)
//-----------------------------------------------------------------------------
var SEQ []byte
var rand_gen = rand.New(rand.NewSource(time.Now().UnixNano()))

type position uint32

type Index struct{
	SA []position 						// suffix array
	C map[byte]position  				// count table
	OCC map[byte][]position 			// occurence table

	END_POS position 					// position of "$" in the text
	SYMBOLS []int  				// sorted symbols
	EP map[byte]position 				// ending row/position of each symbol

	LEN position
	// un-exported variables
	bwt []byte
	freq map[byte]position  // Frequency of each symbol
}
//
//-----------------------------------------------------------------------------
type BySuffix []position

func (s BySuffix) Len() int { return len(s) }
func (s BySuffix) Swap(i, j int) { s[i], s[j] = s[j], s[i] }
func (s BySuffix) Less(i, j int) bool { return (bytes.Compare(SEQ[s[i]:], SEQ[s[j]:]) == -1) }


//-----------------------------------------------------------------------------
func _save(thing interface{}, filename string, error_message string) {
	var out bytes.Buffer
	enc := gob.NewEncoder(&out)
	err := enc.Encode(thing)
	if err != nil {
		log.Fatal(error_message)
	}
	fmt.Println("save", filename)
	ioutil.WriteFile(filename, out.Bytes(), 0600)
}
//-----------------------------------------------------------------------------
func _load(thing interface{}, filename string) {
	fin,err := os.Open(filename)
	decOCC := gob.NewDecoder(fin)
	err = decOCC.Decode(thing)
	if err != nil {
		fmt.Println("Unable to read file ("+filename+"): ",err)
	}
}

//-----------------------------------------------------------------------------
func _load_occ(filename string, Len position) []position {
	thing := make([]position, Len)
	fin,err := os.Open(filename)
	decOCC := gob.NewDecoder(fin)
	err = decOCC.Decode(&thing)
	if err != nil {
		log.Fatal("Error loading occ table:", filename, err)
	}
	return thing
	// fmt.Println(thing[key], key)
}

//-----------------------------------------------------------------------------
// Load FM index
// Usage:  idx := Load(index_file)
func Load (dir string) *Index {
	I := new(Index)
	_load(&I.C, path.Join(dir, "c"))
	_load(&I.SA, path.Join(dir, "sa"))
	_load(&I.END_POS, path.Join(dir, "end_pos"))
	_load(&I.SYMBOLS, path.Join(dir, "symbols"))
	_load(&I.EP, path.Join(dir, "ep"))
	_load(&I.LEN, path.Join(dir, "len"))

	I.OCC = make(map[byte][]position)
	for _,symb := range I.SYMBOLS {
		I.OCC[byte(symb)] = _load_occ(path.Join(dir, "occ."+string(symb)), I.LEN)
	}
	return I
}

//-----------------------------------------------------------------------------
func (I *Index) Save(file string) {
	// save(I, file+".fm", "Fail to save fm index")
	dir := file + ".index"
	os.Mkdir(dir, 0777)

	for symb := range I.OCC {
		_save(I.OCC[symb], path.Join(dir, "occ." + string(symb)),"Fail to save to occ."+string(symb))
	}
	_save(I.SA, path.Join(dir,"sa"), "Fail to save suffix array")
	_save(I.C, path.Join(dir,"c"), "Fail to save count")
	_save(I.END_POS, path.Join(dir,"end_pos"), "Fail to save end_pos")
	_save(I.SYMBOLS, path.Join(dir,"symbols"), "Fail to save symbols")
	_save(I.EP, path.Join(dir,"ep"), "Fail to save ep")
	_save(I.LEN, path.Join(dir,"len"), "Fail to save len")
}
//-----------------------------------------------------------------------------
// BWT is saved into a separate file
func (I *Index) BuildSA_BWT() {
	I.LEN = position(len(SEQ))
	I.freq = make(map[byte]position)
	I.SA = make([]position, I.LEN)
	I.bwt = make([]byte, I.LEN)
	I.C = make(map[byte]position)
	I.OCC = make(map[byte][]position)
	I.EP = make(map[byte]position)
	var i position
	for i = 0; i < I.LEN; i++ {
		I.SA[i] = i
		I.freq[SEQ[i]]++
	}
	for c := range I.freq {
		I.SYMBOLS = append(I.SYMBOLS, int(c))
		I.OCC[c] = make([]position, I.LEN)
		I.C[c] = 0
	}
	sort.Ints(I.SYMBOLS)
	sort.Sort(BySuffix(I.SA))

	for i = 0; i < I.LEN; i++ {
		I.bwt[i] = SEQ[(I.LEN+I.SA[i]-1)%I.LEN]
		if I.bwt[i] == '$' {
			I.END_POS = i
		}
	}
}

//-----------------------------------------------------------------------------
func (I *Index) BuildIndex() {
	for i := 1; i < len(I.SYMBOLS); i++ {
		curr_c, prev_c := byte(I.SYMBOLS[i]), byte(I.SYMBOLS[i-1])
		I.C[curr_c] = I.C[prev_c] + I.freq[prev_c]
		I.EP[curr_c] = I.C[curr_c] + I.freq[curr_c] - 1
	}

	for i := 0; i < len(I.bwt); i++ {
		I.OCC[I.bwt[i]][i] = 1
		if i > 0 {
			for symbol := range I.OCC {
				I.OCC[symbol][i] += I.OCC[symbol][i-1]
			}
		}
	}
	I.SYMBOLS = I.SYMBOLS[1:]
	delete(I.OCC, '$')
	delete(I.C, '$')
}

//-----------------------------------------------------------------------------
// Search for all occurences of SEQ[j:j+read_len] in SEQ
//-----------------------------------------------------------------------------

func (I *Index) Search(j position, read_len position) []position {
	var sp, ep, offset position
	var ok bool

	c := SEQ[j+read_len-1]
	sp, ok = I.C[c]
	if ! ok {
		return make([]position, 0)
	}
	ep = I.EP[c]
	// if Debug { fmt.Println("pattern: ", string(pattern), "\n\t", string(c), sp, ep) }
	for i:=int(read_len-2); sp <= ep && i >= 0; i-- {
  		c = SEQ[j+position(i)]
  		offset, ok = I.C[c]
  		if ok {
			sp = offset + I.OCC[c][sp - 1]
			ep = offset + I.OCC[c][ep] - 1
		} else {
			return make([]position, 0)
		}
  		// if Debug { fmt.Println("\t", string(c), sp, ep) }
	}
	res := make([]position, ep-sp+1)
	for k:=sp; k<=ep; k++ {
		res[k-sp] = I.SA[k]
	}
 	return res
}

//-----------------------------------------------------------------------------
func (I *Index) show() {
	fmt.Printf(" %8s  OCC\n", "C")
	for i := 0; i < len(I.SYMBOLS); i++ {
		c := byte(I.SYMBOLS[i])
		fmt.Printf("%c%8d  %d\n", c, I.C[c], I.OCC[c])
	}
	fmt.Println(I.SYMBOLS)
}


//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------
// Build FM index given the file storing the text.
// Usage:	idx := Build(text_file)
func Build (file string) *Index {
	I := new(Index)
	ReadSequence(file)
	I.BuildSA_BWT()
	I.BuildIndex()
	return I
}

//-----------------------------------------------------------------------------
func print_byte_array(a []byte) {
	for i := 0; i < len(a); i++ {
		fmt.Printf("%c", a[i])
	}
	fmt.Println()
}
//-----------------------------------------------------------------------------

func random_error(base byte) byte {
	not_A, not_T, not_C, not_G := []byte{'C','G','T'}, []byte{'C','G','A'}, []byte{'A','G','T'}, []byte{'C','A','T'}
	var c byte
	switch (base){
		case 'A': c = not_A[rand_gen.Intn(3)]
		case 'C': c = not_C[rand_gen.Intn(3)]
		case 'G': c = not_G[rand_gen.Intn(3)]
		case 'T': c = not_T[rand_gen.Intn(3)]
	}
	if ! Debug {
		return c
	}
	r := []byte{c}
	return bytes.ToLower(r)[0]
}

//-----------------------------------------------------------------------------
func main() {
	var seq_file = flag.String("s", "", "Specify a file containing the sequence.")
	var rl = flag.Int("l", 100, "Read length.")
   var coverage = flag.Float64("c", 2.0, "Coverage")
	var error_rate = flag.Float64("e", 0.01, "Error rate.")
	flag.BoolVar(&Debug, "debug", false, "Turn on debug mode.")
	flag.Parse()
	read_len := position(*rl)
	if *seq_file != "" {
		if *coverage > 0 && read_len > 0 {
			// result := make(chan []position, 100000)
			idx := Build(*seq_file)
         num_of_reads := int(*coverage * float64(idx.LEN) / float64(read_len))
			read_indices := make([][]position, num_of_reads)
			the_read := make([]byte, read_len)

			for i:=0; i<num_of_reads; i++ {
				read_indices[i] = idx.Search(position(rand_gen.Intn(int(idx.LEN - read_len))), read_len)
			}
			for i:=0; i<num_of_reads; i++ {
            var errors []int
            copy(the_read, SEQ[read_indices[i][0]: read_indices[i][0] + read_len])
            for k:=0; k<len(the_read); k++ {
               if rand_gen.Float64() < *error_rate {
                  the_read[k] = random_error(the_read[k])
                  errors = append(errors, k)
               }
            }
            if Debug {
	            for j:=0; j<int(read_indices[i][0]); j++ {
   	            fmt.Printf(" ")
      	      }
      	   }
				fmt.Printf("%s %d ", the_read, len(read_indices[i]))
				for j:=0; j<len(read_indices[i]); j++ {
					fmt.Printf("%d ", read_indices[i][j])
				}
            fmt.Printf("%d ", len(errors))
            for j:=0; j<len(errors); j++ {
               fmt.Printf("%d ", errors[j])
            }
				fmt.Println()
			}
		} else {
			idx := Build(*seq_file)
			idx.Save(*seq_file)
		}
	} else {
		fmt.Println("Must provide sequence file")
	}
}