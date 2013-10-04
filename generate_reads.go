// Copyright 2013 Vinhthuy Phan
// Usage:  go run generate_reads.go --help
// Expected input format: sequence file should have no newlines.  Remove newlines from Fasta files.
// Example output format:
/*
	TAGGA	25 13
	ATCAA	26
	...

	This means TAGGA matches approximately with SEQ[25:30] and SEQ[13:18],  note: read length=5
	The matching is approximate because there are errors.
	It is possible that TAGGA matches approximately with other substrings aside from these.
	Algorithm:
		Let j be a random index
		s = SEQ[j : j+read_len] (e.g. s=ATCAA)
		report all occurence positions of s in SEQ, e.g. 25, 13.
		insert errors randomly in s
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
	// "bufio"
	"path"
	"runtime"
	"math/rand"
	"time"
)

var Debug bool

//-----------------------------------------------------------------------------
// Global variables: sequence (SEQ), suffix array (SA), BWT, FM index (C, OCC)
//-----------------------------------------------------------------------------
var SEQ []byte
var rand_gen = rand.New(rand.NewSource(time.Now().UnixNano()))

/*
type List struct {
	data []int
}
SA = new(List)
OCC = make([]List, len(SYMBOLS))
type BWT struct {
	data []byte
}
*/
type Index struct{
	SA []int 						// suffix array
	C map[byte]int  				// count table
	OCC map[byte][]int 			// occurence table

	END_POS int 					// position of "$" in the text
	SYMBOLS []int  				// sorted symbols
	EP map[byte]int 				// ending row/position of each symbol

	LEN int
	// un-exported variables
	bwt []byte
	freq map[byte]int  // Frequency of each symbol
}
//
//-----------------------------------------------------------------------------
type BySuffix []int

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
		log.Fatal("Load error:", filename, err)
	}
}

//-----------------------------------------------------------------------------
func _load_occ(filename string, Len int) []int {
	thing := make([]int, Len)
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

	I.OCC = make(map[byte][]int)
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
	I.freq = make(map[byte]int)
	I.SA = make([]int, len(SEQ))
	for i := 0; i < len(SEQ); i++ {
		I.SA[i] = i
		I.freq[SEQ[i]]++
	}
	sort.Sort(BySuffix(I.SA))

	I.bwt = make([]byte, len(SEQ))
	for i := 0; i < len(I.SA); i++ {
		I.bwt[i] = SEQ[(len(SEQ)+I.SA[i]-1)%len(SEQ)]
		if I.bwt[i] == '$' {
			I.END_POS = i
		}
	}
}

//-----------------------------------------------------------------------------
func (I *Index) BuildIndex() {
	I.C = make(map[byte]int)
	I.OCC = make(map[byte][]int)
	I.EP = make(map[byte]int)
	I.LEN = len(SEQ)
	for c := range I.freq {
		I.SYMBOLS = append(I.SYMBOLS, int(c))
		I.OCC[c] = make([]int, len(SEQ))
		I.C[c] = 0
	}
	sort.Ints(I.SYMBOLS)
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

func (I *Index) Search(j int, read_len int, result chan []int) {
	var sp, ep, offset int
	var ok bool

	c := SEQ[j+read_len-1]
	sp, ok = I.C[c]
	if ! ok {
		result <- make([]int, 0)
		return
	}
	ep = I.EP[c]
	// if Debug { fmt.Println("pattern: ", string(pattern), "\n\t", string(c), sp, ep) }
	for i:= read_len-2; sp <= ep && i >= 0; i-- {
  		c = SEQ[j+i]
  		offset, ok = I.C[c]
  		if ok {
			sp = offset + I.OCC[c][sp - 1]
			ep = offset + I.OCC[c][ep] - 1
		} else {
			result <- make([]int, 0)
			return
		}
  		// if Debug { fmt.Println("\t", string(c), sp, ep) }
	}
	res := make([]int, ep-sp+1)
	for i:=sp; i<=ep; i++ {
		res[i-sp] = I.SA[i]
	}
 	result <- res
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
	byte_array, err := ioutil.ReadFile(file)
	if err != nil {
		panic(err)
	}
	SEQ = append(bytes.Trim(byte_array, "\n\r "), byte('$'))
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
	not_A := []byte{'C','G','T'}
	not_T := []byte{'C','G','A'}
	not_C := []byte{'A','G','T'}
	not_G := []byte{'C','A','T'}
	if base == 'A' {
		return not_A[rand_gen.Intn(3)]
	} else if base == 'C' {
		return not_C[rand_gen.Intn(3)]
	} else if base == 'G' {
		return not_G[rand_gen.Intn(3)]
	} else {
		return not_T[rand_gen.Intn(3)]
	}
}

//-----------------------------------------------------------------------------
func main() {
	var seq_file = flag.String("seq", "", "Specify a file containing the sequence.")
	var read_len = flag.Int("len", 0, "Read length.")
	var reads = flag.Int("reads", 0, "Number of reads.")
	var error_rate = flag.Float64("erate", 0.01, "Error rate (uniformly distributed).")
	var workers = flag.Int("w", 1, "number of workers")
	flag.BoolVar(&Debug, "debug", false, "Turn on debug mode.")
	flag.Parse()

	if *seq_file != "" {
		if *reads > 0 && *read_len > 0 {
			runtime.GOMAXPROCS(*workers)
			result := make(chan []int, 100000)
			ReadSequence(*seq_file)
			idx := Load(*seq_file + ".index")
			for i:=0; i<*reads; i++ {
				// fmt.Println(j, *read_len, idx.LEN)
				// fmt.Printf("%d\t%s\n", j, string(SEQ[j:j+*read_len]))
				go idx.Search(rand_gen.Intn(idx.LEN - *read_len), *read_len, result)
			}
			for i:=0; i<*reads; i++ {
				indices := <- result
				the_read := SEQ[indices[0]: indices[0] + *read_len]
				var err_pos int
				// fmt.Printf("%s\t", the_read)
				for e:=0; e < int(*error_rate * float64(*read_len)); e++ {
					err_pos = rand_gen.Intn(len(the_read))
					the_read[err_pos] = random_error(the_read[err_pos])
					// fmt.Printf("%d ", err_pos)
				}
				fmt.Printf("%s\t", the_read)
				for j:=0; j<len(indices); j++ {
					fmt.Printf("%d ", indices[j])
				}
				fmt.Println()
			}
		} else {
			idx := Build(*seq_file)
			idx.Save(*seq_file)
		}
	}
}