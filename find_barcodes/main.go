package main

import (
  "os"
  "github.com/biogo/biogo/seq/linear"
  "github.com/biogo/biogo/alphabet"
  "github.com/biogo/biogo/io/seqio"
  "fmt"
  "github.com/biogo/biogo/io/seqio/fasta"
  "github.com/biogo/biogo/io/seqio/fastq"
)

// Hamming computes the hamming distance between the primer string and the read
// starting at the specified index for the read
func Hamming(primer *linear.Seq, read *linear.QSeq, index int) int {
  distance := 0
	for i, x := range primer.Seq {
		if x != read.Seq[i + index].L {
		  distance++
		}
	}
	return distance
}


func main() {
  // Open the fastq file specified on the command line
  // for reading:
  fh, err := os.Open(os.Args[1])
  // Check for open errors and abort:
  if err != nil {
	panic(err)
  }

  // Create a template sequence for the reader:
  template := linear.NewQSeq("", alphabet.QLetters{}, alphabet.DNA, alphabet.Sanger)
  // Create a fastq reader:
  read_scanner := seqio.NewScanner(fastq.NewReader(fh, template))

  // Open the primer fasta file specified on the command line
  // for reading:
  fh, err = os.Open(os.Args[2])
  // Check for open errors and abort:
  if err != nil {
	panic(err)
  }

  // Create a template sequence for the reader:
  primer_template := linear.NewSeq("", alphabet.Letters{}, alphabet.DNA)
  // Create a fastq reader:
  primer_scanner := seqio.NewScanner(fasta.NewReader(fh, primer_template))

  primer_seqs := make([]*linear.Seq, 0)

  barcode_start := 7
  barcode_end := 13

  // Open the output file for writing:
  fho, err := os.Create(os.Args[3])
  // Close the file after we finished writing:
  defer fho.Close()
  if err != nil {
	panic(err)
  }
  // Create a fastq writer for the new sequences
  //writer := fastq.NewWriter(fho)

  for primer_scanner.Next(){
    seq := primer_scanner.Seq().(*linear.Seq)
	primer_seqs = append(primer_seqs, seq)
  }

  for read_scanner.Next(){
	read := read_scanner.Seq().(*linear.QSeq)
	annotation := read.CloneAnnotation()

	for i := range(primer_seqs) {
	  primer := primer_seqs[i]
	  for j := barcode_start ; j < barcode_end ; j++ {
		distance := Hamming(primer, read, j)
		if distance <= 1 {
		  fmt.Println("Hamming distance for", annotation.ID, primer.Seq, ":", distance)
		  //  barcode := read.Seq[0:(j-1)].L
		//  addition := strings.Join("BARCODE:"
		//  // Something like this
		//  id_parts := []string{annotation.ID, replaced}
		//  new_id := strings.Join(id_parts, " ")
		//  seq_replaced := linear.NewQSeq(new_id, read.Seq, read.Alphabet(), read.Encode)
		//  writer.Write(seq_replaced)
		}
	  }
	}

  }







  // Create a fasta writer with width 80:
  //writer := fastq.NewWriter(fho)

  //var i uint64 = 0
  var sequence_description string


  err = read_scanner.Error()
  // handle errors
  if err != nil {
    fmt.Println("last sequence:", sequence_description)
	panic(err)
  }

}
