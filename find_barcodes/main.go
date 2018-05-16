package main

import (
  "os"
  "log"
  "github.com/biogo/biogo/seq/linear"
  "github.com/biogo/biogo/alphabet"
  "github.com/biogo/biogo/io/seqio"
  "fmt"
  "github.com/biogo/biogo/io/seqio/fasta"
  "github.com/biogo/biogo/io/seqio/fastq"
  "strings"
  "sort"
  "encoding/csv"
  "strconv"
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

  barcode_start := 6
  barcode_end := 15

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

    // Strip all Ns from left and right
    seqString := seq.Seq.String()
    stripped := strings.Trim(seqString, "N")
    letters := []alphabet.Letter(stripped)
    newSeq := linear.NewSeq(seq.ID, letters, seq.Alphabet())
	primer_seqs = append(primer_seqs, newSeq)
	//fmt.Fprintf(os.Stderr,"seq.Seq: %s\t newSeq.Seq: %s\n", seq.Seq, newSeq.Seq)
  }


  // Create a fasta writer with width 80:
  //writer := fastq.NewWriter(fho)

  //var i uint64 = 0
  var readID string
  barcodeLengthCount := make(map[int]int)
  primerCount := make(map[string]int)


  for read_scanner.Next(){
	read := read_scanner.Seq().(*linear.QSeq)
	annotation := read.CloneAnnotation()
	readID = annotation.ID

	bestDistance := read.Len()
	var bestPrimer linear.Seq
	var bestJ int

	for i := range(primer_seqs) {
	  primer := primer_seqs[i]
	  for j := barcode_start ; j < barcode_end ; j++ {
		distance := Hamming(primer, read, j)
		if distance <= bestDistance {
		  bestDistance = distance
		  bestPrimer = *primer
		  bestJ = j
		}
	  }
	}
	//fmt.Printf("%s\t%s\t%s\thamming: %d\tindex (j): %d\n", annotation.ID, bestPrimer.Seq, bestPrimer.ID, bestDistance, bestJ)

	barcodeLengthCount[bestJ]++
	primerCount[bestPrimer.ID]++
	//fmt.Println(annotation.ID, "+", bestPrimer.Seq, bestPrimer.ID, "hamming distance:", bestDistance, "index (j):", bestJ)

	writeBarcodeInfo(os.Args[1] + ".barcodeInfo.csv", barcodeLengthCount)
	writePrimerInfo(os.Args[1] + ".primerInfo.csv", primerCount)

	//  barcode := read.Seq[0:(j-1)].L
	//  addition := strings.Join("BARCODE:"
	//  // Something like this
	//  id_parts := []string{annotation.ID, replaced}
	//  new_id := strings.Join(id_parts, " ")
	//  seq_replaced := linear.NewQSeq(new_id, read.Seq, read.Alphabet(), read.Encode)
	//  writer.Write(seq_replaced)

  }



  err = read_scanner.Error()
  // handle errors
  if err != nil {
    fmt.Fprintln(os.Stderr,"last read before error:", readID)
	panic(err)
  }

}

func writeBarcodeInfo(csvFileName string, barcodeLengthCount map[int]int){
  // Create a separate, sorted list of barcode counts so the order is always the same
  var barcodeLengths []int
  for barcodeLength := range barcodeLengthCount {
	barcodeLengths = append(barcodeLengths, barcodeLength)
  }
  sort.Ints(barcodeLengths)

  // --- Write barcode info to file ---
  csvFile, csvErr := os.Create(csvFileName)
  if csvErr != nil {
	panic(csvErr)
  }
  defer csvFile.Close()

  csvWriter := csv.NewWriter(csvFile)
  defer csvWriter.Flush()
  header := []string{"barcode_length", "count"}
  err := csvWriter.Write(header)
  checkError("Cannot write to file", err)

  for _, barcodeLength := range barcodeLengths {
	//fmt.Fprintln(os.Stderr,"Barcode Length:", barcodeLength, "\tCount:", barcodeLengthCount[barcodeLength])
	data := []string{strconv.Itoa(barcodeLength),
	  strconv.Itoa(barcodeLengthCount[barcodeLength])}
	//line := stringify(data)
	err := csvWriter.Write(data)
	checkError("Cannot write to file", err)
  }
}

func writePrimerInfo(csvFileName string, primerCount map[string]int){

  // --- Write primer info to file ---
  csvFile, csvErr := os.Create(csvFileName)
  if csvErr != nil {
	panic(csvErr)
  }
  defer csvFile.Close()

  csvWriter := csv.NewWriter(csvFile)
  defer csvWriter.Flush()
  header := []string{"barcode_length", "count"}
  err := csvWriter.Write(header)
  checkError("Cannot write to file", err)

  for primerID, count := range primerCount {
	//fmt.Fprintf(os.Stderr,"%s\t%d\n", primerID, count)
	data := []string{primerID, strconv.Itoa(count)}
	//line := stringify(data)
	err := csvWriter.Write(data)
	checkError("Cannot write to file", err)
  }
}

//func stringify(intArray []int) []string {
//  stringArray := make([]string, len(intArray))
//  for i, v := range intArray {
//	stringArray[i] += strconv.Itoa(v)
//  }
//
//  return stringArray
//}

func checkError(message string, err error) {
  if err != nil {
	log.Fatal(message, err)
  }
}
