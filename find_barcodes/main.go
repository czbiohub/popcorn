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
  "github.com/urfave/cli"
)

// Hamming computes the hamming distance between the primer string and the read
// starting at the specified index for the read
func Hamming(primer *linear.Seq, read *linear.QSeq, index int) int {
  distance := 0
  for i, x := range primer.Seq {
	if x != read.Seq[i+index].L {
	  distance++
	}
  }
  return distance
}

func main() {
  app := cli.NewApp()
  app.Name = "find_barcodes"
  app.Usage = "Extract barcodes from sequence given known primers and expected positions.\n\n   Example: " +
    "./find_barcodes --primers primers_R1.fasta reads_R1.fastq"
  app.EnableBashCompletion = true
  app.Authors = []cli.Author{{"Olga Botvinnik", "olga.botvinnik@czbiohub.org"},
  {"David Dynerman", "david.dynerman@czbiohub.org"}}

  var primerFilename string
  var primerStart int
  var primerEnd int
  var keyValueSep string
  var pairSep string
  var maxHammingDistance int

  app.Flags = []cli.Flag{
	cli.StringFlag{
	  Name:        "primers",
	  Usage:       "Fasta-formatted file of the primers to use",
	  Destination: &primerFilename,
	},
	//cli.StringFlag{
	//  Name: "output",
	//  Usage: "File to output the barcode-annotated reads to, default [inputfile].barcoded. if '-', will output " +
	//	"to standard output stream (stdout)",
	//},
	cli.IntFlag{
	  Name:  "primerStart",
	  Value: 6,
	  Usage: "0-based index of left side range where to start looking for the primer in the read. To account for " +
		"sequencing errors and indels, a good value is 2 fewer than your shortest expected barcode length, e.g. if " +
		"your shortest barcodes are 8, this should be 6",
	  Destination: &primerStart,
	},
	cli.IntFlag{
	  Name:  "primerEnd",
	  Value: 15,
	  Usage: "0-based index of right side range where to start looking for the primer in the read. To account for " +
		"seuqencing errors and indels, a good value is 3 more (since this is 0-based, have to add one) than your " +
	    "longest expected barcode length, e.g. if " +
		"your longest barcodes are 12, this should be 15",
	  Destination: &primerEnd,
	},
	cli.StringFlag{
	  Name:        "keyValueSep",
	  Value:       ":",
	  Usage:       "Character to use to separate keys and values, default ':', e.g. BARCODELENGTH:11",
	  Destination: &keyValueSep,
	},
	cli.StringFlag{
	  Name:  "pairSep",
	  Value: "_",
	  Usage: "Character to use to separate pairs of keys and values, default '_', e.g. " +
		"BARCODELENGTH:11_BARCODE:GCAGCCAACGG",
	  Destination: &pairSep,
	},
	cli.IntFlag{
	  Name:  "maxHammingDistance",
	  Value: 1,
	  Usage: "Maximum hamming distance from primers to reads. If the distance is greater than this, then the read will " +
		"be reported as NOPRIMERMATCH",
	  Destination: &maxHammingDistance,
	},
  }
  sort.Sort(cli.FlagsByName(app.Flags))
  //fmt.Println("os.Args:", os.Args)

  app.Action = func(c *cli.Context) error {
	fastqFilename := c.Args().First()

	// Open the fastq file specified on the command line
	// for reading:
	fh, err := os.Open(fastqFilename)
	// Check for open errors and abort:
	if err != nil {
	  panic(err)
	}

	// Create a template sequence for the reader:
	template := linear.NewQSeq("", alphabet.QLetters{}, alphabet.DNA, alphabet.Sanger)
	// Create a fastq reader:
	read_scanner := seqio.NewScanner(fastq.NewReader(fh, template))



	//var i uint64 = 0
	var readID string
	barcodeLengthCount := make(map[int]int)
	primerCount := make(map[string]int)

	barcodeLengthCountNoMatch := make(map[int]int)
	primerCountNoMatch := make(map[string]int)

	// Open the primer fasta file specified on the command line
	// for reading:
	primerSeqs := formatPrimers(primerFilename)
	//fmt.Println(primerSeqs)

	barcodedFastq := fastqFilename + ".barcoded"
	fmt.Println("Writing output barcoded fastq:", barcodedFastq)
	fho, err := os.Create(barcodedFastq)
	defer fho.Close()
	if err != nil {
	  panic(err)
	}
	barcodedFastqWriter := fastq.NewWriter(fho)

	noPrimerMatchFastq := fastqFilename + ".noPrimerMatch"
	fmt.Println("Writing output no primer match fastq:", noPrimerMatchFastq)
	fho, err = os.Create(noPrimerMatchFastq)
	defer fho.Close()
	if err != nil {
	  panic(err)
	}
	noPrimerMatchFastqWriter := fastq.NewWriter(fho)


	for read_scanner.Next() {
	  read := read_scanner.Seq().(*linear.QSeq)
	  annotation := read.CloneAnnotation()
	  readID = annotation.ID
	  //fmt.Println(readID)

	  bestDistance := read.Len()
	  var bestPrimer linear.Seq
	  var bestJ int

	  for i := range primerSeqs {
		primer := primerSeqs[i]
		for j := primerStart; j < primerEnd; j++ {
		  distance := Hamming(primer, read, j)
		  if distance <= bestDistance {
			bestDistance = distance
			bestPrimer = *primer
			bestJ = j
		  }
		}
	  }
	  //fmt.Printf("%s\t%s\t%s\thamming: %d\tindex (j): %d\n", annotation.ID, bestPrimer.Seq, bestPrimer.ID, bestDistance, bestJ)

	  barcode := read.String()[:bestJ]
	  //fmt.Println(annotation.ID, "+", bestPrimer.Seq, bestPrimer.ID, "hamming distance:", bestDistance, "index (j):", bestJ)

	  info := map[string]string{"PRIMER": bestPrimer.ID, "PRIMERHAMMINGDISTANCE": strconv.Itoa(bestDistance),
		"BARCODELENGTH": strconv.Itoa(bestJ), "BARCODE": barcode}
	  keyOrder := []string{"PRIMER", "PRIMERHAMMINGDISTANCE", "BARCODELENGTH", "BARCODE"}
	  var infoList []string
	  for _, key := range keyOrder {
		value := info[key]
		infoList = append(infoList, strings.Join([]string{key, value}, keyValueSep))
	  }
	  infoStr := strings.Join(infoList, pairSep)

	  // Don't report a primer if the best distance isn't too big (default 1 or less)
	  //var infoStr string
	  var writer *fastq.Writer

	  if bestDistance <= maxHammingDistance {
		barcodeLengthCount[bestJ]++
		primerCount[bestPrimer.ID]++

		writer = barcodedFastqWriter
	  } else {
		// Count no primer matching within distance as barcode length 0 and primer name "NA"
		barcodeLengthCountNoMatch[bestJ]++
		primerCountNoMatch[bestPrimer.ID]++

		//infoStr = fmt.Sprintf("NOPRIMERMATCH%sPRIMERHAMMINGDISTANCE%s%d", pairSep, keyValueSep, bestDistance)
		writer = noPrimerMatchFastqWriter

	  }

	  //  barcode := read.Seq[0:(j-1)].L
	  //info := fmt.Sprintf("PRIMER%[1]%s%s%[2]sPRIMERHAMMINGDISTANCE%[1]%s%d%[2]sBARCODELENGTH%[1]%s%d%[2]sBARCODE%[1]%s%s",
	  //  keyValueSep, pairSep, bestPrimer.ID, bestDistance, bestJ, barcode)
	  id_parts := []string{annotation.ID, infoStr}
	  newID := strings.Join(id_parts, "_")
	  newRead := linear.NewQSeq(newID, read.Seq, read.Alphabet(), read.Encode)
	  writer.Write(newRead)
	}

	writeBarcodeInfo(barcodedFastq+".barcodeLengths.csv", barcodeLengthCount)
	writePrimerInfo(barcodedFastq+".primers.csv", primerCount)

	writeBarcodeInfo(noPrimerMatchFastq+".barcodeLengths.csv", barcodeLengthCountNoMatch)
	writePrimerInfo(noPrimerMatchFastq+".primers.csv", primerCountNoMatch)



	err = read_scanner.Error()
	// handle errors
	if err != nil {
	  fmt.Fprintln(os.Stderr, "last read before error:", readID)
	  panic(err)
	}
	return nil
  }

  err := app.Run(os.Args)

  checkError("Could not run CLI App", err)

}

func formatPrimers(primerFilename string) []*linear.Seq {
  fh, err := os.Open(primerFilename)
  // Check for open errors and abort:
  if err != nil {
	panic(err)
  }

  // Create a template sequence for the reader:
  primerTemplate := linear.NewSeq("", alphabet.Letters{}, alphabet.DNA)
  // Create a fastq reader:
  primerScanner := seqio.NewScanner(fasta.NewReader(fh, primerTemplate))

  primerSeqs := make([]*linear.Seq, 0)

  for primerScanner.Next() {
	seq := primerScanner.Seq().(*linear.Seq)

	// Strip all Ns from left and right
	seqString := seq.Seq.String()
	stripped := strings.Trim(seqString, "N")
	letters := []alphabet.Letter(stripped)
	newSeq := linear.NewSeq(seq.ID, letters, seq.Alphabet())
	primerSeqs = append(primerSeqs, newSeq)
	//fmt.Fprintf(os.Stderr,"seq.Seq: %s\t newSeq.Seq: %s\n", seq.Seq, newSeq.Seq)
  }
  return primerSeqs
}

func writeBarcodeInfo(csvFileName string, barcodeLengthCount map[int]int) {
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

func writePrimerInfo(csvFileName string, primerCount map[string]int) {

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
