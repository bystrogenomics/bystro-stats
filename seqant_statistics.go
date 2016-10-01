// estab exports elasticsearch fields as tab separated values
package main

import (
	"bufio"
	// "bytes"
	// "encoding/csv"
	// "encoding/json"
	"flag"
	// "fmt"
	"io"
	// "io/ioutil"
	"log"
	"math"
	"os"
	// "sort"
	"strconv"
	"strings"
)

type jsonFloat float64

func (value jsonFloat) MarshalJSON() ([]byte, error) {
	if math.IsNaN(float64(value)) {
		// Can't use "NA", get json: error calling MarshalJSON for type main.jsonFloat: invalid character 'N' looking for beginning of value
		return []byte("null"), nil
	}

	return []byte(strconv.FormatFloat(float64(value), 'f', -1, 64)), nil
}

// NOTE: For now this only supports \n end of line characters
// If we want to support CLRF or whatever, use either csv package, or set a different delimiter
func main() {
	inputFilePath := flag.String("inPath", "", "The input file path (optional: default is stdin)")
	// outputJsonPath := flag.String("outputJsonPath", "", "The output path for the JSON output (optional)")
	// outputTabPath := flag.String("outputTabPath", "", "The output path for tab-delimited file")
	// outputQcTabPath := flag.String("outputQcTabPath", "", "The output path for tab-delimited quality control file")
	referenceColumnIdx := flag.Int("referenceColumnIdx", -9, "The reference base column number")
	alleleColumnIdx := flag.Int("alleleColumnIdx", -9, "The allele column number")
	homozygotesColumnIdx := flag.Int("homozygotesColumnIdx", -9,
		"The homozygous sample column index")
	heterozygotesColumnIdx := flag.Int("heterozygotesColumnIdx", -9,
		"The heterozygous sample column index")
	siteTypeColumnIdx := flag.Int("siteTypeColumnIdx", -9, "The site type column index")
	dbSNPnameColumnIdx := flag.Int("dbSNPnameColumnIdx", -9, "The dbSNP name column index")
	exonicAlleleFunctionColumnIdx := flag.Int("exonicAlleleFunctionColumnIdx", -9,
		"The exonicAlleleFunction column index")
	fieldSeparator := flag.String("fieldSeparator", "\t", "What is used to delimit fields (',', '\t', etc)")
	primaryDelimiter := flag.String("primaryDelimiter", ";",
		"The primary delmiter (1D array string representation separator in input file)")
	emptyFieldString := flag.String("emptyFieldString", "NA",
		"What is used to denoted an empty field (NA by default)")
	secondaryDelimiter := flag.String("secondaryDelimiter", "|",
		"The secondary delmiter (2D array string representation outer separator in input file)")
	numberInputHeaderLines := flag.Int("numberInputHeaderLines", 0, "How many header lines does your input file have (0 is possible)")
	flag.Parse()

	inFh := (*os.File)(nil)

	if *inputFilePath != "" {
		var err error
		inFh, err = os.Open(*inputFilePath)

		if err != nil {
			log.Fatal(err)
		}
	} else {
		inFh = os.Stdin
	}

	// make sure it gets closed
	defer inFh.Close()

	if *referenceColumnIdx == -9 {
		log.Fatal("referenceColumnIdx required")
		os.Exit(1)
	}

	if *alleleColumnIdx == -9 {
		log.Fatal("alleleColumnIdx required")
		os.Exit(1)
	}

	if *homozygotesColumnIdx == -9 {
		log.Fatal("homozygotesColumnIdx required")
		os.Exit(1)
	}

	if *heterozygotesColumnIdx == -9 {
		log.Fatal("heterozygotesColumnIdx required")
		os.Exit(1)
	}

	if *siteTypeColumnIdx == -9 {
		log.Fatal("siteTypeColumnIdx required")
		os.Exit(1)
	}

	if *exonicAlleleFunctionColumnIdx == -9 {
		log.Fatal("exonicAlleleFunctionColumnIdx required")
		os.Exit(1)
	}

	// Beacuse I don't know how to pass a rune...
	if *fieldSeparator != "\t" && *fieldSeparator != "," {
		log.Fatal("fieldSeparator must be a comma or tab")
		os.Exit(1)
	}

	transitionsMap := map[string]map[string]bool{
		"A": map[string]bool{"G": true},
		"G": map[string]bool{"A": true},
		"T": map[string]bool{"C": true},
		"C": map[string]bool{"T": true},
	}

	transversionsMap := map[string]map[string]bool{
		"A": map[string]bool{
			"T": true,
			"C": true,
		},
		"G": map[string]bool{
			"T": true,
			"C": true,
		},
		"T": map[string]bool{
			"A": true,
			"G": true,
		},
		"C": map[string]bool{
			"A": true,
			"G": true,
		},
	}

	// var isTransition bool
	// var isTransversion bool

	// discordantCount := 0
	// weirdSites := 0
	// insCount := 0
	// delCount := 0
	// multiCount := 0
	// snpCount := 0
	// fields := []int{*referenceColumnIdx, *alleleColumnIdx}

	// Expects to place:
	// total => transitions => n
	// total => transversions => n
	// total => trTvRatio => n.y
	// sampleId => siteType transitions => n
	// sampleId => siteType transversions => n
	// total => siteType transitions => ...
	// trCountMap := make(map[string]map[string]int)
	//  tvCountMap := make(map[string]map[string]int)

	// sampleId|total : siteType|total|exonAlleleFunc = N
	trMap := make(map[string]map[string]int)
	tvMap := make(map[string]map[string]int)

	// sampleId|total : siteType|total|exonAlleleFunc = Y
	// ratioMap := make(map[string]map[string]jsonFloat)

	// trTvArray will hold all of the ratios for total trTv
	// trTvRatioArray := []float64{}

	// this one contains both counts and ratios, and is what we put into the return json
	// sampleId|total : "siteType|total|exonAlleleFunc transitions|transversions|ratio" = Z
	// allMap := make(map[string]map[string]interface{})

	totalKey := "total"
	// trKey := "transitions"
	// tvKey := "transversions"
	// trTvRatioKey := "transitions:transversions ratio"
	// trTvRatioMeanKey := "transitions:transversions ratio mean"
	// trTvRatioMedianKey := "transitions:transversions ratio median"
	// trTvRatioModeKey := "transitions:transversions ratio mode"
	// trTvRatioStdDevKey := "transitions:transversions ratio standard deviation"

	// dbSnpKey := "_in_dbSNP"

	// if user gives us dbSNP
	hasDbSnpColumn := *dbSNPnameColumnIdx != -9
	hasDbSnp := false

	rowCount := 0

	reader := bufio.NewReader(inFh)

	trMap[totalKey] = map[string]int{totalKey: 0}
	tvMap[totalKey] = map[string]int{totalKey: 0}

	fillArrayFunc := makeFillArrayFunc(*emptyFieldString, *secondaryDelimiter, *primaryDelimiter)

	nonNullFunc := makeHasNonEmptyRecordFunc(*emptyFieldString, *secondaryDelimiter, *primaryDelimiter)
	for {
		// http://stackoverflow.com/questions/8757389/reading-file-line-by-line-in-go
		// http://www.jeffduckett.com/blog/551119d6c6b86364cef12da7/golang---read-a-file-line-by-line.html
		// Scanner doesn't work well, has buffer restrictions that we need to manually get around
		// and we don't expect any newline characters in a Seqant output body
		row, err := reader.ReadString('\n') // 0x0A separator = newline

		if err == io.EOF {
			// do something here
			break
		} else if err != nil {
			log.Fatal(err)
		}

		if row == "NA" {
			println("blah")
		}

		// // remove the trailing \n
		// // equivalent of chomp https://groups.google.com/forum/#!topic/golang-nuts/smFU8TytFr4
		row = row[:len(row)-1]

		rowCount++

		if rowCount <= *numberInputHeaderLines {
			continue
		}

		record := strings.Split(row, *fieldSeparator)

		if len(record[*referenceColumnIdx]) > 1 || len(record[*alleleColumnIdx]) > 1 {
			continue
		}

		// isTransition = false
		// isTransversion = false

		if transitionsMap[record[*referenceColumnIdx]][record[*alleleColumnIdx]] == true {
			// isTransition = true
			// if record[2] != "SNP" {
			//  weirdSites++
			// }
		} else if transversionsMap[record[*referenceColumnIdx]][record[*alleleColumnIdx]] == true {
			// isTransversion = true
			// if record[2] != "SNP" {
			//  weirdSites++
			// }
		}
		//   else {
		//  if record[2] == "SNP" {
		//    discordantCount++
		//  }
		// }

		if hasDbSnpColumn {
			hasDbSnp = nonNullFunc(record[*dbSNPnameColumnIdx])
		}

		seenSamples := make(map[string]struct{})
		samples := []string{totalKey}

		hetSamples := fillArrayFunc(record[*heterozygotesColumnIdx], seenSamples)

		if hetSamples != nil {
			samples = append(samples, hetSamples...)
		}

		homSamples := fillArrayFunc(record[*homozygotesColumnIdx], seenSamples)

		if homSamples != nil {
			samples = append(samples, homSamples...)
		}

		// // if record[2] == "MULTIALLELIC" {
		// //   multiCount++
		// // } else if record[2] == "DEL" {
		// //   delCount++
		// // } else if record[2] == "INS" {
		// //   insCount++
		// // } else if record[2] == "SNP" {
		// //   snpCount++
		// // }
		// // Make siteTypes, exonicTypes, and samples unique; struct{} takes up no additional space
		// // http://stackoverflow.com/questions/9251234/go-append-if-unique
		// Samples don't really need to be made this way, convenience out of re-using
		// function for siteTypes, exonicAlleleFunctions, and sampels
		seenValues := make(map[string]struct{})
		siteTypes := fillArrayFunc(record[*siteTypeColumnIdx], seenValues)

		if siteTypes == nil {
			siteTypes = fillArrayFunc(record[*exonicAlleleFunctionColumnIdx], seenValues)
		} else {
			siteTypes = append(siteTypes, fillArrayFunc(record[*exonicAlleleFunctionColumnIdx], seenValues)...)
		}

		for _, sample := range samples {
			// if _, exists := trMap[sample]; !exists {
			// 	trMap[sample] = make(map[string]int)
			// 	tvMap[sample] = make(map[string]int)
			// }

			// if isTransition == true {
			// 	trMap[sample][totalKey]++
			// } else if isTransversion == true {
			// 	tvMap[sample][totalKey]++
			// }
			if sample != "" {

			}

			for _, siteType := range siteTypes {
				// if _, exists := trMap[sample][siteType]; !exists {
				// 	trMap[sample][siteType] = 0
				// 	tvMap[sample][siteType] = 0
				// }

				if hasDbSnp {

				}

				if siteType != "" {

				}
				// if hasDbSnp == true {
				// 	var name bytes.Buffer
				// 	// siteType_in_dbSNP
				// 	name.WriteString(siteType)
				// 	name.WriteString(dbSnpKey)

				// 	siteName := name.String()

				// 	if _, exists := trMap[sample][siteName]; !exists {
				// 		trMap[sample][siteName] = 0
				// 		tvMap[sample][siteName] = 0
				// 	}

				// 	if isTransition == true {
				// 		trMap[sample][siteName]++
				// 	} else if isTransversion == true {
				// 		tvMap[sample][siteName]++
				// 	}
				// }

				// if isTransition == true {
				// 	trMap[sample][siteType]++
				// } else if isTransversion == true {
				// 	tvMap[sample][siteType]++
				// }
			}
		}
	}

	// for sampleId := range trMap {
	//  if _, exists := ratioMap[sampleId]; !exists {
	//    ratioMap[sampleId] = make(map[string]jsonFloat)
	//  }

	//  for siteType := range trMap[sampleId] {
	//    // If denominator is 0, NaN will result, which we will store as "NA"
	//    // https://github.com/raintank/metrictank/commit/5de7d6e3751901a23501e5fcd95f0b2d0604e8f4
	//    ratioMap[sampleId][siteType] = jsonFloat(trMap[sampleId][siteType]) / jsonFloat(tvMap[sampleId][siteType])
	//  }
	// }

	// for sampleId := range trMap {
	//  if _, exists := allMap[sampleId]; !exists {
	//    allMap[sampleId] = make(map[string]interface{})
	//  }

	//  for siteType := range trMap[sampleId] {
	//    var trName bytes.Buffer
	//    var tvName bytes.Buffer
	//    var ratioName bytes.Buffer

	//    trName.WriteString(siteType)
	//    trName.WriteString(" ")
	//    trName.WriteString(trKey)

	//    tvName.WriteString(siteType)
	//    tvName.WriteString(" ")
	//    tvName.WriteString(tvKey)

	//    ratioName.WriteString(siteType)
	//    ratioName.WriteString(" ")
	//    ratioName.WriteString(trTvRatioKey)

	//    allMap[sampleId][trName.String()] = trMap[sampleId][siteType]
	//    allMap[sampleId][tvName.String()] = tvMap[sampleId][siteType]
	//    allMap[sampleId][ratioName.String()] = ratioMap[sampleId][siteType]
	//  }
	// }

	// // conduct QC

	// // total names
	// var totalTrName bytes.Buffer
	// var totalTvName bytes.Buffer
	// var totalRatioName bytes.Buffer

	// totalTrName.WriteString(totalKey)
	// totalTrName.WriteString(" ")
	// totalTrName.WriteString(trKey)

	// totalTvName.WriteString(totalKey)
	// totalTvName.WriteString(" ")
	// totalTvName.WriteString(tvKey)

	// totalRatioName.WriteString(totalKey)
	// totalRatioName.WriteString(" ")
	// totalRatioName.WriteString(trTvRatioKey)

	// // numberOfSamples := len(allMap)

	// var sampleNames []string
	// var siteTypes []string
	// var totalSiteTypes []string

	// uniqueSites := make(map[string]struct{})
	// for sampleName, sampleHash := range allMap {
	//  if sampleName != totalKey {
	//    sampleNames = append(sampleNames, sampleName)
	//    trTvRatioArray = append(trTvRatioArray, float64(ratioMap[sampleName][totalKey]))
	//  }

	//  // We don't do this for the "total" sampleName because the totalKey
	//  // must contain only keys found in the sampleName
	//  for siteType := range sampleHash {
	//    if _, exists := uniqueSites[siteType]; exists {
	//      continue
	//    }

	//    // To skip anything we've seen
	//    uniqueSites[siteType] = struct{}{}

	//    // Handle totals differently, we want them at the end
	//    if strings.Contains(siteType, totalKey) {
	//      totalSiteTypes = append(totalSiteTypes, siteType)
	//    } else {
	//      siteTypes = append(siteTypes, siteType)
	//    }
	//  }
	// }

	// // sampleNames =
	// sort.Strings(sampleNames)
	// // We skipped the totalKey above, so that we may put it first
	// sampleNames = append([]string{totalKey}, sampleNames...)

	// sort.Strings(totalSiteTypes)
	// sort.Strings(siteTypes)

	// siteTypes = append(siteTypes, totalSiteTypes...)

	// trTvMean := mean(trTvRatioArray)
	// trTvMedian := median(trTvRatioArray)
	// trTvSd := stdDev(trTvRatioArray, trTvMean)

	// // fmt.Printf("Number of SNP %d\n", snpCount)
	// // fmt.Printf("Number of DEL %d\n", delCount)
	// // fmt.Printf("Number of INS %d\n", insCount)
	// // fmt.Printf("Number of MULTIALLELIC %d\n", multiCount)
	// // fmt.Printf("Number of Discordant %d\n", discordantCount)
	// // fmt.Printf("Number of Weird %d\n", weirdSites)
	// // fmt.Printf("Number of Lines %d\n", rowCount)

	// // fmt.Printf("Transition:Transversion Ratio mean %f\n", trTvMean)
	// // fmt.Printf("Transition:Transversion Ratio median %f\n", trTvMedian)
	// // fmt.Printf("Transition:Transversion Ratio standard deviation %f\n", trTvSd)

	// allMap["qc"] = map[string]interface{}{
	//  trTvRatioMeanKey: jsonFloat(trTvMean), trTvRatioMedianKey: jsonFloat(trTvMedian),
	//  trTvRatioStdDevKey: jsonFloat(trTvSd),
	// }

	// if *outputJsonPath != "" {
	//  json, err := json.Marshal(allMap)

	//  if err != nil {
	//    log.Fatal(err)
	//  }

	//  err = ioutil.WriteFile(*outputJsonPath, json, os.FileMode(0644))

	//  if err != nil {
	//    log.Fatal(err)
	//  }
	// }

	// // Write output as a tabbed file
	// outFh := (*os.File)(nil)

	// if *outputTabPath != "" {
	//  var err error
	//  outFh, err = os.OpenFile(*outputTabPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0600)

	//  if err != nil {
	//    log.Fatal(err)
	//  }
	// } else {
	//  outFh = os.Stdout
	// }

	// defer outFh.Close()
	// writer := csv.NewWriter(outFh)

	// if *fieldSeparator == "\t" {
	//  writer.Comma = '\t'
	// }

	// // first column is for sample names
	// outLines := [][]string{append([]string{"\t"}, siteTypes...)}

	// for _, sampleName := range sampleNames {
	//  // First column is for the sample name
	//  line := []string{sampleName}

	//  for _, siteType := range siteTypes {
	//    switch value := allMap[sampleName][siteType].(type) {
	//    case int:
	//      line = append(line, strconv.Itoa(value))
	//    case jsonFloat:
	//      line = append(line, strconv.FormatFloat(float64(value), 'f', -1, 64))
	//    }
	//  }

	//  outLines = append(outLines, line)
	// }
	// writer.WriteAll(outLines)

	// // Write output as a tabbed file
	// outQcFh := (*os.File)(nil)
	// defer outQcFh.Close()

	// if *outputQcTabPath != "" {
	//  var err error
	//  outFh, err = os.OpenFile(*outputQcTabPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0600)

	//  if err != nil {
	//    log.Fatal(err)
	//  }
	// } else {
	//  outQcFh = os.Stdout
	// }

	// writer = csv.NewWriter(outQcFh)

	// if *fieldSeparator == "\t" {
	//  writer.Comma = '\t'
	// }

	// outQcLines := [][]string{
	//  []string{"Transition:Transversion Ratio Mean", strconv.FormatFloat(trTvMean, 'f', -1, 64)},
	//  []string{"Transition:Transversion Ratio Median", strconv.FormatFloat(trTvMedian, 'f', -1, 64)},
	//  []string{"Transition:Transversion Ratio Median", strconv.FormatFloat(trTvSd, 'f', -1, 64)},
	// }

	// // outQcLines = append(outQcLines, []string{ "Transition:Transversion Ratio Mean", trTvMean})
	// // outQcLines = append(outQcLines, "Transition:Transversion Ratio Median", trTvMedian)
	// // outQcLines = append(outQcLines, "Transition:Transversion Ratio Median", trTvSd)

	// writer.WriteAll(outQcLines)
}

// //https://github.com/dropbox/godropbox/blob/master/sort2/sort.go
// func mean(numbers []float64) (meanVal float64) {
//  return sum(numbers) / float64(len(numbers))
// }

// func sum(numbers []float64) (total float64) {
//  for _, x := range numbers {
//    total += x
//  }
//  return total
// }

// func median(numbers []float64) float64 {
//  middle := len(numbers) / 2
//  result := numbers[middle]
//  if len(numbers)%2 == 0 {
//    result = (result + numbers[middle-1]) / 2
//  }
//  return result
// }

// func mode(numbers []float64) (modes []float64) {
//  frequencies := make(map[float64]int, len(numbers))
//  highestFrequency := 0
//  for _, x := range numbers {
//    frequencies[x]++
//    if frequencies[x] > highestFrequency {
//      highestFrequency = frequencies[x]
//    }
//  }
//  for x, frequency := range frequencies {
//    if frequency == highestFrequency {
//      modes = append(modes, x)
//    }
//  }
//  if highestFrequency == 1 || len(modes) == len(numbers) {
//    modes = modes[:0] // Or: modes = []float64{}
//  }
//  sort.Float64s(modes)
//  return modes
// }

// func stdDev(numbers []float64, mean float64) float64 {
//  total := 0.0
//  for _, number := range numbers {
//    total += math.Pow(number-mean, 2)
//  }
//  variance := total / float64(len(numbers)-1)
//  return math.Sqrt(variance)
// }

func makeFillArrayFunc(emptyField string, secondaryDelim string, primaryDelim string) func(string, map[string]struct{}) []string {
	return func(record string, seen map[string]struct{}) []string {
		var out []string

		if record == emptyField {
			return out
		}

		if strings.Contains(record, secondaryDelim) {
			for _, val := range strings.Split(record, secondaryDelim) {

				if !strings.Contains(val, primaryDelim) {
					if val == emptyField {
						continue
					}

					if _, found := seen[val]; !found {
						seen[val] = struct{}{}
						out = append(out, val)
						continue
					}

					continue
				}

				for _, innerVal := range strings.Split(val, primaryDelim) {
					if innerVal == emptyField {
						continue
					}

					if _, found := seen[innerVal]; !found {
						seen[innerVal] = struct{}{}
						out = append(out, innerVal)
					}
				}
			}

			return out
		}

		if !strings.Contains(record, primaryDelim) {
			if record == emptyField {
				return out
			}

			if _, found := seen[record]; !found {
				seen[record] = struct{}{}

				// avoid an append
				return []string{record}
			}

			return out
		}

		for _, innerVal := range strings.Split(record, primaryDelim) {
			if innerVal != emptyField {
				if _, found := seen[innerVal]; !found {
					seen[innerVal] = struct{}{}
					out = append(out, innerVal)
				}
			}
		}

		return out
	}
}

func makeHasNonEmptyRecordFunc(emptyField string, secondaryDelim string, primaryDelim string) func(string) bool {
	return func(record string) bool {
		if record == emptyField {
			return false
		}

		if strings.Contains(record, secondaryDelim) {
			for _, val := range strings.Split(record, secondaryDelim) {

				if !strings.Contains(val, primaryDelim) {
					if val != emptyField {
						return true
					}
				}

				for _, innerVal := range strings.Split(val, primaryDelim) {
					if innerVal != emptyField {
						return true
					}
				}
			}

			return false
		}

		if !strings.Contains(record, primaryDelim) {
			if record != emptyField {
				return true
			}

			return false
		}

		for _, innerVal := range strings.Split(record, primaryDelim) {
			if innerVal != emptyField {
				return true
			}
		}

		return false
	}
}
