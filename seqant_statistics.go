// estab exports elasticsearch fields as tab separated values
package main

import (
	"bufio"
	"bytes"
	"encoding/csv"
	"encoding/json"
	"flag"
	// "fmt"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	// "github.com/davecgh/go-spew/spew"
	"runtime/pprof"
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
	outputJSONPath := flag.String("outputJSONPath", "", "The output path for the JSON output (optional)")
	outputTabPath := flag.String("outputTabPath", "", "The output path for tab-delimited file")
	outputQcTabPath := flag.String("outputQcTabPath", "", "The output path for tab-delimited quality control file")
	referenceColumnIdx := flag.Int("referenceColumnIdx", -9, "The reference base column index")
	referenceColumnName := flag.String("referenceColumnName", "",
		"The reference base column name. This is usually the name of the assembly")
	alleleColumnIdx := flag.Int("alleleColumnIdx", -9, "The allele column index")
	alleleColumnName := flag.String("alleleColumnName", "minorAlleles", "The alleles column name")
	homozygotesColumnIdx := flag.Int("homozygotesColumnIdx", -9,
		"The homozygous sample column index")
	homozygotesColumnName := flag.String("homozygotesColumnName", "homozygotes",
		"The homozygous sample column name")
	heterozygotesColumnIdx := flag.Int("heterozygotesColumnIdx", -9,
		"The heterozygous sample column index")
	heterozygotesColumnName := flag.String("heterozygotesColumnName", "heterozygotes",
		"The homozygous sample column name")
	siteTypeColumnIdx := flag.Int("siteTypeColumnIdx", -9, "The site type column index")
	siteTypeColumnName := flag.String("siteTypeColumnName", "refSeq.siteType", "The site type column name")
	dbSNPnameColumnIdx := flag.Int("dbSNPnameColumnIdx", -9, "Optional. The dbSNP name column index")
	dbSNPnameColumnName := flag.String("dbSNPnameColumnName", "dbSNP146.name", "Optional. The snp name column name")
	exonicAlleleFunctionColumnIdx := flag.Int("exonicAlleleFunctionColumnIdx", -9,
		"Optional. The exonicAlleleFunction column index")
	exonicAlleleFunctionColumnName := flag.String("exonicAlleleFunctionColumnName",
		"refSeq.codonEffect", `Optional. The name of the column that has exonicAlleleFunction, aka the one that has
    nonSynonymous, synonymous, etc values`)
	fieldSeparator := flag.String("fieldSeparator", "\t", "What is used to delimit fields (',', '\t', etc)")
	primaryDelimiter := flag.String("primaryDelimiter", ";",
		"The primary delmiter (1D array string representation separator in input file)")
	emptyFieldString := flag.String("emptyFieldString", "NA",
		"What is used to denoted an empty field (NA by default)")
	secondaryDelimiter := flag.String("secondaryDelimiter", "|",
		"The secondary delmiter (2D array string representation outer separator in input file)")
	numberInputHeaderLines := flag.Int("numberInputHeaderLines", 0, "How many header lines does your input file have (0 is possible)")
	cpuprofile := flag.String("cpuProfile", "", "write cpu profile to file")
	flag.Parse()

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

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

	if *referenceColumnIdx == -9 && (*referenceColumnName == "" || *numberInputHeaderLines != 1) {
		log.Fatal(`If referenceColumnIdx not provided, referenceColumnName must be, and
      the numberInputHeaderLines must equal 1`)
		os.Exit(1)
	}

	if *alleleColumnIdx == -9 && (*alleleColumnName == "" || *numberInputHeaderLines != 1) {
		log.Fatal(`If alleleColumnIdx not provided, alleleColumnName must be, and
      the numberInputHeaderLines must equal 1`)
		os.Exit(1)
	}

	if *homozygotesColumnIdx == -9 && (*homozygotesColumnName == "" || *numberInputHeaderLines != 1) {
		log.Fatal(`If homozygotesColumnIdx not provided, homozygotesColumnName must be, and
      the numberInputHeaderLines must equal 1`)
		os.Exit(1)
	}

	if *heterozygotesColumnIdx == -9 && (*heterozygotesColumnName == "" || *numberInputHeaderLines != 1) {
		log.Fatal(`If heterozygotesColumnIdx not provided, heterozygotesColumnName must be, and
      the numberInputHeaderLines must equal 1`)
		os.Exit(1)
	}

	if *siteTypeColumnIdx == -9 && (*siteTypeColumnName == "" || *numberInputHeaderLines != 1) {
		log.Fatal(`If siteTypeColumnIdx not provided, siteTypeColumnName must be, and
      the numberInputHeaderLines must equal 1`)
		os.Exit(1)
	}

	// if *exonicAlleleFunctionColumnIdx == -9 && (*exonicAlleleFunctionColumnName == "" || *numberInputHeaderLines != 1) {
	// 	log.Fatal(`If exonicAlleleFunctionColumnIdx not provided, exonicAlleleFunctionColumnName must be,
	//     and numberInputHeaderLines must be equal 1`)
	// 	os.Exit(1)
	// }

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

	var isTransition bool
	var isTransversion bool

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

	// sampleId|total : siteType|total|exonAlleleFunc = N
	trMap := make(map[string]map[string]int, 1000)
	tvMap := make(map[string]map[string]int, 1000)

	// sampleId|total : siteType|total|exonAlleleFunc = Y
	ratioMap := make(map[string]map[string]jsonFloat, 1000)

	totalKey := "total"
	trKey := "transitions"
	tvKey := "transversions"
	trTvRatioKey := "transitions:transversions ratio"
	trTvRatioMeanKey := "transitions:transversions ratio mean"
	trTvRatioMedianKey := "transitions:transversions ratio median"
	trTvRatioStdDevKey := "transitions:transversions ratio standard deviation"

	dbSnpKey := "_in_dbSNP"

	// if user gives us dbSNP
	hasDbSnpColumn := *dbSNPnameColumnIdx != -9
	hasDbSnp := false

	hasExonicColumn := *exonicAlleleFunctionColumnIdx != -9

	rowCount := 0

	reader := bufio.NewReader(inFh)

	trMap[totalKey] = make(map[string]int, 200)
	tvMap[totalKey] = make(map[string]int, 200)
	trMap[totalKey][totalKey] = 0
	tvMap[totalKey][totalKey] = 0

	dbSNPfeatureMap := make(map[string]string, 200)

	fillArrayFunc := makeFillArrayFunc(*emptyFieldString, *secondaryDelimiter, *primaryDelimiter)

	nonNullFunc := makeHasNonEmptyRecordFunc(*emptyFieldString, *secondaryDelimiter, *primaryDelimiter)

	// samples that are variant in a single row, capacity
	samples := make([]string, 1, 1000)
	samples[0] = totalKey

	siteTypes := make([]string, 0, 200)

	var record []string
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

		// // remove the trailing \n
		// // equivalent of chomp https://groups.google.com/forum/#!topic/golang-nuts/smFU8TytFr4
		row = row[:len(row)-1]

		if rowCount < *numberInputHeaderLines {
			rowCount++

			if rowCount == 1 && *numberInputHeaderLines == 1 {
				record = strings.Split(row, *fieldSeparator)

				if *referenceColumnName != "" && *referenceColumnIdx == -9 {
					*referenceColumnIdx = findIndex(record, *referenceColumnName)

					if *referenceColumnIdx == -9 {
						log.Fatal("referenceColumnName not found")
					}
				}

				if *alleleColumnName != "" && *alleleColumnIdx == -9 {
					*alleleColumnIdx = findIndex(record, *alleleColumnName)

					if *alleleColumnIdx == -9 {
						log.Fatal("alleleColumnName not found")
					}
				}

				if *heterozygotesColumnName != "" && *heterozygotesColumnIdx == -9 {
					*heterozygotesColumnIdx = findIndex(record, *heterozygotesColumnName)

					if *heterozygotesColumnIdx == -9 {
						log.Fatal("heterozygotesColumnName not found")
					}
				}

				if *homozygotesColumnName != "" && *homozygotesColumnIdx == -9 {
					*homozygotesColumnIdx = findIndex(record, *homozygotesColumnName)

					if *homozygotesColumnIdx == -9 {
						log.Fatal("homozygotesColumnName not found")
					}
				}

				if *siteTypeColumnName != "" && *siteTypeColumnIdx == -9 {
					*siteTypeColumnIdx = findIndex(record, *siteTypeColumnName)

					if *siteTypeColumnIdx == -9 {
						log.Fatal("siteTypeColumnNme not found")
					}
				}

				if *dbSNPnameColumnName != "" && *dbSNPnameColumnIdx == -9 {
					*dbSNPnameColumnIdx = findIndex(record, *dbSNPnameColumnName)

					if *dbSNPnameColumnIdx == -9 {
						hasDbSnpColumn = true
					}
				}

				if *exonicAlleleFunctionColumnName != "" && *exonicAlleleFunctionColumnIdx == -9 {
					*exonicAlleleFunctionColumnIdx = findIndex(record, *exonicAlleleFunctionColumnName)

					if *exonicAlleleFunctionColumnIdx != -9 {
						hasExonicColumn = true
					}
				}
			}

			continue
		}

		record = strings.Split(row, *fieldSeparator)

		// spew.Dump(record)
		if len(record[*alleleColumnIdx]) > 1 {
			continue
		}

		isTransition = false
		isTransversion = false

		if transitionsMap[record[*referenceColumnIdx]][record[*alleleColumnIdx]] == true {
			isTransition = true
			// if record[2] != "SNP" {
			//  weirdSites++
			// }
		} else if transversionsMap[record[*referenceColumnIdx]][record[*alleleColumnIdx]] == true {
			isTransversion = true
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

		// remove everything but the total key
		samples = samples[:1]

		// var allSamples bytes.Buffer
		// allSamples.WriteString(record[*heterozygotesColumnIdx])
		// allSamples.WriteString(*secondaryDelimiter)
		// allSamples.WriteString(record[*homozygotesColumnIdx])
		// samples = append(samples, fillArrayFunc(allSamples.String(), make(map[string]struct{}))...)
		samples = append(samples, fillArrayFunc(record[*heterozygotesColumnIdx], false, false)...)

		samples = append(samples, fillArrayFunc(record[*homozygotesColumnIdx], false, false)...)

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
		// function for siteTypes, exonicAlleleFunctions, and samples

		// var allSiteTypes bytes.Buffer
		// allSiteTypes.WriteString(record[*siteTypeColumnIdx])
		// allSiteTypes.WriteString(*secondaryDelimiter)
		// allSiteTypes.WriteString(record[*exonicAlleleFunctionColumnIdx])

		// seenValues := make(map[string]struct{})
		siteTypes = fillArrayFunc( /*allSiteTypes.String()*/ record[*siteTypeColumnIdx], true, true)

		if hasExonicColumn == true {
			if siteTypes == nil {
				siteTypes = fillArrayFunc(record[*exonicAlleleFunctionColumnIdx], true, true)
			} else {
				siteTypes = append(siteTypes, fillArrayFunc(record[*exonicAlleleFunctionColumnIdx], true, true)...)
			}
		}

		for _, sample := range samples {
			if _, exists := trMap[sample]; !exists {
				trMap[sample] = make(map[string]int, 200)
				tvMap[sample] = make(map[string]int, 200)
				// fmt.Printf("making map %s %s\n", sample, totalKey)
				// trMap[sample][totalKey] = 0
				// tvMap[sample][totalKey] = 0
			}

			if isTransition == true {
				trMap[sample][totalKey]++
			} else if isTransversion == true {
				tvMap[sample][totalKey]++
			}

			for _, siteType := range siteTypes {
				// if _, exists := trMap[sample][siteType]; !exists {
				//  // fmt.Printf("making map %s %s\n", sample, siteType)
				//  trMap[sample][siteType] = 0
				//  tvMap[sample][siteType] = 0
				// }

				if isTransition == true {
					trMap[sample][siteType]++
				} else if isTransversion == true {
					tvMap[sample][siteType]++
				}

				if hasDbSnp == true {
					if dbSNPfeatureMap[siteType] == "" {
						var name bytes.Buffer
						// siteType_in_dbSNP
						name.WriteString(siteType)
						name.WriteString(dbSnpKey)

						dbSNPfeatureMap[siteType] = name.String()

						trMap[sample][dbSNPfeatureMap[name.String()]] = 0
						tvMap[sample][dbSNPfeatureMap[name.String()]] = 0
					}

					if isTransition == true {
						trMap[sample][dbSNPfeatureMap[siteType]]++
					} else if isTransversion == true {
						tvMap[sample][dbSNPfeatureMap[siteType]]++
					}
				}
			}
		}
	}

	for sampleID := range trMap {
		if _, exists := ratioMap[sampleID]; !exists {
			ratioMap[sampleID] = make(map[string]jsonFloat, 100)
		}

		for siteType := range trMap[sampleID] {
			// If denominator is 0, NaN will result, which we will store as "NA"
			// https://github.com/raintank/metrictank/commit/5de7d6e3751901a23501e5fcd95f0b2d0604e8f4
			ratioMap[sampleID][siteType] = jsonFloat(trMap[sampleID][siteType]) / jsonFloat(tvMap[sampleID][siteType])
		}
	}

	trSiteTypeMap := make(map[string]string, 200)
	tvSiteTypeMap := make(map[string]string, 200)
	ratioSiteTypeMap := make(map[string]string, 200)

	samplesMap := make(map[string]map[string]interface{}, 1000)

	for sampleID := range trMap {
		if _, exists := samplesMap[sampleID]; !exists {
			samplesMap[sampleID] = make(map[string]interface{}, 200)
		}

		for siteType := range trMap[sampleID] {
			if _, exists := trSiteTypeMap[siteType]; !exists {
				var trName bytes.Buffer
				var tvName bytes.Buffer
				var ratioName bytes.Buffer

				trName.WriteString(siteType)
				trName.WriteString(" ")
				trName.WriteString(trKey)

				trSiteTypeMap[siteType] = trName.String()

				tvName.WriteString(siteType)
				tvName.WriteString(" ")
				tvName.WriteString(tvKey)

				tvSiteTypeMap[siteType] = tvName.String()

				ratioName.WriteString(siteType)
				ratioName.WriteString(" ")
				ratioName.WriteString(trTvRatioKey)

				ratioSiteTypeMap[siteType] = ratioName.String()
			}

			samplesMap[sampleID][trSiteTypeMap[siteType]] = trMap[sampleID][siteType]
			samplesMap[sampleID][tvSiteTypeMap[siteType]] = tvMap[sampleID][siteType]
			samplesMap[sampleID][ratioSiteTypeMap[siteType]] = ratioMap[sampleID][siteType]
		}
	}

	// // conduct QC
	// trTvArray will hold all of the ratios for total trTv
	trTvRatioArray := make([]float64, 0, 1000)

	// total names
	var totalTrName bytes.Buffer
	var totalTvName bytes.Buffer
	var totalRatioName bytes.Buffer

	totalTrName.WriteString(totalKey)
	totalTrName.WriteString(" ")
	totalTrName.WriteString(trKey)

	totalTvName.WriteString(totalKey)
	totalTvName.WriteString(" ")
	totalTvName.WriteString(tvKey)

	totalRatioName.WriteString(totalKey)
	totalRatioName.WriteString(" ")
	totalRatioName.WriteString(trTvRatioKey)

	// numberOfSamples := len(allMap)

	var sampleNames []string

	siteTypes = nil

	var totalSiteTypes []string

	uniqueSites := make(map[string]struct{})
	for sampleName, sampleHash := range samplesMap {
		if sampleName != totalKey {
			sampleNames = append(sampleNames, sampleName)
			trTvRatioArray = append(trTvRatioArray, float64(ratioMap[sampleName][totalKey]))
		}

		// We don't do this for the "total" sampleName because the totalKey
		// must contain only keys found in the sampleName
		for siteType := range sampleHash {
			if _, exists := uniqueSites[siteType]; exists {
				continue
			}

			// To skip anything we've seen
			uniqueSites[siteType] = struct{}{}

			// Handle totals differently, we want them at the end
			if strings.Contains(siteType, totalKey) {
				totalSiteTypes = append(totalSiteTypes, siteType)
			} else {
				siteTypes = append(siteTypes, siteType)
			}
		}
	}

	// sampleNames =
	sort.Strings(sampleNames)
	// We skipped the totalKey above, so that we may put it first
	sampleNames = append([]string{totalKey}, sampleNames...)

	sort.Strings(totalSiteTypes)
	sort.Strings(siteTypes)

	siteTypes = append(siteTypes, totalSiteTypes...)

	trTvMean := mean(trTvRatioArray)
	trTvMedian := median(trTvRatioArray)
	trTvSd := stdDev(trTvRatioArray, trTvMean)

	// fmt.Printf("Number of SNP %d\n", snpCount)
	// fmt.Printf("Number of DEL %d\n", delCount)
	// fmt.Printf("Number of INS %d\n", insCount)
	// fmt.Printf("Number of MULTIALLELIC %d\n", multiCount)
	// fmt.Printf("Number of Discordant %d\n", discordantCount)
	// fmt.Printf("Number of Weird %d\n", weirdSites)
	// fmt.Printf("Number of Lines %d\n", rowCount)

	// fmt.Printf("Transition:Transversion Ratio mean %f\n", trTvMean)
	// fmt.Printf("Transition:Transversion Ratio median %f\n", trTvMedian)
	// fmt.Printf("Transition:Transversion Ratio standard deviation %f\n", trTvSd)
	// this one contains both counts and ratios, and is what we put into the return json
	// sampleId|total : "siteType|total|exonAlleleFunc transitions|transversions|ratio" = Z
	allMap := make(map[string]map[string]interface{}, 2)

	//later we will have failedSamples
	allMap["stats"] = map[string]interface{}{
		trTvRatioMeanKey: jsonFloat(trTvMean), trTvRatioMedianKey: jsonFloat(trTvMedian),
		trTvRatioStdDevKey: jsonFloat(trTvSd),
	}

	allMap["results"] = map[string]interface{}{
		"order":   siteTypes,
		"samples": samplesMap,
	}

	if *outputJSONPath != "" {
		json, err := json.Marshal(allMap)

		if err != nil {
			log.Fatal(err)
		}

		err = ioutil.WriteFile(*outputJSONPath, json, os.FileMode(0644))

		if err != nil {
			log.Fatal(err)
		}
	}

	// Write output as a tabbed file
	outFh := (*os.File)(nil)

	if *outputTabPath != "" {
		var err error
		outFh, err = os.OpenFile(*outputTabPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0600)

		if err != nil {
			log.Fatal(err)
		}
	} else {
		outFh = os.Stdout
	}

	defer outFh.Close()
	writer := csv.NewWriter(outFh)

	if *fieldSeparator == "\t" {
		writer.Comma = '\t'
	}

	// first column is for sample names
	outLines := [][]string{append([]string{"\t"}, siteTypes...)}

	for _, sampleName := range sampleNames {
		// First column is for the sample name
		line := []string{sampleName}

		for _, siteType := range siteTypes {
			switch value := allMap[sampleName][siteType].(type) {
			case int:
				line = append(line, strconv.Itoa(value))
			case jsonFloat:
				line = append(line, strconv.FormatFloat(float64(value), 'f', -1, 64))
			}
		}

		outLines = append(outLines, line)
	}
	writer.WriteAll(outLines)

	// Write output as a tabbed file
	outQcFh := (*os.File)(nil)
	defer outQcFh.Close()

	if *outputQcTabPath != "" {
		var err error
		outFh, err = os.OpenFile(*outputQcTabPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0600)

		if err != nil {
			log.Fatal(err)
		}
	} else {
		outQcFh = os.Stdout
	}

	writer = csv.NewWriter(outQcFh)

	if *fieldSeparator == "\t" {
		writer.Comma = '\t'
	}

	outQcLines := [][]string{
		[]string{"Transition:Transversion Ratio Mean", strconv.FormatFloat(trTvMean, 'f', -1, 64)},
		[]string{"Transition:Transversion Ratio Median", strconv.FormatFloat(trTvMedian, 'f', -1, 64)},
		[]string{"Transition:Transversion Ratio Standard Deviation", strconv.FormatFloat(trTvSd, 'f', -1, 64)},
	}

	// outQcLines = append(outQcLines, []string{ "Transition:Transversion Ratio Mean", trTvMean})
	// outQcLines = append(outQcLines, "Transition:Transversion Ratio Median", trTvMedian)
	// outQcLines = append(outQcLines, "Transition:Transversion Ratio Median", trTvSd)

	writer.WriteAll(outQcLines)

	pprof.StopCPUProfile()
}

// //https://github.com/dropbox/godropbox/blob/master/sort2/sort.go
func mean(numbers []float64) (meanVal float64) {
	return sum(numbers) / float64(len(numbers))
}

func sum(numbers []float64) (total float64) {
	for _, x := range numbers {
		total += x
	}
	return total
}

func median(numbers []float64) float64 {
	middle := len(numbers) / 2
	result := numbers[middle]
	if len(numbers)%2 == 0 {
		result = (result + numbers[middle-1]) / 2
	}
	return result
}

func mode(numbers []float64) (modes []float64) {
	frequencies := make(map[float64]int, len(numbers))
	highestFrequency := 0
	for _, x := range numbers {
		frequencies[x]++
		if frequencies[x] > highestFrequency {
			highestFrequency = frequencies[x]
		}
	}
	for x, frequency := range frequencies {
		if frequency == highestFrequency {
			modes = append(modes, x)
		}
	}
	if highestFrequency == 1 || len(modes) == len(numbers) {
		modes = modes[:0] // Or: modes = []float64{}
	}
	sort.Float64s(modes)
	return modes
}

func stdDev(numbers []float64, mean float64) float64 {
	total := 0.0
	for _, number := range numbers {
		total += math.Pow(number-mean, 2)
	}
	variance := total / float64(len(numbers)-1)
	return math.Sqrt(variance)
}

func makeFillArrayFunc(emptyField string, secondaryDelim string, primaryDelim string) func(string, bool, bool) []string {
	return func(record string, checkSecondary bool, checkDuplicates bool) []string {
		if record == emptyField {
			return nil
		}

		if checkSecondary && strings.Contains(record, secondaryDelim) {
			out := make([]string, 0, 50)

		OUTER:
			for _, val := range strings.Split(record, secondaryDelim) {

				if !strings.Contains(val, primaryDelim) {
					if val == emptyField {
						continue
					}

					if checkDuplicates {
						for _, haveVal := range out {
							if haveVal == val {
								continue OUTER
							}
						}
					}

					continue
				}

			INNER:
				for _, innerVal := range strings.Split(val, primaryDelim) {
					if innerVal == emptyField {
						continue
					}

					if checkDuplicates {
						for _, haveVal := range out {
							if haveVal == innerVal {
								continue INNER
							}
						}
					}

					out = append(out, innerVal)
				}
			}

			return out
		}

		if !strings.Contains(record, primaryDelim) {
			return []string{record}
		}

		out := make([]string, 0, 50)

	PRIMARY:
		for _, innerVal := range strings.Split(record, primaryDelim) {
			if innerVal != emptyField {
				if checkDuplicates {
					for _, haveVal := range out {
						if haveVal == innerVal {
							continue PRIMARY
						}
					}
				}
				out = append(out, innerVal)
			}
		}

		return out
	}
}

func makeHasNonEmptyRecordFunc(emptyField string, secondaryDelim string, primaryDelim string) func(string) bool {
	return func(record string) bool {
		if strings.Contains(record, secondaryDelim) {
			for _, val := range strings.Split(record, secondaryDelim) {

				if !strings.Contains(val, primaryDelim) {
					if val != emptyField {
						return true
					}

					continue
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

func findIndex(record []string, field string) int {
	for idx, val := range record {
		if val == field {
			return idx
		}
	}

	return -9
}
