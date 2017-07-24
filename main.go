// estab exports elasticsearch fields as tab separated values
package main

import (
	"bufio"
	"bytes"
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
	"github.com/akotlar/sequtils/parse"

	// "github.com/davecgh/go-spew/spew"
	// "math/big"
	"runtime/pprof"
)

type jsonFloat float64

func (value jsonFloat) MarshalJSON() ([]byte, error) {
	if math.IsInf(float64(value), 1) {
		return []byte("null"), nil
	}

	if math.IsNaN(float64(value)) {
		// Can't use "NA", get json: error calling MarshalJSON for type main.jsonFloat: invalid character 'N' looking for beginning of value
		return []byte("null"), nil
	}

	return []byte(fmt.Sprintf("%.3f", value)), nil
}


// NOTE: For now this only supports \n end of line characters
// If we want to support CLRF or whatever, use either csv package, or set a different delimiter
func main() {
	var inputFilePath string
	flag.StringVar(&inputFilePath, "inPath", "", "The input file path (optional: default is stdin)")
	
	var outputJSONPath string
	flag.StringVar(&outputJSONPath, "outputJSONPath", "", "The output path for the JSON output (optional)")
	
	var outputTabPath string
	flag.StringVar(&outputTabPath, "outputTabPath", "", "The output path for tab-delimited file")
	
	var outputQcTabPath string
	flag.StringVar(&outputQcTabPath, "outputQcTabPath", "", "The output path for tab-delimited quality control file")
	
	var trTvColumnName string
	flag.StringVar(&trTvColumnName, "trTvColumnName", "trTv", "The trTv column name.")
	
	var referenceColumnName string
	flag.StringVar(&referenceColumnName, "referenceColumnName", "",
		"The reference base column name. This is usually the name of the assembly")

	var alleleColumnName string
	flag.StringVar(&alleleColumnName, "alleleColumnName", "minorAlleles", "The alleles column name")
		
	var homozygotesColumnName string
	flag.StringVar(&homozygotesColumnName, "homozygotesColumnName", "homozygotes",
		"The homozygous sample column name")
	
	var heterozygotesColumnName string 
	flag.StringVar(&heterozygotesColumnName, "heterozygotesColumnName", "heterozygotes",
		"The homozygous sample column name")
	
	var siteTypeColumnName string
	flag.StringVar(&siteTypeColumnName, "siteTypeColumnName", "refSeq.siteType", "The site type column name")

	var dbSNPnameColumnName string
	flag.StringVar(&dbSNPnameColumnName, "dbSNPnameColumnName", "dbSNP.name", "Optional. The snp name column name")
	
	var exonicAlleleFunctionColumnName string
	flag.StringVar(&exonicAlleleFunctionColumnName, "exonicAlleleFunctionColumnName",
		"refSeq.codonEffect", `Optional. The name of the column that has exonicAlleleFunction, aka the one that has
    nonSynonymous, synonymous, etc values`)
	
	var fieldSeparator string
	flag.StringVar(&fieldSeparator, "fieldSeparator", "\t", "What is used to delimit fields (',', '\t', etc)")
	
	var primaryDelimiter string
	flag.StringVar(&primaryDelimiter, "primaryDelimiter", ";",
		"The primary delmiter (1D array string representation separator in input file)")
	
	var emptyFieldString string
	flag.StringVar(&emptyFieldString, "emptyFieldString", "!",
		"What is used to denoted an empty field (\\ by default)")
	
	var countSNPmulti bool
	flag.BoolVar(&countSNPmulti, "countSNPmulti", false, "Count SNP sites that have 2 non-reference alleles")

	var cpuProfile string
	flag.StringVar(&cpuProfile, "cpuProfile", "", "write cpu profile to file")
	flag.Parse()

	if cpuProfile != "" {
		f, err := os.Create(cpuProfile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	inFh := (*os.File)(nil)

	if inputFilePath != "" {
		var err error
		inFh, err = os.Open(inputFilePath)

		if err != nil {
			log.Fatal(err)
		}
	} else {
		inFh = os.Stdin
	}

	// make sure it gets closed
	defer inFh.Close()

	trTvColumnIdx := -9
	referenceColumnIdx := -9
	alleleColumnIdx := -9
	homozygotesColumnIdx := -9
	heterozygotesColumnIdx := -9
	siteTypeColumnIdx := -9
	dbSNPnameColumnIdx := -9
	exonicAlleleFunctionColumnIdx := -9

	numberInputHeaderLines := 1

	if referenceColumnName == ""{
		log.Fatal(`If referenceColumnIdx not provided, referenceColumnName must be, and
      the numberInputHeaderLines must equal 1`)
	}

	if alleleColumnName == "" {
		log.Fatal(`If alleleColumnIdx not provided, alleleColumnName must be, and
      the numberInputHeaderLines must equal 1`)
	}

	if homozygotesColumnName == "" {
		log.Fatal(`If homozygotesColumnIdx not provided, homozygotesColumnName must be, and
      the numberInputHeaderLines must equal 1`)
	}

	if heterozygotesColumnName == "" {
		log.Fatal(`If heterozygotesColumnIdx not provided, heterozygotesColumnName must be, and
      the numberInputHeaderLines must equal 1`)
	}

	if siteTypeColumnName == ""{
		log.Fatal(`If siteTypeColumnIdx not provided, siteTypeColumnName must be, and
      the numberInputHeaderLines must equal 1`)
	}

	var trTv string

	discordantCount := 0

	//Form: sampleId|total : siteType|total|exonAlleleFunc = N
	trMap := make(map[string]map[string]int, 1000)
	tvMap := make(map[string]map[string]int, 1000)

	//Form: sampleId|total : siteType|total|exonAlleleFunc = Y
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
	hasDbSnpColumn := dbSNPnameColumnIdx != -9
	hasDbSnp := false

	hasExonicColumn := exonicAlleleFunctionColumnIdx != -9

	rowCount := 0

	reader := bufio.NewReader(inFh)

	trMap[totalKey] = make(map[string]int, 200)
	tvMap[totalKey] = make(map[string]int, 200)
	trMap[totalKey][totalKey] = 0
	tvMap[totalKey][totalKey] = 0

	dbSNPfeatureMap := make(map[string]string, 200)

	fillArrayFunc := makeFillArrayFunc(emptyFieldString, primaryDelimiter)

	nonNullFunc := makeHasNonEmptyRecordFunc(emptyFieldString, primaryDelimiter)

	// samples that are variant in a single row, capacity
	samples := make([]string, 1, 1000)
	samples[0] = totalKey

	siteTypes := make([]string, 0, 200)

	discordantRows := make([][]string, 0, 400)

	var record []string

	var simpleTrTv bool
	// Initialize a default row lenght; this is just long enough to contain
	// the default ref field index
	// We will update this in our header reader
	rowLength := 9
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

		if rowCount < numberInputHeaderLines {
			rowCount++

			if rowCount == 1 && numberInputHeaderLines == 1 {
				record = strings.Split(row, fieldSeparator)

				rowLength = len(record)

				if trTvColumnIdx == -9 && trTvColumnName != "" {
					trTvColumnIdx = findIndex(record, trTvColumnName)

					simpleTrTv = trTvColumnIdx != -9
				}

				if referenceColumnName != "" {
					referenceColumnIdx = findIndex(record, referenceColumnName)

					if referenceColumnIdx == -9 {
						log.Fatal("referenceColumnName not found")
					}
				}

				if alleleColumnName != "" {
					alleleColumnIdx = findIndex(record, alleleColumnName)

					if alleleColumnIdx == -9 {
						log.Fatal("alleleColumnName not found")
					}
				}

				if heterozygotesColumnName != "" {
					heterozygotesColumnIdx = findIndex(record, heterozygotesColumnName)

					if heterozygotesColumnIdx == -9 {
						log.Fatal("heterozygotesColumnName not found")
					}
				}

				if homozygotesColumnName != "" {
					homozygotesColumnIdx = findIndex(record, homozygotesColumnName)

					if homozygotesColumnIdx == -9 {
						log.Fatal("homozygotesColumnName not found")
					}
				}

				if siteTypeColumnName != "" {
					siteTypeColumnIdx = findIndex(record, siteTypeColumnName)

					if siteTypeColumnIdx == -9 {
						log.Fatal("siteTypeColumnNme not found")
					}
				}

				if dbSNPnameColumnName != "" {
					dbSNPnameColumnIdx = findIndex(record, dbSNPnameColumnName)

					if dbSNPnameColumnIdx != -9 {
						hasDbSnpColumn = true
					}
				}

				if exonicAlleleFunctionColumnName != "" {
					exonicAlleleFunctionColumnIdx = findIndex(record, exonicAlleleFunctionColumnName)

					if exonicAlleleFunctionColumnIdx != -9 {
						hasExonicColumn = true
					}
				}
			}

			continue
		}

		record = strings.Split(row, fieldSeparator)

		// Skip very short lines
		if len(record) < rowLength {
			continue
		}

		if record[2] != "SNP" {
			continue
		}

		if !simpleTrTv {
			trTv = string(parse.TrTv(record[referenceColumnIdx], record[alleleColumnIdx]))
		} else {
			trTv = record[trTvColumnIdx]
		}

		if trTv != "0" {
			continue;
		}

		if hasDbSnpColumn {
			hasDbSnp = nonNullFunc(record[dbSNPnameColumnIdx])
		}

		// remove everything but the total key
		samples = samples[:1]

		// Make siteTypes, exonicTypes, and samples unique; struct{} takes up no additional space
		// http://stackoverflow.com/questions/9251234/go-append-if-unique
		// Samples are assumed to be unique already, so uniqueness test flags disabled
		samples = append(samples, fillArrayFunc(record[heterozygotesColumnIdx], false)...)
		samples = append(samples, fillArrayFunc(record[homozygotesColumnIdx], false)...)
		siteTypes = fillArrayFunc(record[siteTypeColumnIdx], true)

		if hasExonicColumn == true {
			if siteTypes == nil {
				siteTypes = fillArrayFunc(record[exonicAlleleFunctionColumnIdx], true)
			} else {
				siteTypes = append(siteTypes, fillArrayFunc(record[exonicAlleleFunctionColumnIdx], true)...)
			}
		}

		for _, sample := range samples {
			if _, exists := trMap[sample]; !exists {
				trMap[sample] = make(map[string]int, 200)
				tvMap[sample] = make(map[string]int, 200)
			}

			if trTv == "1" {
				trMap[sample][totalKey]++
			} else {
				tvMap[sample][totalKey]++
			}

			for _, siteType := range siteTypes {
				if trTv == "1" {
					trMap[sample][siteType]++
				} else {
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

					if trTv == "1" {
						trMap[sample][dbSNPfeatureMap[siteType]]++
					} else {
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

	numSamples := len(sampleNames)

	// sampleNames =
	sort.Strings(sampleNames)
	// We skipped the totalKey above, so that we may put it first
	sampleNames = append([]string{totalKey}, sampleNames...)

	sort.Strings(totalSiteTypes)
	sort.Strings(siteTypes)

	siteTypes = append(siteTypes, totalSiteTypes...)

	var trTvMean float64
	var trTvMedian float64
	var trTvSd float64

	if len(trTvRatioArray) > 0 {
		trTvMean = mean(trTvRatioArray)
		trTvMedian = median(trTvRatioArray)

		if trTvMean != 0.0 {
			trTvSd = stdDev(trTvRatioArray, trTvMean)
		}
	}

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
		"samples":          numSamples,
	}

	allMap["results"] = map[string]interface{}{
		"order":   siteTypes,
		"samples": samplesMap,
	}

	if outputJSONPath != "" {
		json, err := json.Marshal(allMap)

		if err != nil {
			log.Fatal(err)
		}

		err = ioutil.WriteFile(outputJSONPath, json, os.FileMode(0644))

		if err != nil {
			log.Fatal(err)
		}
	}

	// Write Tab output
	outFh := (*os.File)(nil)

	if outputTabPath != "" {
		var err error
		outFh, err = os.OpenFile(outputTabPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0600)

		if err != nil {
			log.Fatal(err)
		}
	} else {
		outFh = os.Stdout
	}

	defer outFh.Close()
	writer := csv.NewWriter(outFh)

	writer.Comma = rune(fieldSeparator[0])

	// first column is for sample names
	outLines := [][]string{append([]string{fieldSeparator}, siteTypes...)}

	for _, sampleName := range sampleNames {
		// First column is for the sample name

		line := []string{sampleName}

		for _, siteType := range siteTypes {
			switch value := samplesMap[sampleName][siteType].(type) {
			case int:
				line = append(line, strconv.Itoa(value))
			case jsonFloat:
				line = append(line, strconv.FormatFloat(float64(value), 'f', -1, 64))
			default:
				line = append(line, emptyFieldString)
			}
		}

		outLines = append(outLines, line)
	}

	writer.WriteAll(outLines)

	// Write output as a tabbed file
	outQcFh := (*os.File)(nil)
	defer outQcFh.Close()

	if countSNPmulti == true {
		fmt.Printf("\n\nFound %d SNP sites that look multi-allelic\n\n", discordantCount)

		if discordantCount != 0 {
			discordantWriter := csv.NewWriter(os.Stdout)
			discordantWriter.Comma = '\t'
			discordantWriter.WriteAll(discordantRows)
		}
	}

	if outputQcTabPath != "" {
		var err error
		outFh, err = os.OpenFile(outputQcTabPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0600)

		if err != nil {
			log.Fatal(err)
		}
	} else {
		outQcFh = os.Stdout
	}

	writer = csv.NewWriter(outQcFh)

	writer.Comma = rune(fieldSeparator[0])

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

func makeFillArrayFunc(emptyField string, primaryDelim string) func(string, bool) []string {
	return func(record string, checkDuplicates bool) []string {
		if record == emptyField {
			return nil
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

func makeHasNonEmptyRecordFunc(emptyField string, primaryDelim string) func(string) bool {
	return func(record string) bool {
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
