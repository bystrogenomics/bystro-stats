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
	"github.com/akotlar/bystro-utils/parse"

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

type Config struct {
  inPath string
  outTabPath string
  outQcTabPath  string
  outJSONPath string
  trTvColumnName string
  refColumnName string
  altColumnName string
  typeColumnName string
  homozygotesColumnName string
  heterozygotesColumnName string
  siteTypeColumnName string
  dbSNPnameColumnName string
  exonicAlleleFunctionColumnName string
  fieldSeparator string
  primaryDelimiter string
  cpuProfile string
  emptyField string
  fieldDelimiter string
  keepId bool
  keepInfo bool
  cpuProfile string
  allowedFilters map[string]bool
}

// NOTE: For now this only supports \n end of line characters
// If we want to support CLRF or whatever, use either csv package, or set a different delimiter
func setup(args []string) {
	config := &Config{}
	flag.StringVar(&config.inPath, "inPath", "", "The input file path (default: stdin)")
	flag.StringVar(&config.outTabPath, "outTabPath", "", "The output path for tab-delimited file (default: stdout)")
	flag.StringVar(&config.outQcTabPath, "outQcTabPath", "", "The output path for tab-delimited quality control file (default: stdout)")
	flag.StringVar(&config.outJSONPath, "outJSONPath", "", "The output path for JSON output if you wish for it (default: '')")
	flag.StringVar(&config.trTvColumnName, "typeColumnName", "type", "The type column name (default: type)")
  flag.StringVar(&config.trTvColumnName, "trTvColumnName", "trTv", "The trTv column name (default: trTv)")
	flag.StringVar(&config.refColumnName, "refColumnName", "ref",
		"The reference base column name. This is usually the name of the assembly (default: ref)")
	flag.StringVar(&config.altColumnName, "altColumnName", "alt", "The alleles column name (default: alt)")
	flag.StringVar(&config.homozygotesColumnName, "homozygotesColumnName", "homozygotes",
		"The homozygous sample column name (default: homozygotes)")
	flag.StringVar(&config.heterozygotesColumnName, "heterozygotesColumnName", "heterozygotes",
		"The homozygous sample column name (default: heterozygotes)")
	flag.StringVar(&config.siteTypeColumnName, "siteTypeColumnName", "refSeq.siteType", "The site type column name (default: refSeq.siteType)")
	flag.StringVar(&config.dbSNPnameColumnName, "dbSNPnameColumnName", "dbSNP.name", "Optional. The snp name column name (default: dbSNP.name)")
	flag.StringVar(&config.exonicAlleleFunctionColumnName, "exonicAlleleFunctionColumnName",
		"refSeq.exonicAlleleFunction", `The name of the column that has nonSynonymous, synonymous, etc values (default: refSeq.exonicAlleleFunction)`)
	flag.StringVar(&config.fieldSeparator, "fieldSeparator", "\t", "What is used to delimit fields (deault '\\t')")
	flag.StringVar(&config.primaryDelimiter, "primaryDelimiter", ";",
		"The value delimiter (default ';')")
	flag.StringVar(&config.emptyField, "emptyField", "!",
		"What is used to denoted an empty field (default: '!')")
	flag.StringVar(&config.cpuProfile, "cpuProfile", "", "write cpu profile to file")
	
  // allows args to be mocked https://github.com/nwjlyons/email/blob/master/inputs.go
  // can only run 1 such test, else, redefined flags error
  a := os.Args[1:]
  if args != nil {
    a = args
  }
  flag.CommandLine.Parse(a)

  return config
}

func init() {
  log.SetFlags(0)
}

func main() {
  config := setup(nil)

	if cpuProfile != "" {
		f, err := os.Create(cpuProfile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	inFh := (*os.File)(nil)

	if inPath != "" {
		var err error
		inFh, err = os.Open(inPath)

		if err != nil {
			log.Fatal(err)
		}
	} else {
		inFh = os.Stdin
	}

	// make sure it gets closed
	defer inFh.Close()
}

func readAnnotation (config *Config, reader *bufio.Reader, resultFunc func(row string)) {
  foundHeader := false

  var header []string

  // Read buffer
  workQueue := make(chan string, 100)
  complete := make(chan bool)
  // Write buffer
  results := make(chan string, 100)
  var wg sync.WaitGroup

  endOfLineByte, numChars, header, err := parse.FindEndOfLine(reader, "")

  headerFields := strings.Split(header[:len(row) - numChars], config.fieldSeparator)

  trTvIdx, typeIx, refIdx, altIdx, hetIdx, homIdx, siteTypeIdx, dbSnpNameIdx, exonicAlleleFunctionIdx := findFeatures(headerFields, config)

  // Read the lines into the work queue.
  go func() {
    for {
      row, err := reader.ReadString(endOfLineByte) // 0x0A separator = newline

      if err == io.EOF {
        break
      } else if err != nil {
        log.Fatal(err)
      } else if row == "" {
        // We may have not closed the pipe, but not have any more information to send
        // Wait for EOF
        continue
      }

      workQueue <- row[:len(row) - numChars];
    }

    // Close the channel so everyone reading from it knows we're done.
    close(workQueue)
  }()

  wg.Add(1)
  go func() {
    defer wg.Done()
    for line := range results {
      resultFunc(line)
    }
  }()

  // Now read them all off, concurrently.
  for i := 0; i < concurrency; i++ {
    go processLines(trTvIdx, refIdx, altIdx, hetIdx, homIdx, siteTypeIdx,
      dbSnpNameIdx, exonicAlleleFunctionIdx, workQueue, results, complete)
  }

  // Wait for everyone to finish.
  for i := 0; i < concurrency; i++ {
    <-complete
  }

  close(results)

  wg.Wait()
}

func findFeatures (record []string, config *Config) (int, int, int, int, int, int, int, int, int) {
  trTvIdx = -9
  typeIdx = -9
  refIdx = -9
  altIdx = -9
  hetIdx = -9
  homIdx = -9
  siteTypeIdx = -9
  dbSnpNameIdx = -9
  exonicAlleleFunctionIdx = -9

  if config.trTvColumnName != "" {
    trTvColumnIdx = findIndex(record, config.trTvColumnName)
  }

  if config.typeColumnName != "" {
    typeIdx = findIndex(record, config.typeColumnNamee)
  }

  if config.refColumnName != "" {
    refIdx = findIndex(record, config.refColumnName)
  }

  if config.altColumnName != "" {
    altIdx = findIndex(record, config.altColumnName)
  }

  if config.heterozygotesColumnName != "" {
    hetIdx = findIndex(record, config.heterozygotesColumnName)
  }

  if config.homozygotesColumnName != "" {
    homIdx = findIndex(record, config.homozygotesColumnName)
  }

  if config.siteTypeColumnName != "" {
    siteTypeIdx = findIndex(record, config.siteTypeColumnName)
  }

  if config.dbSNPnameColumnName != "" {
    dbSnpNameIdx = findIndex(record, config.dbSNPnameColumnName)
  }

  if config.exonicAlleleFunctionColumnName != "" {
    exonicAlleleFunctionIdx = findIndex(record, config.exonicAlleleFunctionColumnName)
  }

  return trTvIdx, refIdx, altIdx, hetIdx, homIdx, siteTypeIdx, dbSnpNameIdx, exonicAlleleFunctionIdx
}

func processLines (trTvIdx int, typeIdx int, refIdx int, altIdx int, hetIdx int, homIdx int,
  siteTypeIdx int, dbSnpNameIdx int, exonicAlleleFunctionIdx int, config,
  queue chan string, results chan string, complete chan bool) {
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

	fillArrayFunc := makeFillArrayFunc(emptyField, primaryDelimiter)

	// samples that are variant in a single row, capacity
	samples := make([]string, 1, 1000)
	samples[0] = totalKey

	siteTypes := make([]string, 0, 200)

	discordantRows := make([][]string, 0, 400)

	var record []string

	simpleTrTv = trTvIdx > -9
  hasDbSnp = dbSnpNameIdx > -9

  for line := range queue {
    record := strings.Split(line, config.fieldSeparator)

		// Skip very short lines
		if len(record) < 10 {
			continue
		}

		if record[typeIdx] != "SNP" {
			continue
		}

		if !simpleTrTv {
			trTv = parse.GetTrTv(record[refIdx], record[altIdx])
		} else {
			trTv = record[trTvColumnIdx]
		}

		if trTv == parse.NotTrTv {
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

			if trTv == parse.Tr {
				trMap[sample][totalKey]++
			} else {
				tvMap[sample][totalKey]++
			}

			for _, siteType := range siteTypes {
				if trTv == parse.Tr {
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

					if trTv == parse.Tr {
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
		sort.Slice(trTvRatioArray, func(a, b int) bool {
			return trTvRatioArray[a] < trTvRatioArray[b];
		});

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

	if outJSONPath != "" {
		json, err := json.Marshal(allMap)

		if err != nil {
			log.Fatal(err)
		}

		err = ioutil.WriteFile(outJSONPath, json, os.FileMode(0644))

		if err != nil {
			log.Fatal(err)
		}
	}

	// Write Tab output
	outFh := (*os.File)(nil)

	if outTabPath != "" {
		var err error
		outFh, err = os.OpenFile(outTabPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0600)

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
				line = append(line, emptyField)
			}
		}

		outLines = append(outLines, line)
	}

	writer.WriteAll(outLines)

	// Write output as a tabbed file
	outQcFh := (*os.File)(nil)
	defer outQcFh.Close()

	if outQcTabPath != "" {
		var err error
		outFh, err = os.OpenFile(outQcTabPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0600)

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

func findIndex(record []string, field string) int {
	for idx, val := range record {
		if val == field {
			return idx
		}
	}

	return -9
}
