// estab exports elasticsearch fields as tab separated values
package main

import (
	"bufio"
	"bytes"
	"encoding/csv"
	// "encoding/json"
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
	"github.com/akotlar/bystro-utils/parse"
  "github.com/pquerna/ffjson/ffjson"

	// "github.com/davecgh/go-spew/spew"
	// "math/big"
  "sync"
	"runtime/pprof"
)

type jsonFloat float64

func (value jsonFloat) MarshalJSON() ([]byte, error) {
  if float64(value) == float64(int64(value)) {
    return []byte(strconv.Itoa(int(value))), nil
  }

  //inf will not work
	if math.IsInf(float64(value), 1) {
		return []byte("null"), nil
	}

  // Can't use "NA", get json: error calling MarshalJSON for type main.jsonFloat: invalid character 'N' looking for beginning of value
	if math.IsNaN(float64(value)) {
		return []byte("null"), nil
	}

	return []byte(strconv.FormatFloat(float64(value), 'f', 4, 64)), nil
}

type Config struct {
  inPath string
  outTabPath string
  outQcTabPath  string
  outJsonPath string
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
  emptyField string
  cpuProfile string
  maxThreads int
}

const totalKey string = "total"
const trKey string = "transitions"
const tvKey string = "transversions"
const trTvRatioKey string = "transitions:transversions ratio"
const trTvRatioMeanKey string = "transitions:transversions ratio mean"
const trTvRatioMedianKey string = "transitions:transversions ratio median"
const trTvRatioStdDevKey = "transitions:transversions ratio standard deviation"
const dbSnpKey string = "_in_dbSNP"
// NOTE: For now this only supports \n end of line characters
// If we want to support CLRF or whatever, use either csv package, or set a different delimiter
func setup(args []string) *Config {
	config := &Config{}
	flag.StringVar(&config.inPath, "inPath", "", "The input file path (default: stdin)")
	flag.StringVar(&config.outTabPath, "outTabPath", "", "The output path for tab-delimited file (default: stdout)")
	flag.StringVar(&config.outQcTabPath, "outQcTabPath", "", "The output path for tab-delimited quality control file (default: stdout)")
	flag.StringVar(&config.outJsonPath, "outJsonPath", "", "The output path for JSON output if you wish for it (default: '')")
	flag.StringVar(&config.typeColumnName, "typeColumn", "type", "The type column name (default: type)")
  flag.StringVar(&config.trTvColumnName, "trTvColumn", "trTv", "The trTv column name (default: trTv)")
	flag.StringVar(&config.refColumnName, "refColumnName", "ref",
		"The reference base column name. This is usually the name of the assembly (default: ref)")
	flag.StringVar(&config.altColumnName, "altColumnName", "alt", "The alleles column name (default: alt)")
	flag.StringVar(&config.homozygotesColumnName, "homozygotesColumn", "homozygotes",
		"The homozygous sample column name (default: homozygotes)")
	flag.StringVar(&config.heterozygotesColumnName, "heterozygotesColumn", "heterozygotes",
		"The homozygous sample column name (default: heterozygotes)")
	flag.StringVar(&config.siteTypeColumnName, "siteTypeColumn", "refSeq.siteType", "The site type column name (default: refSeq.siteType)")
	flag.StringVar(&config.dbSNPnameColumnName, "dbSnpNameColumn", "dbSNP.name", "Optional. The snp name column name (default: dbSNP.name)")
	flag.StringVar(&config.exonicAlleleFunctionColumnName, "exonicAlleleFunctionColumn",
		"refSeq.exonicAlleleFunction", `The name of the column that has nonSynonymous, synonymous, etc values (default: refSeq.exonicAlleleFunction)`)
	flag.StringVar(&config.fieldSeparator, "fieldSeparator", "\t", "What is used to delimit fields (deault '\\t')")
	flag.StringVar(&config.primaryDelimiter, "primaryDelimiter", ";",
		"The value delimiter (default ';')")
	flag.StringVar(&config.emptyField, "emptyField", "!",
		"What is used to denoted an empty field (default: '!')")
	flag.StringVar(&config.cpuProfile, "cpuProfile", "", "write cpu profile to file")
  flag.IntVar(&config.maxThreads, "maxThreads", 8, "Number of goroutines to use")
	
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

	if config.cpuProfile != "" {
		f, err := os.Create(config.cpuProfile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	inFh := (*os.File)(nil)

	if config.inPath != "" {
		var err error
		inFh, err = os.Open(config.inPath)

		if err != nil {
			log.Fatal(err)
		}
	} else {
		inFh = os.Stdin
	}

	// make sure it gets closed
	defer inFh.Close()

  reader := bufio.NewReader(inFh)

  processAnnotation(config, reader)
}

func processAnnotation(config *Config, reader *bufio.Reader) {
  // Read buffer
  workQueue := make(chan string, 100)
  complete := make(chan bool)
  // Write buffer
  trResults := make(chan map[string]map[string]int, 100)
  tvResults := make(chan map[string]map[string]int, 100)

  trOverallMap := make(map[string]map[string]int)
  tvOverallMap := make(map[string]map[string]int)

  var wg sync.WaitGroup

  endOfLineByte, numChars, header, err := parse.FindEndOfLine(reader, "")

  if err != nil {
    log.Fatal(err)
  }

  headerFields := strings.Split(header[:len(header) - numChars], config.fieldSeparator)

  trTvIdx, typeIdx, refIdx, altIdx, hetIdx, homIdx, siteTypeIdx, dbSnpNameIdx, exonicAlleleFunctionIdx := findFeatures(headerFields, config)

  // TODO check if any -9

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

  // Now read them all off, concurrently.
  for i := 0; i < config.maxThreads; i++ {
    go processLines(trTvIdx, typeIdx, refIdx, altIdx, hetIdx, homIdx, siteTypeIdx,
      dbSnpNameIdx, exonicAlleleFunctionIdx, config, workQueue, trResults, tvResults, complete)
  }

  wg.Add(1)
  go func() {
    defer wg.Done()

    for trMap := range trResults {
      if trMap != nil {
        for k, _ := range trMap {
          if trOverallMap[k] == nil {
            trOverallMap[k] = make(map[string]int)
          }

          for iK, iV := range trMap[k] {
            trOverallMap[k][iK] = iV
          }
        }
      }
    }
  }()

  wg.Add(1)
  go func() {
    defer wg.Done()

    for tvMap := range tvResults {
      if tvMap != nil {
        for k, _ := range tvMap {
          if tvOverallMap[k] == nil {
            tvOverallMap[k] = make(map[string]int)
          }

          for iK, iV := range tvMap[k] {
            tvOverallMap[k][iK] = iV
          }
        }
      }
    }
  }()

  // Wait for everyone to finish.
  for i := 0; i < config.maxThreads; i++ {
    <-complete
  }

  close(trResults)
  close(tvResults)

  wg.Wait()

  trSiteTypeMap := make(map[string]string, 200)
  tvSiteTypeMap := make(map[string]string, 200)
  ratioSiteTypeMap := make(map[string]string, 200)

  samplesMap := make(map[string]map[string]jsonFloat, 1000)

  for sampleID := range trOverallMap {
    if samplesMap[sampleID] == nil{
      samplesMap[sampleID] = make(map[string]jsonFloat, 100)
      // ratioMap[sampleID] = make(map[string]jsonFloat, 100)
    }

    for siteType := range trOverallMap[sampleID] {
      if trSiteTypeMap[siteType] == "" {
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

      samplesMap[sampleID][trSiteTypeMap[siteType]] = jsonFloat(trOverallMap[sampleID][siteType])
      samplesMap[sampleID][tvSiteTypeMap[siteType]] = jsonFloat(tvOverallMap[sampleID][siteType])

      // If denominator is 0, NaN will result, which we will store as config.emptyField
      // https://github.com/raintank/metrictank/commit/5de7d6e3751901a23501e5fcd95f0b2d0604e8f4
      samplesMap[sampleID][ratioSiteTypeMap[siteType]] = jsonFloat(trOverallMap[sampleID][siteType]) / jsonFloat(tvOverallMap[sampleID][siteType])
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

  var totalSiteTypes []string
  var siteTypes []string

  uniqueSites := make(map[string]bool)
  for sampleName, sampleHash := range samplesMap {
    if sampleName != totalKey {
      sampleNames = append(sampleNames, sampleName)
      trTvRatioArray = append(trTvRatioArray, float64(samplesMap[sampleName][totalKey]))
    }

    // We don't do this for the "total" sampleName because the totalKey
    // must contain only keys found in the sampleName
    for siteType := range sampleHash {
      if uniqueSites[siteType] {
        continue
      }

      // To skip anything we've seen
      uniqueSites[siteType] = true

      // Handle totals differently, we want them at the end
      if strings.Contains(siteType, totalKey) {
        totalSiteTypes = append(totalSiteTypes, siteType)
      } else {
        siteTypes = append(siteTypes, siteType)
      }
    }
  }

  numSamples := float64(len(sampleNames))

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

  // this one contains both counts and ratios, and is what we put into the return json
  // sampleId|total : "siteType|total|exonAlleleFunc transitions|transversions|ratio" = Z
  allMap := make(map[string]map[string]interface{}, 2)

  //later we will have failedSamples
  allMap["stats"] = map[string]interface{} {
    trTvRatioMeanKey: jsonFloat(trTvMean),
    trTvRatioMedianKey: jsonFloat(trTvMedian),
    trTvRatioStdDevKey: jsonFloat(trTvSd),
    "samples": numSamples,
  }

  allMap["results"] = map[string]interface{} {
    "samples": samplesMap,
  }

  if config.outJsonPath != "" {
    json, err := ffjson.Marshal(allMap)

    if err != nil {
      log.Fatal(err)
    }

    err = ioutil.WriteFile(config.outJsonPath, json, os.FileMode(0644))

    if err != nil {
      log.Fatal(err)
    }
  }

  // Write Tab output
  outFh := (*os.File)(nil)

  if config.outTabPath != "" {
    var err error
    outFh, err = os.OpenFile(config.outTabPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0600)

    if err != nil {
      log.Fatal(err)
    }
  } else {
    outFh = os.Stdout
  }

  defer outFh.Close()
  writer := csv.NewWriter(outFh)

  writer.Comma = rune(config.fieldSeparator[0])

  // first column is for sample names
  outLines := [][]string{append([]string{config.fieldSeparator}, siteTypes...)}

  for _, sampleName := range sampleNames {
    // First column is for the sample name

    line := []string{sampleName}

    for _, siteType := range siteTypes {
      val := samplesMap[sampleName][siteType]

      if float64(val) == float64(int64(val)) {
        line = append(line, strconv.Itoa(int(val)))
      } else {
        line = append(line, strconv.FormatFloat(float64(val), 'f', 4, 64))
      }
    }

    outLines = append(outLines, line)
  }

  writer.WriteAll(outLines)

  // // Write output as a tabbed file
  // outQcFh := (*os.File)(nil)
  // defer outQcFh.Close()

  // if config.outQcTabPath != "" {
  //   var err error
  //   outFh, err = os.OpenFile(config.outQcTabPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0600)

  //   if err != nil {
  //     log.Fatal(err)
  //   }
  // } else {
  //   outQcFh = os.Stdout
  // }

  // writer = csv.NewWriter(outQcFh)

  // writer.Comma = rune(config.fieldSeparator[0])

  // outQcLines := [][]string{
  //   []string{"Transition:Transversion Ratio Mean", strconv.FormatFloat(trTvMean, 'f', -1, 64)},
  //   []string{"Transition:Transversion Ratio Median", strconv.FormatFloat(trTvMedian, 'f', -1, 64)},
  //   []string{"Transition:Transversion Ratio Standard Deviation", strconv.FormatFloat(trTvSd, 'f', -1, 64)},
  // }

  // // outQcLines = append(outQcLines, []string{ "Transition:Transversion Ratio Mean", trTvMean})
  // // outQcLines = append(outQcLines, "Transition:Transversion Ratio Median", trTvMedian)
  // // outQcLines = append(outQcLines, "Transition:Transversion Ratio Median", trTvSd)

  // writer.WriteAll(outQcLines)
}

func findFeatures (record []string, config *Config) (int, int, int, int, int, int, int, int, int) {
  trTvIdx := -9
  typeIdx := -9
  refIdx := -9
  altIdx := -9
  hetIdx := -9
  homIdx := -9
  siteTypeIdx := -9
  dbSnpNameIdx := -9
  exonicAlleleFunctionIdx := -9

  if config.trTvColumnName != "" {
    trTvIdx = findIndex(record, config.trTvColumnName)
  }

  if config.typeColumnName != "" {
    typeIdx = findIndex(record, config.typeColumnName)
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

  return trTvIdx, typeIdx, refIdx, altIdx, hetIdx, homIdx, siteTypeIdx, dbSnpNameIdx, exonicAlleleFunctionIdx
}

func processLines (trTvIdx int, typeIdx int, refIdx int, altIdx int, hetIdx int, homIdx int,
siteTypeIdx int, dbSnpNameIdx int, exonicAlleleFunctionIdx int, config *Config,
queue chan string, trResults chan map[string]map[string]int, tvResults chan map[string]map[string]int,
complete chan bool) {
	//Form: sampleId|total : siteType|total|exonAlleleFunc = N
	trMap := make(map[string]map[string]int, 1000)
	tvMap := make(map[string]map[string]int, 1000)

	//Form: sampleId|total : siteType|total|exonAlleleFunc = Y
	// ratioMap := make(map[string]map[string]jsonFloat, 1000)

	trMap[totalKey] = make(map[string]int, 20)
	tvMap[totalKey] = make(map[string]int, 20)

	featureCache := make(map[string]string, 20)

  var name bytes.Buffer
  name.WriteString(totalKey)
  name.WriteString(dbSnpKey)

  featureCache[totalKey] = name.String()

	siteTypes := make([]string, 0, 10)
  exonicTypes := make([]string, 0, 10)

	simpleTrTv := trTvIdx > -9
  hasDbSnpColumn := dbSnpNameIdx > -9
  hasExonicColumn := exonicAlleleFunctionIdx != -9

  inDbSnp := false
  isTr := false

  var trTv string

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
			trTv = record[trTvIdx]
		}

		if trTv == parse.NotTrTv {
			continue;
		}

    isTr = trTv == parse.Tr

		if hasDbSnpColumn {
			inDbSnp = strings.Contains(record[dbSnpNameIdx], "rs")
		}

		siteTypes = uniqSlice(record[siteTypeIdx], config.emptyField, config.primaryDelimiter)

		if hasExonicColumn {
			exonicTypes = uniqSlice(record[exonicAlleleFunctionIdx], config.emptyField, config.primaryDelimiter)
    }

    if isTr {
      trMap[totalKey][totalKey]++
    } else {
      tvMap[totalKey][totalKey]++
    }

    if inDbSnp {
      if isTr {
        trMap[totalKey][featureCache[totalKey]]++
      } else {
        tvMap[totalKey][featureCache[totalKey]]++
      }
    }

    fillType(strings.Split(record[hetIdx],config.primaryDelimiter), siteTypes, exonicTypes,
      trMap, tvMap, featureCache, isTr, inDbSnp, config.emptyField)

    fillType(strings.Split(record[homIdx], config.primaryDelimiter), siteTypes, exonicTypes,
      trMap, tvMap, featureCache, isTr, inDbSnp, config.emptyField)
	}

  trResults <- trMap
  tvResults <- tvMap

  complete <- true
}

func fillType(samples []string, siteTypes []string, exonicTypes []string,
trMap map[string]map[string]int, tvMap map[string]map[string]int,
featureCache map[string]string, isTr bool, inDbSnp bool, emptyField string) {
  for _, sample := range samples {
    if sample == emptyField {
      continue
    }

    if trMap[sample] == nil {
      trMap[sample] = make(map[string]int)
      tvMap[sample] = make(map[string]int)
    }

    if isTr {
      trMap[sample][totalKey]++
    } else {
      tvMap[sample][totalKey]++
    }

    if inDbSnp {
      if isTr {
        trMap[sample][featureCache[totalKey]]++
      } else {
        tvMap[sample][featureCache[totalKey]]++
      }
    }

    for _, siteType := range siteTypes {
      if isTr {
        trMap[sample][siteType]++
      } else {
        tvMap[sample][siteType]++
      }

      if featureCache[siteType] == "" {
        var name bytes.Buffer
        name.WriteString(siteType)
        name.WriteString(dbSnpKey)

        featureCache[siteType] = name.String()
      }

      if inDbSnp {
        if isTr {
          trMap[sample][featureCache[siteType]]++
        } else {
          tvMap[sample][featureCache[siteType]]++
        }
      }
    }

    for _, siteType := range exonicTypes {
      if isTr {
        trMap[sample][siteType]++
      } else {
        tvMap[sample][siteType]++
      }

      if featureCache[siteType] == "" {
        var name bytes.Buffer
        name.WriteString(siteType)
        name.WriteString(dbSnpKey)

        featureCache[siteType] = name.String()
      }

      if inDbSnp {
        if isTr {
          trMap[sample][featureCache[siteType]]++
        } else {
          tvMap[sample][featureCache[siteType]]++
        }
      }
    }
  }
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

func uniqSlice(record string, emptyField string, primaryDelimiter string) []string {
  var out []string

  if record == emptyField {
    return out
  }

  dup := make(map[string]bool)

  for _, val := range strings.Split(record, primaryDelimiter) {
    if val != emptyField {
      if !dup[val] {
        out = append(out, val)
        dup[val] = true
      }
    }
  }

  return out
}

func findIndex(record []string, field string) int {
	for idx, val := range record {
		if val == field {
			return idx
		}
	}

	return -9
}
