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

  // "github.com/kr/pretty"
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
  trTvColumn string
  refColumn string
  altColumn string
  typeColumn string
  homozygotesColumn string
  heterozygotesColumn string
  siteTypeColumn string
  dbSNPnameColumn string
  exonicAlleleFunctionColumn string
  silentSite string
  replacementSite string
  intergenicSite string
  fieldSeparator string
  primaryDelimiter string
  emptyField string
  cpuProfile string
  maxThreads int
}

type Indices struct {
  trTv int
  typeIdx int
  ref int
  alt int
  het int
  hom int
  siteType int
  dbSNPname int
  exonicAlleleFunction int
}

const totalKey string = "total"
const trKey string = "transitions"
const tvKey string = "transversions"
const trTvRatioKey string = "transitions:transversions ratio"
const dbSnpKey string = "_in_dbSNP"
const hetKey string = "heterozygotes"
const homKey string = "homozygotes"
const hetHomRatioKey string = "heterozygotes:homozygotes ratio"
const silentKey string = "silent"
const replacementKey string = "replacement"
const silentRepRatioKey string = "silent:replacement ratio"
const sitesKey string = "sites"
const thetaKey string = "theta"

// NOTE: For now this only supports \n end of line characters
// If we want to support CLRF or whatever, use either csv package, or set a different delimiter
func setup(args []string) *Config {
  config := &Config{}
  flag.StringVar(&config.inPath, "inPath", "", "The input file path (default: stdin)")
  flag.StringVar(&config.outTabPath, "outTabPath", "", "The output path for tab-delimited file (default: stdout)")
  flag.StringVar(&config.outQcTabPath, "outQcTabPath", "", "The output path for tab-delimited quality control file (default: stdout)")
  flag.StringVar(&config.outJsonPath, "outJsonPath", "", "The output path for JSON output if you wish for it (default: '')")
  flag.StringVar(&config.typeColumn, "typeColumn", "type", "The type column name (default: type)")
  flag.StringVar(&config.trTvColumn, "trTvColumn", "trTv", "The trTv column name (default: trTv)")
  flag.StringVar(&config.refColumn, "refColumn", "ref",
    "The reference base column name. This is usually the name of the assembly (default: ref)")
  flag.StringVar(&config.altColumn, "altColumn", "alt", "The alleles column name (default: alt)")
  flag.StringVar(&config.homozygotesColumn, "homozygotesColumn", "homozygotes",
    "The homozygous sample column name (default: homozygotes)")
  flag.StringVar(&config.heterozygotesColumn, "heterozygotesColumn", "heterozygotes",
    "The homozygous sample column name (default: heterozygotes)")
  flag.StringVar(&config.siteTypeColumn, "siteTypeColumn", "refSeq.siteType", "The site type column name (default: refSeq.siteType)")
  flag.StringVar(&config.dbSNPnameColumn, "dbSnpNameColumn", "dbSNP.name", "Optional. The snp name column name (default: dbSNP.name)")
  flag.StringVar(&config.exonicAlleleFunctionColumn, "exonicAlleleFunctionColumn",
    "refSeq.exonicAlleleFunction", `The name of the column that has nonSynonymous, synonymous, etc values (default: refSeq.exonicAlleleFunction)`)
  flag.StringVar(&config.silentSite, "silentName",
    "synonymous", `What do we call silent/synonymous sites`)
  flag.StringVar(&config.replacementSite, "replacementName",
    "nonSynonymous", `What do we call replacement/nonSynonymous sites`)
  flag.StringVar(&config.intergenicSite, "intergenicName",
    "intergenic", `What do call sites not in any transcript`)
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
      log.Println("Couldn't create cpuProfile file")
      return
    }

    pprof.StartCPUProfile(f)
    defer pprof.StopCPUProfile()
  }

  inFh := (*os.File)(nil)

  if config.inPath != "" {
    var err error
    inFh, err = os.Open(config.inPath)

    if err != nil {
      log.Println("Couldn't read annotation file")
      return
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
  workQueue := make(chan string, 100)

  // Stores total results
  complete := make(chan int)
  // Stores ts & tv info from each thread
  trTvChan := make(chan map[string]map[string][]int, 100)
  hetHomChan := make(chan map[string][]int, 100)

  // Accumulate results
  trTvResults := make(map[string]map[string][]int)
  hetHomResults := make(map[string][]int)
  var totalVariants int

  var wg sync.WaitGroup

  endOfLineByte, numChars, header, err := parse.FindEndOfLine(reader, "")

  if err != nil {
    log.Println("Couldn't read annotation file")
    return
  }

  headerFields := strings.Split(header[:len(header) - numChars], config.fieldSeparator)

  indices := findFeatures(headerFields, config)

  // TODO check if any -9
  if indices.trTv == -9 {
    log.Println("trTv field is required")
    return
  }

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
    go processLines(indices, config, workQueue, trTvChan, hetHomChan, complete)
  }

  wg.Add(1)
  go func() {
    defer wg.Done()

    for trTvMap := range trTvChan {
      if trTvMap != nil {
        for k, _ := range trTvMap {
          if trTvResults[k] == nil {
            trTvResults[k] = make(map[string][]int)
          }

          for iK, iV := range trTvMap[k] {
            if trTvResults[k][iK] == nil {
              trTvResults[k][iK] = []int{0, 0}
            }

            trTvResults[k][iK][0] += iV[0]
            trTvResults[k][iK][1] += iV[1]
          }
        }
      }
    }
  }()

  wg.Add(1)
  go func() {
    defer wg.Done()

    for hhMap := range hetHomChan {
      if hhMap != nil {
        for k, v := range hhMap {
          if hetHomResults[k] == nil {
            hetHomResults[k] = []int{0, 0}
          }

          hetHomResults[k][0] += v[0]
          hetHomResults[k][1] += v[1]
        }
      }
    }
  }()

  // Wait for everyone to finish.
  for i := 0; i < config.maxThreads; i++ {
    totalVariants += <-complete
  }

  close(trTvChan)
  close(hetHomChan)

  wg.Wait()

  trNames := make(map[string]string, 50)
  tvNames := make(map[string]string, 50)
  ratioNames := make(map[string]string, 50)

  var sampleNames []string
  var siteTypes []string

  samplesMap := make(map[string]map[string]jsonFloat, 1000)
  // // conduct QC
  // trTvArray will hold all of the ratios for total trTv
  trTvRatios := make([]float64, 0, 1000)
  silentRepRatios := make([]float64, 0, 1000)
  hetHomRatios := make([]float64, 0, 1000)

  var tr jsonFloat
  var tv jsonFloat
  var trTv jsonFloat

  for sampleID := range trTvResults {
    if samplesMap[sampleID] == nil{
      samplesMap[sampleID] = make(map[string]jsonFloat, 100)
      // ratioMap[sampleID] = make(map[string]jsonFloat, 100)
    }

    if sampleID != totalKey {
      sampleNames = append(sampleNames, sampleID)
    }

    var silent jsonFloat
    var repl jsonFloat
    var silentRepl jsonFloat
    for siteType := range trTvResults[sampleID] {
      if trNames[siteType] == "" {
        var trName bytes.Buffer
        var tvName bytes.Buffer
        var ratioName bytes.Buffer

        trName.WriteString(siteType)
        trName.WriteString(" ")
        trName.WriteString(trKey)

        trNames[siteType] = trName.String()

        tvName.WriteString(siteType)
        tvName.WriteString(" ")
        tvName.WriteString(tvKey)

        tvNames[siteType] = tvName.String()

        ratioName.WriteString(siteType)
        ratioName.WriteString(" ")
        ratioName.WriteString(trTvRatioKey)

        ratioNames[siteType] = ratioName.String()

        if siteType != totalKey {
          siteTypes = append(siteTypes, siteType)
        }
      }

      tr = jsonFloat(trTvResults[sampleID][siteType][0])
      tv = jsonFloat(trTvResults[sampleID][siteType][1])

      samplesMap[sampleID][trNames[siteType]] = tr
      samplesMap[sampleID][tvNames[siteType]] = tv

      // If denominator is 0, NaN will result, which we will store as config.emptyField
      // https://github.com/raintank/metrictank/commit/5de7d6e3751901a23501e5fcd95f0b2d0604e8f4
      trTv = tr / tv

      samplesMap[sampleID][ratioNames[siteType]] = trTv

      if sampleID != totalKey && siteType == totalKey {
        trTvRatios = append(trTvRatios, float64(trTv))
      }
    }

    silent = jsonFloat(trTvResults[sampleID][config.silentSite][0] + trTvResults[sampleID][config.silentSite][1])
    repl = jsonFloat(trTvResults[sampleID][config.replacementSite][0] + trTvResults[sampleID][config.replacementSite][1])
    silentRepl = silent / repl

    samplesMap[sampleID][silentKey] = silent
    samplesMap[sampleID][replacementKey] = repl
    samplesMap[sampleID][silentRepRatioKey] = silentRepl

    // Calculate the sample mean
    if sampleID != totalKey {
      silentRepRatios = append(silentRepRatios, float64(silentRepl))
    }
  }

  for sampleID, val := range hetHomResults {
    if samplesMap[sampleID] == nil{
      samplesMap[sampleID] = make(map[string]jsonFloat, 100)
      // ratioMap[sampleID] = make(map[string]jsonFloat, 100)
    }

    hets := float64(val[0])
    homs := float64(val[1])
    hetHom := hets / homs

    samplesMap[sampleID][hetKey] = jsonFloat(hets)
    samplesMap[sampleID][homKey] = jsonFloat(homs)

    samplesMap[sampleID][hetHomRatioKey] = jsonFloat(hetHom)

    hetHomRatios = append(hetHomRatios, hetHom)
  }

  numSamples := float64(len(sampleNames))

  sort.Strings(sampleNames)
  // We skipped the totalKey above, so that we may put it first
  sampleNames = append([]string{totalKey}, sampleNames...)

  sort.Strings(siteTypes)

  siteTypes = append([]string{totalKey}, siteTypes...)

  /************************ Write Tab Delimited Output ***********************/
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

  var siteOrder []string

  for _, siteType := range siteTypes {
    siteOrder = append(siteOrder, trNames[siteType], tvNames[siteType], ratioNames[siteType])
  }

  siteOrder = append(siteOrder, silentKey, replacementKey, silentRepRatioKey,
    hetKey, homKey, hetHomRatioKey)

  outLines := [][]string{siteOrder}

  for _, sampleName := range sampleNames {
    // First column is for the sample name

    line := []string{sampleName}

    for _, siteType := range siteTypes {
      tr = samplesMap[sampleName][trNames[siteType]]
      tv = samplesMap[sampleName][tvNames[siteType]]
      trTv = samplesMap[sampleName][ratioNames[siteType]]

      line = append(line, roundVal(tr), roundVal(tv), roundVal(trTv))
    }

    line = append(
      line,
      roundVal(samplesMap[sampleName][silentKey]),
      roundVal(samplesMap[sampleName][replacementKey]),
      roundVal(samplesMap[sampleName][silentRepRatioKey]),
    )

    line = append(
      line,
      roundVal(samplesMap[sampleName][hetKey]),
      roundVal(samplesMap[sampleName][homKey]),
      roundVal(samplesMap[sampleName][hetHomRatioKey]),
    )

    outLines = append(outLines, line)
  }

  fmt.Fprint(outFh, config.fieldSeparator)
  writer.WriteAll(outLines)

  /*************************** Calc Statistics && Write JSON **********************************/
  if config.outJsonPath != "" {
    // this one contains both counts and ratios, and is what we put into the return json
    // sampleId|total : "siteType|total|exonAlleleFunc transitions|transversions|ratio" = Z
    allMap := make(map[string]map[string]interface{}, 2)

    allMap["stats"] = map[string]interface{} {
      "samples": numSamples,
      "variants": totalVariants,
    }

    for key, val := range stats(trTvRatios, trTvRatioKey) {
      allMap["stats"][key] = jsonFloat(val)
    }

    for key, val := range stats(silentRepRatios, silentRepRatioKey) {
      allMap["stats"][key] = jsonFloat(val)
    }

    for key, val := range stats(hetHomRatios, hetHomRatioKey) {
      allMap["stats"][key] = jsonFloat(val)
    }

    allMap["results"] = map[string]interface{} {
      "samples": samplesMap,
      "order": siteOrder,
    }

    json, err := json.Marshal(allMap)

    if err != nil {
      log.Fatal(err)
    }

    err = ioutil.WriteFile(config.outJsonPath, json, os.FileMode(0644))

    if err != nil {
      log.Fatal(err)
    }
  }
}

func findFeatures (record []string, config *Config) (Indices) {
  trTvIdx := -9
  typeIdx := -9
  refIdx := -9
  altIdx := -9
  hetIdx := -9
  homIdx := -9
  siteTypeIdx := -9
  dbSNPnameIdx := -9
  exonicAlleleFunctionIdx := -9

  if config.trTvColumn != "" {
    trTvIdx = findIndex(record, config.trTvColumn)
  }

  if config.typeColumn != "" {
    typeIdx = findIndex(record, config.typeColumn)
  }

  if config.refColumn != "" {
    refIdx = findIndex(record, config.refColumn)
  }

  if config.altColumn != "" {
    altIdx = findIndex(record, config.altColumn)
  }

  if config.heterozygotesColumn != "" {
    hetIdx = findIndex(record, config.heterozygotesColumn)
  }

  if config.homozygotesColumn != "" {
    homIdx = findIndex(record, config.homozygotesColumn)
  }

  if config.siteTypeColumn != "" {
    siteTypeIdx = findIndex(record, config.siteTypeColumn)
  }

  if config.dbSNPnameColumn != "" {
    dbSNPnameIdx = findIndex(record, config.dbSNPnameColumn)
  }

  if config.exonicAlleleFunctionColumn != "" {
    exonicAlleleFunctionIdx = findIndex(record, config.exonicAlleleFunctionColumn)
  }

  indices := Indices{
    trTv: trTvIdx, typeIdx: typeIdx, ref: refIdx, alt: altIdx,
    het: hetIdx, hom: homIdx, siteType: siteTypeIdx, dbSNPname: dbSNPnameIdx,
    exonicAlleleFunction: exonicAlleleFunctionIdx,
  }

  return indices
}

func processLines (indices Indices, config *Config,
queue chan string, trTvChan chan map[string]map[string][]int,
hetHomChan chan map[string][]int, complete chan int) {
  //Form: sampleId|total : siteType|total|exonAlleleFunc = N
  trTvMap := make(map[string]map[string][]int, 1000)

  hetHomMap := make(map[string][]int, 1000)
  hetHomMap[totalKey] = []int{0, 0}

  hets := make([]string, 0, 1000)
  homs := make([]string, 0, 100)

  trTvMap[totalKey] = make(map[string][]int, 20)
  featureCache := make(map[string]string, 20)

  var name bytes.Buffer
  name.WriteString(totalKey)
  name.WriteString(dbSnpKey)

  featureCache[totalKey] = name.String()

  siteTypes := make([]string, 0, 10)
  exonicTypes := make([]string, 0, 10)

  // No more support for files that were created before trTv field
  // simpleTrTv := trTvIdx > -9
  hasDbSnpColumn := indices.dbSNPname > -9
  hasExonicColumn := indices.exonicAlleleFunction > -9

  inDbSnp := false
  isTr := false

  trTvIdx := indices.trTv
  hetIdx := indices.het
  homIdx := indices.hom
  funcIdx := indices.exonicAlleleFunction
  dbSNPnameIdx := indices.dbSNPname
  siteTypeIdx := indices.siteType

  fakeTotal := []string{totalKey}

  emptyField := config.emptyField
  pDelim := config.primaryDelimiter

  // var trTv string
  var totalSites int
  // var hetCount int
  // var homCount int

  empty := []string{}
  intergenic := []string{config.intergenicSite}

  var hasHet bool
  var hasHom bool
  for line := range queue {
    record := strings.Split(line, config.fieldSeparator)

    // Skip very short lines
    if len(record) < 10 {
      continue
    }

    totalSites++

    hasHet = record[hetIdx] != emptyField
    hasHom = record[homIdx] != emptyField

    if hasHet {
      hets = strings.Split(record[hetIdx], pDelim)

      for _, sample := range hets {
        if hetHomMap[sample] == nil {
          hetHomMap[sample] = []int{0, 0}
        }

        hetHomMap[sample][0]++
        hetHomMap[totalKey][0]++
      }
    }

    if hasHom {
      homs = strings.Split(record[homIdx], pDelim)

      for _, sample := range homs {
        if hetHomMap[sample] == nil {
          hetHomMap[sample] = []int{0, 0}
        }

        hetHomMap[sample][1]++
        hetHomMap[totalKey][1]++
      }
    }

    // If not transition or transversion we're in an indel
    // which uses different delimiters, and which we're not as interested in
    if record[trTvIdx] == parse.NotTrTv {
      continue
    }

    isTr = record[trTvIdx] == parse.Tr

    if hasDbSnpColumn {
      inDbSnp = record[dbSNPnameIdx] != emptyField
    }

    // If siteType doesn't have the delimiter, it's a single site type
    if !strings.Contains(record[siteTypeIdx], pDelim) {
      // no site type if intergenic (b11.0.0 on)
      if record[siteTypeIdx] == emptyField {
        siteTypes = intergenic
      } else {
        siteTypes = []string{record[siteTypeIdx]}
      }
    } else {
      siteTypes = uniqSlice(record[siteTypeIdx], emptyField, pDelim)
    }

    if hasExonicColumn {
      if !strings.Contains(record[funcIdx], pDelim) {
        if record[funcIdx] == emptyField {
          exonicTypes = empty
        } else {
          exonicTypes = []string{record[funcIdx]}
        }
      } else {
        exonicTypes = uniqSlice(record[funcIdx], emptyField, pDelim)
      }
    }

    fillType(fakeTotal, siteTypes, exonicTypes,
      trTvMap, featureCache, isTr, inDbSnp, emptyField)

    if hasHet {
      fillType(hets, siteTypes, exonicTypes, trTvMap, featureCache, isTr, inDbSnp, emptyField)
    }

    if hasHom{
      fillType(homs, siteTypes, exonicTypes, trTvMap, featureCache, isTr, inDbSnp, emptyField)
    }
  }

  trTvChan <- trTvMap
  hetHomChan <- hetHomMap

  complete <- totalSites
}

func fillType(samples []string, siteTypes []string, exonicTypes []string,
trTvMap map[string]map[string][]int, featureCache map[string]string,
isTr bool, inDbSnp bool, emptyField string) {
  for _, sample := range samples {
    if trTvMap[sample] == nil {
      trTvMap[sample] = make(map[string][]int)
    }
    
    if trTvMap[sample][totalKey] == nil {
      trTvMap[sample][totalKey] = []int{0,0}
    }
    
    if isTr {
      trTvMap[sample][totalKey][0]++
    } else {
      trTvMap[sample][totalKey][1]++
    }

    fillInner(siteTypes, trTvMap[sample], featureCache, isTr, inDbSnp)
    fillInner(exonicTypes, trTvMap[sample], featureCache, isTr, inDbSnp)
  }
}

func fillInner(siteTypes []string, sMap map[string][]int, featureCache map[string]string,
isTr bool, inDbSnp bool) {
  for _, siteType := range siteTypes {
    if sMap[siteType] == nil {
      sMap[siteType] = []int{0,0}
    }

    if isTr {
      sMap[siteType][0]++
    } else {
      sMap[siteType][1]++
    }

    if featureCache[siteType] == "" {
      var name bytes.Buffer
      name.WriteString(siteType)
      name.WriteString(dbSnpKey)

      featureCache[siteType] = name.String()
    }

    if inDbSnp {
      if sMap[featureCache[siteType]] == nil {
        sMap[featureCache[siteType]] = []int{0,0}
      }

      if isTr {
        sMap[featureCache[siteType]][0]++
      } else {
        sMap[featureCache[siteType]][1]++
      }
    }
  }
}

// Will mutate ratios by sorting, but who cares
func stats(ratios []float64, baseKey string) (map[string]float64) {
  meanKey := baseKey + " mean"
  medKey := baseKey + " median"
  sdKey := baseKey + " sd"

  data := map[string]float64 {
    meanKey: 0,
    medKey: 0,
    sdKey: 0,
  }

  if len(ratios) > 0 {
    sort.Slice(ratios, func(a, b int) bool {
      return ratios[a] < ratios[b];
    });

    data[meanKey] = mean(ratios)
    data[medKey] = median(ratios)

    if data[meanKey] != 0 {
      data[sdKey] = stdDev(ratios, data[meanKey])
    }
  }

  return data
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

func roundVal (val jsonFloat) string {
  if float64(val) == float64(int64(val)) {
    return strconv.Itoa(int(val))
  }

  return strconv.FormatFloat(float64(val), 'f', 3, 64)
}