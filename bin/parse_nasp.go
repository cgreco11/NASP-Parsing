package main

import (
    "fmt"
    "os"
    "bufio"
    "strings"
    "strconv"
    "sort"
    "flag"
)

func readLines(path string) ([]string, error) {
  file, err := os.Open(path)
  if err != nil {
    return nil, err
  }
  defer file.Close()

  var lines []string
  scanner := bufio.NewScanner(file)
  for scanner.Scan() {
    lines = append(lines, scanner.Text())
  }
  return lines, scanner.Err()
}

//------------------------------------------------------------------------------------//
// IO Functions
// Parsing TSV file
func parseTSV(tsv_file string) (map[int]string, map[string]string){
  lines, err := readLines(tsv_file)
  if err != nil {
    fmt.Println("Problem reading in TSV File")
  }
  tsvHeader := make(map[int]string)
  tsvMatrix := make(map[string]string)
  for line_num, line := range lines{
    if line_num == 0{
    for i, data := range strings.Split(line, "\t"){
      tsvHeader[i] = data
    }} else{
        key := strings.Split(line, "\t")[0]
        value := strings.Join(strings.Split(line, "\t")[1:], "\t")
        tsvMatrix[key] = value
  }}
  return tsvHeader, tsvMatrix
}

//Parsing Group Files
func parseGroupFile(group_file string) (map[string]string){
  group_dict := make(map[string][]string) // Group : Genomes in Group
  genome_dict := make(map[string]string) // Genome : Group
  lines, err := readLines(group_file)
  if err != nil{
    fmt.Println("Problem reading in Group File")
  }
  for _, line := range lines{
    group_data := strings.Split(line, "\t")
    group_dict[group_data[1]] = append(group_dict[group_data[1]], group_data[0])
  }
  for k, v := range(group_dict){
      for _,b := range(v){
        genome_dict[b] = k
      }
  }

  return genome_dict
}
//Parsing ATT Files
func parseAttFile(att_file string)(map[int]map[string]string){
  att_dict := make(map[int]map[string]string) //Contig : Locus Start End Protein_Product
  lines, err := readLines(att_file)
  if err != nil{
    fmt.Println(err)
    os.Exit(2)
  }

  for i, line := range lines{
    att_data := strings.Split(line, "\t")
    contig := att_data[0]
    locus := att_data[1]
    start,_ := strconv.Atoi(att_data[2])
    end,_ := strconv.Atoi(att_data[3])
    product := att_data[4]
    if att_dict[i] == nil{
      att_dict[i] = make(map[string]string)
    }
    if start >= end{
        start, end := end, start
        out_start := strconv.Itoa(start)
        out_end := strconv.Itoa(end)
        att_dict[i]["Contig"] = contig
        att_dict[i]["Start"] = out_start
        att_dict[i]["End"] = out_end
        att_dict[i]["Locus"] = locus
        att_dict[i]["Product"] = product
    } else {
    out_start := strconv.Itoa(start)
    out_end := strconv.Itoa(end)
    att_dict[i]["Contig"] = contig
    att_dict[i]["Locus"] = locus
    att_dict[i]["Start"] = out_start
    att_dict[i]["End"] = out_end
    att_dict[i]["Product"] = product
  }}
  return att_dict
}
//-------------------------------------------------------------------------------------//
//Add Annotation, Length, Locus to TSV Matrix
func annotateTSV(tsvMatrix map[string]string, attMatrix map[int]map[string]string, tsvHeader map[int]string) (map[string]string, map[int]string) {
  tsv_dict := make(map[string]string)
  header_dict := make(map[int]string)
  for k,v := range tsvMatrix{
    split_v := strings.Split(v, "\t")
    genome_loc,_ := strconv.Atoi(strings.Split(k, "::")[1])
    contig_with_version := strings.Split(k, "::")[0]
    contig := strings.Split(contig_with_version, ".")[0]
    var locus []string
    var product []string
    var length []string
    var locations []string
    var snpPosition []string

    for i,_ := range attMatrix{ // Account for overlapping contig positions? ATT File designation
      start,_ := strconv.Atoi(attMatrix[i]["Start"])
      end,_ := strconv.Atoi(attMatrix[i]["End"])
      att_contig := attMatrix[i]["Contig"]
      if (start <= genome_loc) && (genome_loc <= end) && (contig == att_contig){
        product = append(product, attMatrix[i]["Product"])
        locus = append(locus, attMatrix[i]["Locus"])
        length = append(length, strconv.Itoa((end - start)))
        snpPosition = append(snpPosition, strconv.Itoa((genome_loc - start)))
        locations = append(locations, (strconv.Itoa(start) + "-" + strconv.Itoa(end)))
      }
    }
    new_locus := strings.Join(locus, "/")
    new_product := strings.Join(product, "/")
    new_length := strings.Join(length, "/")
    new_location := strings.Join(locations, "/")
    new_snpPosition := strings.Join(snpPosition, "/")
    new_GenomeLoc := strings.Split(k, "::")[1]
    split_v = append(split_v, new_GenomeLoc, new_locus, new_product, new_length, new_location, new_snpPosition)
    tsv_dict[k] = strings.Join(split_v, "\t")
  }
  header_dict = tsvHeader
  headerDictLength := len(header_dict)
  header_dict[headerDictLength + 1] = "SNP's Genome Location"
  header_dict[headerDictLength + 2] = "Locus"
  header_dict[headerDictLength + 3] = "Protein Product"
  header_dict[headerDictLength + 4] = "Locus Length"
  header_dict[headerDictLength + 5] = "Locus Locations"
  header_dict[headerDictLength + 6] = "SNP Location within Locus"
  return tsv_dict, header_dict
}

func groupToIndex(tsvHeaders map[int]string, groupDict map[string]string) (map[string][]int){
  outDict := make(map[string][]int)
  for _,v := range groupDict{
    var indices []int
    for index, info := range tsvHeaders{
      if groupDict[info] == v{
          indices = append(indices, index)
        }
      }
      outDict[v] = indices
    }
  return outDict
}
//Convert this function to take in a row from Annotated TSV
func groupStatsPerLocus(groupIndexDict map[string][]int, tsvMatrix map[string]string, thresholdVal float64, tsvHeader map[int]string) (map[string]string, string) {

  var keys []string
  for k := range groupIndexDict{
    keys = append(keys, k)
  }
  sort.Strings(keys)
  for _, val := range keys{
    headerLength := len(tsvHeader)
    tsvHeader[headerLength + 1] = "Group " + val + " SNP Counts"
  }
  for k,_ := range tsvMatrix{
    rowInfo := tsvMatrix[k]
    splitRow := strings.Split(rowInfo, "\t")
    refAllele := splitRow[0]
    var groupInfo []string
    var groupPerc []string
    //fmt.Println(refAllele, k)
    for _,group := range keys{
      alleleCounts := make(map[string]string)
      groupIndices := groupIndexDict[group]
      var count float64 = 0.0
      for _,individualIndex := range groupIndices{
        //fmt.Println(individualIndex, tsvHeader[individualIndex])
        //fmt.Println(splitRow[individualIndex])
        allele := splitRow[individualIndex - 1] // Header is 1 based (Don't want to keep index column name), tsvHeader is Zero Based -- Need to subtract 1 to only get Reference => Last Genome (Not #SNP Call)
        //Reference Genome is Zero, would be 1 if we didn't subtract 1
        if allele != refAllele {
          count = count + 1.0
          if val, ok := alleleCounts[allele]; ok {
              new_val,_ := strconv.Atoi(val)
              alleleCounts[allele] = strconv.Itoa(new_val + 1)
            } else {
              alleleCounts[allele] = strconv.Itoa(1)
            }
          }

        }
        floatGroupIndices := float64(len(groupIndices))
        snpPercentage := (count / floatGroupIndices) * 100
        if snpPercentage > thresholdVal{
          groupInfo = append(groupInfo, group)
        } else {
          groupInfo = append(groupInfo, "Not " + group)
        }
        s := strconv.FormatFloat(snpPercentage, 'f', 2, 64)
        outPerc := "Group " + group + "-- " + s
        groupPerc = append(groupPerc, outPerc)
        //fmt.Println(alleleCounts, group)
        var groupCounts []string
        if len(alleleCounts) == 0{
          groupCounts = append(groupCounts, "")
        } else {
            for key, value := range alleleCounts {
            str := key + ":" + value
            groupCounts = append(groupCounts,  str)
          }
        }
      outGroupCounts := strings.Join(groupCounts, ",")
      if len(outGroupCounts) == 0 {
      splitRow = append(splitRow, "")
    } else {
      splitRow = append(splitRow, outGroupCounts)
      }
    }
    splitRow = append(splitRow, strings.Join(groupInfo, ","), strings.Join(groupPerc, " "))
    tsvMatrix[k] = strings.Join(splitRow, "\t")
  }
  tsvFinalIndex := len(tsvHeader)
  tsvHeader[tsvFinalIndex + 2] = "Group SNP Percentages"
  tsvHeader[tsvFinalIndex + 1] = "Groups Passing Threshold"

  var intKeys []int
  var values []string
  for k := range tsvHeader{
    intKeys = append(intKeys, k)
  }
  sort.Ints(intKeys)
  for _,v := range intKeys{
    values = append(values, tsvHeader[v])
  }
  outHeader := strings.Join(values, "\t")
  return tsvMatrix, outHeader
}

//------------------------------------------------------------------------------------//
func main(){

  tsvFilePtr := flag.String("tsvFile", "", "The bestsnp.tsv file from NASP results")
  thresholdPtr := flag.Float64("threshold", 90, "Threshold at which to cut off SNP reporting-- Default 90 (90 Percent of genomes must have a SNP for reporting to occur)")
  attFilePtr := flag.String("attFile", "", "The ATT file for Annotation-- Can be obtained from the GenBank downloader/parser")
  groupFilePtr := flag.String("groupFile", "", "The Group File for Annotation")
  flag.Parse()

  tsvFile := *tsvFilePtr
  threshold := *thresholdPtr
  attFile := *attFilePtr
  groupFile := *groupFilePtr
  //fmt.Println(tsvFile, threshold)

  if tsvFile == "" {
    fmt.Println()
    fmt.Println()
    fmt.Println("ERROR: Must include bestsnp.tsv for the command line argument -tsvFile")
    fmt.Println()
    os.Exit(1)
  }

  if attFile == "" {
    fmt.Println()
    fmt.Println()
    fmt.Println("ERROR: Must include .att file for the command line argument -attFile")
    fmt.Println()
    os.Exit(1)
  }

  if groupFile == "" {
    fmt.Println()
    fmt.Println()
    fmt.Println("ERROR: Must include group file for the command line argument -attFile")
    fmt.Println()
    os.Exit(1)
  }
  tsvHeaderNames, tsvMatrixValues := parseTSV(tsvFile)
  //groupToGenome, genomeToGroup := parseGroupFile("group_subclade.txt")
  genomeToGroup := parseGroupFile(groupFile)
  attInfo := parseAttFile(attFile)
  annotatedTsvMatrix, annotatedTsvHeader := annotateTSV(tsvMatrixValues, attInfo, tsvHeaderNames)
  groupToIndexDict := groupToIndex(annotatedTsvHeader, genomeToGroup)
  body, header := groupStatsPerLocus(groupToIndexDict, annotatedTsvMatrix, threshold, annotatedTsvHeader)
  fmt.Println(header)

  for k,v := range body{
    fmt.Printf("%s\t%s\n", k, v)
  }
}
