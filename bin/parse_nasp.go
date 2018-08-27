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

func parseSnpEff(snpEffFile string) (map[string]map[string]string){
  outMap := make(map[string]map[string]string)
  lines, err := readLines(snpEffFile)
  if err != nil {
    fmt.Println("Problem reading snpEff output file.")
  }
  for _, line := range lines {
    internalDict := make(map[string]string)
    if !(strings.HasPrefix(line, "#")){
      splitStrings := strings.Split(line, "\t")
      outKey := splitStrings[0] + ":" + splitStrings[1]
      info := splitStrings[7]
      if strings.Contains(info, "("){
        realInfo := strings.Split(info, "(")[1]
        infoOfInterest := strings.Split(realInfo, "|")
        internalDict["Effect Impact"] = infoOfInterest[0]
        internalDict["SNP Effect"] = infoOfInterest[1]
        internalDict["Nucleotide Change"] = infoOfInterest[2]
        internalDict["Peptide Change"] = infoOfInterest[3]
        outMap[outKey] = internalDict
      }
    }
  }
  return outMap
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
func parseAttFile(att_file string)(map[string][]map[string]string){
  contig_array := make(map[string][]map[string]string)
  lines, err := readLines(att_file)
  if err != nil{
    fmt.Println(err)
    os.Exit(1)
  }
  for _, line := range lines{
    att_dict := make(map[string]string) //Contig : Locus Start End Protein_Product
    att_data := strings.Split(line, "\t")
    contig := att_data[0]
    locus := att_data[1]
    start,_ := strconv.Atoi(att_data[2])
    end,_ := strconv.Atoi(att_data[3])
    product := att_data[4]
    if start >= end{
        start, end := end, start
        out_start := strconv.Itoa(start)
        out_end := strconv.Itoa(end)
        att_dict["Orientation"] = "-"
        att_dict["Start"] = out_start
        att_dict["End"] = out_end
    } else {
    out_start := strconv.Itoa(start)
    out_end := strconv.Itoa(end)
    att_dict["Orientation"] = "+"
    att_dict["Start"] = out_start
    att_dict["End"] = out_end
  }
    att_dict["Product"] = product
    att_dict["Locus"] = locus

    contig_array[contig] = append(contig_array[contig], att_dict)
  }
  return contig_array
}
//-------------------------------------------------------------------------------------//
//Add Annotation, Length, Locus to TSV Matrix
func annotateTSV(tsvMatrix map[string]string, attMatrix map[string][]map[string]string, tsvHeader map[int]string, snpEffDict map[string]map[string]string) (map[string]string, map[int]string) {
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
    var orientation []string
    for _, dict := range attMatrix[contig] { // Account for overlapping contig positions? ATT File designation
      start,_ := strconv.Atoi(dict["Start"])
      end,_ := strconv.Atoi(dict["End"])
      if (start <= genome_loc) && (genome_loc <= end){
        product = append(product, dict["Product"])
        locus = append(locus, dict["Locus"])
        length = append(length, strconv.Itoa((end - start)))
        snpPosition = append(snpPosition, strconv.Itoa((genome_loc - start)))
        locations = append(locations, (strconv.Itoa(start) + "-" + strconv.Itoa(end)))
        orientation = append(orientation, dict["Orientation"])
      }
    }
    snpKeySplit := strings.Split(k, "::")
    snpKey := snpKeySplit[0] + ":" + snpKeySplit[1]

    new_locus := strings.Join(locus, "/")
    new_product := strings.Join(product, "/")
    new_length := strings.Join(length, "/")
    new_location := strings.Join(locations, "/")
    new_snpPosition := strings.Join(snpPosition, "/")
    new_GenomeLoc := strings.Split(k, "::")[1]
    new_orientation := strings.Join(orientation, "/")
    split_v = append(split_v, new_GenomeLoc, new_locus, new_product, new_length, new_location, new_snpPosition, new_orientation)
    if _, ok := snpEffDict[snpKey]; ok {
      snpInfo := snpEffDict[snpKey]
      new_effectImpact := snpInfo["Effect Impact"]
      new_snpEffect := snpInfo["SNP Effect"]
      new_nucleotideChange := snpInfo["Nucleotide Change"]
      new_peptideChange := snpInfo["Peptide Change"]
      split_v := append(split_v, new_effectImpact, new_snpEffect, new_nucleotideChange, new_peptideChange)
      tsv_dict[k] = strings.Join(split_v, "\t")
    } else {
      new_effectImpact := ""
      new_snpEffect := ""
      new_nucleotideChange := ""
      new_peptideChange := ""
      split_v := append(split_v, new_effectImpact, new_snpEffect, new_nucleotideChange, new_peptideChange)
      tsv_dict[k] = strings.Join(split_v, "\t")
    }

  }
  header_dict = tsvHeader
  headerDictLength := len(header_dict)
  header_dict[headerDictLength + 1] = "SNP's Genome Location"
  header_dict[headerDictLength + 2] = "Locus"
  header_dict[headerDictLength + 3] = "Protein Product"
  header_dict[headerDictLength + 4] = "Locus Length"
  header_dict[headerDictLength + 5] = "Locus Locations"
  header_dict[headerDictLength + 6] = "SNP Location within Locus"
  header_dict[headerDictLength + 7] = "DNA Strand Orientation"
  header_dict[headerDictLength + 8] = "Effect Impact"
  header_dict[headerDictLength + 9] = "SNP Effect"
  header_dict[headerDictLength + 10] = "Nucleotide Change"
  header_dict[headerDictLength + 11] = "Peptide Change"
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
  tsvHeader[tsvFinalIndex + 1] = "Groups Passing Threshold"
  tsvHeader[tsvFinalIndex + 2] = "Group SNP Percentages"


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

func writeOutput(header string, body map[string]string){
  f, err := os.Create("annotated_bestsnp.tsv")
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }
  defer f.Close()
  _, err = f.WriteString(header)
  if err != nil {
    fmt.Println("Error writing annotated_bestsnp.tsv")
  }
  for k,v := range(body){
  _, err = f.WriteString(fmt.Sprintf("%s\t%s\n", k, v))
  if err != nil {
    fmt.Println("Error writing annotated_bestsnp.tsv")
  }
  }
}
//------------------------------------------------------------------------------------//
func main(){

  tsvFilePtr := flag.String("tsvFile", "", "The bestsnp.tsv file from NASP results")
  thresholdPtr := flag.Float64("threshold", 90, "Threshold at which to cut off SNP reporting-- Default 90 (90 Percent of genomes must have a SNP for reporting to occur)")
  //attFilePtr := flag.String("attFile", "", "The ATT file for Annotation-- Can be obtained from the GenBank downloader/parser")
  groupFilePtr := flag.String("groupFile", "", "The Group File for Annotation")
  strainNamePtr := flag.String("strainName", "", "The RefSeq GenBank file prefix of the reference.")
  flag.Parse()

  tsvFile := *tsvFilePtr
  threshold := *thresholdPtr
  //attFile := *attFilePtr
  strainName := *strainNamePtr
  groupFile := *groupFilePtr
  //fmt.Println(tsvFile, threshold)

  if tsvFile == "" {
    fmt.Println()
    fmt.Println()
    fmt.Println("ERROR: Must include bestsnp.tsv for the command line argument -tsvFile")
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
  //scriptPath := getScriptPath()
  //makeGenomeList(genbankFile)
  //parseGenBankFile(genbankFile, scriptPath)
  attFile := strainName + ".att"
  snpEffFile := "bestsnpAnnotated.vcf"
  snpEffDict := parseSnpEff(snpEffFile)
  tsvHeaderNames, tsvMatrixValues := parseTSV(tsvFile)
  genomeToGroup := parseGroupFile(groupFile)
  attInfo := parseAttFile(attFile)
  annotatedTsvMatrix, annotatedTsvHeader := annotateTSV(tsvMatrixValues, attInfo, tsvHeaderNames, snpEffDict)
  groupToIndexDict := groupToIndex(annotatedTsvHeader, genomeToGroup)
  body, header := groupStatsPerLocus(groupToIndexDict, annotatedTsvMatrix, threshold, annotatedTsvHeader)
  writeOutput(header, body)
}
