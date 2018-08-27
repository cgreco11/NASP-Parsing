import pandas as pd
import sys, subprocess, argparse, shlex, io, os, shutil
from Bio import SeqIO

def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def getSNPeff():
	subprocess.call(["wget", "-q", "http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip"])
	subprocess.call(["unzip", "-qq", "snpEff_latest_core.zip"])
	subprocess.call(["rm", "snpEff_latest_core.zip"])
	subprocess.call(['rm', "-r", "clinEff/"])

def create_group_mapping(cluster_file):
	genome_to_group_dict = {} # Genome : Group
	with open(cluster_file, 'r') as f:
		for row in f:
			genome, group = row.split("\t")[0].strip(), row.split("\t")[1].strip()
			genome_to_group_dict[genome] = group
	return genome_to_group_dict

def readGoStdOut(group_file, threshold, tsv_file, strain_name):
    p = subprocess.call([get_script_path() + '/bin/parse_nasp', '-groupFile', group_file, '-threshold', threshold, '-tsvFile', tsv_file, "-strainName", strain_name])
    dataframe = pd.read_csv("annotated_bestsnp.tsv", sep = "\t", header = 0, low_memory = False)
    return dataframe

def formatDataframe(dataframe, groupMap):
	columns_to_drop = ["#Indelcall", "#CallWasMade", "#PassedDepthFilter", "#PassedProportionFilter", "#Indel", "#NXdegen", "Contig", "InDupRegion", "SampleConsensus", "CallWasMade", "PassedDepthFilter", "PassedProportionFilter", "Pattern", "Position"]
	dataframe =dataframe.drop(columns_to_drop, axis = 1)
	col_names = list(dataframe)
	gene_groups = []
	for i in col_names:
		if i in groupMap:
			gene_groups.append(groupMap[i])
		else:
			gene_groups.append("")

	dataframe.columns = pd.MultiIndex.from_tuples(zip(col_names, gene_groups))
	return dataframe

def createATTfile(gb_file):
    with open('genome.list', 'w') as f:
        f.write(gb_file)
    subprocess.call(['perl', get_script_path() + '/bin/parse_genbank_files.pl', '--file_list', "genome.list", '--no_dos2unix', '--nuc', '--no_check'], shell = False)

def createReferenceFasta(gb_file, strain_name):
	headers = []
	seqs = []
	for record in SeqIO.parse(gb_file, 'gb'):
		header = ">" + record.id + " " + record.description
		sequence = str(record.seq)
		headers.append(header)
		seqs.append(sequence)

	with open(strain_name + ".fa", 'w') as seq_file:
		for i in range(0, len(headers)):
			seq_file.write(headers[i] + "\n")
			seq_file.write(seqs[i] + "\n")

def createSNPeffDirectories(strain_name, genbank_file):
	current_directory = os.getcwd()
	if not os.path.exists(current_directory + "/snpEff/data/"):
		os.makedirs("snpEff/data/")
		os.makedirs("snpEff/data/genomes/")
		os.makedirs("snpEff/data/" + strain_name)
		shutil.copy(genbank_file, "snpEff/data/" + strain_name + "/" + "genes.gbk")
		shutil.copy2(strain_name + ".fa", "snpEff/data/genomes/")
	else:
		subprocess.call(["rm", "-r", "snpEff/data/"])
		createSNPeffDirectories(strain_name, genbank_file)

def editGBfile(genbank_file, chromosome_dict, strain_name):
	records = []
	for record in SeqIO.parse(genbank_file, 'gb'):
		record.name = chromosome_dict[record.name]
		records.append(record)
	SeqIO.write(records, "snpEff/data/" + strain_name + "/genes.gbk", "genbank")

def getChromosomeNames(strain_name):
	chrom_dict = {}
	chromosomes = []
	for record in SeqIO.parse(strain_name + ".fa", 'fasta'):
		chromosomes.append(record.id)
		chrom_dict[record.id.rsplit(".", 1)[0]] = record.id
	return chromosomes, chrom_dict

def createConfigString(chrom_list, strain_name):
	str1 = strain_name + ".genome : " + strain_name
	str2 = "\t" + strain_name + ".chromosomes : " + ", ".join(chrom_list)
	str3 = []
	for i in chrom_list:
		str3.append("\t" + strain_name + "." + i + ".codonTable : Bacterial_and_Plant_Plastid")

	str3 = "\n".join(str3)
	outStr = str1 + "\n" + str2 + "\n" + str3 + "\n"
	return outStr

def editConfigFile(config_file, config_string):
	with open(config_file, 'a') as f:
		f.write(config_string)

def buildSNPeffDatabase(config_file, strain_name):
	subprocess.call(["java", "-jar", "snpEff/snpEff.jar", "build", "-genbank", "-v", strain_name])

def runSNPeff(config_file, strain_name, vcf_file):
	p = subprocess.Popen(["java", "-Xmx4g", "-jar", "snpEff/snpEff.jar", "ann", "-classic", "-no-downstream", "-no-upstream", "-no-utr", "-no-intergenic", "-c", "snpEff/snpEff.config", strain_name, vcf_file], stdout=subprocess.PIPE)
	out = p.stdout.read().decode()
	with open("bestsnpAnnotated.vcf", 'w') as out_vcf:
		out_vcf.write(out)


def main():
    parser = argparse.ArgumentParser(description = "Parse NASP Output.")
    parser.add_argument('--strain_name', help = 'Strain Name of Reference Genome', required = True)
    parser.add_argument('-g', '--genbank_file', help = 'Reference RefSeq file from NCBI', required = True)
    parser.add_argument('-m', '--matrix_file', help = 'nasp_results/matrices/bestsnp.tsv', required = True)
    parser.add_argument("-b", '--vcf_file', help = 'nasp_results/matrices/bestsnp.vcf', required = True)
    parser.add_argument("-c", '--cluster_file', help = 'Tab separated grouping file with 2 columns -- Genome Name and then whatever group that genome belongs to', required = True)
    parser.add_argument("-t", '--threshold', help = 'Threshold at which SNP Groups are reported (100 means that all genomes in group must have SNP for it to be reported) -- Default is 90', default = 90)
    parser.add_argument("-o", '--out_file', help = 'Out File Name', default = 'nasp_parse_results.txt')

    args = parser.parse_args()
    strain_name = args.strain_name
    genbank_file = args.genbank_file
    tsv_file = args.matrix_file
    vcf_file = args.vcf_file
    group_file = args.cluster_file
    out_file = args.out_file
    threshold = str(float(args.threshold))
    config_file = "snpEff/snpEff.config"

    #Install SNPeff directory locally
    try:
    	getSNPeff()
    except:
    	print
    	print "ERROR: Unable to install snpEff"
    	print
    	sys.exit()

	#Set Up SNPEff run
    createReferenceFasta(genbank_file, strain_name)
    createSNPeffDirectories(strain_name, genbank_file)
    chrom_list, chrom_dict = getChromosomeNames(strain_name)
    groupMappingDict = create_group_mapping(group_file)
    print "Setting up snpEff run."
    editGBfile('snpEff/data/' + strain_name + '/genes.gbk', chrom_dict, strain_name)
    config_string = createConfigString(chrom_list, strain_name)
    editConfigFile(config_file, config_string)
    buildSNPeffDatabase(config_file, strain_name)
    print "Running snpEff."
    runSNPeff(config_file, strain_name, vcf_file)

    #Parse RefSeq File
    gb_prefix = genbank_file.split(".")[0]

    createATTfile(genbank_file)
    att_file = gb_prefix + ".att"

    #Create DataFrame
    print "Parsing NASP results and adding snpEff results."
    df = readGoStdOut(group_file, threshold, tsv_file, strain_name)
    df = formatDataframe(df, groupMappingDict)
    df.fillna("").sort_values(by = ["Groups Passing Threshold", "SNP's Genome Location"]).to_csv(out_file, sep = "\t", index = False)


if __name__ == '__main__':
	main()
