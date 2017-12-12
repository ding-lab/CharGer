#!/usr/bin/env python
# analyze_charger_summary.py
# Count membership numbers of pathways in a gene list file

import sys
import getopt
import glob
from scipy import stats
#from autovivification import autovivification as AV

path2count = {}

def main( argv ):

	helpText = "analyze_charger_summary.py\n\n\n"
	helpText += "USAGE:\npython analyze_charger_summary.py -s <charger summary>\n"

	# default option
	sumFile = "/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/MSigDB/c2.cp.v5.1.symbols.gmt" 
	# headLine = delim.join( ["HUGO_Symbol" , "Chromosome" , "Start" , \
	# 		"Stop" , "Reference" , "Alternate" , "Strand" , "Assembly" , \
	# 		"Variant_Type" , "Variant_Classification" , \
	# 		"Sample" , "Transcript" , "Codon_Position" , "HGVSg", "Protein" , \
	# 		"Peptide_Reference" , "Peptide_Position" , "Peptide_Alternate" , \
	# 		"HGVSp","Allele_Frequency","VEP_Most_Severe_Consequence" , "ClinVar_Pathogenicity" , \
	# 		"Positive_Evidence" , "Negative_Evidence" , \
	# 		"Positive_CharGer_Score" , "Negative_CharGer_Score" , \
	# 		"CharGer_Classification" , "ACMG_Classification" , \
	# 		"PubMed_Link" , "ClinVar_Traits"] )
	posEvColNum = 
	negEvColNum =

	try:
		opts, args = getopt.getopt( argv , "hg:l:" )
		for opt, arg in opts:
			if opt in ( "-h" , "--help" ):
				print( helpText )
				sys.exit()
			elif opt in ( "-s" , "--gmt" ):
				sumFile = arg
			elif opt in ( "-l" , "--list" ):
				geneListFile = arg
	except getopt.GetoptError:
		print "Command not recognized"
		print( helpText ) 
		sys.exit(2)

	
	geneList = readGeneList( geneListFile )
	#print "Probe_Id	Symbol	SNPID	pvalue	beta	R.squared	t.test	n.samples	FDR	relativePosition	TAGID"
	gene2path = readGMT( gmtFile )

	pathCount = {}
	for gene in gene2path:
		if gene in geneList:
			paths = gene2path[gene]
			for path in paths:
				if path in pathCount:
					pathCount[path] += 1
				else:
					pathCount[path] = 1

	path_list = []
	d_view = [ (v,k) for k,v in pathCount.iteritems() ]
	d_view.sort(reverse=True) # natively sort tuples by first element
	for v,k in d_view:
		path_list.append(k)
	# 	print "%s: %d" % (k,v)


	for gene in geneList:
		if gene in gene2path:
			first_path = ""
			for path in path_list[]:
				if path in gene2path[gene]:
					first_path = path
					next

			#print gene + "\t" + ', '.join(gene2path[gene])
			print gene + "\t" + first_path)
		else: 
			print gene

	# inFile = open( inputFile , 'r' )

	
	# header = inFile.readline()
	# print header.rstrip() + "\texpression_quantile" 
	# next(inFile)

	# try:
	# 	for line in inFile:
	# 		fields = line.split( "\t" )
	# 		gene = fields[geneColumn]
	# 		sample = fields[sampleColumn]
	# 		var_exp = "NA"
	# 		if expression[sample][gene]:
	# 			var_exp = expression[sample][gene]
	# except: 
	# 	print "can't annotate file" + inFile


def readGMT( gmtFile ): # expect sample(col)-gene(row) matrixes

	gene2path = {}
	inFile = open( gmtFile , 'r' )
	#header = inFile.readline()
	#print header

	for line in inFile:
		line = line.strip()
		fields = line.split( "\t" )		
		path = fields[0]
		genes = fields[2:len(fields)]
		for gene in genes:
			if gene in gene2path:
				gene2path[gene].append(path)
			else: 
				gene2path[gene] = [path]

	inFile.close()

	return(gene2path)

def readGeneList( geneList ): # expect sample(col)-gene(row) matrixes

	genes = {}
	inFile = open( geneList , 'r' )
	#header = inFile.readline()
	#print header

	for line in inFile:
		line = line.strip()
		fields = line.split( "\t" )		
		gene = fields[0]
		genes[gene] = 1

	inFile.close()

	return(genes)


if __name__ == "__main__":
	main( sys.argv[1:] )
