#!/usr/bin/python
# CharGer - Characterization of Germline variants
# author: Adam D Scott (ascott@genome.wustl.edu)
# version: v0.0 - 2015*12

import sys
import getopt
from entrezAPI import entrezAPI
import variant

def parseArgs( argv ):
	helpText = "python main.py" + " "
	helpText += "-i \"input (.maf)\" "
	helpText += "(-o \"output\")\n"
	inputFile = ""
	output = ""
	database = ""
	query = ""
	try:
		opts, args = getopt.getopt( argv , "h:i:o:" , ["input=" , "output="] )
	except getopt.GetoptError:
		print "ADSERROR: Command not recognized"
		print( helpText ) 
		sys.exit(2)
	if not opts:
		print "ADSERROR: Expected flagged input"
		print( helpText ) 
		sys.exit(2)
	for opt, arg in opts:
		#print opt + " " + arg
		if opt in ( "-h" , "--help" ):
			print( helpText )
			sys.exit()
		elif opt in ( "-i" , "--inputFile" ):
			inputFile = arg
		elif opt in ( "-o" , "--output" ):
			output = arg
	return { "input" : inputFile , "output" : output }
	
def checkConnection():
	print "\tChecking Connection!"
	entrezInstance = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=asthma[mesh]+AND+leukotrienes[mesh]+AND+2009[pdat]"
	summaryTest = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=20113659,20074456"
	res = requests.get( entrezInstance )
	if res:
		print "have response"
	else:
		print res.status_code

def splitByVariantType( inputFile ):
	print "\tSplitting .maf by variant type!"
	variants = {}
	if inputFile:
		inFile = open( inputFile , 'r' )
		next(inFile)
		gv = {}
		for line in inFile:
			fields = line.split( "\t" )
			var = variant.variant()
			var.mafLine2Variant( line )
			#var.printVariant("\t")
			sample = fields[15]
			variantClass = var.variantClass
			genVar = var.genomicVar()
			vc = {}
			s = {}
			gv = {}
			if variantClass in variants:
				#print variantClass
				vc = variants[variantClass] #dict of variants by variantClass=>samples dict
			if sample in vc:
				#print sample
				s = vc[sample] #dict of samples by sample=>genVars dict
			if genVar in s:
				#print genVar
				gv = s[genVar]
			tempgv = { "line" : line , "variant" : var }
			gv.update( tempgv )
			temps = { genVar : gv }
			s.update( temps )
			tempvc = { sample : s }
			vc.update( tempvc )
			#print genVar + "=>" + s[genVar]
			tempVariants = { variantClass : tempvc }
			variants.update( tempVariants )
	return variants
	#nonsense, frameshift, canonical 1 or 2 splice sites, initiation codon, single exon or multiexon deletion

def PVS1( variants ):
	return None
def PM1( entrezAPI , variants ):
	for uid in variants:
		var = variants[uid]
		var.printVariant(' ')
	
def main( argv ):
	values = parseArgs( argv )
	inputFile = values["input"]
	outputFormat = values["output"]

	ent = entrezAPI()	

	variants = splitByVariantType( inputFile )
	for variantClass in variants:
		#print "\tVariant Class: " + variantClass
		for sample in variants[variantClass]:
			#print "\t\tSample: " + sample
			#variants[variantClass][sample][
			for genVar in variants[variantClass][sample]:
				#print "\t\t\tGenomic Variant: " + genVar
				#print "\t\t\t\tVariant: " , 
				#variants[variantClass][sample][genVar]["variant"].printVariant("  ")
				thisGroup = sample+genVar
				var = variants[variantClass][sample][genVar]["variant"]
				ent.addQuery( var.gene , field="gene" , group=thisGroup )
				ent.addQuery( var.chromosome , field="chr" , group=thisGroup )
				ent.addQuery( var.start + ":" + var.stop , field="chrpos37" , group=thisGroup )
				ent.addQuery( "human" , field="orgn" , group=thisGroup )
				#ent.addQuery( var.variantClass , "vartype" )
				#ent.addQuery( var.referencePeptide + var.positionPeptide + var.mutantPeptide , "Variant name" )

				#var.referencePeptide , var.positionPeptide , var.mutantPeptide

	ent.database = entrezAPI.clinvar
	clinvarEntries = ent.doBatch( 5 )
	ClinVarVariants = clinvarEntries["variants"]
	ClinVarTraits = clinvarEntries["traits"]
	PM1( ent , ClinVarVariants )

if __name__ == "__main__":
	main( sys.argv[1:] )
