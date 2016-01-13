#!/usr/bin/python
# CharGer - Characterization of Germline variants
# author: Adam D Scott (ascott@genome.wustl.edu) & Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2015*12

import sys
import getopt
import math
from entrezAPI import entrezAPI
from exacAPI import exacAPI
from variant import variant
from variant import MAFVariant
import autovivification

class CharGerVariant(object):
# TODO

def parseArgs( argv ):
	helpText = "python main.py" + " "
	helpText += "-i \"input (.maf)\" "
	helpText += "(-c suppress ClinVar, "
	helpText += "-x suppress ExAC, "
	helpText += "-o \"output\")\n"
	inputFile = ""
	output = ""
	clinvar = True
	exac = True
	try:
		opts, args = getopt.getopt( argv , "cxh:i:o:" , ["input=" , "output="] )
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
		elif opt in ( "-c" , "--help" ):
			clinvar = False
		elif opt in ( "-x" , "--help" ):
			exac = False
	return { "input" : inputFile , "output" : output , "clinvar" : clinvar , "exac" : exac }
	
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
	variants = []
	if inputFile:
		inFile = open( inputFile , 'r' )
		next(inFile)
		for line in inFile:
			fields = line.split( "\t" )
			var = MAFVariant()
			var.mafLine2Variant( line )
			variants.append( var )
	return variants
	#nonsense, frameshift, canonical 1 or 2 splice sites, initiation codon, single exon or multiexon deletion

### Evidence levels ### 
def PVS1( userVariants, geneListFile , expressionFile, expThres = 0.05 ):
	truncations = ["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Nonstop_Mutation","Splice_Site"]
	geneList = readGeneList( inputFile ) 
	expression = getExpression( expressionFile )
	calls = autovivification.autovivification({})
	
	if geneList: #gene, disease, mode of inheritance
		for var in userVariants:
			call = False # default
			varGene = var.gene
			varDisease = var.disease	
			varSample = var.sample
			varType = var.variantClass
			if varType in truncations:
				if varGene in geneList: # check if in gene list
					if ( "dominant" in geneList[varGene][varDisease] or "dominant" in geneList[varGene]["all"]):
						call = True # if call is true then check expression effect
						if expression: # consider expression data only if the user has supplied an expression matrix
							if expression[varSample][varGene] >= expThres:
								call = False 
				
			calls[var] = call

	else: 
		raise Exception("No gene list file supplied.")

	return calls

def PS1( inputVariants , clinvarVariants , clinvarClinical ):
	print "CharGer module PS1"
	print "- same peptide change that is pathogenic and is a different genomic variant of the same reference peptide"
	return peptideChange( inputVariants , clinvarVariants , clinvarClinical , "PS1" )

def PS4( inputVariants ):
	calls = autovivification.autovivification({})
	
	for var in userVariants:
		call = False # default

		caseVarFreq = "NEED UPDATE" # may take input from current MAF
		controlVarFreq = "NEED UPDATE" # may take input from ExAC
		OR = (caseVarFreq/controlVarFreq) / (1-caseVarFreq)/(1-controlVarFreq)
		 # Adam will update
		if OR >= 5:
			CIlower = math.log(OR) - math.sqrt( 1/caseVarFreq + 1/controlVarFreq + 1/caseVarFreq + 1/controlVarFreq)
			if (CIlower > 1):
				call = True
		
		calls[var] = call

	return calls

def PM2( inputVariants , clinvarVariants , clinvarClinical ):
	calls = autovivification.autovivification({})
	
	for var in userVariants:
		call = False # default
		varMAF = var.getExACasdf # Adam will update
		if isFrequentAllele(varMAF):
			call = True
		calls[var] = call

	return calls

def PM4( inputVariants ):
	lenShift = ["In_Frame_Del","In_Frame_Ins","Nonstop_Mutation"]
	calls = autovivification.autovivification({})
	varType = var.variantClass
	
	for var in userVariants:
		call = False # default
		if varType in lenShift:
			call = True
		calls[var] = call

	return calls

def PM5( inputVariants , clinvarVariants , clinvarClinical ):
	print "CharGer module PM5"
	print "- different peptide change of a pathogenic variant at the same reference peptide"
	return peptideChange( inputVariants , clinvarVariants , clinvarClinical , "PM5" )

### helper function of evidence levels ###
def peptideChange( inputVariants , clinvarVariants , clinvarClinical , mod ):
	calls = autovivification.autovivification({})
	for inVar in inputVariants:
		uniVar = inVar.uniqueVar()
		#print "\tInput variant: " + genVar , 
		canBePM1 = True
		canBePM5 = True
		pm1Call = False
		pm5Call = False
		call = calls[uniVar]
		#print "Call: " + genVar ,
		#print " => " + str(call)
		if not call: #is already true
			#print "checking"
			call = False
			for uid in clinvarVariants:
				var = clinvarVariants[uid]
				if uid in clinvarClinical:
					clin = clinvarClinical[uid]
					if inVar.chromosome == var.chromosome and \
						inVar.start == var.start and \
						inVar.stop == var.stop and \
						inVar.reference == var.reference and \
						inVar.referencePeptide == var.referencePeptide and \
						inVar.positionPeptide == var.positionPeptide: #same genomic position & reference
						if inVar.alternatePeptide == var.alternatePeptide: #same amino acid change
							if clin["description"] == "Pathogenic":
								#print "Already called pathogenic: " ,
								#var.printVariant(' ')
								canBePM1 = False
								canBePM5 = False
							else:
								#print "This is NOT called as pathogenic: " ,
								#var.printVariant(' ')
								if mod == "PM1":
									pm1Call = True
						else: #different amino acid change ( CAN BE USED FOR PM5 )
							if clin["description"] == "Pathogenic":
								#print "Alternate peptide change called pathogenic: " ,
								#var.printVariant(' ')
								if mod == "PM5":
									pm5Call = True
							else:
								print "" , 
								#print "Alternate peptide change NOT called as pathogenic: " ,
								#var.printVariant(' ')
				else:
					print "" , 
					#print "Not given a clinical call: " ,
					#var.printVariant(' ')
			if mod == "PM1":
				if canBePM1:
					call = pm1Call
			if mod == "PM5":
				if canBePM5:
					call = pm5Call
		calls[uniVar] = call
	return calls

def getExpression( inputFile ): # expect a sample(col)-gene(row) matrix
	expression = autovivification.autovivification({})
	
	if inputFile:
		inFile = open( inputFile , 'r' )
		header = inFile.readline() # for future fetch header to get other field
		samples = header.split( "\t" )

		for line in inFile:
			fields = line.split( "\t" )
			gene = fields[0]
			for i in range(1,len(fields)):
				expression[samples[i]][gene] = fields[i]

	return expression

def isFrequentAllele( freq , threshold ):
	if freq > threshold:
		return True
	return False

def readGeneList( inputFile, diseaseSpecific = True ): # gene list formatted "gene", "disease", "mode of inheritance"
	geneList = autovivification.autovivification({})
	if inputFile:
		inFile = open( inputFile , 'r' )
		#header = inFile.readline() # for future fetch header to get other fields
		#gv = {}
		for line in inFile:
			fields = line.split( "\t" )
			gene = fields[0]
			if diseaseSpecific:
				disease = fields[1]
			else: #set the gene to match all disease
				disease = "all" 
			mode_inheritance = fields[2]
			geneList[gene][disease] = mode_inheritance
	return geneList

def prepQuery( inputFile , ent , userVariants ):
	var = MAFVariant()
	for var in userVariants:
		thisGroup = var.sample+var.genomicVar()
		ent.addQuery( var.gene , field="gene" , group=thisGroup )
		ent.addQuery( var.chromosome , field="chr" , group=thisGroup )
		ent.addQuery( var.start + ":" + var.stop , field="chrpos37" , group=thisGroup )
		ent.addQuery( "human" , field="orgn" , group=thisGroup )
		#ent.addQuery( var.variantClass , "vartype" )
		#ent.addQuery( var.referencePeptide + var.positionPeptide + var.alternatePeptide , "Variant name" )
		#var.referencePeptide , var.positionPeptide , var.alternatePeptide
	return ent

### main ### 
def main( argv ):
	values = parseArgs( argv )
	inputFile = values["input"]
	outputFormat = values["output"]
	doClinVar = values["clinvar"]
	doExAC = values["exac"]

	userVariants = splitByVariantType( inputFile )
	
	calls = autovivification.autovivification({})
	if doClinVar:
		ent = entrezAPI()	
		ent = prepQuery( inputFile , ent , userVariants )
		ent.database = entrezAPI.clinvar
		clinvarEntries = ent.doBatch( 5 )
		clinvarVariants = clinvarEntries["variants"]
		clinvarTraits = clinvarEntries["traits"]
		clinvarClinical = clinvarEntries["clinical"]

		calls["PVS1"] = PVS1( userVariants , None, None , None )
		calls["PS1"] = PS1( userVariants , clinvarVariants , clinvarClinical )
		calls["PM4"] = PM4( userVariants )
		calls["PM5"] = PM5( userVariants , clinvarVariants , clinvarClinical )

	for module in calls:
		print module
		for uniVar in calls[module]:
			print uniVar ,
			if calls[module][uniVar]:
				print "\tis " + module
			else:
				print "\tis NOT " + module

	if doExAC:
		exac = exacAPI(harvard=True)
		exacEntries = exac.getAlleleFrequencies( userVariants )
		thresh = 0
		for genVar in exacEntries:
			alleleFrequency = exacEntries[genVar]
			if isFrequentAllele( alleleFrequency , thresh ):
				print genVar + " is NOT rare(" + str(thresh) + "): " + str(alleleFrequency)

if __name__ == "__main__":
	main( sys.argv[1:] )
