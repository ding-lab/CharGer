#!/usr/bin/python
# CharGer - Characterization of Germline variants
# author: Adam D Scott (ascott@genome.wustl.edu) & Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2015*12

import sys
import getopt
import charger

def parseArgs( argv ):
	helpText = "python main.py" + " "
	helpText += "-m \"maf\" "
	helpText += "(-l suppress ClinVar, "
	helpText += "-x suppress ExAC, "
	helpText += "-e expression, "
	helpText += "-g gene list, "
	helpText += "-n de novo, "
	helpText += "-a assumed de novo, "
	helpText += "-c co-segregation, "
	helpText += "-o \"output\")\n"
	mafFile = ""
	expressionFile = ""
	geneListFile = ""
	deNovoFile = ""
	assumedDeNovoFile = ""
	coSegregationFile = ""
	output = ""
	clinvar = True
	exac = True
	try:
		opts, args = getopt.getopt( argv , "lxh:m:o:g:e:n:a:c:" , \
		["maf=" , "output=" , "geneList=" , "expression=" , "deNovo=" , "assumedDeNovo=" , "coSegregation="] )
	except getopt.GetoptError:
		print "CharGer ERROR: Command not recognized"
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
		elif opt in ( "-m" , "--maf" ):
			mafFile = arg
		elif opt in ( "-o" , "--output" ):
			output = arg
		elif opt in ( "-g" , "--geneList" ):
			geneListFile = arg
		elif opt in ( "-e" , "--expression" ):
			expressionFile = arg
		elif opt in ( "-n" , "--deNovo" ):
			deNovoFile = arg
		elif opt in ( "-a" , "--assumedDeNovo" ):
			assumedDeNovoFile = arg
		elif opt in ( "-c" , "--coSegregation" ):
			coSegregationFile = arg
		elif opt in ( "-l" , "--help" ):
			clinvar = False
		elif opt in ( "-x" , "--help" ):
			exac = False
	return { "maf" : mafFile , \
	"output" : output , \
	"clinvar" : clinvar , \
	"exac" : exac , \
	"expression" : expressionFile , \
	"deNovo" : deNovoFile , \
	"assumedDeNovo" : assumedDeNovoFile , \
	"coSegregation" : coSegregationFile , \
	"geneList" : geneListFile }

### main ### 
def main( argv ):
	values = parseArgs( argv )
	mafFile = values["maf"]
	expressionFile = values["expression"]
	deNovoFile = values["deNovo"]
	assumedDeNovoFile = values["assumedDeNovo"]
	coSegregationFile = values["coSegregation"]
	geneListFile = values["geneList"]
	outputFormat = values["output"]
	doClinVar = values["clinvar"]
	doExAC = values["exac"]

	CharGer = charger.charger()
	CharGer.getInputData( maf=mafFile , \
	geneList=geneListFile , \
	expression=geneListFile , \
	deNovo=deNovoFile , \
	assumedDeNovo=assumedDeNovoFile , \
	coSegregation=coSegregationFile )

	CharGer.getExternalData( clinvar=doClinVar , exac=doExAC , batchSize=200 )

	CharGer.PVS1( )
	CharGer.PS1( )
	CharGer.PS2( )
	CharGer.PS3( )
	CharGer.PS4( )
	CharGer.PM1( )
	CharGer.PM2( threshold=0.05 )
	CharGer.PM3( )
	CharGer.PM4( )
	CharGer.PM5( )
	CharGer.PM6( )
	CharGer.PP1( )
	CharGer.PP2( )
	CharGer.PP3( )
	CharGer.PP4( )
	CharGer.PP5( )
	#CharGer.printResult( )
	CharGer.classify()
	CharGer.printClassifications( )

if __name__ == "__main__":
	main( sys.argv[1:] )
