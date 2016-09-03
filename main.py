#!/usr/bin/python
# CharGer - Characterization of Germline variants
# author: Adam D Scott (adamscott@wustl.edu) & Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2015*12

import sys
import getopt
import charger
import time

def parseArgs( argv ):
	helpText = "python main.py" + " "
	helpText += "Accepted input data files:\n"
	helpText += "  -m Standard .maf\n"
	helpText += "  -f Standard .vcf\n"
	helpText += "  -T Custom .tsv\n\n"
	helpText += "Output:\n"
	helpText += "  -o output file\n"
	helpText += "  -w output as HTML (flag)\n"
	helpText += "Suppress accession of data:\n"
	helpText += "  -l suppress ClinVar (flag)\n"
	helpText += "  -x suppress ExAC (flag)\n"
	helpText += "  -E suppress VEP (flag)\n"
	helpText += "  -t suppress TCGA cancer types (flag)\n"
	helpText += "  -O override with ClinVar description (flag)\n"
	helpText += "  -D disease specific (flag)\n"
	helpText += "Cross-reference data files:\n"
	helpText += "  -z pathogenic variants .vcf\n"
	helpText += "  -e expression matrix file .tsv\n"
	helpText += "  -g gene list file .txt\n"
	helpText += "  -d diseases file (format: gene\\tdisease) .tsv\n"
	helpText += "  -n de novo file .?\n"
	helpText += "  -a assumed de novo file .?\n"
	helpText += "  -c co-segregation file .?\n"
	helpText += "  -H HotSpot3D clusters file .clusters\n"
	helpText += "ReST batch sizes:\n"
	helpText += "  -v VEP\n"
	helpText += "  -b ClinVar summary\n"
	helpText += "  -B ClinVar searchsize\n"
	helpText += "  -r recurrence threshold\n"
	helpText += "Custom columns (0-based)\n"
	helpText += "  -G HUGO gene symbol\n"
	helpText += "  -X chromosome\n"
	helpText += "  -S start position\n"
	helpText += "  -P stop position\n"
	helpText += "  -R reference allele\n"
	helpText += "  -A alternate allele\n"
	helpText += "  -s strand\n"
	helpText += "  -M sample name\n"
	helpText += "  -C codon\n"
	helpText += "  -p peptide change\n"
	helpText += "  -F allele frequency\n"
	helpText += "\n"
	helpText += "  -h this message\n"

	mafFile = ""
	vcfFile = ""
	tsvFile = ""
	expressionFile = ""
	geneListFile = ""
	deNovoFile = ""
	assumedDeNovoFile = ""
	coSegregationFile = ""
	diseasesFile = ""
	output = "charger_summary.tsv"
	clinvarSummaryBatchSize = 100
	clinvarSearchBatchSize = 100
	vepBatchSize = 400
	chrColumn = 0
	startColumn = 1
	stopColumn = 2
	refColumn = 3
	altColumn = 4
	geneColumn = 6
	strandColumn = 11
	codonColumn = 14
	peptideChangeColumn = 15
	sampleColumn = 21
	alleleFrequencyColumn= 33
	specific = True
	tcga = True
	clinvar = True
	exac = True
	vep = True
	asHTML = False
	override = False
	recurrenceThreshold = 2
	clustersFile = ""
	pathogenicVariantsFile = ""
	try:
		#opts, args = getopt.getopt( argv , "DEtlxhwOX:s:A:R:S:P:M:G:m:f:T:o:v:b:B:p:C:F:g:d:e:n:a:c:r:H:z:" , \
		opts, args = getopt.getopt( argv , "DEtlxhwOX:s:A:R:S:P:M:G:m:f:T:o:v:b:B:p:C:F:g:d:e:n:a:c:r:H:z:" , \
		["maf=" , "vcf=" , "tsv=" , "output=" , "vepBatchSize=" , "summaryBatchSize=" , "searchBatchSize=" , \
		"peptideChange=" , "codon=" , "alleleFrequency=" , "geneList=" , "diseases=" , \
		"expression=" , "deNovo=" , "assumedDeNovo=" , "coSegregation=" , \
		"recurrence=" , "hotspot3d=" , "pathogenicVariants=" ] )
	except getopt.GetoptError:
		print "CharGer ERROR: Command not recognized"
		print( helpText ) 
		sys.exit(2)
	if not opts:
		print "CharGer ERROR: Expected flagged input"
		print( helpText ) 
		sys.exit(2)
	for opt, arg in opts:
		if opt in ( "-h" , "--help" ):
			print( helpText )
			sys.exit()
		elif opt in ( "-m" , "--maf" ):
			mafFile = arg
			codonColumn = 48
			peptideChangeColumn = 49
		elif opt in ( "-f" , "--vcf" ):
			vcfFile = arg
		elif opt in ( "-T" , "--tsv" ):
			tsvFile = arg
		elif opt in ( "-z" , "--pathogenicVariants" ):
			pathogenicVariantsFile = arg
		elif opt in ( "-o" , "--output" ):
			output = arg
		elif opt in ( "-v" , "--vepBatchSize" ):
			vepBatchSize = int( arg )
		elif opt in ( "-b" , "--summaryBatchSize" ):
			clinvarSummaryBatchSize = int( arg )
		elif opt in ( "-B" , "--searchBatchSize" ):
			clinvarSearchBatchSize = int( arg )
		elif opt in ( "-p" , "--peptideChange" ):
			peptideChangeColumn = arg
		# all customized .tsv options are in caps for now
		elif opt in ( "-X" , "--chromosome" ):
			chrColumn = arg
		elif opt in ( "-s" , "--strand" ):
			strandColumn = arg
		elif opt in ( "-A" , "--alt" ):
			altColumn = arg
		elif opt in ( "-R" , "--ref" ):
			refColumn = arg
		elif opt in ( "-S" , "--start" ):
			output = arg
		elif opt in ( "-v" , "--vepBatchSize" ):
			vepBatchSize = int( arg )
		elif opt in ( "-b" , "--summaryBatchSize" ):
			clinvarSummaryBatchSize = int( arg )
		elif opt in ( "-B" , "--searchBatchSize" ):
			clinvarSearchBatchSize = int( arg )
		elif opt in ( "-p" , "--peptideChange" ):
			peptideChangeColumn = arg
		# all customized .tsv options are in caps for now
		elif opt in ( "-X" , "--chromosome" ):
			chrColumn = arg
		elif opt in ( "-s" , "--strand" ):
			strandColumn = arg
		elif opt in ( "-A" , "--alt" ):
			altColumn = arg
		elif opt in ( "-R" , "--ref" ):
			refColumn = arg
		elif opt in ( "-S" , "--start" ):
			startColumn = arg
		elif opt in ( "-P" , "--stop" ):
			stopColumn = arg
		elif opt in ( "-G" , "--gene" ):
			geneColumn = arg
		elif opt in ( "-M" , "--sample" ):
			sampleColumn = arg
		elif opt in ( "-C" , "--codon" ):
			codonColumn = arg
		elif opt in ( "-F" , "--alleleFrequency" ):
			alleleFrequencyColumn = arg
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
		elif opt in ( "-d" , "--diseases" ):
			diseasesFile = arg
		elif opt in ( "-D" , "--diseaseSpecific" ):
			specific = False
		elif opt in ( "-t" , "--notcga" ):
			tcga = False
		elif opt in ( "-l" , "--noclinvar" ):
			clinvar = False
		elif opt in ( "-E" , "--noVEP" ):
			vep = False
		elif opt in ( "-x" , "--noexac" ):
			exac = False
		elif opt in ( "-w" , "--html" ):
			asHTML = True
		elif opt in ( "-O" , "--override" ):
			override = True
		elif opt in ( "-r" , "--recurrenc" ):
			recurrenceThreshold = float( arg )
		elif opt in ( "-H" , "--hotspot3d" ):
			clustersFile = arg
	return { "maf" : mafFile , \
	"vcf" : vcfFile , \
	"tsv" : tsvFile , \
	"pathogenicVariantsFile" : pathogenicVariantsFile , \
	"output" : output , \
	"specific" : specific , \
	"tcga" : tcga , \
	"clinvar" : clinvar , \
	"vep" : vep , \
	"exac" : exac , \
	"html" : asHTML , \
	"override" : override , \
	"vepBatchSize" : vepBatchSize , \
	"clinvarSummaryBatchSize" : clinvarSummaryBatchSize , \
	"clinvarSearchBatchSize" : clinvarSearchBatchSize , \
	"peptideChangeColumn" : peptideChangeColumn , \
	"codonColumn" : codonColumn , \
	"expression" : expressionFile , \
	"deNovo" : deNovoFile , \
	"assumedDeNovo" : assumedDeNovoFile , \
	"coSegregation" : coSegregationFile , \
	"diseases" : diseasesFile , \
	"geneList" : geneListFile , \
	"chrColumn" : chrColumn, \
	"strandColumn" : strandColumn, \
	"startColumn" : startColumn, \
	"stopColumn" : stopColumn, \
	"refColumn" : refColumn, \
	"altColumn" : altColumn, \
	"geneColumn" : geneColumn, \
	"sampleColumn" : sampleColumn, \
	"alleleFrequencyColumn" : alleleFrequencyColumn, \
	"recurrenceThreshold" : recurrenceThreshold , \
	"clustersFile" : clustersFile , \
	}

### main ### 
def main( argv ):
	t0 = time.time()
	values = parseArgs( argv )
	mafFile = values["maf"]
	vcfFile = values["vcf"]
	tsvFile = values["tsv"]
	expressionFile = values["expression"]
	deNovoFile = values["deNovo"]
	assumedDeNovoFile = values["assumedDeNovo"]
	coSegregationFile = values["coSegregation"]
	geneListFile = values["geneList"]
	diseasesFile = values["diseases"]
	outputFile = values["output"]
	diseaseSpecific = values["specific"]
	doTCGA = values["tcga"]
	doClinVar = values["clinvar"]
	doExAC = values["exac"]
	doVEP = values["vep"]
	vepBatchSize = values["vepBatchSize"]
	clinvarSummaryBatchSize = values["clinvarSummaryBatchSize"]
	clinvarSearchBatchSize = values["clinvarSearchBatchSize"]
	peptideChangeColumn = values["peptideChangeColumn"]
	codonColumn = values["codonColumn"]
	chrColumn = values["chrColumn"]
	strandColumn = values["strandColumn"]
	startColumn = values["startColumn"]
	stopColumn = values["stopColumn"]
	refColumn = values["refColumn"]
	altColumn = values["altColumn"]
	sampleColumn = values["sampleColumn"]
	alleleFrequencyColumn = values["alleleFrequencyColumn"]
	asHTML = values["html"]
	override = values["override"]
	recurrenceThreshold = values["recurrenceThreshold"]
	clustersFile = values["clustersFile"]
	pathogenicVariantsFile = values["pathogenicVariantsFile"]
	
	t1 = time.time()

	CharGer = charger.charger()

	[ vepDone , preVEP , exacDone ] = CharGer.getInputData( maf=mafFile , \
	vcf=vcfFile , \
	tsv=tsvFile , \
	pathogenicVariants=pathogenicVariantsFile , \
	specific=diseaseSpecific , \
	tcga=doTCGA , \
	geneList=geneListFile , \
	expression=expressionFile , \
	deNovo=deNovoFile , \
	assumedDeNovo=assumedDeNovoFile , \
	coSegregation=coSegregationFile , \
	diseases=diseasesFile , \
	peptideChange=peptideChangeColumn , \
	codon=codonColumn , \
	chr=chrColumn , \
	strand=strandColumn , \
	start=startColumn , \
	stop=stopColumn , \
	ref=refColumn , \
	alt=altColumn , \
	sample=sampleColumn , \
	alleleFrequency=alleleFrequencyColumn \
	)

	if doVEP:
		if vepDone:
			doVEP = False
	if doExAC:
		if exacDone:
			doExac = False

	t2 = time.time() 

	CharGer.getExternalData( clinvar=doClinVar , \
	exac=doExAC , \
	vep=doVEP , \
	prevep=preVEP , \
	summaryBatchSize=clinvarSummaryBatchSize , \
	searchBatchSize=clinvarSearchBatchSize , \
	allOptions=False , \
	maxPost=vepBatchSize , \
	#timeout=(20,20) , \
	)

	t3 = time.time() 

	rareThreshold = 0.0005
	commonThreshold = 0.05
	minimumEvidence = 2

	CharGer.PVS1( )
	CharGer.PS1( )
	CharGer.PS2( )
	CharGer.PS3( )
	CharGer.PS4( )
	CharGer.PM1( recurrenceThreshold , hotspot3d=clustersFile )
	CharGer.PM2( rareThreshold )
	CharGer.PM3( )
	CharGer.PM4( )
	CharGer.PM5( )
	CharGer.PM6( )
	CharGer.PP1( )
	CharGer.PP2( )
	CharGer.PP3( minimumEvidence )
	CharGer.PP4( )
	CharGer.PP5( )

	CharGer.BA1( commonThreshold )
	CharGer.BS1( )
	CharGer.BS2( )
	CharGer.BS3( )
	CharGer.BS4( )
	CharGer.BP1( )
	CharGer.BP2( )
	CharGer.BP3( )
	CharGer.BP4( minimumEvidence )
	CharGer.BP5( )
	CharGer.BP6( )
	CharGer.BP7( )

	t4 = time.time() 

	CharGer.classify( system="ACMG" )
	CharGer.classify( system="CharGer" , override=override )

	t5 = time.time() 

	CharGer.printClassifications( )

	CharGer.writeSummary( outputFile , delim='\t' , html=asHTML )


	#CharGer.pdfSummary( outputFile )

	print "\nCharGer run Times:"
	dt1_0 = t1-t0
	print "input parse time (s): " + str(dt1_0)
	dt2_1 = t2-t1
	print "get input data time (s): " + str(dt2_1)
	dt3_2 = t3-t2
	print "get external data time (s): " + str(dt3_2)
	dt4_3 = t4-t3
	print "modules run time (s): " + str(dt4_3)
	dt5_4 = t5-t4
	print "classification time (s): " + str(dt5_4)
	dt5_0 = t5-t0
	print "CharGer full run time (s): " + str(dt5_0)

if __name__ == "__main__":
	main( sys.argv[1:] )
