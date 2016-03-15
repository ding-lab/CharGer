#!/usr/bin/python
# CharGer - Characterization of Germline variants
# author: Adam D Scott (ascott@genome.wustl.edu) & Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2015*12

import time
import math
import re
from WebAPI.Ensembl.ensemblAPI import ensemblAPI
import glob
from scipy import stats
from WebAPI.Entrez.entrezAPI import entrezAPI
from WebAPI.ExAC.exacAPI import exacAPI
from WebAPI.Variant.clinvarVariant import clinvarVariant
from WebAPI.Variant.vepVariant import vepVariant
from WebAPI.Variant.MAFVariant import MAFVariant
from WebAPI.Variant.variant import variant
from chargerVariant import chargerVariant
from autovivification import autovivification as AV
import vcf

class charger(object):
	''' Example usage:
			CharGer = charger()
			CharGer.getInput( maf=mafFile , expression=expressionFile , geneList=geneListFile )
			CharGer.getWebData( clinvar=doClinVar , exac=doExAC )
			CharGer.characterize( )
			CharGer.PVS1( )
			CharGer.PS1( )
			CharGer.PM4( )
			CharGer.PM5( )
	'''
	allDiseases = "all"
	def __init__( self , **kwargs ):
		self.userVariants = kwargs.get( 'variants' , [] )
		self.userExpression = kwargs.get( 'expressions' , AV({}) )
		self.userGeneList = kwargs.get( 'geneList' , AV({}) )
		self.userDeNovoVariants = kwargs.get( 'deNovo' , {} )
		self.userAssumedDeNovoVariants = kwargs.get( 'assumedDeNovo' , {} )
		self.userCoSegregateVariants = kwargs.get( 'coSegregate' , {} )
		#self.clinvarVariants = kwargs.get( 'clinvarVariants' , {} )
		#self.vepVariants = kwargs.get( 'vepVariants' , [] )
		self.diseases = kwargs.get( 'diseases' , {} )

### Retrieve input data from user ###
	def getInputData( self  , **kwargs ):
		mafFile = kwargs.get( 'maf' , "" )
		vcfFile = kwargs.get( 'vcf' , "" )
		tsvFile = kwargs.get( 'tsv' , "" )
		expressionFile = kwargs.get( 'expression' , "" )
		geneListFile = kwargs.get( 'geneList' , "" )
		deNovoFile = kwargs.get( 'deNovo' , "" )
		assumedDeNovoFile = kwargs.get( 'assumedDeNovo' , "" )
		coSegregateFile = kwargs.get( 'coSegregate' , "" )
		tcga = kwargs.get( 'tcga' , True )
		diseaseFile = kwargs.get( 'diseases' , "" )
		specific = kwargs.get( 'specific' , False )
		self.getDiseases( diseaseFile , **kwargs )
		if mafFile:
			self.readMAF( mafFile , **kwargs )
		if vcfFile:
			self.readVCF( vcfFile , **kwargs )
		if tsvFile:
			self.readTSV( tsvFile , **kwargs )
		if expressionFile:
			self.readExpression( expressionFile )
		else: 
			print "No expression file uploaded. CharGer will allow all passed truncations without expression data in PVS1."
		if geneListFile:
			self.readGeneList( geneListFile , specific=specific )
		else:
			print "No gene list file uploaded. CharGer will not make PVS1 calls."
		for var in self.userVariants:
			if str(var.reference) == "0" or not var.reference:
				var.reference = "-"
			if str(var.alternate) == "0" or not var.alternate:
				var.alternate = "-"
	def readMAF( self , inputFile , **kwargs ):
		inFile = self.safeOpen( inputFile , 'r' )
		codonColumn = kwargs.get( 'codon' , 48 )
		peptideChangeColumn = kwargs.get( 'peptideChange' , 49 )
		tcga = kwargs.get( 'tcga' , True )
		specific = kwargs.get( 'specific' , True )
		try:
			next(inFile)
			for line in inFile:
				var = chargerVariant()
				var.mafLine2Variant( line , peptideChange=peptideChangeColumn , codon=codonColumn )
				if specific:
					if tcga:
						match = re.match( "TCGA\-(\w\w)" , var.sample )
						if match:
							var.disease = self.diseases[match.groups()[0]]
				else:
					var.disease = charger.allDiseases
				self.userVariants.append( var )
		except:
			raise Exception( "CharGer Error: bad .maf file" )
	def readVCF( self , inputFile , **kwargs ):
		inFile = vcf.Reader( open( inputFile , 'r' ) )
		for record in inFile:
			reference = record.REF
			alternates = record.ALT
			start = record.POS #1-base beginning of ref
			stop = record.end+1 #0-base ending of ref
			for alternate in alternates:
				if alternate == "None":
					alternate = None
				var = chargerVariant( \
					chromosome = record.CHROM , \
					start = start , \
					stop = stop , \
					dbsnp = record.ID , \
					reference = reference , \
					alternate = str(alternate) , \
				)
				#quality = record.QUAL
				#vcfFilter = record.FILTER
				#info = record.INFO
				#if ( vcfFilter == "PASS" or quality >= 10 ) or \
				#	( not vcfFilter and not quality ):
				#	self.userVariants.append( var )
				self.userVariants.append( var )

	def readTSV( self , inputFile , **kwargs ):
		print "\tReading .tsv!"
		inFile = self.safeOpen( inputFile , 'r' )
		chrColumn = kwargs.get( 'chr' , 0 )
		startColumn = kwargs.get( 'start' , 1 )
		stopColumn = kwargs.get( 'stop' , 2 )
		refColumn = kwargs.get( 'ref' , 3 )
		altColumn = kwargs.get( 'alt' , 4 )
		geneColumn = kwargs.get( 'gene' , 6 )
		strandColumn = kwargs.get( 'strand' , 11 )
		peptideColumn = kwargs.get( 'peptideChange' , 14 )
		codonColumn = kwargs.get( 'codon' , 15 )
		sampleColumn = kwargs.get( 'sample' , 21 )
		tcga = kwargs.get( 'tcga' , True )
		specific = kwargs.get( 'specific' , True )
		headerLine = next(inFile).split( "\t" ) 
		try:
			for line in inFile:
				fields = line.split( "\t" )
				chrom = fields[int(chrColumn)]
				alt = fields[int(altColumn)]
				ref = fields[int(refColumn)]
				start = fields[int(startColumn)]
				stop = fields[int(stopColumn)]
				#strand = fields[int(strandColumn)]
				#gene = fields[int(geneColumn)]
				sample = fields[int(sampleColumn)]
				
				var = chargerVariant( \
					chromosome = chrom , \
					start = start , \
					stop = stop , \
					reference = ref , \
					alternate = alt , \
					#strand = strand , \
					#gene = gene , \
					sample = sample , \
					#codon = codonChange , \
					#peptide = peptideChange , \
				)
				var.splitHGVSc( fields[int(codonColumn)] )
				var.splitHGVSp( fields[int(peptideColumn)] )

				if specific:
					if tcga:
						match = re.match( "TCGA\-(\w\w)" , var.sample )
						if match:
							var.disease = self.diseases[match.groups()[0]]
				else:
					var.disease = charger.allDiseases
				self.userVariants.append( var )
		except:
			raise Exception( "CharGer Error: bad .tsv file" )
			print "Hint: correct columns?"
			print headerLine
			
	def readExpression( self , inputFile ): # expect sample(col)-gene(row) matrixes
		try:
			fileNames = glob.glob(inputFile + "*")
			for fileName in fileNames:
				print "Processing expression from ", fileName
				inFile = self.safeOpen( fileName , 'r' , warning=True )
				header = inFile.readline()
				samples = header.split( "\t" )

				# make sure the sample names are matching up, assuming
				# sample name in RSEM: TCGA-OR-A5J1-01A-11R-A29S-07
				# sample name in MAF: TCGA-H4-A2HQ
				for idx, sample in enumerate(samples):
					sample_s = sample.split( "-" )
					sample_only = "-".join(sample_s[0:3])
					samples[idx] = sample_only
					
				for line in inFile:
					fields = line.split( "\t" )		
					gene = fields[0] 
					gene_exp = [ self.toFloat(x) for x in fields[1:]] # manage potential NA strings that may cause errors
					gene_exp_p = (stats.rankdata(gene_exp, 'min')-1)/len(gene_exp) # convert to percentile
					for i in range(1,len(gene_exp_p)):
						self.userExpression[samples[i+1]][gene] = gene_exp_p[i]

		except:
			print "CharGer Error: bad expression file"

	def readGeneList( self , inputFile , **kwargs ): # gene list formatted "gene", "disease", "mode of inheritance"
		specific = kwargs.get( 'specific', True )
		try:
			inFile = self.safeOpen( inputFile , 'r' , warning=True )
			for line in inFile:
				fields = line.split( "\t" )
				gene = fields[0]
				if specific:
					disease = fields[1]
				else: #set the gene to match all disease
					disease = charger.allDiseases
				mode_inheritance = fields[2]
				self.userGeneList[gene][disease] = mode_inheritance
		except:
			print "CharGer Error: bad gene list file"
	def readDeNovo( self , inputFile ):
		self.readOtherMAF( inputFile , varDict = self.deNovoVariants )
	def readCoSegregate( self , inputFile ):
		self.readOtherMAF( inputFile , varDict = self.coSegregateVariants )
	def readAssumedDeNovo( self , inputFile ):
		self.readOtherMAF( inputFile , varDict = self.assumedDeNovoVariants )
	def readOtherMAF( self , inputFile, varDict ):
		try:
			inFile = self.safeOpen( inputFile , 'r' , warning=True )
			next(inFile)
			for line in inFile:
				fields = line.split( "\t" )
				var = MAFVariant()
				var.mafLine2Variant( line )
				varDict[var.uniqueVar()] = 1
		except:
			print "CharGer Warning: bad .maf for " + inputFile

### Retrieve external reference data ###
	def getExternalData( self , **kwargs ):
		t = time.time()
		self.getClinVar( **kwargs )
		self.printRunTime( "ClinVar" , self.runTime( t ) )
		t = time.time()
		self.getExAC( **kwargs )
		self.printRunTime( "ExAC" , self.runTime( t ) )
		t = time.time()
		self.getVEP( **kwargs )
		self.printRunTime( "VEP" , self.runTime( t ) )
		self.fillMissingVariantInfo()
	def fillMissingVariantInfo( self ):
		for var in self.userVariants:
			var.fillMissingInfo()
	def getClinVar( self , **kwargs ):
		doClinVar = kwargs.get( 'clinvar' , True )
		summaryBatchSize = kwargs.get( 'summaryBatchSize' , 500 )
		searchBatchSize = kwargs.get( 'searchBatchSize' , 50 )
		if doClinVar:
			ent = entrezAPI()
			i = 0
			for varsStart in range( 0 , len( self.userVariants ) , int(searchBatchSize) ):
				varsEnd = varsStart + int(searchBatchSize)
				varsSet = self.userVariants[varsStart:varsEnd]
				ent.prepQuery( varsSet )
				ent.subset = entrezAPI.esearch
				ent.database = entrezAPI.clinvar
				clinvarsSet = ent.doBatch( summaryBatchSize )
				varsBoth = self.matchClinVar( varsSet , clinvarsSet )
				self.userVariants[varsStart:varsEnd] = varsBoth["userVariants"]
				#self.clinvarVariants.update( varsBoth["clinvarVariants"] )
				#self.userVariants[varsStart:varsEnd] = self.matchClinVar( varsSet , clinvarsSet )
	def getExAC( self , **kwargs ):
		doExAC = kwargs.get( 'exac' , True )
		useHarvard = kwargs.get( 'harvard' , True )
		threshold = kwargs.get( 'threshold' , 0 )
		if doExAC:
			common = 0
			rare = 0
			totalVars = len( self.userVariants )
			exac = exacAPI(harvard=useHarvard)
#entries by genomivVar
			exacIn = self.getUniqueGenomicVariantList( self.userVariants )
			#entries = exac.getAlleleFrequencies( self.userVariants )
			entries = exac.getAlleleFrequencies( exacIn )
			for var in self.userVariants:
				if var.genomicVar() in entries:
					alleleFrequency = entries[var.genomicVar()]
				else:
					alleleFrequency = None
				var.alleleFrequency = alleleFrequency
				if var.isFrequentAllele( threshold ):
					common += 1
				else:
					rare += 1
			elen = len(entries.keys())
			print "ExAC found " + str(common) + "common & " + str(rare) + "rare variants out of " + str(totalVars) + "total variants and " + str(elen) + "unique variants"
	def getVEP( self , **kwargs ):
		doVEP = kwargs.get( 'vep' , True )
		if doVEP:
			vep = ensemblAPI()
			luv = len(self.userVariants)
			vepVariants = vep.annotateVariantsPost( self.userVariants )
			self.matchVEP( vepVariants )
			aluv = 0
			for var in self.userVariants:
				print var.proteogenomicVar()
				if var.vepVariant:
					aluv += 1
			luvafter = len(self.userVariants)
			print "\nVEP annotated userVariants " + str( luvafter )
			print "VEP annotated " + str(aluv) + " from the original set of " + str(luv)

#### Helper methods for data retrieval ####
	def matchClinVar( self , userVariants , clinvarVariants ):
		for var in userVariants:
			for uid in clinvarVariants:
				cvar = clinvarVariants[uid]
				if var.sameGenomicVariant( cvar ):
					var.clinical = cvar.clinical
					var.clinvarVariant = cvar
		return { "userVariants" : userVariants , "clinvarVariants" : clinvarVariants }
	def matchVEP( self , vepVariants ):
		for var in self.userVariants:
			genVar = var.vcf()
			vepVar = vepVariant()
			if genVar in vepVariants:
				vepVar = vepVariants[genVar]
			print ""
			print var.proteogenomicVar() + "\t" ,
			if var.sameGenomicVariant( vepVar ):
				var.vepVariant = vepVar
				var.copyMostSevereConsequence()
			if var.vepVariant:
				print "VEP+" ,
				if var.vepVariant.consequences:
					print "CONS\t" ,
			print ""
	def getDiseases( self , diseasesFile , **kwargs ):
		tcga = kwargs.get( 'tcga' , True )
		try:
			if tcga:
				conversion = open( diseasesFile , "r" )
				for line in conversion:
					fields = line.split('\t')
					self.diseases[fields[0]] = fields[2].strip()
				return self.diseases
			return
		except:
			print "CharGer Warning: No diseases file provided"
			return


### Evidence levels ### 
##### Very Strong #####
	def PVS1( self , expressionThreshold = 0.2 ):
		print "CharGer module PVS1"
		print "- truncations in genes where LOF is a known mechanism of the disease"
		print "- require the mode of inheritance to be dominant (assuming heterzygosity) and co-occurence with reduced gene expression"
		maf_truncations = ["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Splice_Site"] #,"Nonstop_Mutation"
		vep_truncations = ["transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained",\
							"frameshift_variant","start_lost"]
		if self.userGeneList: #gene, disease, mode of inheritance
			for var in self.userVariants:
				varGene = var.gene
				varDisease = var.disease # no disease field in MAF; may require user input	
				varSample = var.sample
				varClass = var.variantClass
				varVEPClass = ""
				if var.vepVariant:
					varVEPClass = var.vepVariant.mostSevereConsequence
				altPeptide = var.alternatePeptide
				if (varClass in maf_truncations) or \
					(varVEPClass in vep_truncations) or \
					altPeptide == "*" or \
					altPeptide == "fs":
					if varGene in self.userGeneList: # check if in gene list
						if ( "dominant" in self.userGeneList[varGene][varDisease] or \
							"dominant" in self.userGeneList[varGene][charger.allDiseases]):
							var.PVS1 = True # if call is true then check expression effect
							if self.userExpression: # consider expression data only if the user has supplied an expression matrix
								if self.userExpression[varSample][varGene] >= expressionThreshold:
									var.PVS1 = False 
		else: 
			print "CharGer Error: Cannot evaluate PVS1: No gene list supplied."

##### Strong #####
	def PS1( self ):
		print "CharGer module PS1"
		print "- same peptide change as a previously established pathogenic variant"
		self.peptideChange( "PS1" )
	def PS2( self ):
		print "CharGer module PS2"
		print "- de novo with maternity and paternity confirmation and no family history"
		for var in self.userVariants:
			if var.uniqueVar() in self.userDeNovoVariants:
				var.PS2 = True
	def PS3( self ):
		print "CharGer module PS3: not yet implemented"
#		print "- "
	def PS4( self ): # not relevant in rare variants, such big effect size common variants are so rare may as well just take a input variant list
		print "CharGer module PS4: not yet implemented"
#		print "- variant prevalence in cases significantly greater than controls"
		for var in self.userVariants:
			return
			caseVarFreq = "NEED UPDATE" # may take input from current MAF
			controlVarFreq = "NEED UPDATE" # may take input from ExAC
			if ( caseVarFreq != 0 and controlVarFreq != 0):
				OR = (caseVarFreq/controlVarFreq) / ( (1-caseVarFreq)/(1-controlVarFreq) )
				# Adam will update
				if OR >= 5:
					CIlower = math.log(OR) - math.sqrt( 1/caseVarFreq + 1/controlVarFreq + 1/caseVarFreq + 1/controlVarFreq)
					if (CIlower > 1):
						var.PS4 = True

##### Moderate #####
	def PM1( self ):
		print "CharGer module PM1: not yet implemented"
#		print "- "
	def PM2( self , threshold ):
		print "CharGer module PM2"
		print "- absent or extremely low frequency in controls"
		self.checkAlleleFrequencies( "PM2" , threshold )
	def PM3( self ):
		print "CharGer module PM3: not yet implemented"
#		print "- "
	def PM4( self ):
		print "CharGer module PM4"
		print "- protein length changes due to inframe indels or nonstop variant"
		lenShift = ["In_Frame_Del","In_Frame_Ins","Nonstop_Mutation"]
		for var in self.userVariants:
			varClass = var.variantClass
			if varClass in lenShift:
				var.PM4 = True
	def PM5( self ):
		print "CharGer module PM5"
		print "- different peptide change of a pathogenic variant at the same reference peptide"
		self.peptideChange( "PM5" )
	def PM6( self ):
		print "CharGer module PM6"
		print "- assumed de novo without maternity and paternity confirmation"
		for var in self.userVariants:
			if var.uniqueVar() in self.userAssumedDeNovoVariants:
				var.PM6 = True

##### Supporing #####
	def PP1( self ):
		print "CharGer module PP1"
		print "- cosegregation with disease in family members in a known disease gene"
		for var in self.userVariants:
			if var.uniqueVar() in self.userCoSegregateVariants:
				var.PP1 = True
	def PP2( self ):
		print "CharGer module PP2: not yet implemented"
#		print "- "
	def PP3( self , minimumEvidence ):
		print "CharGer module PP3"
		print "- multiple lines of in silico evidence of deliterous effect"
		callSIFTdam = "damaging"
		callSIFTdel = "deleterious"
		callPolyphen = "probably damaging"
		callBlosum62 = -2
		callCompara = 2
		callImpact = "high"
		fracMaxEntScan = 0.8
		callGeneSplicer = ""
		for var in self.userVariants:
			if var.vepVariant:
				if var.vepVariant.consequences:
					for vcVar in var.vepVariant.consequences:
						if not var.PP3:
							evidence = 0
							if vcVar.blosum:
								if vcVar.blosum < callBlosum62:
									evidence += 1
							if vcVar.predictionSIFT:
								if vcVar.predictionSIFT.lower() == callSIFTdam or \
								vcVar.predictionSIFT.lower() == callSIFTdel:
									evidence += 1
							if vcVar.predictionPolyphen:
								if vcVar.predictionPolyphen.lower() == callPolyphen:
									evidence += 1
							if vcVar.compara:
								if vcVar.compara > callCompara:
									evidence += 1
							if vcVar.impact:
								if vcVar.impact.lower() == callImpact:
									evidence += 1
							if vcVar.maxentscan:
								callMaxEntScan = vcVar.maxentscan[0]*fracMaxEntScan
								if vcVar.maxentscan[1] <= callMaxEntScan:
									evidence += 1
							if vcVar.genesplicer:
								if vcVar.genesplicer.lower() == callGeneSplicer:
									evidence += 1
							if evidence >= minimumEvidence:
								var.PP3 = True

	def PP4( self ):
		print "CharGer module PP4: not yet implemented"
#		print "- "
	def PP5( self ):
		print "CharGer module PP5: not yet implemented"
#		print "- "

### helper functions of evidence levels ###
	def peptideChange( self , mod ):
		called = 0
		for var in self.userVariants:
			uniVar = var.uniqueVar()
			ps1Call = False
			pm5Call = False
			call = False
			if mod == "PS1":
				call = var.PS1
			if mod == "PM5":
				call = var.PM5
			if not call: #is already true
				#for genVar in self.clinvarVariants:
					#clinvarVar = self.clinvarVariants[genVar]
				if var.clinvarVariant:
					clinvarVar = var.clinvarVariant
					clin = clinvarVar.clinical
					if var.vepVariant:
						if var.vepVariant.consequences:
							for consequence in var.vepVariant.consequences:
								if consequence.sameGenomicVariant( clinvarVar ):
								#if genomic change is the same, then PS1
									if clin["description"] == clinvarVariant.pathogenic:
										if mod == "PS1":
											var.PS1 = True # already pathogenic still suffices to be PS1
											called += 1
								elif consequence.sameGenomicReference( clinvarVar ):
								#if genomic change is different, but the peptide change is the same, then PS1
									if clinvarVar.alternatePeptide == consequence.alternatePeptide: #same amino acid change
										if clin["description"] == clinvarVariant.pathogenic:
											if mod == "PS1":
												var.PS1 = True
												called += 1
								if consequence.samePeptideReference( clinvarVar ):
									if not consequence.samePeptideChange( clinvarVar ):
									#if peptide change is different, but the peptide reference is the same, then PM5
										if consequence.plausibleCodonFrame( clinvarVar ):
											if clin["description"] == clinvarVariant.pathogenic:
												if mod == "PM5":
													var.PM5 = True # already pathogenic still suffices to be PS1
													called += 1
		print mod + " found " + str(called) + " pathogenic variants"
	def printResult( self ):
		for var in self.userVariants:
			for module in var.modules():
				print var.uniqueVar() ,
				if var.check( module ):
					print "\tis " + module
				else:
					print "\tis NOT " + module
	def checkAlleleFrequencies( self , mod , threshold ):
		for var in self.userVariants:
			if mod == "PM2":
				if not var.isFrequentAllele( threshold ):
					var.PM2 = True
			if mod == "BA1":
				if var.isFrequentAllele( threshold ):
					var.BA1 = True

### Benign Modules ###
#### Stand-alone ####
	def BA1( self , threshold ):
		print "CharGer module BA1"
		print "- allele frequency >5%"
		self.checkAlleleFrequencies( "BA1" , threshold )
#### Strong ####
	def BS1( self ):
		print "CharGer module BS1: not yet implemented"
	def BS2( self ):
		print "CharGer module BS2: not yet implemented"
	def BS3( self ):
		print "CharGer module BS3: not yet implemented"
		#print " - in vitro or in vivo functional studies with no damaging effect on protein function or splicing"
	def BS4( self ):
		print "CharGer module BS4: not yet implemented"
#### Supporting ####
	def BP1( self ):
		print "CharGer module BP1: not yet implemented"
	def BP2( self ):
		print "CharGer module BP2: not yet implemented"
	def BP3( self ):
		print "CharGer module BP3: not yet implemented"
	def BP4( self , minimumEvidence ):
		print "CharGer module BP4"
		print " - in silico evidence of no damage"
		callSIFTdam = "damaging"
		callSIFTdel = "deleterious"
		callPolyphen = "probably damaging"
		callBlosum62 = -2
		callCompara = 2
		callImpact = "high"
		fracMaxEntScan = 0.8
		callGeneSplicer = ""
		for var in self.userVariants:
			if var.vepVariant:
				if var.vepVariant.consequences:
					for vcVar in var.vepVariant.consequences:
						if not var.PP3:
							evidence = 0
							if vcVar.blosum:
								if vcVar.blosum > callBlosum62:
									evidence += 1
							if vcVar.predictionSIFT:
								if vcVar.predictionSIFT.lower() != callSIFTdam and \
								vcVar.predictionSIFT.lower() != callSIFTdel:
									evidence += 1
							if vcVar.predictionPolyphen:
								if vcVar.predictionPolyphen.lower() != callPolyphen:
									evidence += 1
							if vcVar.compara:
								if vcVar.compara <= callCompara:
									evidence += 1
							if vcVar.impact:
								if vcVar.impact.lower() != callImpact:
									evidence += 1
							if vcVar.maxentscan:
								callMaxEntScan = vcVar.maxentscan[0]*fracMaxEntScan
								if vcVar.maxentscan[1] > callMaxEntScan:
									evidence += 1
							if vcVar.genesplicer:
								if vcVar.genesplicer.lower() != callGeneSplicer:
									evidence += 1
							if evidence >= minimumEvidence:
								var.BP4 = True

	def BP5( self ):
		print "CharGer module BP5: not yet implemented"
		#print" - multiple lines of computational evidence suggesting no impact on gene or product"
	def BP6( self ):
		print "CharGer module BP6: not yet implemented"
	def BP7( self ):
		print "CharGer module BP7: not yet implemented"

### Classifier ###
	def classify( self ):
		for var in self.userVariants:
			var.isPathogenic( )
			var.isLikelyPathogenic( )
			var.isLikelyBenign( )
			var.isBenign( )
			var.isUncertainSignificance( )
	def printClassifications( self ):
		headLine = '\t'.join( ["Variant" , "PositiveEvidence" , \
			"CharGerClassification" , "ClinVarAnnoation"] )
		print headLine
		i = 0
		for var in self.userVariants:
			i += 1
			print '\t'.join( [ str(i) , var.uniqueProteogenomicVar() , \
				var.positiveEvidence() , var.pathogenicity , \
				var.clinical["description"] ] )
	def writeSummary( self , outFile , **kwargs ):
		delim = kwargs.get( 'delim' , '\t' )
		outFH = self.safeOpen( outFile , 'w' , warning=True )
		headLine = delim.join( ["HUGO_Symbol" , "Chromosome" , "Start" , \
			"Stop" , "Reference" , "Alternate" , "Strand" , "Assembly" , \
			"Variant_Type" , "Variant_Classification" , \
			"Sample" , "Transcript" , "Codon_Position" , "Protein" , \
			"Peptide_Reference" , "Peptide_Position" , "Peptide_Alternate" , \
			"VEP_Most_Severe_Consequence" , "ClinVar_Pathogenicity" , \
			"PositiveEvidence" , "NegativeEvidence" , "CharGerClassification"] )
		try:
			outFH.write( headLine )
			outFH.write( "\n" )
			for var in self.userVariants:
				fields = []

				self.appendStr( fields,var.gene)
				self.appendStr( fields,var.chromosome)
				self.appendStr( fields,var.start)
				self.appendStr( fields,var.stop)
				self.appendStr( fields,var.reference)
				self.appendStr( fields,var.alternate)
				self.appendStr( fields,var.strand)
				self.appendStr( fields,var.assembly)
				self.appendStr( fields,var.variantType)
				self.appendStr( fields,var.variantClass)
				self.appendStr( fields,var.sample)
				self.appendStr( fields,var.transcriptCodon)
				self.appendStr( fields,var.positionCodon)
				self.appendStr( fields,var.transcriptPeptide)
				self.appendStr( fields,var.referencePeptide)
				self.appendStr( fields,var.positionPeptide)
				self.appendStr( fields,var.alternatePeptide)
				#self.appendStr( fields,var.vepVariant.mostSevereConsequence) ## this line will fail you on insertions regardless of all the checks in appendStr
				#self.appendStr( fields,var.clinvarVariant.trait)
				self.appendStr( fields,var.clinical["description"])
				self.appendStr( fields,var.positiveEvidence())
				self.appendStr( fields,var.negativeEvidence())
				self.appendStr( fields,var.pathogenicity)

				outFH.write( delim.join( fields ) )
				outFH.write( "\n" )
		except:
			print "CharGer Error: Could not write output summary"
			pass

	@staticmethod
	def appendStr( array, value ):
		try:
			array.append( str( value ) )
		except:
			print "failed to append a value\n"
			array.append( str( "NA" ) )
			pass

	@staticmethod
	def safeOpen( inputFile , rw , **kwargs ):
		errwar = kwargs.get( 'warning' , False )
		try:
			return open( inputFile , rw )
		except:
			if errwar:
				print "CharGer Warning: could not open " + inputFile
			else:
				print "CharGer Error: could not open " + inputFile
			pass

	@staticmethod
	def getUniqueGenomicVariantList( aVarList ):
		uniqueVarList = []
		uniqueVarDict = AV({})
		for var in aVarList:
			generalVar = variant()
			generalVar.copyInfo( var )
			if generalVar.genomicVar() not in uniqueVarDict:
				uniqueVarDict[generalVar.genomicVar()] = 1
				uniqueVarList.append( generalVar )
		return uniqueVarList
	@staticmethod
	def toFloat(x):
		try:
			float(x)
			return float(x)
		except:
			return "nan"
	@staticmethod
	def runTime( initialTime ):
	   return time.time() - initialTime
	@staticmethod
	def printRunTime( step , interval ):
		print "Running " + step + " took " + str( interval ) + "seconds"
	@staticmethod
	def printVariants( variants ):
		for var in variants:
			print var.proteogenomicVar()
