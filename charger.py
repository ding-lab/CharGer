#!/usr/bin/python
# CharGer - Characterization of Germline variants
# author: Adam D Scott (ascott@genome.wustl.edu) & Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2015*12

import math
import re
from WebAPI.Ensembl.ensemblAPI import ensemblAPI
from WebAPI.Entrez.entrezAPI import entrezAPI
from WebAPI.ExAC.exacAPI import exacAPI
from WebAPI.Variant.clinvarVariant import clinvarVariant
from WebAPI.Variant.variant import variant
from chargerVariant import chargerVariant
from autovivification import autovivification as AV

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
		self.clinvarVariants = kwargs.get( 'clinvarVariants' , {} )
		self.vepVariants = kwargs.get( 'vepVariants' , [] )
		self.diseases = kwargs.get( 'diseases' , {} )

### Retrieve input data from user ###
	def getInputData( self  , **kwargs ):
		mafFile = kwargs.get( 'maf' , "" )
		expressionFile = kwargs.get( 'expression' , "" )
		geneListFile = kwargs.get( 'geneList' , "" )
		deNovoFile = kwargs.get( 'deNovo' , "" )
		assumedDeNovoFile = kwargs.get( 'assumedDeNovo' , "" )
		coSegregateFile = kwargs.get( 'coSegregate' , "" )
		tcga = kwargs.get( 'tcga' , True )
		diseaseFile = kwargs.get( 'diseases' , "" )
		specific = kwargs.get( 'specific' , False )
		self.getDiseases( diseaseFile , **kwargs )
		self.readMAF( mafFile , **kwargs )
		self.readExpression( expressionFile )
		self.readGeneList( geneListFile , specific=specific )
	def readMAF( self , inputFile , **kwargs ):
#		print "\tReading .maf!"
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
				print var.genomicVar() + "\t" + var.disease
				self.userVariants.append( var )
		except:
			raise Exception( "CharGer Error: bad .maf file" )
			
	def readExpression( self , inputFile ): # expect a sample(col)-gene(row) matrix
		try:
			inFile = self.safeOpen( inputFile , 'r' , warning=True )
			header = inFile.readline() # for future fetch header to get other field
			samples = header.split( "\t" )
			for line in inFile:
				fields = line.split( "\t" )
				gene = fields[0]
				for i in range(1,len(fields)):
					self.userExpression[samples[i]][gene] = fields[i]
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
		self.getClinVar( **kwargs )
		self.getExAC( **kwargs )
		self.getVEP( **kwargs )
	def getClinVar( self , **kwargs ):
#		print "charger - getClinVar"
		doClinVar = kwargs.get( 'clinvar' , True )
		summaryBatchSize = kwargs.get( 'summaryBatchSize' , 500 )
		searchBatchSize = kwargs.get( 'searchBatchSize' , 50 )
		if doClinVar:
			ent = entrezAPI()
			i = 0
			for varsStart in range( 0 , len( self.userVariants ) , int(searchBatchSize) ):
				varsEnd = varsStart + int(searchBatchSize)
				varsSet = self.userVariants[varsStart:varsEnd]
#				for var in varsSet:
#					i += 1
#					print str(i) + "\t" + var.genomicVar()
#				print str(varsStart) + ":" + str(varsEnd)
				ent.prepQuery( varsSet )
				ent.subset = entrezAPI.esearch
				ent.database = entrezAPI.clinvar
				clinvarsSet = ent.doBatch( summaryBatchSize )
				varsBoth = self.matchClinVar( varsSet , clinvarsSet )
				self.userVariants[varsStart:varsEnd] = varsBoth["userVariants"]
				self.clinvarVariants.update( varsBoth["clinvarVariants"] )
	def getExAC( self , **kwargs ):
#		print "charger - getExac"
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
#			print "Genomic_Variant\tAllele_Frequency"
			for var in self.userVariants:
				if var.genomicVar() in entries:
					alleleFrequency = entries[var.genomicVar()]
				else:
					alleleFrequency = None
				var.alleleFrequency = alleleFrequency
				if var.isFrequentAllele( threshold ):
#					print var.uniqueVariant() + " is NOT rare(" + str(threshold) + "): " + str(var.alleleFrequency)
					common += 1
				else:
					rare += 1
#				print var.genomicVar() + "\t" + str(var.alleleFrequency)
			elen = len(entries.keys())
			print "ExAC found " + str(common) + "common & " + str(rare) + "rare variants out of " + str(totalVars) + "total variants and " + str(elen) + "unique variants"
	def getVEP( self , **kwargs ):
#		print "charger - getVEP"
		doVEP = kwargs.get( 'vep' , True )
		if doVEP:
			vep = ensemblAPI()
			luv = len(self.userVariants)
			vepVariants = vep.annotateVariantsPost( self.userVariants )
			self.vepVariants = self.matchVEP( vepVariants )
			aluv = 0
#			print "\nvepVariants"
			if self.vepVariants:
				aluv = len(self.vepVariants)
#				for genVar in self.vepVariants:
#					print self.vepVariants[genVar].uniqueProteogenomicVar()
			luvafter = len(self.userVariants)
			print "\nVEP annotated userVariants " + str( luvafter )
#			if self.userVariants:
#				for var in self.userVariants:
#					for consequence in var.consequences:
#						print consequence.uniqueProteogenomicVar()
			print "VEP annotated " + str(aluv) + " from the original set of " + str(luv)

#### Helper methods for data retrieval ####
	def matchClinVar( self , userVariants , clinvarVariants ):
#		print "charger::matchClinVar - "
		for var in userVariants:
#			print "userVariant: " ,
#			var.printVariant(',')
			for uid in clinvarVariants:
				cvar = clinvarVariants[uid]
#				print "clinvarVariant:" ,
#				cvar.printVariant(',')
				if var.sameGenomicVariant( cvar ):
					var.clinical = cvar.clinical
		return { "userVariants" : userVariants , "clinvarVariants" : clinvarVariants }
	def matchVEP( self , vepVariants ):
#		print "charger::matchVEP - "
		for var in self.userVariants:
#			print ""
			for genVar in vepVariants:
				vepVar = vepVariants[genVar]
#				print "userVariant: " + var.genomicVar() ,
#				print " =? vepVariant: " + vepVar.genomicVar() + " with N=" + str( len( vepVar.consequences ) ) + " consequences"
				if var.sameGenomicVariant( vepVar ):
#					print "\tsame! consequences N=" ,
#					print str( len( vepVar.consequences ) ) ,
#					print " transfered to userVariant=" ,
					var.consequences = vepVar.consequences
#					print str( len( var.consequences ) ) ,
					var.colocatedVariants = vepVar.colocatedVariants
					if vepVar.strand: #MGI .maf annotator puts on positive strand
#						print str(vepVar.strand) + " from " + str(var.strand)
						var.strand = vepVar.strand
		return vepVariants
	def getDiseases( self , diseasesFile , **kwargs ):
#		print "charger::getDiseases - " ,
		tcga = kwargs.get( 'tcga' , True )
#		print tcga
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
	def PVS1( self , expressionThreshold = 0.05 ):
		print "CharGer module PVS1"
		print "- truncations in genes where LOF is a known mechanism of the disease"
		print "- require the mode of inheritance to be dominant (assuming heterzygosity) and co-occurence with reduced gene expression"
		truncations = ["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Nonstop_Mutation","Splice_Site"]
#		for gene in sorted(self.userGeneList.keys()):
#			for disease in sorted(self.userGeneList[gene].keys()):
#				print '\t'.join( [ gene , disease , self.userGeneList[gene][disease] ] )
		if self.userGeneList: #gene, disease, mode of inheritance
			for var in self.userVariants:
				varGene = var.gene
				varDisease = var.disease # no disease field in MAF; may require user input	
				varSample = var.sample
				varClass = var.variantClass
#				print '\t'.join( [ varGene , str(varDisease) , str(varClass) ] )
				if varClass in truncations:
					if varGene in self.userGeneList: # check if in gene list
						#var.PVS1 = True
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
		for var in self.userVariants:
			#varMAF = var.getExACasdf # Adam will update use alleleFrequency method
#			print var.genomicVar() + "\t" + str(var.alleleFrequency) + "\t" + str(threshold)
			if not var.isFrequentAllele( threshold ):
				var.PM2 = True
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
	def PP3( self ):
		print "CharGer module PP3: not yet implemented"
		print "- multiple lines of in silico evidence of deliterous effect"
		for var in self.userVariants:
#			print var.genomicVar()
#			print " consequences N=" ,
#			print str( len( var.consequences ) )
			for vcVar in var.consequences:
#				print var.codingHGVS()
				if not var.PP3:
					evidence = 0
					if vcVar.blosum:
						evidence += 1
					if vcVar.predictionSIFT:
						evidence += 1
					if vcVar.predictionPolyphen:
						evidence += 1
					if vcVar.compara:
						evidence += 1
					if vcVar.impact:
						evidence += 1
					if vcVar.maxentscan:
						evidence += 1
					if vcVar.genesplicer:
						evidence += 1
					if evidence > 2:
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
#			print "\tInput variant: " ,
#			var.printVariant(',')
#			print var.proteogenomicVar()
			ps1Call = False
			pm5Call = False
			call = False
			if mod == "PS1":
				call = var.PS1
			if mod == "PM5":
				call = var.PM5
			if not call: #is already true
				for genVar in self.clinvarVariants:
					cvar = self.clinvarVariants[genVar]
					clin = cvar.clinical
#					print str(cvar.uid) + ": " ,
#					print cvar.proteogenomicVar()
					if var.sameGenomicVariant( cvar ):
					#if genomic change is the same, then PS1
						if clin["description"] == clinvarVariant.pathogenic:
							if mod == "PS1":
								var.PS1 = True # already pathogenic still suffices to be PS1
								called += 1
					elif var.sameGenomicReference( cvar ):
					#if genomic change is different, but the peptide change is the same, then PS1
						if cvar.alternatePeptide == var.alternatePeptide: #same amino acid change
							if clin["description"] == clinvarVariant.pathogenic:
								if mod == "PS1":
									var.PS1 = True
									called += 1
					if var.samePeptideReference( cvar ):
#						print var.genomicVar() ,
#						print ',' ,
#						print var.codingHGVS() ,
#						print "\t" ,
#						print cvar.genomicVar() ,
#						print ',' ,
#						print cvar.codingHGVS()
						if not var.samePeptideChange( cvar ):
						#if peptide change is different, but the peptide reference is the same, then PM5
							if var.plausibleCodonFrame( cvar ):
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

### Classifier ###
	def classify( self ):
		for var in self.userVariants:
			var.isPathogenic( )
			var.isLikelyPathogenic( )
			var.isLikelyBenign( )
			var.isBenign( )
			var.isUncertainSignificance( )
	def printClassifications( self ):
		print '\t'.join( ["Variant" , "PositiveEvidence" , "CharGerClassification" , "ClinVarAnnoation"] )
		i = 0
		for var in self.userVariants:
			i += 1
			print '\t'.join( [ str(i) , var.uniqueVar() , var.positiveEvidence() , var.pathogenicity , var.clinical["description"] ] )

	@staticmethod
	def safeOpen( inputFile , rw , **kwargs ):
		errwar = kwargs.get( 'warning' , False )
		try:
			return open( inputFile , rw )
		except:
			if errwar:
				return "CharGer Warning: could not open " + inputFile
			else:
				return "CharGer Error: could not open " + inputFile
	@staticmethod
	def getUniqueGenomicVariantList( aVarList ):
#		print "charger::getUniqueGenomicVariantList"
		uniqueVarList = []
		uniqueVarDict = AV({})
		for var in aVarList:
			generalVar = variant()
			generalVar.copyInfo( var )
#			var.printVariant('_')
			if generalVar.genomicVar() not in uniqueVarDict:
#				generalVar.printVariant('__')
				uniqueVarDict[generalVar.genomicVar()] = 1
				uniqueVarList.append( generalVar )
		return uniqueVarList
