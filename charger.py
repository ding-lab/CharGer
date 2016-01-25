#!/usr/bin/python
# CharGer - Characterization of Germline variants
# author: Adam D Scott (ascott@genome.wustl.edu) & Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2015*12

import math
from WebAPI.Entrez.entrezAPI import entrezAPI
from WebAPI.ExAC.exacAPI import exacAPI
from WebAPI.Variant.clinvarVariant import clinvarVariant
from chargerVariant import chargerVariant
from autovivification import autovivification

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
	def __init__( self , **kwargs ):
		self.userVariants = kwargs.get( 'variants' , [] )
		self.userExpression = kwargs.get( 'expressions' , autovivification({}) )
		self.userGeneList = kwargs.get( 'geneList' , autovivification({}) )
		self.userDeNovoVariants = kwargs.get( 'deNovo' , {} )
		self.userAssumedDeNovoVariants = kwargs.get( 'assumedDeNovo' , {} )
		self.userCoSegregateVariants = kwargs.get( 'coSegregate' , {} )
		self.clinvarVariants = kwargs.get( 'clinvarVariants' , {} )

### Retrieve input data from user ###
	def getInputData( self  , **kwargs ):
		mafFile = kwargs.get( 'maf' , "" )
		expressionFile = kwargs.get( 'expression' , "" )
		geneListFile = kwargs.get( 'geneList' , "" )
		deNovoFile = kwargs.get( 'deNovo' , "" )
		assumedDeNovoFile = kwargs.get( 'assumedDeNovo' , "" )
		coSegregateFile = kwargs.get( 'coSegregate' , "" )
		self.readMAF( mafFile , **kwargs )
		self.readExpression( expressionFile )
		self.readGeneList( geneListFile )
	def readMAF( self , inputFile , **kwargs ):
#		print "\tReading .maf!"
		inFile = self.safeOpen( inputFile , 'r' )
		codonColumn = kwargs.get( 'codon' , 48 )
		peptideChangeColumn = kwargs.get( 'peptideChange' , 49 )
		try:
			next(inFile)
			for line in inFile:
				var = chargerVariant()
				var.mafLine2Variant( line , peptideChange=peptideChangeColumn , codon=codonColumn )
				var.genomicVar()
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
	def readGeneList( self , inputFile , diseaseSpecific = True ): # gene list formatted "gene", "disease", "mode of inheritance"
		try:
			inFile = self.safeOpen( inputFile , 'r' , warning=True )
			for line in inFile:
				fields = line.split( "\t" )
				gene = fields[0]
				if diseaseSpecific:
					disease = fields[1]
				else: #set the gene to match all disease
					disease = "all" 
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
		doClinVar = kwargs.get( 'clinvar' , True )
		doExAC = kwargs.get( 'exac' , True )
		if doClinVar:
			self.getClinVar( **kwargs )
		if doExAC:
			self.getExAC( **kwargs )
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
				for var in varsSet:
					i += 1
			#		print str(i) + "\t" + var.genomicVar()
				#print str(varsStart) + ":" + str(varsEnd)
				ent.prepQuery( varsSet )
				ent.subset = entrezAPI.esearch
				ent.database = entrezAPI.clinvar
				clinvarsSet = ent.doBatch( summaryBatchSize )
				varsBoth = self.matchClinVar( varsSet , clinvarsSet )
				self.userVariants[varsStart:varsEnd] = varsBoth["userVariants"]
				self.clinvarVariants.update( varsBoth["clinvarVariants"] )
	def getExAC( self , **kwargs ):
		print "charger - getExac"
		doExAC = kwargs.get( 'exac' , True )
		useHarvard = kwargs.get( 'harvard' , True )
		threshold = kwargs.get( 'threshold' , 0 )
		if doExAC:
			common = 0
			rare = 0
			totalVars = len( self.userVariants )
			exac = exacAPI(harvard=useHarvard)
			entries = exac.getAlleleFrequencies( self.userVariants )
			print "Genomic_Variant\tAllele_Frequency"
			alleleFrequency = None
			for var in self.userVariants:
				if var.genomicVar() in entries:
					alleleFrequency = entries[var.genomicVar()]
				else:
					alleleFrequency = None
				var.alleleFrequency = alleleFrequency
				if var.isFrequentAllele( threshold ):
			#		print var.uniqueVariant() + " is NOT rare(" + str(threshold) + "): " + str(var.alleleFrequency)
					common += 1
				else:
					rare += 1
				print var.genomicVar() + "\t" + str(var.alleleFrequency)
			print "ExAC found " + str(common) + "common & " + str(rare) + "rare variants out of " + str(totalVars) + "total variants"
	def getVEP( self , **kwargs ):
#		print "charger - getVEP"
		doVEP = kwargs.get( 'vep' , True )
		if doVEP:
			vep = ensemblAPI( )
			vep.annotateVariants( self.userVariants )

#### Helper methods for data retrieval ####
	def matchClinVar( self , userVariants , clinvarVariants ):
		#print "charger::matchClinVar - "
		for var in userVariants:
			#print "userVariant: " ,
			#var.printVariant(',')
			for uid in clinvarVariants:
				cvar = clinvarVariants[uid]
				#print "clinvarVariant:" ,
				#cvar.printVariant(',')
				if var.sameGenomicVariant( cvar ):
					var.clinical = cvar.clinical
		return { "userVariants" : userVariants , "clinvarVariants" : clinvarVariants }

### Evidence levels ### 
##### Very Strong #####
	def PVS1( self , expressionThreshold = 0.05 ):
		print "CharGer module PVS1"
		print "- truncations in genes where LOF is a known mechanism of the disease"
		print "- require the mode of inheritance to be dominant (assuming heterzygosity) and co-occurence with reduced gene expression"
		truncations = ["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Nonstop_Mutation","Splice_Site"]
		if self.userGeneList: #gene, disease, mode of inheritance
			for var in self.userVariants:
				varGene = var.gene
				varDisease = var.disease # no disease field in MAF; may require user input	
				varSample = var.sample
				varClass = var.variantClass
				if varClass in truncations:
					if varGene in self.userGeneList: # check if in gene list
						if ( "dominant" in self.userGeneList[varGene][varDisease] or \
							"dominant" in self.userGeneList[varGene]["all"]):
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
	def PS4( self ): # not relevant in rare variants, such big effect size common variants are so rare may as well just take a input variant list
		print "CharGer module PS4: not yet implemented"
		#print "- variant prevalence in cases significantly greater than controls"
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
		#print "- "
	def PM2( self , threshold ):
		print "CharGer module PM2"
		print "- absent or extremely low frequency in controls"
		for var in self.userVariants:
			#varMAF = var.getExACasdf # Adam will update use alleleFrequency method
			print var.genomicVar() + "\t" + str(var.alleleFrequency) + "\t" + str(threshold)
			if not var.isFrequentAllele( threshold ):
				var.PM2 = True
	def PM3( self ):
		print "CharGer module PM3: not yet implemented"
		#print "- "
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
		#print "- "
	def PP3( self ):
		print "CharGer module PP3: not yet implemented"
		#print "- "
	def PP4( self ):
		print "CharGer module PP4: not yet implemented"
		#print "- "
	def PP5( self ):
		print "CharGer module PP5: not yet implemented"
		#print "- "

### helper functions of evidence levels ###
	def peptideChange( self , mod ):
		called = 0
		for var in self.userVariants:
			uniVar = var.uniqueVar()
	#		print "\tInput variant: " ,
			#var.printVariant(',')
			print var.proteogenomicVar()
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
