#!/usr/bin/python
# CharGer - Characterization of Germline variants
# author: Adam D Scott (ascott@genome.wustl.edu) & Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2015*12

import math
from entrezAPI import entrezAPI
from exacAPI import exacAPI
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
		self.userExpression = kwargs.get( 'expressions' , [] )
		self.userGeneList = kwargs.get( 'geneList' , autovivification({}) )
		self.clinvarVariants = kwargs.get( 'clinvarVariants' , {} )

### Retrieve input data from user ###
	def getInputData( self  , **kwargs ):
		mafFile = kwargs.get( 'maf' , "" )
		expressionFile = kwargs.get( 'expression' , "" )
		geneListFile = kwargs.get( 'geneList' , "" )
		self.readMAF( mafFile )
		self.readExpression( expressionFile )
		self.readGeneList( geneListFile )
	def readMAF( self , inputFile ):
		print "\tSplitting .maf by variant type!"
		if inputFile:
			inFile = open( inputFile , 'r' )
			next(inFile)
			for line in inFile:
				fields = line.split( "\t" )
				var = chargerVariant()
				var.mafLine2Variant( line )
				self.userVariants.append( var )
	def readExpression( self , inputFile ): # expect a sample(col)-gene(row) matrix
		if inputFile:
			inFile = open( inputFile , 'r' )
			header = inFile.readline() # for future fetch header to get other field
			samples = header.split( "\t" )
			for line in inFile:
				fields = line.split( "\t" )
				gene = fields[0]
				for i in range(1,len(fields)):
					self.userExpression[samples[i]][gene] = fields[i]
	def readGeneList( self , inputFile , diseaseSpecific = True ): # gene list formatted "gene", "disease", "mode of inheritance"
		if inputFile:
			inFile = open( inputFile , 'r' )
			for line in inFile:
				fields = line.split( "\t" )
				gene = fields[0]
				if diseaseSpecific:
					disease = fields[1]
				else: #set the gene to match all disease
					disease = "all" 
				mode_inheritance = fields[2]
				self.userGeneList[gene][disease] = mode_inheritance

### Retrieve external reference data ###
	def getExternalData( self , **kwargs ):
		doClinVar = kwargs.get( 'clinvar' , True )
		doExAC = kwargs.get( 'exac' , True )
		if doClinVar:
			self.getClinVar( **kwargs )
		if doExAC:
			self.getExAC( **kwargs )
	def getClinVar( self , **kwargs ):
		doClinVar = kwargs.get( 'clinvar' , True )
		if doClinVar:
			ent = entrezAPI()
			ent.prepQuery( self.userVariants )
			ent.database = entrezAPI.clinvar
			self.clinvarVariants = ent.doBatch( 5 )
	def getExAC( self , **kwargs ):
		doExAC = kwargs.get( 'exac' , True )
		useHarvard = kwargs.get( 'harvard' , True )
		threshold = kwargs.get( 'threshold' , 0 )
		if doExAC:
			exac = exacAPI(harvard=useHarvard)
			exac.getAlleleFrequencies( self.userVariants )
			for var in self.userVariants:
				if var.isFrequentAllele( threshold ):
					print var.uniqueVariant() + " is NOT rare(" + str(threshold) + "): " + str(var.alleleFrequency)

### Evidence levels ### 
##### Very Strong #####
	def PVS1( self , expressionThreshold = 0.05 ):
		truncations = ["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Nonstop_Mutation","Splice_Site"]
		if self.userGeneList: #gene, disease, mode of inheritance
			for var in self.userVariants:
				varGene = var.gene
				varDisease = var.disease	
				varSample = var.sample
				varClass = var.variantClass
				if varClass in truncations:
					if varGene in geneList: # check if in gene list
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
		print "- same peptide change that is pathogenic and is a different genomic variant of the same reference peptide"
		self.peptideChange( "PS1" )
	def PS2( self ):
		NotImplemented
	def PS3( self ):
		NotImplemented
	def PS4( self ):
		for var in self.userVariants:
			caseVarFreq = "NEED UPDATE" # may take input from current MAF
			controlVarFreq = "NEED UPDATE" # may take input from ExAC
			OR = (caseVarFreq/controlVarFreq) / (1-caseVarFreq)/(1-controlVarFreq)
			# Adam will update
			if OR >= 5:
				CIlower = math.log(OR) - math.sqrt( 1/caseVarFreq + 1/controlVarFreq + 1/caseVarFreq + 1/controlVarFreq)
				if (CIlower > 1):
					var.PS4 = True

##### Moderate #####
	def PM1( self ):
		NotImplemented
	def PM2( self ):
		for var in self.userVariants:
			#varMAF = var.getExACasdf # Adam will update use alleleFrequency method
			if isFrequentAllele(var):
				var.PM2 = True
	def PM3( self ):
		NotImplemented
	def PM4( self ):
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
		NotImplemented

##### Supporing #####
	def PP1( self ):
		NotImplemented
	def PP2( self ):
		NotImplemented
	def PP3( self ):
		NotImplemented
	def PP4( self ):
		NotImplemented
	def PP5( self ):
		NotImplemented

### helper functions of evidence levels ###
	def peptideChange( self , mod ):
		for var in self.userVariants:
			uniVar = var.uniqueVar()
			#print "\tInput variant: " + genVar , 
			canBePS1 = True
			canBePM5 = True
			pm1Call = False
			pm5Call = False
			call = var.PS1
			#print "Call: " + genVar ,
			#print " => " + str(call)
			if not call: #is already true
				#print "checking"
				call = False
				for genVar in self.clinvarVariants:
					cvar = self.clinvarVariants[genVar]
					clin = cvar.clinical
					if cvar.chromosome == var.chromosome and \
						cvar.start == var.start and \
						cvar.stop == var.stop and \
						cvar.reference == var.reference and \
						cvar.referencePeptide == var.referencePeptide and \
						cvar.positionPeptide == var.positionPeptide: #same genomic position & reference
						if cvar.alternatePeptide == var.alternatePeptide: #same amino acid change
							if clin["description"] == "Pathogenic":
								#print "Already called pathogenic: " ,
								#var.printVariant(' ')
								canBePS1 = False
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
					if canBePS1:
						call = pm1Call
				if mod == "PM5":
					if canBePM5:
						call = pm5Call
			if mod == "PS1":
				var.PS1 = call
			if mod == "PM5":
				var.PM5 = call
		print ""
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
