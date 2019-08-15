#!/usr/bin/env python
# chargervariant - CharGer annotated variants
# author: Adam D Scott (ascott@genome.wustl.edu) & Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2016*01*13

import pdb
from biomine.variant.clinvarvariant import clinvarvariant
from biomine.variant.vepvariant import vepvariant
from biomine.variant.mafvariant import mafvariant
from autovivification import autovivification

class chargervariant(mafvariant):
	# (rjm) module scores and pathogenicity threshold values now set in __main__
	pathogenic = "Pathogenic"
	likelyPathogenic = "Likely Pathogenic"
	likelyBenign = "Likely Benign"
	benign = "Benign"
	uncertain = "Uncertain Significance"
	def __init__( self , **kwargs ):
		super(chargervariant,self).__init__(**kwargs)
		self.PVS1 = kwargs.get( 'PVS1' , False )
		self.PS1 = kwargs.get( 'PS1' , False )
		self.PS2 = kwargs.get( 'PS2' , False )
		self.PS3 = kwargs.get( 'PS3' , False )
		self.PS4 = kwargs.get( 'PS4' , False )
		self.PM1 = kwargs.get( 'PM1' , False )
		self.PM2 = kwargs.get( 'PM2' , False )
		self.PM3 = kwargs.get( 'PM3' , False )
		self.PM4 = kwargs.get( 'PM4' , False )
		self.PM5 = kwargs.get( 'PM5' , False )
		self.PM6 = kwargs.get( 'PM6' , False )
		self.PP1 = kwargs.get( 'PP1' , False )
		self.PP2 = kwargs.get( 'PP2' , False )
		self.PP3 = kwargs.get( 'PP3' , False )
		self.PP4 = kwargs.get( 'PP4' , False )
		self.PP5 = kwargs.get( 'PP5' , False )

		self.BA1 = kwargs.get( 'BA1' , False )
		self.BS1 = kwargs.get( 'BS1' , False )
		self.BS2 = kwargs.get( 'BS2' , False )
		self.BS3 = kwargs.get( 'BS3' , False )
		self.BS4 = kwargs.get( 'BS4' , False )
		self.BP1 = kwargs.get( 'BP1' , False )
		self.BP2 = kwargs.get( 'BP2' , False )
		self.BP3 = kwargs.get( 'BP3' , False )
		self.BP4 = kwargs.get( 'BP4' , False )
		self.BP5 = kwargs.get( 'BP5' , False )
		self.BP6 = kwargs.get( 'BP6' , False )
		self.BP7 = kwargs.get( 'BP7' , False )

		self.PSC1 = kwargs.get( 'PSC1' , False )
		self.PMC1 = kwargs.get( 'PMC1' , False )
		self.PPC1 = kwargs.get( 'PPC1' , False )
		self.PPC2 = kwargs.get( 'PPC2' , False )

		self.BSC1 = kwargs.get( 'BSC1' , False )
		self.BMC1 = kwargs.get( 'BMC1' , False )

		self.otherTranscripts = kwargs.get( 'otherTranscripts' , {} )
		self.alleleFrequency = kwargs.get( 'alleleFrequency' , None )
		self.vepAnnotations = kwargs.get( 'VEP' , None )
		self.vcfHeaders = kwargs.get( 'headers' , None )
		self.vcfInfo = kwargs.get( 'INFO' , [] )
		self.pathogenicity = kwargs.get( 'pathogenicity' , \
			{ "CharGer" : chargervariant.uncertain , \
			  "ACMG" : chargervariant.uncertain 
			}
		)
		self.clinical = kwargs.get( 'clinical' , { "description" : chargervariant.uncertain , "review_status" : "" } )
		self.pathogenicScore = kwargs.get( 'pathogenicScore' , 0 )
		self.benignScore = kwargs.get( 'benignScore' , 0 )
		self.chargerScore = kwargs.get( 'chargerScore' , 0 )
		self.callSummary = kwargs.get( 'summary' , [] )
		aParentVariant = kwargs.get( 'parentVariant' , None )
		if aParentVariant:
			super( chargervariant , self ).copyInfo( aParentVariant )
		# make vep and clinvar variants attributes of chargervariant
		self.vepVariant = kwargs.get( 'vepvariant' , vepvariant() )
		self.clinvarVariant = kwargs.get( 'clinvarvariant' , clinvarvariant() )
		self.transvarVariant = kwargs.get( 'transvarvariant' , None )
		self.scoresMap = kwargs.get( 'scoresMap' , dict() )

	def annotations( self , delim = "\t" ):
		line = delim.join( [ self.gene , \
							 self.chromosome , \
							 self.start , \
							 self.stop , \
							 self.reference , \
							 self.alternate , \
							 self.variantClass , \
							 self.HGVSg() , \
							 self.HGVSc() , \
							 self.HGVSp() , \
							 self.alleleFrequency , \
							 self.returnMostSevereConsequence() , \
							 self.positiveEvidence() , \
							 self.negativeEvidence() , \
							 self.pathogenicScore , \
							 self.benignScore , \
							 self.chargerScore , \
							 self.clinicalDescription() , \
							 self.pathogenicity( system = "ACMG" ) , \
							 self.pathogenicity( system = "CharGer" ) , \
							 self.linkPubMed() , \
							 self.traits( delim = '|' ) , \
							 self.callSummary( delim = ' -- ' ) ] )
		return line

	def callSummary( self , delim = " -- " ):
		return delim.join( self.callSummary )

	def pathogenicity( self , system = "ACMG" ):
		try:
			return self.pathogenicity[system]
		except:
			return "NA" 

	def clinicalDescription( self ):
		try:
			return self.clinical["description"]
		except:
			return "NA"

	def traits( self , delim = '|' ):
		try:
			return self.clinvarVariant.getTraits( delim )
		except:
			return "NA"

	def linkPubMed( self ):
		try:
			return self.clinvarVariant.linkPubMed()
		except:
			return "NA"

	def returnMostSevereConsequence( self ):
		try:
			return self.vepVariant.mostSevereConsequence
		except:
			return "NA"
		
	@classmethod
	def annotationsHeader( cls , delim = "\t" ):
		delim = "\tCharGer_"
		line = "CharGer_"
		line += delim.join( [ "gene" , \
							  "chromosome" , \
							  "start" , \
							  "stop" , \
							  "reference" , \
							  "alternate" , \
							  "variant_class" , \
							  "HGVSg" , \
							  "HGVSc" , \
							  "HGVSp" , \
							  "allele_frequency" , \
							  "VEP_most_severe_consequence" , \
							  "positive_evidence" , \
							  "negative_evidence" , \
							  "pathogenic_score" , \
							  "benign_score" , \
							  "cumulative_score" , \
							  "ClinVar_clinical_description" , \
							  "pathogenicity_call_ACMG" , \
							  "pathogenicity_call_CharGer" , \
							  "PubMed_Link" , \
							  "ClinVar_traits" , \
							  "call_summary" ] )
		return line

	def copyInfo( self , other ):
		self.vepAnnotations = other.vepAnnotations
		self.vcfHeaders = other.vcfHeaders
		self.vcfInfo = other.vcfInfo
		self.pathogenicity = other.pathogenicity
		self.clinical = other.clinical
		self.otherTranscripts = other.otherTranscripts
		self.alleleFrequency = other.alleleFrequency
	def printVariant( self , delim , **kwargs ):
		onlyThisVariant = kwargs.get( 'minimal' , False )
		print( "chargervariant { " )
		if not onlyThisVariant:
			super( mafvariant , self ).printVariant( delim , **kwargs )
			if self.vepVariant:
				self.vepVariant.printVariant( delim , **kwargs )
			if self.clinvarVariant:
				self.clinvarVariant.printVariant( delim , **kwargs )
		print( "}" )
	def __nonzero__( self ):
		pathogenicitySet = True
		clinicalSet = True
		for k , v in self.__dict__.iteritems():
			if ( k == "pathogenicity" ):
				if ( v["CharGer"] == chargervariant.uncertain \
				and v["ACMG"] == chargervariant.uncertain ):
					pathogenicitySet = False
			elif ( k == "clinical" ):
				if ( v["description"] == chargervariant.uncertain \
				and v["review_status"] == "" ):
					clinicalSet = False
			elif ( self.checkIfRefAltStrand( k ) ):
				if ( self.nonzeroRefAltStrand( k ) ):
					return True
			else:
				if ( bool( v ) ):
					return True
		if ( pathogenicitySet and clinicalSet ):
			return True
		return False
	def fillMissingInfo( self , copy ):
		#print( "FMI: before" )
		#self.printVariant( ", " )
		super( chargervariant , self ).fillMissingInfo( copy )
		#self.printVariant( ", " )
		if isinstance( copy , chargervariant ):
			if not self.vepVariant:
				self.vepVariant = copy.vepVariant
				self.vepVariant.fillMissingInfo( copy.vepVariant )
			#self.printVariant( ", " )
			if not self.clinvarVariant:
				self.clinvarVariant = copy.clinvarVariant
				self.clinvarVariant.fillMissingInfo( copy.clinvarVariant )
			#self.printVariant( ", " )
		elif isinstance( copy , vepvariant ):
			if not self.vepVariant:
				self.vepVariant = copy
				self.vepVariant.fillMissingInfo( copy )
			#self.printVariant( ", " )
		elif isinstance( copy , clinvarvariant ):
			if not self.clinvarVariant:
				self.clinvarVariant = copy
				self.clinvarVariant.fillMissingInfo( copy )
			#self.printVariant( ", " )

	def copyMostSevereConsequence( self ):
		for consequence in self.vepVariant.consequences:
			if self.vepVariant.mostSevereConsequence in consequence.terms:
				if consequence.gene:
					self.vepVariant.gene = consequence.gene
					self.gene = consequence.gene
				if consequence.referencePeptide:
					self.vepVariant.referencePeptide = consequence.referencePeptide
					self.referencePeptide = consequence.referencePeptide
				if consequence.positionPeptide:
					self.vepVariant.positionPeptide = consequence.positionPeptide
					self.positionPeptide = consequence.positionPeptide
				if consequence.alternatePeptide:
					self.vepVariant.alternatePeptide = consequence.alternatePeptide
					self.alternatePeptide = consequence.alternatePeptide
				if consequence.transcriptPeptide:
					self.vepVariant.transcriptPeptide = consequence.transcriptPeptide
					self.transcriptPeptide = consequence.transcriptPeptide
				if consequence.positionCodon:
					self.vepVariant.positionCodon = consequence.positionCodon
					self.positionCodon = consequence.positionCodon
				if consequence.transcriptCodon:
					self.vepVariant.transcriptCodon = consequence.transcriptCodon
					self.transcriptCodon = consequence.transcriptCodon
				#print self.proteogenomicVar()
				#print self.vepVariant.proteogenomicVar()
				#print consequence.proteogenomicVar()
				#print ""
				if consequence.canonical:
					return
	def check( self , mod ):
		checks = self.checks()
		return checks[mod]
	def checks( self , **kwargs ):
		doPositive = kwargs.get( 'positive' , True )
		doNegative = kwargs.get( 'negative' , True )
		checks = autovivification({})
		mods = self.modules()
		if doPositive:
			checks[mods[0]] = self.PVS1
			checks[mods[1]] = self.PS1
			checks[mods[2]] = self.PS2
			checks[mods[3]] = self.PS3
			checks[mods[4]] = self.PS4
			checks[mods[5]] = self.PM1
			checks[mods[6]] = self.PM2
			checks[mods[7]] = self.PM3
			checks[mods[8]] = self.PM4
			checks[mods[9]] = self.PM5
			checks[mods[10]] = self.PM6
			checks[mods[11]] = self.PP1
			checks[mods[12]] = self.PP2
			checks[mods[13]] = self.PP3
			checks[mods[14]] = self.PP4
			checks[mods[15]] = self.PP5

			checks[mods[28]] = self.PSC1
			checks[mods[29]] = self.PMC1
			checks[mods[30]] = self.PPC1
			checks[mods[31]] = self.PPC2
		if doNegative:
			checks[mods[16]] = self.BA1
			checks[mods[17]] = self.BS1
			checks[mods[18]] = self.BS2
			checks[mods[19]] = self.BS3
			checks[mods[20]] = self.BS4
			checks[mods[21]] = self.BP1
			checks[mods[22]] = self.BP2
			checks[mods[23]] = self.BP3
			checks[mods[24]] = self.BP4
			checks[mods[25]] = self.BP5
			checks[mods[26]] = self.BP6
			checks[mods[27]] = self.BP7

			checks[mods[32]] = self.BSC1
			checks[mods[33]] = self.BMC1
		return checks
	def modules( self ):
		return [ 'PVS1' , \
		'PS1' , 'PS2' , 'PS3' , 'PS4' ,  \
		'PM1' , 'PM2' , 'PM3' , 'PM4' , 'PM5' , 'PM6' , \
		'PP1' , 'PP2' , 'PP3' , 'PP4' , 'PP5' , \
		'BA1' , \
		'BS1' , 'BS2' , 'BS3' , 'BS4' , \
		'BP1' , 'BP2' , 'BP3' , 'BP4' , 'BP5' , 'BP6' , 'BP7', \
		'PSC1' , 'PMC1' , 'PPC1' , 'PPC2' , \
		'BSC1' , 'BMC1' ]
	def readVCFInfo( self , **kwargs ):
		for info in self.vcfInfo:
			print info
	def hasAlleleFrequency( self ):
		if var.alleleFrequency is None:
			return False
		return True
	def isFrequentAllele( self , threshold ):
		#print( "isFrequencyAllele? " + str( self.alleleFrequency ) + " > " + str( threshold ) )
		try:
			if self.alleleFrequency is None:
				return False
			if float( self.alleleFrequency ) > float( threshold ):
				return True
		except:
			print( "CharGer chargervariant::isFrequentAllele Warning: alleleFrequency not floatable (" \
				+ str( self.alleleFrequency ) + ") for " + self.genomicVar() )
			pass
		return False
	def countModule( self , module , weighted , scoresMap ):
		if ( self.__dict__.get( module ) ):
			if ( weighted ):
				return scoresMap.get( module , 0 )
			else:
				return 1
		return 0
	def countPathogenicVeryStrong( self , weighted , **kwargs ):
		count = 0
		if self.PVS1:
			count += self.countModule( 'PVS1' , weighted , kwargs[ 'scoresMap' ] )
		return count
	def countPathogenicStrong( self , weighted , **kwargs ):
		count = 0
		count += self.countModule( 'PS1' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PS2' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PS3' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PS4' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PSC1' , weighted , kwargs[ 'scoresMap' ] )
		return count
	def countPathogenicModerate( self , weighted , **kwargs ):
		count = 0
		count += self.countModule( 'PM1' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PM2' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PM3' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PM4' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PM5' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PM6' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PMC1' , weighted , kwargs[ 'scoresMap' ] )
		return count
	def countPathogenicSupport( self , weighted , **kwargs ):
		count = 0
		count += self.countModule( 'PP1' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PP2' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PP3' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PP4' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PP5' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PPC1' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'PPC2' , weighted , kwargs[ 'scoresMap' ] )
		return count
	def isPathogenic( self , **kwargs ):
		weighted = False
		numPathogenicStrong = self.countPathogenicStrong( weighted , **kwargs )
		numPathogenicModerate = self.countPathogenicModerate( weighted , **kwargs )
		numPathogenicSupport = self.countPathogenicSupport( weighted , **kwargs )
		if self.PVS1:
			if numPathogenicStrong >= 1 or \
			numPathogenicModerate >= 2 or \
			(numPathogenicModerate+numPathogenicSupport) >= 2 or \
			numPathogenicSupport >= 2:
				self.setAsPathogenic( **kwargs )
				return True
		elif numPathogenicStrong >= 2:
			self.setAsPathogenic( **kwargs )
			return True
		elif numPathogenicStrong >= 1:
			if numPathogenicModerate >= 3 or \
			(numPathogenicModerate == 2 and numPathogenicSupport >= 2) or \
			(numPathogenicModerate == 1 and numPathogenicSupport >= 4):
				self.setAsPathogenic( **kwargs )
				return True
	def isLikelyPathogenic( self , **kwargs ):
		weighted = False
		numPathogenicStrong = self.countPathogenicStrong( weighted , **kwargs )
		numPathogenicModerate = self.countPathogenicModerate( weighted , **kwargs )
		numPathogenicSupport = self.countPathogenicSupport( weighted , **kwargs )
		if numPathogenicStrong and numPathogenicModerate == 1:
			self.setAsLikelyPathogenic( **kwargs )
			return True
		if numPathogenicStrong == 1 and ( numPathogenicModerate == 1 or numPathogenicModerate == 2 ):
			self.setAsLikelyPathogenic( **kwargs )
			return True
		if numPathogenicStrong == 1 and numPathogenicSupport >= 2:
			self.setAsLikelyPathogenic( **kwargs )
			return True
		if numPathogenicModerate >= 3:
			self.setAsLikelyPathogenic( **kwargs )
			return True
		if numPathogenicModerate == 2 and numPathogenicSupport >= 2:
			self.setAsLikelyPathogenic( **kwargs )
			return True
		if numPathogenicModerate == 1 and numPathogenicSupport >= 4:
			self.setAsLikelyPathogenic( **kwargs )
			return True
	def countBenignStandAlone( self , weighted , **kwargs ):
		count = 0
		count += self.countModule( 'BA1' , weighted , kwargs[ 'scoresMap' ] )
		return count
	def countBenignStrong( self , weighted , **kwargs ):
		count = 0
		count += self.countModule( 'BS1' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'BS2' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'BS3' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'BS4' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'BSC1' , weighted , kwargs[ 'scoresMap' ] )
		return count
	def countBenignModerate( self , weighted , **kwargs ):
		count = 0
		count += self.countModule( 'BMC1' , weighted , kwargs[ 'scoresMap' ] )
		return count
	def countBenignSupport( self , weighted , **kwargs ):
		count = 0
		count += self.countModule( 'BP1' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'BP2' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'BP3' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'BP4' , weighted , kwargs[ 'scoresMap' ] )
		count += self.countModule( 'BP5' , weighted , kwargs[ 'scoresMap' ] )
		return count
	def isLikelyBenign( self , **kwargs ):
		weighted = False
		numBenignStrong = self.countBenignStrong( weighted , **kwargs )
		numBenignModerate = self.countBenignModerate( weighted , **kwargs )
		numBenignSupport = self.countBenignSupport( weighted , **kwargs )
		if numBenignStrong == 1 and \
		numBenignSupport == 1:
			return True
		if numBenignSupport >= 2:
			return True
		return False
	def isBenign( self , **kwargs ):
		weighted = False
		numBenignStandAlone = self.countBenignStandAlone( weighted , **kwargs )
		numBenignStrong = self.countBenignStrong( weighted , **kwargs )
		numBenignModerate = self.countBenignModerate( weighted , **kwargs )
		numBenignSupport = self.countBenignSupport( weighted , **kwargs )
		if numBenignStandAlone:
			return True
		if numBenignStrong >= 2:
			return True
		if numBenignStrong == 1 \
		and numBenignModerate == 1 \
		and numBenignSupport >= 2:
			return True
		return False

	def isUncertainSignificance( self , **kwargs ):
		if ( self.isPathogenic( **kwargs ) or self.isLikelyPathogenic( **kwargs ) )and \
		( self.isBenign( **kwargs ) or self.isLikelyBenign( **kwargs ) ):
			self.setAsUncertainSignificance( **kwargs )
			return True
	def setAsPathogenic( self , **kwargs ):
		scoreSystem = kwargs.get( "system" , "CharGer" )
		self.pathogenicity[scoreSystem] = chargervariant.pathogenic
	def setAsLikelyPathogenic( self , **kwargs ):
		scoreSystem = kwargs.get( "system" , "CharGer" )
		self.pathogenicity[scoreSystem] = chargervariant.likelyPathogenic
	def setAsUncertainSignificance( self , **kwargs ):
		scoreSystem = kwargs.get( "system" , "CharGer" )
		self.pathogenicity[scoreSystem] = chargervariant.uncertain
	def setAsLikelyBenign( self , **kwargs ):
		scoreSystem = kwargs.get( "system" , "CharGer" )
		self.pathogenicity[scoreSystem] = chargervariant.likelyBenign
	def setAsBenign( self , **kwargs ):
		scoreSystem = kwargs.get( "system" , "CharGer" )
		self.pathogenicity[scoreSystem] = chargervariant.benign
	def positiveEvidence( self ):
		positive = []
		checks = self.checks( negative=False )
		for k in sorted(checks.keys()):
			if checks[k]:
				positive.append(k)
		return ",".join(positive)
	def negativeEvidence( self ):
		negative = []
		checks = self.checks( positive=False )
		for k in sorted(checks.keys()):
			if checks[k]:
				negative.append(k)
		return ",".join(negative)

	def compareDescription( self , checkBetter , **kwargs ):
		desc = self.clinvarVariant.clinical["description"]
		if desc.lower() == checkBetter.lower():
			return True
		else:
			return False
	def tallyScore( self , **kwargs ):
		scoreSystem = kwargs.get( "system" , "CharGer" )
		if scoreSystem != "CharGer":
			return
		weighted = True
		self.pathogenicScore = 0
		self.pathogenicScore += self.countPathogenicSupport( weighted , **kwargs )
		self.pathogenicScore += self.countPathogenicModerate( weighted , **kwargs )
		self.pathogenicScore += self.countPathogenicStrong( weighted , **kwargs )
		self.pathogenicScore += self.countPathogenicVeryStrong( weighted , **kwargs )
		self.benignScore = 0
		self.benignScore += self.countBenignSupport( weighted , **kwargs )
		self.benignScore += self.countBenignModerate( weighted , **kwargs )
		self.benignScore += self.countBenignStrong( weighted , **kwargs )
		self.benignScore += self.countBenignStandAlone( weighted , **kwargs )
		self.chargerScore = self.pathogenicScore + self.benignScore
		scoreSystem = kwargs.get( "system" , "CharGer" )
		self.compositeScore( **kwargs )
	def compositeScore( self , **kwargs ):
		scoreSystem = kwargs.get( "system" , "CharGer" )
		override = kwargs.get( 'override' , False )
		minPathogenicScore = kwargs[ 'scoresMap' ][ 'minPathogenicScore' ]
		minLikelyPathogenicScore = kwargs[ 'scoresMap' ][ 'minLikelyPathogenicScore' ]
		maxBenignScore = kwargs[ 'scoresMap' ][ 'maxBenignScore' ]
		maxLikelyBenignScore = kwargs[ 'scoresMap' ][ 'maxLikelyBenignScore' ]
		if self.chargerScore >= minPathogenicScore:
			self.setAsPathogenic( **kwargs )
		elif self.chargerScore >= minLikelyPathogenicScore:
			self.setAsLikelyPathogenic( **kwargs )
		elif self.benignScore <= maxLikelyBenignScore:
			if self.pathogenicity[scoreSystem] == chargervariant.pathogenic \
			or self.pathogenicity[scoreSystem] == chargervariant.likelyPathogenic:
				self.setAsUncertainSignificance( **kwargs ) 
			else:
				if self.benignScore <= maxBenignScore:
					self.setAsBenign( **kwargs )
				elif self.benignScore <= maxLikelyBenignScore:
					self.setAsLikelyBenign( **kwargs )
		else:
			self.setAsUncertainSignificance( **kwargs ) 
		if override:
			self.clinvarOverride( **kwargs )
	def clinvarOverride( self , **kwargs ):
		scoreSystem = kwargs.get( "system" , "CharGer" )
		if not self.clinvarVariant:
			return
		desc = ""
		try:
			desc = repr( self.clinvarVariant.clinical["description"].lower() )
		except:
			print( "CharGer::chargervariant::clinvarOverride Warning: " \
				+ "- no description for clinical information for " \
				+ self.genomicVar() )
		call = repr( self.pathogenicity[scoreSystem].lower() )
		path = repr( chargervariant.pathogenic.lower() )
		lPath = repr( chargervariant.likelyPathogenic.lower() )
		lBen = repr( chargervariant.likelyBenign.lower() )
		ben = repr( chargervariant.benign.lower() )
		if call == lPath:
			if desc == path:
				self.pathogenicity[scoreSystem] = chargervariant.pathogenic
		elif call == lBen:
			if desc == path:
				self.pathogenicity[scoreSystem] = chargervariant.pathogenic
			elif desc == lPath:
				self.pathogenicity[scoreSystem] = chargervariant.likelyPathogenic
		elif call == ben:
			if desc == path:
				self.pathogenicity[scoreSystem] = chargervariant.pathogenic
			elif desc == lPath:
				self.pathogenicity[scoreSystem] = chargervariant.likelyPathogenic
			elif desc == lBen:
				self.pathogenicity[scoreSystem] = chargervariant.likelyBenign
	def addSummary( self , string , **kwargs ):
		self.callSummary.append( string )
