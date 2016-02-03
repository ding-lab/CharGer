#!/usr/bin/python
# chargerVariant - CharGer annotated variants
# author: Adam D Scott (ascott@genome.wustl.edu) & Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2016*01*13

from WebAPI.Variant.clinvarVariant import clinvarVariant
from WebAPI.Variant.vepVariant import vepVariant
from autovivification import autovivification

class chargerVariant(clinvarVariant,vepVariant):
	pathogenic = "Pathogenic"
	likelyPathogenic = "Likely Pathogenic"
	likelyBenign = "Likely Benign"
	benign = "Benign"
	uncertain = "Uncertain Significance"
	def __init__( self , **kwargs ):
		super(chargerVariant,self).__init__(**kwargs)
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
		self.otherTranscripts = kwargs.get( 'otherTranscripts' , {} )
		self.alleleFrequency = kwargs.get( 'alleleFrequency' , None )
		self.pathogenicity = kwargs.get( 'pathogenicity' , chargerVariant.uncertain )
		self.clinical = kwargs.get( 'clinical' , { "description" : chargerVariant.uncertain , "review_status" : "" } )
		aParentVariant = kwargs.get( 'parentVariant' , None )
		if aParentVariant:
			super( chargerVariant , self ).copyInfo( aParentVariant )
	def check( self , mod ):
		checks = self.checks()
		return checks[mod]
	def checks( self ):
		#print "chargerVariant::checks"
		checks = autovivification({})
		mods = self.modules()
		#print mods
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
		#print checks
		return checks
	def modules( self ):
		return ['PVS1' , \
		'PS1' , 'PS2' , 'PS3' , 'PS4' , \
		'PM1' , 'PM2' , 'PM3' , 'PM4' , 'PM5' , 'PM6' , \
		'PP1' , 'PP2' , 'PP3' , 'PP4' , 'PP5' , \
		'BA1' , \
		'BS1' , 'BS2' , 'BS3' , 'BS4' , \
		'BP1' , 'BP2' , 'BP3' , 'BP4' , 'BP5' , 'BP6' , 'BP7' \
		]
	def hasAlleleFrequency( self ):
		if var.alleleFrequency == None:
			return False
		return True
	def isFrequentAllele( self , threshold ):
		if self.alleleFrequency > threshold:
			return True
		return False
	def countPathogenicStrong( self ):
		count = 0
		if self.PS1:
			count += 1
		if self.PS2:
			count += 1
		if self.PS3:
			count += 1
		if self.PS4:
			count += 1
		return count
	def countPathogenicModerate( self ):
		count = 0
		if self.PM1:
			count += 1
		if self.PM2:
			count += 1
		if self.PM3:
			count += 1
		if self.PM4:
			count += 1
		if self.PM5:
			count += 1
		if self.PM6:
			count += 1
		return count
	def countPathogenicSupport( self ):
		count = 0
		if self.PP1:
			count += 1
		if self.PP2:
			count += 1
		if self.PP3:
			count += 1
		if self.PP4:
			count += 1
		if self.PP5:
			count += 1
		return count
	def isPathogenic( self ):
		numPathogenicStrong = self.countPathogenicStrong()
		numPathogenicModerate = self.countPathogenicModerate()
		numPathogenicSupport = self.countPathogenicSupport()
		if self.PVS1:
			if numPathogenicStrong >= 1 or \
			numPathogenicModerate >= 2 or \
			(numPathogenicModerate+numPathogenicSupport) >= 2 or \
			numPathogenicSupport >= 2:
				self.setAsPathogenic()
				return True
		elif numPathogenicStrong >= 2:
			self.setAsPathogenic()
			return True
		elif numPathogenicStrong >= 1:
			if numPathogenicModerate >= 3 or \
			(numPathogenicModerate == 2 and numPathogenicSupport >= 2) or \
			(numPathogenicModerate == 1 and numPathogenicSupport >= 4):
				self.setAsPathogenic()
				return True
	def isLikelyPathogenic( self ):
		numPathogenicStrong = self.countPathogenicStrong()
		numPathogenicModerate = self.countPathogenicModerate()
		numPathogenicSupport = self.countPathogenicSupport()
		if numPathogenicStrong and numPathogenicModerate == 1:
			self.setAsLikelyPathogenic()
			return True
		if numPathogenicStrong == 1 and ( numPathogenicModerate == 1 or numPathogenicModerate == 2 ):
			self.setAsLikelyPathogenic()
			return True
		if numPathogenicStrong == 1 and numPathogenicSupport >= 2:
			self.setAsLikelyPathogenic()
			return True
		if numPathogenicModerate >= 3:
			self.setAsLikelyPathogenic()
			return True
		if numPathogenicModerate == 2 and numPathogenicSupport >= 2:
			self.setAsLikelyPathogenic()
			return True
		if numPathogenicModerate == 1 and numPathogenicSupport >= 4:
			self.setAsLikelyPathogenic()
			return True
	def countBenignStrong( self ):
		count = 0
		if self.BS1:
			count += 1
		if self.BS2:
			count += 1
		if self.BS3:
			count += 1
		if self.BS4:
			count += 1
		return count
	def countBenignSupport( self ):
		count = 0
		if self.BP1:
			count += 1
		if self.BP2:
			count += 1
		if self.BP3:
			count += 1
		if self.BP4:
			count += 1
		if self.BP5:
			count += 1
		return count
	def isLikelyBenign( self ):
		numBenignStrong = self.countBenignStrong()
		numBenignSupport = self.countBenignSupport()
		if numBenignStrong == 1 and \
		numBenignSupport == 1:
			return True
		if numBenignSupport >= 2:
			return True
		return False
	def isBenign( self ):
		numBenignStrong = self.countBenignStrong()
		if self.BA1:
			return True
		if numBenignStrong >= 2:
			return True
		return False

	def isUncertainSignificance( self ):
		if ( self.isPathogenic() or self.isLikelyPathogenic() )and \
		( self.isBenign() or self.isLikelyBenign() ):
			self.setAsUncertainSignificance()
			return True
	def setAsPathogenic( self ):
		self.pathogenicity = chargerVariant.pathogenic
	def setAsLikelyPathogenic( self ):
		self.pathogenicity = chargerVariant.likelyPathogenic
	def setAsLikelyBenign( self ):
		self.pathogenicity = chargerVariant.likelyBenign
	def setAsBenign( self ):
		self.pathogenicity = chargerVariant.benign
	def positiveEvidence( self ):
		positive = []
		checks = self.checks()
		for k in sorted(checks.keys()):
#			print k + "\t" + str(checks[k])
			if checks[k]:
				positive.append(k)
		return ",".join(positive)

