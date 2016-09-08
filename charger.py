#!/usr/bin/python
# CharGer - Characterization of Germline variants
# author: Adam D Scott (ascott@genome.wustl.edu) & Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2015*12

import time
import math
import re
from biomine.webapi.ensembl.ensemblapi import ensemblapi
import glob
from scipy import stats
from biomine.webapi.entrez.entrezapi import entrezapi
from biomine.webapi.exac.exacapi import exacapi
from biomine.variant.clinvarvariant import clinvarvariant
from biomine.variant.vepvariant import vepvariant
from biomine.variant.vepconsequencevariant import vepconsequencevariant
from biomine.variant.mafvariant import mafvariant
from biomine.variant.variant import variant
from chargervariant import chargervariant
from autovivification import autovivification as AV
import vcf
from collections import OrderedDict as OD

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
		self.pathogenicVariants = kwargs.get( 'pathogenic' , AV({}) )
		self.userExpression = kwargs.get( 'expressions' , AV({}) )
		self.userGeneList = kwargs.get( 'geneList' , AV({}) )
		self.userDeNovoVariants = kwargs.get( 'deNovo' , {} )
		self.userAssumedDeNovoVariants = kwargs.get( 'assumedDeNovo' , {} )
		self.userCoSegregateVariants = kwargs.get( 'coSegregate' , {} )
		#self.clinvarVariants = kwargs.get( 'clinvarVariants' , {} )
		#self.vepVariants = kwargs.get( 'vepVariants' , [] )
		self.diseases = kwargs.get( 'diseases' , {} )
		self.vcfHeaderInfo = kwargs.get( 'vcfHeaderInfo' , [] )
		self.vcfKeyIndex = kwargs.get( 'vcfKeyIndex' , {} )

		#### ADD key index here to access vcf info for individual variant

### Retrieve input data from user ###
	def getInputData( self  , **kwargs ):
		mafFile = kwargs.get( 'maf' , "" )
		vcfFile = kwargs.get( 'vcf' , "" )
		tsvFile = kwargs.get( 'tsv' , "" )
		pathogenicVariantsFile = kwargs.get( 'pathogenicVariants' , "" )
		expressionFile = kwargs.get( 'expression' , "" )
		geneListFile = kwargs.get( 'geneList' , "" )
		deNovoFile = kwargs.get( 'deNovo' , "" )
		assumedDeNovoFile = kwargs.get( 'assumedDeNovo' , "" )
		coSegregateFile = kwargs.get( 'coSegregate' , "" )
		tcga = kwargs.get( 'tcga' , True )
		diseaseFile = kwargs.get( 'diseases' , "" )
		specific = kwargs.get( 'specific' , False )
		self.getDiseases( diseaseFile , **kwargs )
		preVEP = []
		vepDone = False
		exacDone = False
		if mafFile:
			self.readMAF( mafFile , **kwargs )
		if vcfFile:
			[ vepDone , preVEP , exacDone ] = self.readVCF( vcfFile , appendTo="user" , **kwargs )
			#exacDone=False # currently only has 1000G
		if tsvFile:
			exacDone = self.readTSV( tsvFile , **kwargs )
		if pathogenicVariantsFile:
			self.readVCF( pathogenicVariantsFile , appendTo="pathogenic" , **kwargs )
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
		return [ vepDone , preVEP , exacDone ]
	def readMAF( self , inputFile , **kwargs ):
		inFile = self.safeOpen( inputFile , 'r' )
		codonColumn = kwargs.get( 'codon' , 48 )
		peptideChangeColumn = kwargs.get( 'peptideChange' , 49 )
		alleleFrequencyColumn = kwargs.get( 'alleleFrequency' , None )
		tcga = kwargs.get( 'tcga' , True )
		specific = kwargs.get( 'specific' , True )
		try:
			next(inFile)
			for line in inFile:
				var = chargervariant()
				var.mafLine2Variant( line , peptideChange=peptideChangeColumn , codon=codonColumn )
				#print var.proteogenomicVar()
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
		inFile = None
		if ( re.match( "\.gz" , inputFile ) ):
			inFile = vcf.Reader( open( inputFile , 'r' ) , compressed=True )    
		else:
			inFile = vcf.Reader( open( inputFile , 'r' ) )
		appendTo = kwargs.get( "appendTo" , "user" )
		preVEP = []
		vepDone = False
		exacDone = False
		vepInfo = OD()
		self.vcfHeaderInfo = []
		metadata = inFile.metadata
		#print( str( metadata ) )
		for pairs in metadata:
			if pairs == 'VEP':
				print "This .vcf has VEP annotations!"
				infos = inFile.infos
				for info_ID in infos.items():
					if info_ID[0] == "CSQ": #CSQ tag the VEP annotation, probably means consequence
						csq = info_ID
						Info = csq[1] #Info(...)
						if Info:
							desc = Info[3] #Consequence type...Format: Allele|Gene|...
							keysString = desc.split( "Format: " )[1]
							self.vcfHeaderInfo = keysString.split( "|" )
							self.vcfKeyIndex = {}
							i = 0
							for key in self.vcfHeaderInfo:
								vepInfo[key] = None
								self.vcfKeyIndex[key] = i
								#print str(i) + " => " + key
								i = i + 1
			if pairs == 'AF':
				print "This .vcf has AF!"
				exacDone = True
		for record in inFile:
			chrom = record.CHROM
			reference = record.REF
			alternates = record.ALT
			start = record.start + 1 #1-base beginning of ref
			stop = record.end #0-base ending of ref
			info = record.INFO
			alti = -1
			for alternate in alternates:
				alti += 1
				alt = str( alternate )
				if alt == "None":
					alt = None
				if record.is_indel and not record.is_deletion: #insertion
					reference = "-"
					alt = alt[1:len(alt)]
					stop = stop + 1
				elif record.is_deletion:
					reference = reference[1:len(reference)] #assumes only one base overlap
					alt = "-"
					start = start + 1
					stop = stop

				parentVar = mafvariant( \
					chromosome = chrom , \
					start = start , \
					stop = stop , \
					dbsnp = record.ID , \
					reference = reference , \
					alternate = alt , \
				)

				var = chargervariant( \
					parentVariant=parentVar
				)

				hasAF = False
				var.alleleFrequency = info.get( 'AF' , "noAF" )
				if ( not var.alleleFrequency == "noAF" ):
					afs = var.alleleFrequency[alti]
					var.alleleFrequency = afs
					hasAF = True
				#print( str( var.alleleFrequency ) )

				csq = info.get( 'CSQ' , "noCSQ" )
				if not csq == "noCSQ":
					vepDone = True
					#exacDone = True
					var.vepVariant = vepvariant()
					for thisCSQ in csq:
						values = thisCSQ.split( "|" )
						var.vcfInfo = values
						aas = [None , None] 
						if self.getVCFKeyIndex( values , "Amino_acids" ): #8 => Amino_acids
							aas = self.getVCFKeyIndex( values , "Amino_acids" ).split("/") 
							if len( aas ) > 1:
								aas[0] = mafvariant().convertAA( aas[0] )
								aas[1] = mafvariant().convertAA( aas[1] )
							else:
								#28 => HGVSc
								#29 => HGVSp
								hgvsp = self.getVCFKeyIndex( values , "HGVSp" ).split( ":" )
								changep = None
								if len( hgvsp ) > 1:
									changep = re.match( "p\." , hgvsp[1] )
								if changep:
									aas = mafvariant().splitHGVSp( hgvsp[1] )
									aas[0] = mafvariant().convertAA( aas[0] )
									aas[2] = mafvariant().convertAA( aas[2] )
								else:
									aas.append( None )
									needVEP = True
									preVEP.append( var )
						exons = [None , None]
						if self.getVCFKeyIndex( values , "EXON" ): #25 => EXON
							exons = self.getVCFKeyIndex( values , "EXON" ).split( "/" )
							if len( exons ) == 1:
								exons.append(None)
						introns = [None , None]
						if self.getVCFKeyIndex( values , "INTRON" ): #26 => INTRON
							introns = self.getVCFKeyIndex( values , "INTRON" ).split( "/" )
							if len( introns ) == 1:
								introns.append(None)
						siftStuff = [None , None]
						if self.getVCFKeyIndex( values , "SIFT" ):
							siftStuff = self.getVCFKeyIndex( values , "SIFT" ).split( "(" ) 
							if len( siftStuff ) == 1:
								siftStuff.append( None )
							else:
								siftStuff[1] = siftStuff[1].rstrip( ")" )
						polyPhenStuff = [None , None]
						if self.getVCFKeyIndex( values , "PolyPhen" ):
							polyPhenStuff = self.getVCFKeyIndex( values , "PolyPhen" ).split( "(" ) 
							if len( polyPhenStuff ) == 1:
								polyPhenStuff.append( None )
							else:
								polyPhenStuff[1] = polyPhenStuff[1].rstrip( ")" )
						consequence_terms = self.getVCFKeyIndex( values , "Consequence" )
						csq_terms = []
						if consequence_terms:
							csq_terms = self.getVCFKeyIndex( values , "Consequence" ).split( "&" )
						vcv = vepconsequencevariant( \
							#parentVariant=var
							chromosome = chrom , \
							start = start , \
							stop = stop , \
							dbsnp = record.ID , \
							reference = reference , \
							alternate = alt , \
							#1 => Gene
							gene_id=self.getVCFKeyIndex( values , "Gene" ) , \
							#2 => Feature
							transcriptCodon=self.getVCFKeyIndex( values , "Feature" ) , \
							#4 => Consequence
							consequence_terms=csq_terms , \
							#5 => cDNA_position
							positionCodon=self.getVCFKeyIndex( values , "cDNA_position" ) , \
							#7 => Protein_position
							positionPeptide=self.getVCFKeyIndex( values , "Protein_position" ) , \
							referencePeptide=aas[0] , \
							alternatePeptide=aas[1] , \
							#12 => STRAND
							strand=self.getVCFKeyIndex( values , "STRAND" ) , \
							#13 => SYMBOL
							gene=self.getVCFKeyIndex( values , "SYMBOL" ) , \
							#14 => SYMBOL_SOURCE
							gene_symbol_source=self.getVCFKeyIndex( values , "SYMBOL_SOURCE" ) , \
							#15 => HGNC_ID
							hgnc_id=self.getVCFKeyIndex( values , "HGNC_ID" ) , \
							#16 => BIOTYPE
							biotype=self.getVCFKeyIndex( values , "BIOTYPE" ) , \
							#17 => CANONICAL

							canonical=self.getVCFKeyIndex( values , "CANONICAL" ) , \
							#18 => CCDS
							ccds=self.getVCFKeyIndex( values , "CCDS" ) , \
							#19 => ENSP
							transcriptPeptide=self.getVCFKeyIndex( values , "ENSP" ) , \
							#23 => SIFT
							predictionSIFT=siftStuff[0] , \
							scoreSIFT=siftStuff[1] , \
							#24 => POLYPHEN
							predictionPolyphen=polyPhenStuff[0] , \
							scorePolyphen=polyPhenStuff[1] , \
							exon=exons[0] , \
							totalExons=exons[1] , \
							intron=introns[0] , \
							totalIntrons=introns[1] , \
						)
#6 => CDS_position
#10 => Existing_variation
#11 => DISTANCE
#20 => SWISSPROT
#21 => TREMBL
#22 => UNIPARC
#27 => DOMAINS
#30 => GMAF
						if ( not hasAF ):
							var.alleleFrequency = self.getVCFKeyIndex( values , "GMAF" )
#31 => AFR_MAF
#32 => AMR_MAF
#33 => ASN_MAF
#34 => EUR_MAF
#35 => AA_MAF
#36 => EA_MAF
#37 => CLIN_SIG
#38 => SOMATIC
#39 => PUBMED
#40 => MOTIF_NAME
#41 => MOTIF_POS
#42 => HIGH_INF_POS
#43 => MOTIF_SCORE_CHANGE

						var.vepVariant.consequences.append( vcv )
				severeRank = [ 	"transcript_ablation" , \
								"splice_acceptor_variant" , \
								"splice_donor_variant" , \
								"stop_gained" , \
								"frameshift_variant" , \
								"stop_lost" , \
								"start_lost" , \
								"transcript_amplification" , \
								"inframe_insertion" , \
								"inframe_deletion" , \
								"missense_variant" , \
								"protein_altering_variant" , \
								"splice_region_variant" , \
								"incomplete_terminal_codon_variant" , \
								"stop_retained_variant" , \
								"synonymous_variant" , \
								"coding_sequence_variant" , \
								"mature_miRNA_variant" , \
								"5_prime_UTR_variant" , \
								"3_prime_UTR_variant" , \
								"non_coding_transcript_exon_variant" , \
								"intron_variant" , \
								"NMD_transcript_variant" , \
								"non_coding_transcript_variant" , \
								"upstream_gene_variant" , \
								"downstream_gene_variant" , \
								"TFBS_ablation" , \
								"TFBS_amplification" , \
								"TF_binding_site_variant" , \
								"regulatory_region_ablation" , \
								"regulatory_region_amplification" , \
								"feature_elongation" , \
								"regulatory_region_variant" , \
								"feature_truncation" , \
								"intergenic_variant" ]
				mostSevere = None
				rankMostSevere = 10000
				mostSevereCons = severeRank[-1]
				try:
					for cons in var.vepVariant.consequences:
						for term in cons.terms:
							if term in severeRank:
								rank = severeRank.index( term )
							else:
								rank = 10000
							if rank < rankMostSevere:
								mostSevere = cons
								rankMostSevere = rank
								mostSevereCons = term
							elif rank == rankMostSevere:
								if cons.canonical:
									mostSevere = cons
					if mostSevere:
						var.gene = mostSevere.gene
						var.referencePeptide = mostSevere.referencePeptide
						var.positionPeptide = mostSevere.positionPeptide
						var.alternatePeptide = mostSevere.alternatePeptide
						var.transcriptPeptide = mostSevere.transcriptPeptide
						var.transcriptCodon = mostSevere.transcriptCodon
						var.positionCodon = mostSevere.positionCodon
						var.vepVariant.mostSevereConsequence = mostSevereCons
						var.variantClass = mostSevereCons
				except:
					print( "CharGer::readVCF Warning: no consequences" + var.genomicVar() )
					pass
				#print var.proteogenomicVar()

				if appendTo == "user":
					self.userVariants.append( var )
				elif appendTo == "pathogenic":
					pathKey = self.pathogenicKey( var )
					self.pathogenicVariants[pathKey] = var
		return [ vepDone , preVEP , exacDone ]

	def readTSV( self , inputFile , **kwargs ):
		print "\tReading .tsv!"
		inFile = self.safeOpen( inputFile , 'r' )
		chrColumn = kwargs.get( 'chromosome' , 0 )
		startColumn = kwargs.get( 'start' , 1 )
		stopColumn = kwargs.get( 'stop' , 2 )
		refColumn = kwargs.get( 'ref' , 3 )
		altColumn = kwargs.get( 'alt' , 4 )
		geneColumn = kwargs.get( 'gene' , None )
		strandColumn = kwargs.get( 'strand' , None )
		codonColumn = kwargs.get( 'codon' , None )
		peptideColumn = kwargs.get( 'peptideChange' , None )
		variantClassificationColumn = kwargs.get( 'variantClassification' , None )
		sampleColumn = kwargs.get( 'sample' , None )
		alleleFrequencyColumn = kwargs.get( 'alleleFrequency' , None )
		tcga = kwargs.get( 'tcga' , True )
		specific = kwargs.get( 'specific' , True )
		headerLine = next(inFile).split( "\t" ) 
		exacDone = False
		try:
			for line in inFile:
				var = chargervariant()
				fields = line.split( "\t" )
				chro = self.getCustom( fields , chrColumn , desc="chromosome" )
				var.chromosome = self.getChrNum( chro )
				var.reference = self.getCustom( fields , refColumn , desc="reference" )
				var.alternate = self.getCustom( fields , altColumn , desc="alternate" )
				var.start = self.getCustom( fields , startColumn , desc="start" )
				var.stop = self.getCustom( fields , stopColumn , desc="stop" )
				if geneColumn is not None:
					var.gene = self.getCustom( fields , geneColumn , desc="gene" )
				if sampleColumn is not None:
					var.sample = self.getCustom( fields , sampleColumn , desc="sample" )
				if codonColumn is not None:
					codon = self.getCustom( fields , codonColumn , desc="codonChange" )
					var.splitHGVSc( codon )
				if peptideColumn is not None:
					peptide = self.getCustom( fields , peptideColumn , desc="peptideChange" )
					var.splitHGVSp( peptide )
				if variantClassificationColumn is not None:
					var.variantClass = self.getCustom( fields , variantClassificationColumn , desc="variantClassification" )
				if alleleFrequencyColumn is not None:
					var.alleleFrequency = self.getCustom( fields , alleleFrequencyColumn , desc="alleleFrequency" )
					exacDone = True

				if specific:
					if tcga:
						match = re.match( "TCGA\-(\w\w)" , var.sample )
						if match:
							var.disease = self.diseases[match.groups()[0]]
				else:
					var.disease = charger.allDiseases

				self.userVariants.append( var )
		except:
			raise Exception( "CharGer::readTSV Error: bad .tsv file" )
			print "Hint: correct columns?"
			print headerLine
		return exacDone

	@staticmethod
	def getCustom( array , column , desc="" ):
		thing = None
		try:
			thing = array[int( column )]
		except:
			print( "CharGer::getCustom Error: " + desc + " column at " + str( column ) + \
				   " not available in row from input .tsv" )
			pass
		return thing
			
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
			print "CharGer::readExpression Error: bad expression file"

	def readGeneList( self , inputFile , **kwargs ): # gene list formatted "gene", "disease", "mode of inheritance"
		specific = kwargs.get( 'specific', True )
		try:
			inFile = self.safeOpen( inputFile , 'r' , warning=True )
			for line in inFile:
				fields = line.split( "\t" )
				gene = fields[0].rstrip()
				if specific:
					disease = fields[1].rstrip()
				else: #set the gene to match all disease
					disease = charger.allDiseases
				mode_inheritance = fields[2].rstrip()
				self.userGeneList[gene][disease] = mode_inheritance
				#print '__'.join( [ gene , disease , mode_inheritance ] )
		except:
			print "CharGer::readGeneList Error: bad gene list file"
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
				var = mafvariant()
				var.mafLine2Variant( line )
				varDict[var.uniqueVar()] = 1
		except:
			print "CharGer::readOtherMAF Warning: bad .maf for " + inputFile

### Retrieve external reference data ###
	def getExternalData( self , **kwargs ):
		t = time.time()
		self.getClinVar( **kwargs )
		self.printRunTime( "ClinVar" , self.runTime( t ) )
		t = time.time()
		self.getExAC( **kwargs )
		self.printRunTime( "exac" , self.runTime( t ) )
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
			ent = entrezapi()
			i = 0
			for varsStart in range( 0 , len( self.userVariants ) , int(searchBatchSize) ):
				varsEnd = varsStart + int(searchBatchSize)
				varsSet = self.userVariants[varsStart:varsEnd]
				ent.prepQuery( varsSet )
				ent.subset = entrezapi.esearch
				ent.database = entrezapi.clinvar
				clinvarsSet = ent.doBatch( summaryBatchSize )
				varsBoth = self.matchClinVar( varsSet , clinvarsSet )
				self.userVariants[varsStart:varsEnd] = varsBoth["userVariants"]
				#self.clinvarVariants.update( varsBoth["clinvarVariants"] )
				#self.userVariants[varsStart:varsEnd] = self.matchClinVar( varsSet , clinvarsSet )
	def getExAC( self , **kwargs ):
		doExAC = kwargs.get( 'exac' , True )
		useHarvard = kwargs.get( 'harvard' , True )
		threshold = kwargs.get( 'threshold' , 0 )
		alleleFrequencyColumn = kwargs.get( 'alleleFrequency' , None )
		if doExAC and not alleleFrequencyColumn:
			common = 0
			rare = 0
			totalVars = len( self.userVariants )
			exac = exacapi(harvard=useHarvard)
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
		preVEP = kwargs.get( 'prevep' , [] )
		doREST = kwargs.get( 'rest' , False )
		doCMD = kwargs.get( 'cmd' , False )

		if doVEP:
			if doREST or ( not doREST and not doCMD ):
				self.getVEPviaREST( **kwargs )
			if doCMD:
				self.getVEPviaCMD( **kwargs )
		elif len( preVEP ) > 0:
			if doREST or ( not doREST and not doCMD ):
				self.getVEPviaREST( **kwargs )
			if doCMD:
				self.getVEPviaCMD( **kwargs )

	def getVEPviaREST( self , **kwargs ):
		doAllOptions = kwargs.get( 'allOptions' , True )
		doVEP = kwargs.get( 'vep' , [] )
		preVEP = kwargs.get( 'prevep' , [] )
		vep = ensemblapi()
		luv = len(self.userVariants)
		vepVariants = []
		if doVEP:
			vepVariants = vep.annotateVariantsPost( self.userVariants , **kwargs )
		elif len( preVEP ) > 0:
			vepVariants = vep.annotateVariantsPost( preVEP , **kwargs )
		self.matchVEP( vepVariants )
		aluv = 0
		for var in self.userVariants:
			if var.vepVariant:
				aluv += 1
		luvafter = len(self.userVariants)
		print "\nVEP annotated userVariants " + str( luvafter )
		print "VEP annotated " + str(aluv) + " from the original set of " + str(luv)

	def getVEPviaCMD( self , **kwargs ):
		NotImplemented
		#defaultVEPDir = "./"
		#vepDir = kwargs.get( 'vepDir' , defaultVEPDir )
		#defaultVEPOutputDir = '/'.join( [ vepDir , ".vep" ] ) + "/"
		#defaultVEPCache = '/'.join( [ vepDir , ".vep" ] ) + "/"
		#defaultEnsemblVersion = "75"
		#defaultVEPVersion = "80"
		#defaultAssembly = "GRCh37"
		#defaultFasta = '/'.join( [ vepDir , \
		#	".vep" , \
		#	"homo_sapiens" , \
		#	"*" + assembly , \
		#	"Homo_sapiens." + assembly + "." + vepVersion + ".dna.primary_assembly.fa" \
		#] )
		#assembly = kwargs.get( 'vepAssembly' , defaultAssembly )
		#fasta = kwargs.get( 'vepFasta' , defaultFasta )
		#vcfFile = kwargs.get( 'vcf' , "" )
		#outputFile = kwargs.get( 'outputVCF' , "" )
		#if vcfFile:
		#	vep_command = ' '.join( [ "perl" , vepDir + , \
		#		"--everything" , \
		#		"--species homo_sapiens" , \
		#		"--assembly" , assembly , \
		#		"--fasta" , fasta , \
		#		"--input_file" , vcfFile , \
		#		"--output_file" , vcfTemp , \
		#		"--dir" , vepDir + "/.vep" , \
		#		"--dir_cache" , vepDir + "/.vep" , \
		#		"--cache" , \
		#		"--no_progress" , \
		#		"--quiet" , \
		#		"--total_length" , \
		#		"--no_escape" , \
		#		"--xref_refseq" , \
		#		"--force_overwrite" , \
		#		"--format vcf" , \
		#		"--vcf" , \
		#		"--no_stats" , \
		#		"--fork 4" , \
		#	] )
		#	os.system(vep_command)

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
			vepVar = vepvariant()
			if genVar in vepVariants:
				vepVar = vepVariants[genVar]
			if var.sameGenomicVariant( vepVar ):
				var.vepVariant = vepVar
				var.copyMostSevereConsequence()
			if var.vepVariant:
				var.vepAnnotations = var.vepVariant.printVariant( delim="|" , minimal=True )
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
			print "CharGer::getDiseases Warning: No diseases file provided"
			return


### Evidence levels ### 

# any level ending with C (ex. PSC1, PPC1) are CharGer-invented modules #

##### Very Strong #####
	def PVS1( self , expressionThreshold = 0.2 ):
		print "CharGer module PVS1"
		print "- truncations in genes where LOF is a known mechanism of the disease"
		print "- require the mode of inheritance to be dominant (assuming heterzygosity) and co-occurence with reduced gene expression"
		print "- run concurrently with PSC1, PM4, and PPC1 - "
		self.runIndelModules()
		

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
				var.addSummary( "PS2(de novo with parent confirmation and no history)" )
	def PS3( self ):
		print "CharGer module PS3: Well-established in vitro or in vivo functional studies \
		supportive of a damaging effect on the gene or gene product"
#		print "- "
	def PS4( self ): # not relevant in rare variants, such big effect size common variants are so rare may as well just take a input variant list
		print "CharGer module PS4: not yet implemented"
#		print "- variant prevalence in cases significantly greater than controls"
		for var in self.userVariants:
			return
			caseVarFreq = "NEED UPDATE" # may take input from current MAF
			controlVarFreq = "NEED UPDATE" # may take input from exac
			if ( caseVarFreq != 0 and controlVarFreq != 0):
				OR = (caseVarFreq/controlVarFreq) / ( (1-caseVarFreq)/(1-controlVarFreq) )
				# Adam will update
				if OR >= 5:
					CIlower = math.log(OR) - math.sqrt( 1/caseVarFreq + 1/controlVarFreq + 1/caseVarFreq + 1/controlVarFreq)
					if (CIlower > 1):
						var.PS4 = True
						var.addSummary( "PS4(Prevalence significantly greater than controls)" )

	def PSC1( self ):
		print "CharGer module PSC1"
		print "Truncations of other, not-specificied genes"

##### Moderate #####
	def PM1( self , recurrenceThreshold , **kwargs ):
		print "CharGer module PM1:  Located in a mutational hot spot and/or critical and well-established \
		functional domain (e.g., active site of an enzyme) without benign variation"
		clustersFile = kwargs.get( 'hotspot3d' , "" )
		if clustersFile:
			print "Reading HotSpot3D clusters file: " + clustersFile
			with open( clustersFile , 'r' ) as clustersFH:
				clustersFH.next()
				for line in clustersFH:
					line.strip()
					fields = line.split('\t')
					recurrence = fields[6]
					if recurrence >= recurrenceThreshold:
						chromosome = fields[9]
						start = fields[10]
						stop = fields[11]
						reference = fields[12]
						alternate = fields[13]
						hotspot3dVar = mafvariant( 
							chromosome = chromosome , 
							start = start , 
							stop = stop , 
							reference = reference , 
							alternate = alternate 
							)
						for var in self.userVariants:
							if hotspot3dVar.sameGenomicVariant( var ):
								var.PM1 = True
								var.addSummary( "PM1(HotSpot3D: somatic hotspot among 19 TCGA cancer types with " + str( recurrence ) + "samples)" )
		else:
			print "CharGer::PM1 Warning: clustersFile is not supplied. PM1 was not executed. "
							#break #maybe don't break if possible redundant genomic variants
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
		print "- protein length changes due to inframe indels or nonstop variant of selected genes - "
		
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
				var.addSummary( "PM6(Assumed de novo without parent confirmation)" )

##### Supporting #####
	def PP1( self ):
		print "CharGer module PP1"
		print "- cosegregation with disease in family members in a known disease gene"
		for var in self.userVariants:
			if var.uniqueVar() in self.userCoSegregateVariants:
				var.PP1 = True
				var.addSummary( "PP1(Cosegregation with disease in family from known disease gene " + self.gene + ")" )
	def PP2( self ):
		print "CharGer module PP2: Missense variant in a gene of the given gene list"
		if self.userGeneList: #gene, disease, mode of inheritance
			for var in self.userVariants:
				varGene = var.gene
				varDisease = var.disease # no disease field in MAF; may require user input	
				varSample = var.sample
				varClass = var.variantClass
				varVEPClass = "blahblah"
				if var.vepVariant:
					varVEPClass = var.vepVariant.mostSevereConsequence
				if ( "missense" in varClass.lower() ) or \
					( "missense" in varVEPClass.lower() ):
					if varGene in self.userGeneList: # check if in gene list
						var.PP2 = True # if call is true then check expression effect
						var.addSummary( "PP2(Missense variant in gene from gene list)" )
		else: 
			print "CharGer::PP2 Error: Cannot evaluate PP2: No gene list supplied."
#		print "- "
	def PP3( self , minimumEvidence ):
		print "CharGer module PP3"
		print "- multiple lines of in silico evidence of deliterous effect"
		callSIFTdam = "damaging"
		callSIFTdel = "deleterious"
		thresholdSIFT = 0.05
		callPolyphen = "probably damaging"
		thresholdPolyphen = 0.432
		callBlosum62 = -2
		callCompara = 2
		callImpact = "high"
		fracMaxEntScan = 0.8
		callGeneSplicer = ""
		for var in self.userVariants:
			case = []
			if var.vepVariant:
				if var.vepVariant.consequences:
					for vcVar in var.vepVariant.consequences:
						if not var.PP3:
							evidence = 0
							if vcVar.blosum:
								if vcVar.blosum < callBlosum62:
									case.append( "Blosum62:" \
										+ str( vcVar.blosum ) \
										+ "<" + str( callBlosum62 ) )
									evidence += 1
							if vcVar.scoreSIFT and (vcVar.scoreSIFT < thresholdSIFT):
								case.append( "SIFT:" \
									+ str( vcVar.scoreSIFT ) \
									+ "<" + str( thresholdSIFT ) )
								evidence += 1
							# elif vcVar.predictionSIFT:
							# 	if vcVar.predictionSIFT.lower() == callSIFTdam and \
							# 	vcVar.predictionSIFT.lower() == callSIFTdel:
							# 		case.append( "SIFT:" \
							# 			+ str( vcVar.predictionSIFT ) )
							# 		evidence += 1
							if vcVar.scorePolyphen and (vcVar.scorePolyphen > thresholdPolyphen):
								case.append( "PolyPhen:" \
									+ str( vcVar.scorePolyphen ) \
									+ ">" + str( thresholdPolyphen ) )
								evidence += 1
							# elif vcVar.predictionPolyphen:
							# 	if vcVar.predictionPolyphen.lower().replace( "_" , " " ) == callPolyphen:
							# 		case.append( "PolyPhen:" \
							# 			+ str( vcVar.predictionPolyphen ) )
							# 		evidence += 1
							if vcVar.compara:
								if vcVar.compara > callCompara:
									case.append( "Compara:" \
										+ str( vcVar.compara ) \
										+ ">" + str( callCompara ) )
									evidence += 1
							if vcVar.impact:
								if vcVar.impact.lower() == callImpact:
									case.append( "VEP_Impact:" \
										+ str( vcVar.impact ) )
									evidence += 1
							if vcVar.maxentscan:
								callMaxEntScan = vcVar.maxentscan[0]*fracMaxEntScan
								if vcVar.maxentscan[1] <= callMaxEntScan:
									case.append( "MaxEntScan:" \
										+ str( vcVar.maxentscan[1] ) \
										+ "<=" + str( callMaxEntScan ) )
									evidence += 1
							if vcVar.genesplicer:
								if vcVar.genesplicer.lower() == callGeneSplicer:
									case.append( "GeneSplicer:" \
										+ str( vcVar.genesplicer ) )
									evidence += 1
							if evidence >= minimumEvidence:
								var.PP3 = True
				if var.PP3:
					var.addSummary( "PP3(Multiple (>=" + str( minimumEvidence ) \
						+ ") in silico predictions of deliterious effect=" \
						+ "|".join( case ) + ")" )

	def PP4( self ):
		print "CharGer module PP4: not yet implemented"
#		print "- "
	def PP5( self ):
		print "CharGer module PP5: not yet implemented"
#		print "- "
	def PPC1( self ):
		print "CharGer module PPC1"
		print "- protein length changes due to inframe indels or nonstop variant of other, not-specificied genes - "

### helper functions of evidence levels ###
	def runIndelModules( self ):
		maf_truncations = ["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Splice_Site"] #,"Nonstop_Mutation"
		vep_truncations = ["transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained",\
							"frameshift_variant","start_lost"]
		lenShift = ["In_Frame_Del","In_Frame_Ins","Nonstop_Mutation"]
		vep_inframe = [	"inframe_insertion" , "inframe_deletion" , \
						"stop_lost" ]

		if not self.userGeneList:
			print "CharGer::runIndelModules Error: Cannot evaluate PVS1 or PM4: No gene list supplied."
		for var in self.userVariants:
			varClass = var.variantClass
			varGene = var.gene
			#print varClass
			varVEPClass = ""
			if var.vepVariant:
				varVEPClass = var.vepVariant.mostSevereConsequence
			altPeptide = var.alternatePeptide

			if (varClass in maf_truncations) or \
					(varVEPClass in vep_truncations) or \
					altPeptide == "*" or \
					altPeptide == "fs":
				if self.userGeneList:
					if varGene in self.userGeneList: # check if in gene list
						varDisease = var.disease # no disease field in MAF; may require user input	
						if ( "dominant" in self.userGeneList[varGene][varDisease].lower() or \
							"dominant" in self.userGeneList[varGene][charger.allDiseases].lower() ):
							var.PVS1 = True # if call is true then check expression effect
							if self.userExpression: # consider expression data only if the user has supplied an expression matrix
								varSample = var.sample
								if self.userExpression[varSample][varGene] >= expressionThreshold:
									var.PVS1 = False
					else:
						var.PSC1 = True 
						var.addSummary( "PSC1(truncation)" )
				else: 
					# in case of no gene list, make all truncations PSC1
					var.PSC1 = True
			elif (varClass in lenShift \
				or varClass in vep_inframe):
				if self.userGeneList:
					if varGene in self.userGeneList: # check if in gene list
						varDisease = var.disease # no disease field in MAF; may require user input	
						if ( "dominant" in self.userGeneList[varGene][varDisease].lower() or \
							"dominant" in self.userGeneList[varGene][charger.allDiseases].lower() ):
							var.PM4 = True
					else: 
						var.PPC1 = True
				else: 
					# in case of no gene list, make all inframes PSC1
					var.PPC1 = True

			if var.PVS1:
				var.addSummary( "PVS1(" + str( varClass ) + " in susceptible gene " + str( varGene ) + ")" )
			if var.PSC1:
				var.addSummary( "PSC1(" + str( varClass ) + " in gene " + str( varGene ) + ")" )
			if var.PM4:
				var.addSummary( "PM4(" + str( varClass ) + " in susceptible gene " + str( varGene ) + ")" )
			if var.PPC1:
				var.addSummary( "PPC1(" + str( varClass ) + " in gene " + str( varGene ) + ")" )

	def peptideChange( self , mod , **kwargs ):
		called = 0
		for var in self.userVariants:
			call = False
			if mod == "PS1":
				call = var.PS1
			if mod == "PM5":
				call = var.PM5
			if not call: #is already true
				CVchecked = 0
				PVchecked = 0
				if var.clinvarVariant or self.pathogenicVariants:
					if var.vepVariant:
						if var.vepVariant.consequences:
							for consequence in var.vepVariant.consequences:
								CVchecked = self.checkClinVarPC( var , mod , consequence )
								#for pathVar in self.pathogenicVariants:
								PVchecked = self.checkPathogenicVariants( var , mod , consequence )
								if CVchecked or PVchecked:
									called += 1
					else:
						CVchecked = self.checkClinVarPC( var , mode , consequence )
						PVchecked = self.checkPathogenicVariants( var , mod , consequence )
						if CVchecked or PVchecked:
							called += 1
			if var.PS1 and mod == "PS1":
				var.addSummary( "PS1(Peptide change is known pathogenic)" )
			if var.PM5 and mod == "PM5":
				var.addSummary( "PM5(Peptide change at the same location of a known pathogenic change)" )
		print mod + " found " + str(called) + " pathogenic variants"
	def checkClinVarPC( self , var , mod , consequence ):
		called = 0
		if var.clinvarVariant:
			clinvarVar = var.clinvarVariant
			clin = clinvarVar.clinical
			if consequence.sameGenomicVariant( clinvarVar ):
			#if genomic change is the same, then PS1
				if clin["description"] == clinvarvariant.pathogenic:
					if mod == "PS1":
						var.PS1 = True # already pathogenic still suffices to be PS1
						called = 1
			elif consequence.sameGenomicReference( clinvarVar ):
			#if genomic change is different, but the peptide change is the same, then PS1
				if clinvarVar.alternatePeptide == consequence.alternatePeptide: #same amino acid change
					if clin["description"] == clinvarvariant.pathogenic:
						if mod == "PS1":
							var.PS1 = True
							called = 1
			if consequence.samePeptideReference( clinvarVar ):
				if not consequence.samePeptideChange( clinvarVar ):
				#if peptide change is different, but the peptide reference is the same, then PM5
					if consequence.plausibleCodonFrame( clinvarVar ):
						if clin["description"] == clinvarvariant.pathogenic:
							if mod == "PM5":
								var.PM5 = True # already pathogenic still suffices to be PS1
								called = 1
		return called
	def checkPathogenicVariants( self , var , mod , consequence ):
		called = 0
		if self.pathogenicVariants:
			pathKey = self.pathogenicKey( consequence )
			if pathKey in self.pathogenicVariants:
				if self.pathogenicVariants[pathKey].vepVariant:
					if self.pathogenicVariants[pathKey].vepVariant.consequences:
						for pathVar in self.pathogenicVariants[pathKey].vepVariant.consequences:
							if consequence.sameGenomicVariant( pathVar ):
							#if genomic change is the same, then PS1
								if mod == "PS1":
									var.PS1 = True # already pathogenic still suffices to be PS1
									called = 1
									#print var.genomicVar() ,
									#print " matched a pathogenic Variant in PS1 by same genomic change"
							elif consequence.sameGenomicReference( pathVar ):
							#if genomic change is different, but the peptide change is the same, then PS1
								if pathVar.alternatePeptide == consequence.alternatePeptide: #same amino acid change
									if mod == "PS1":
										var.PS1 = True
										called = 1
										#print var.genomicVar() ,
										#print " matched a pathogenic Variant in PS1 by same peptide change"
							if consequence.samePeptideReference( pathVar ):
								if not consequence.samePeptideChange( pathVar ):
								#if peptide change is different, but the peptide reference is the same, then PM5
									if consequence.plausibleCodonFrame( pathVar ):
										if mod == "PM5":
											var.PM5 = True # already pathogenic still suffices to be PS1
											called = 1
											#print var.genomicVar() ,
											#print " matched a pathogenic Variant in PM5"
		return called
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
				#print( "PM2: Checking if " + str( var.alleleFrequency ) + " <= " + str( threshold ) + str( var.isFrequentAllele( threshold ) ) )
				if not var.isFrequentAllele( threshold ):
					var.PM2 = True
					var.addSummary( "PM2(Low/absent allele frequency " \
						+ str( var.alleleFrequency ) \
						+ "<=" + str( threshold ) + ")" )
			if mod == "BA1":
				#print( "BA1: Checking if " + str( var.alleleFrequency ) + " > " + str( threshold ) + str( var.isFrequentAllele( threshold ) ) )
				if var.isFrequentAllele( threshold ):
					var.BA1 = True
					var.addSummary( "BA1(Frequent allele " \
						+ str( var.alleleFrequency ) \
						+ ">" + str( threshold ) + ")" )

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
		thresholdSIFT = 0.05
		callPolyphen = "probably damaging"
		thresholdPolyphen = 0.432
		callBlosum62 = -2
		callCompara = 2
		callImpact = "high"
		fracMaxEntScan = 0.8
		callGeneSplicer = ""
		for var in self.userVariants:
			case = []
			if var.vepVariant:
				if var.vepVariant.consequences:
					for vcVar in var.vepVariant.consequences:
						if not var.BP4:
							evidence = 0
							if vcVar.blosum:
								if vcVar.blosum > callBlosum62:
									case.append( "Blosum62:" \
										+ str( vcVar.blosum ) \
										+ ">" + str( callBlosum62 ) )
									evidence += 1
							if vcVar.scoreSIFT and (vcVar.scoreSIFT >= thresholdSIFT):
								case.append( "SIFT:" \
									+ str( vcVar.scoreSIFT ) \
									+ ">=" + str( thresholdSIFT ) )
								evidence += 1
							# elif vcVar.predictionSIFT:
							# 	if vcVar.predictionSIFT.lower() != callSIFTdam and \
							# 	vcVar.predictionSIFT.lower() != callSIFTdel:
							# 		case.append( "SIFT:" \
							# 			+ str( vcVar.predictionSIFT ) )
							# 		evidence += 1
							if vcVar.scorePolyphen and (vcVar.scorePolyphen <= thresholdPolyphen):
								case.append( "PolyPhen:" \
									+ str( vcVar.scorePolyphen ) \
									+ "<=" + str( thresholdPolyphen ) )
								evidence += 1
							# elif vcVar.predictionPolyphen:
							# 	if vcVar.predictionPolyphen.lower().replace( "_" , " " ) != callPolyphen:
							# 		case.append( "PolyPhen:" \
							# 			+ str( vcVar.predictionPolyphen ) )
							# 		evidence += 1
							if vcVar.compara:
								if vcVar.compara <= callCompara:
									case.append( "Compara:" \
										+ str( vcVar.compara ) )
									evidence += 1
							if vcVar.impact:
								if vcVar.impact.lower() != callImpact:
									case.append( "VEP_Impact:" \
										+ str( vcVar.impact ) )
									evidence += 1
							if vcVar.maxentscan:
								callMaxEntScan = vcVar.maxentscan[0]*fracMaxEntScan
								if vcVar.maxentscan[1] > callMaxEntScan:
									case.append( "MaxEntScan:" \
										+ str( vcVar.maxentscan[1] ) \
										+ ">" + str( callMaxEntScan ) )
									evidence += 1
							if vcVar.genesplicer:
								if vcVar.genesplicer.lower() != callGeneSplicer:
									case.append( "GeneSplicer:" \
										+ str( vcVar.genesplicer ) )
									evidence += 1
							if evidence >= minimumEvidence:
								var.BP4 = True
				if var.BP4:
					var.addSummary( "BP4(Multiple (>=" + str( minimumEvidence ) \
						+ ") in silico predictions of non-deleterious effect=" \
						+ "|".join( case ) + ")" )

	def BP5( self ):
		print "CharGer module BP5: not yet implemented"
		#print" - multiple lines of computational evidence suggesting no impact on gene or product"
	def BP6( self ):
		print "CharGer module BP6: not yet implemented"
	def BP7( self ):
		print "CharGer module BP7: not yet implemented"

### Classifier ###
	def classify( self , **kwargs ):
		scoreSystem = kwargs.get( "system" , "CharGer" )
		for var in self.userVariants:
			var.isPathogenic( **kwargs )
			var.isLikelyPathogenic( **kwargs )
			var.isLikelyBenign( **kwargs )
			var.isBenign( **kwargs )
			var.isUncertainSignificance( **kwargs )
			var.tallyScore( **kwargs )

	def printClassifications( self , **kwargs ):
		scoreSystem = kwargs.get( "system" , "CharGer" )
		headLine = '\t'.join( ["Variant" , "PositiveEvidence" , \
			"CharGerClassification" , "ClinVarAnnoation"] )
		print headLine
		i = 0
		for var in self.userVariants:
			i += 1
			print '\t'.join( [ str(i) , var.uniqueProteogenomicVar() , \
				var.positiveEvidence() , var.pathogenicity[scoreSystem] , \
				var.clinical["description"] ] )
	def writeSummary( self , outFile , **kwargs ):
		delim = kwargs.get( 'delim' , '\t' )
		asHTML = kwargs.get( 'html' , False )
		print( "write to " + outFile )
		outFH = self.safeOpen( outFile , 'w' , warning=True )
		if asHTML:
			print "Writing to .html with delim = </td><td>"
			delim = "</td><td>"
			outFH.write( "<html><head></head><body><table style=\"width: 100%;\"><tr><td>" )
		headLine = delim.join( ["HUGO_Symbol" , "Chromosome" , "Start" , \
			"Stop" , "Reference" , "Alternate" , \
			# "Strand" , "Assembly" ,"Variant_Type" , \
			#"SIFT" , "PolyPhen", \
			"Variant_Classification" , \
			#"Sample" , \
			"HGVSg", "HGVSc", "HGVSp" , \
			#"Sample" , "Transcript" , "Codon_Position" , "HGVSg", "HGVSc", "Protein" , \
			#"Peptide_Reference" , "Peptide_Position" , "Peptide_Alternate" , \
			#"HGVSp","Allele_Frequency","VEP_Most_Severe_Consequence" , "ClinVar_Pathogenicity" , \
			"Allele_Frequency","VEP_Most_Severe_Consequence" ,  \
			"Positive_Evidence" , "Negative_Evidence" , \
			"Positive_CharGer_Score" , "Negative_CharGer_Score" , "CharGer_Score", "ClinVar_Pathogenicity" , \
			"ACMG_Classification" , "CharGer_Classification" , \
			"PubMed_Link" , "ClinVar_Traits" , \
			"CharGer_Summary"] )
			# "VEP_Annotations" , \
			# "VCF_Headers" , "VCF_INFO" , "CharGer_Summary"] )
		try:
			outFH.write( headLine )
			if asHTML:
				outFH.write( "</td></tr>" )
			else:
				outFH.write( "\n" )
			for var in self.userVariants:
				fields = []

				self.appendStr( fields,var.gene)
				self.appendStr( fields,var.chromosome)
				self.appendStr( fields,var.start)
				self.appendStr( fields,var.stop)
				self.appendStr( fields,var.reference)
				self.appendStr( fields,var.alternate)
				# self.appendStr( fields,var.strand)
				# self.appendStr( fields,var.assembly)
				# self.appendStr( fields,var.variantType)
				#elf.appendStr( fields, self.vcfKeyIndex["SIFT"])
				#self.appendStr( fields, var.vcfInfo)
				# if var.vcfInfo:
				# 	if "SIFT" in self.vcfKeyIndex:
				# 		self.appendStr( fields,var.vcfInfo[self.vcfKeyIndex["SIFT"]])
				# 	else:
				# 		self.appendStr( fields , "NA" )
				# 	if "PolyPhen" in self.vcfKeyIndex:
				# 		self.appendStr( fields,var.vcfInfo[self.vcfKeyIndex["PolyPhen"]])
				# 	else:
				# 		self.appendStr( fields , "NA" )
				# else:
				# 	self.appendStr( fields , "NA" )
				# 	self.appendStr( fields , "NA" )
				#self.appendStr( fields,var.polyphen)
				self.appendStr( fields,var.variantClass)
				#self.appendStr( fields,var.sample)
				#self.appendStr( fields,var.transcriptCodon)
				#self.appendStr( fields,var.positionCodon)
				self.appendStr( fields,var.HGVSg() ) # need to be corrected to cDNA change on the most severe peptide
				self.appendStr( fields,var.HGVSc() ) # need to be corrected to cDNA change on the most severe peptide
				#self.appendStr( fields,var.transcriptPeptide)
				#self.appendStr( fields,var.referencePeptide)
				#self.appendStr( fields,var.positionPeptide)
				#self.appendStr( fields,var.alternatePeptide)
				self.appendStr( fields, var.HGVSp() ) # need to be corrected too
				self.appendStr( fields,var.alleleFrequency)
				#self.appendStr( fields,var.vepVariant.mostSevereConsequence) ## this line will fail you on insertions regardless of all the checks in appendStr
				try:
					self.appendStr( fields,var.vepVariant.mostSevereConsequence)
				except:
					self.appendStr( fields , "NA" )
					pass
				#self.appendStr( fields,var.clinvarVariant.trait)
				self.appendStr( fields,var.positiveEvidence())
				self.appendStr( fields,var.negativeEvidence())
				self.appendStr( fields,var.pathogenicScore)
				self.appendStr( fields,var.benignScore)
				self.appendStr( fields,var.chargerScore)
				self.appendStr( fields,var.clinical["description"])
				self.appendStr( fields,var.pathogenicity["ACMG"])
				self.appendStr( fields,var.pathogenicity["CharGer"])
				try:
					if asHTML:
						text = "<a href=\""
						text += var.clinvarVariant.linkPubMed()
						text += "\">uid="
						text += str( var.clinvarVariant.uid )
						text += "</a>"
						self.appendStr( fields, text )
					else:
						self.appendStr( fields,var.clinvarVariant.linkPubMed())
				except:
					self.appendStr( fields , "NA" )
					pass
				try:
					self.appendStr( fields,var.clinvarVariant.getTraits( '|' ) )
				except:
					self.appendStr( fields , "NA" )
					pass
				# TODO: add all the bioinformatic good stuff from VEP
				# self.appendStr( fields , var.vepAnnotations ) #make sure this works
				# self.appendStr( fields , self.vcfHeaderInfo )
				# self.appendStr( fields , var.vcfInfo )
				self.appendStr( fields , ' -- '.join( var.callSummary ) )

				if asHTML:
					outFH.write( "<tr><td>" )
				outFH.write( delim.join( fields ) )
				if asHTML:
					outFH.write( "</td></tr>" )
				else:
					outFH.write( "\n" )
		except:
			print "CharGer Error: Could not write output summary"
			pass
		if asHTML:
			outFH.write( "</table></body></html>" )

	def getVCFKeyIndex( self , values , field ):
		if field in self.vcfKeyIndex:
			return values[self.vcfKeyIndex[field]]
		return None

	@staticmethod
	def appendStr( array, value ):
		try:
			array.append( str( value ) )
		except:
			print "CharGer Warning: failed to append a value\n"
			array.append( "NA" )
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
	@staticmethod
	def getChrNum( chrString ):
		''' Get the chromosome number in case chr or Chr is present'''
		return chrString.upper().replace("CHR", "")
	@staticmethod
	def pathogenicKey( var ):
		return str( var.gene ) + str( var.referencePeptide ) + str( var.positionPeptide )
