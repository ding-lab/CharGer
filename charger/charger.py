#!/usr/bin/env python
# CharGer - Characterization of Germline variants
# author: Adam D Scott (ascott@genome.wustl.edu) & Kuan-lin Huang (khuang@genome.wustl.edu)
# version: v0.0 - 2015*12

import os
import sys
import subprocess
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
from biomine.parsers.exacparser import exacparser
from chargervariant import chargervariant
from autovivification import autovivification as AV
import vcf
from collections import OrderedDict as OD
import gzip
import json
import pdb

class charger(object):
	''' Example usage:
			CharGer = charger()
			CharGer.getInputData( maf=mafFile , expression=expressionFile , geneList=geneListFile )
			CharGer.getExternalData( clinvar=doClinVar , exac=doExAC )
			CharGer.PVS1( )
			CharGer.PS1( )
			CharGer.PM4( )
			CharGer.PM5( )
			CharGer.characterize( )
	'''
	allDiseases = "all"
	DOMINANT = "dominant"
	RECESSIVE = "recessive"
	def __init__( self , **kwargs ):
		self.userVariants = kwargs.get( 'variants' , [] )   # (rjm) list of variants of class chargervariant
		self.pathogenicVariants = kwargs.get( 'pathogenic' , AV({}) )
		self.userExpression = kwargs.get( 'expressions' , AV({}) )
		self.inheritanceGeneList = kwargs.get( 'inheritanceGeneList' , AV({}) )
		self.PP2GeneList = kwargs.get( 'PP2List' , set() )
		self.BP1GeneList = kwargs.get( 'BP1List' , set() )
		self.userDeNovoVariants = kwargs.get( 'deNovo' , {} )
		self.userAssumedDeNovoVariants = kwargs.get( 'assumedDeNovo' , {} )
		self.userCoSegregateVariants = kwargs.get( 'coSegregate' , {} )
		#self.clinvarVariants = kwargs.get( 'clinvarVariants' , {} )
		#self.vepVariants = kwargs.get( 'vepVariants' , [] )
		self.diseases = kwargs.get( 'diseases' , {} )
		self.vcfHeaderInfo = kwargs.get( 'vcfHeaderInfo' , [] )
		self.vcfKeyIndex = kwargs.get( 'vcfKeyIndex' , {} )
		self.mutationTypesFilter = kwargs.get( 'mutationTypes' , [] )
		self.thresholdAF = kwargs.get( 'thresholdAF' , 0.0005 )
		self.keepAF = kwargs.get( 'keepAF' , 1 )
		#self.filtered = []

		#### ADD key index here to access vcf info for individual variant

	def testLocalInputs( self , **kwargs ):
		mafFile = kwargs.get( 'maf' , "" )
		vcfFile = kwargs.get( 'vcf' , "" )
		tsvFile = kwargs.get( 'tsv' , "" )
		if not os.path.exists( mafFile ) and \
		   not os.path.exists( vcfFile ) and \
		   not os.path.exists( tsvFile ):
			sys.exit()

		self.checkInputExistence( 'pathogenicVariants' , **kwargs )
		self.checkInputExistence( 'expression' , **kwargs )
		self.checkInputExistence( 'inheritanceGeneList' , **kwargs )
		self.checkInputExistence( 'PP2GeneList' , **kwargs )
		self.checkInputExistence( 'BP1GeneList' , **kwargs )
		self.checkInputExistence( 'deNovo' , **kwargs )
		self.checkInputExistence( 'assumedDeNovo' , **kwargs )
		self.checkInputExistence( 'coSegregate' , **kwargs )
		self.checkInputExistence( 'diseases' , **kwargs )

		self.checkInputExistence( 'macClinVarTSV' , **kwargs )
		self.checkInputExistence( 'macClinVarVCF' , **kwargs )
		self.checkInputExistence( 'exacVCF' , **kwargs )

		self.checkInputExistence( 'perl' , **kwargs )
		self.checkInputExistence( 'vepScript' , **kwargs )
		self.checkInputExistence( 'vepConfig' , **kwargs )
		#self.checkInputExistence( 'vepDir' , **kwargs )
		self.checkInputExistence( 'vepCache' , **kwargs )
		self.checkInputExistence( 'referenceFasta' , **kwargs )
		vepOutputFile = kwargs.get( 'vepOutput' , "" )
		if vepOutputFile:
			vepOutHandle = self.safeOpen( vepOutputFile , "w" )
			self.checkInputExistence( 'vepOutput' , **kwargs )
			vepOutHandle.close()

		outputFile = kwargs.get( 'output' , "" )
		if outputFile:
			outHandle = self.safeOpen( outputFile , "w" )
			self.checkInputExistence( 'output' , **kwargs )
			outHandle.close()

	def checkInputExistence( self , key , **kwargs ):
		keyFile = kwargs.get( key , "" )
		if keyFile :
			if not os.path.exists( keyFile ):
				print( "CharGer files error: this file does not exist: " + keyFile )
				sys.exit()

### Retrieve input data from user ###
	def getInputData( self  , **kwargs ):
		mafFile = kwargs.get( 'maf' , "" )
		vcfFile = kwargs.get( 'vcf' , "" )
		tsvFile = kwargs.get( 'tsv' , "" )
		pathogenicVariantsFile = kwargs.get( 'pathogenicVariants' , "" )
		expressionFile = kwargs.get( 'expression' , "" )
		inheritanceGeneListFile = kwargs.get( 'inheritanceGeneList' , "" )
		PP2GeneListFile = kwargs.get( 'PP2GeneList' , "" )
		BP1GeneListFile = kwargs.get( 'BP1GeneList' , "" )
		deNovoFile = kwargs.get( 'deNovo' , "" )
		assumedDeNovoFile = kwargs.get( 'assumedDeNovo' , "" )
		coSegregateFile = kwargs.get( 'coSegregate' , "" )
		diseaseFile = kwargs.get( 'diseases' , "" )
		specific = kwargs.get( 'specific' , False )
		self.getDiseases( diseaseFile , **kwargs )
		preVEP = []
		vepDone = False
		exacDone = False
		clinvarDone = False
		if mafFile:
			self.readMAF( mafFile , **kwargs )
		if vcfFile:
			[ vepDone , preVEP , exacDone , clinvarDone , _ ] = self.readVCF( \
				vcfFile , appendTo="user" , **kwargs )
			#exacDone=False # currently only has 1000G
		if tsvFile:
			exacDone = self.readTSV( tsvFile , **kwargs )
		if inheritanceGeneListFile:
			self.readModesGeneList( inheritanceGeneListFile , specific=specific )
		else:
			print "No gene list file uploaded. CharGer will not make PVS1 calls."
		if PP2GeneListFile:
			self.readPP2GeneList( PP2GeneListFile )
		else:
			print "No PP2 gene list file uploaded. CharGer will not make PP2 calls."
		if BP1GeneListFile:
			self.readBP1GeneList( BP1GeneListFile )
		else:
			print "No BP1 gene list file uploaded. CharGer will not make BP1 calls."
		if pathogenicVariantsFile:
			self.readVCF( pathogenicVariantsFile , appendTo="pathogenic" , **kwargs )
		if deNovoFile:
			self.readDeNovo( deNovoFile )
		if assumedDeNovoFile:
			self.readAssumedDeNovo( assumedDeNovoFile )
		if coSegregateFile:
			self.readCoSegregate( coSegregateFile )
		if expressionFile:
			self.readExpression( expressionFile )
		else: 
			print "No expression file uploaded. CharGer will allow all passed truncations without expression data in PVS1."
		for var in self.userVariants:
			if str(var.reference) == "0" or not var.reference:
				var.reference = "-"
			if str(var.alternate) == "0" or not var.alternate:
				var.alternate = "-"
		return [ vepDone , preVEP , exacDone , clinvarDone ]
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
		""" read & parse input .vcf
			Look for VEP CSQ info field & extract all information possible
			Can get allele frequency & clinical annotations after VEP v81 (or so)
			http://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#output
		"""
		inFile = None
		if ( re.match( "\.gz" , inputFile ) ):
			inFile = vcf.Reader( open( inputFile , 'r' ) , compressed=True )    
		else:
			inFile = vcf.Reader( open( inputFile , 'r' ) )
		variantDict = {}
		appendTo = kwargs.get( "appendTo" , "" )
		anyFilter = kwargs.get( "anyFilter" , False )
		self.mutationTypes = kwargs.get( "mutationTypes" , [] )
		includeVCFDetails = kwargs.get( "includeVCFDetails" , False )
		print( "Will capture vcf details for output: " + str( includeVCFDetails ) )
		preVEP = []
		vepDone = False
		exacDone = False
		clinvarDone = False
		vepInfo = OD()
		self.vcfHeaderInfo = []
		metadata = inFile.metadata
		#print( str( metadata ) )
		failedFilter = 0
		failedAF = 0
		failedMT = 0
		[ vepDone , exacDone , clinvarDone ] = self.readMetaData( metadata , inFile.infos , vepInfo )
		for record in inFile:
			try:
				chrom = record.CHROM
				reference = record.REF
				alternates = record.ALT
				if anyFilter and record.FILTER != "PASS":
					failedFilter += 1
					continue
				start = record.start + 1 #1-base beginning of ref
				stop = record.end #0-base ending of ref
				info = record.INFO
				strALT = json.dumps( record.ALT , cls = charger.ALTEncoder , separators = ( ',' , ':' ) )
				strINFO = json.dumps( record.INFO , cls = charger.ALTEncoder , separators = ( ',' , ':' ) ) 
				if ( includeVCFDetails ):
					vcfInfo = '::'.join( [ str( record.CHROM ) , str( record.POS ) , str( record.ID ) , str( record.REF ) , strALT , strINFO ] )
				else:
					vcfInfo = ""
				alti = -1
				for alternate in alternates:
					alti += 1
					alt = str( alternate )
					if alt == "None":
						alt = None
					if record.is_sv:
						if not record.is_sv_precise:
							print( "CharGer::readVCF - SV is not precise" )
							print( record )
					elif record.is_indel and not record.is_deletion: #insertion
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

					var.vcfInfo = vcfInfo
					var.disease = charger.allDiseases

					hasAF = False
					hasAF = self.getAF( info , var , alti )

					self.getVEPConsequences( info , var , preVEP )
					if self.skipIfHighAF( var ):
						#self.filtered.append( var.vcf() )
						failedAF += 1
						continue
					if self.skipIfNotInMutationTypes( var ):
						#self.filtered.append( var.vcf() )
						failedMT += 1
						continue
						
					variantDict = self.appendToList( appendTo , var , \
													 variantDict = variantDict )
			except Exception as E:
				print( "CharGer::readVCF - something went wrong with this variant record:\n\t" )
				print( record )
		totalVars = len( self.userVariants ) + failedFilter + failedAF + failedMT
		print(  "Skipping: " + str( failedFilter ) + " for filters and " + \
				str( failedAF ) + " for AF and " + \
				str( failedMT ) + " for mutation types out of " + \
				str( totalVars ) )
		if totalVars == ( failedFilter + failedAF + failedMT ) and appendTo != "pathogenic":
			print( "All variants failed to pass filters. Stopping CharGer." )
			sys.exit( )
		return [ vepDone , preVEP , exacDone , clinvarDone , variantDict ]

	def getAF( self , info , var , alti ):
		hasAF = False
		withAF = info.get( 'AF' , "noAF" )
		if ( withAF != "noAF" ):
			afs = withAF[alti]
			var.alleleFrequency = afs
			hasAF = True
		return hasAF

	def skipIfHighAF( self , var ):
		if var.alleleFrequency is not None:
			if var.isFrequentAllele( self.keepAF ):
				return True
		return False
		
	def skipIfNotInMutationTypes( self , var ):
		if self.mutationTypes:
			if var.variantClass is None:
				return False
			if var.variantClass not in self.mutationTypes:
				return True
		return False

	def getVEPConsequences( self , info , var , preVEP ):
		csq = info.get( 'CSQ' , "noCSQ" )
		if not csq == "noCSQ":
			vepDone = True
			var.vepVariant = vepvariant()
			for thisCSQ in csq:
				values = thisCSQ.split( "|" )
				aas = self.getRefAltAminoAcids( values , var , preVEP )
				positionPeptide = self.getCodingPosition( values , var , preVEP , "Protein_position" )
				positionCodon = self.getCodingPosition( values , var , preVEP , "cDNA_position" )
				exons = self.getExons( values )
				introns = self.getIntrons( values )
				siftStuff = self.getSIFT( values )
				polyPhenStuff = self.getPolyPhen( values )
				csq_terms = self.getConsequence( values )
				vcv = vepconsequencevariant( \
					chromosome = var.chromosome , \
					start = var.start , \
					stop = var.stop , \
					dbsnp = var.dbsnp , \
					reference = var.reference , \
					alternate = var.alternate , \
					gene_id = self.getVCFKeyIndex( values , "Gene" ) , \
					transcriptCodon = self.getVCFKeyIndex( values , "Feature" ) , \
					consequence_terms = csq_terms , \
					positionCodon = positionCodon , \
					positionPeptide =  positionPeptide , \
					referencePeptide = aas[0] , \
					alternatePeptide = aas[1] , \
					strand = self.getVCFKeyIndex( values , "STRAND" ) , \
					gene = self.getVCFKeyIndex( values , "SYMBOL" ) , \
					gene_symbol_source = self.getVCFKeyIndex( values , "SYMBOL_SOURCE" ) , \
					hgnc_id = self.getVCFKeyIndex( values , "HGNC_ID" ) , \
					biotype = self.getVCFKeyIndex( values , "BIOTYPE" ) , \
					canonical = self.getVCFKeyIndex( values , "CANONICAL" ) , \
					ccds = self.getVCFKeyIndex( values , "CCDS" ) , \
					transcriptPeptide = self.getVCFKeyIndex( values , "ENSP" ) , \
					predictionSIFT = siftStuff[0] , \
					scoreSIFT = siftStuff[1] , \
					predictionPolyphen = polyPhenStuff[0] , \
					scorePolyphen = polyPhenStuff[1] , \
					exon = exons[0] , \
					totalExons = exons[1] , \
					intron = introns[0] , \
					totalIntrons = introns[1] , \
				)
				var.vepVariant.consequences.append( vcv )
			self.getMostSevereConsequence( var )
			hasAF = self.getExAC_MAF( values , var )
			hasAF = self.getGMAF( values , var )
			self.getCLIN_SIG( values , var )

	def getCodingPosition( self , values , var , preVEP , key ):
		pos = None
		if self.getVCFKeyIndex( values , key ):
			position = self.getVCFKeyIndex( values , key )
			poss = position.split( "/" )
			if len( poss ) > 1:
				pos = poss[0]
			else:
				pos = poss[0]
		return pos

	def getRefAltAminoAcids( self , values , var , preVEP ):
		aas = [None , None] 
		if self.getVCFKeyIndex( values , "Amino_acids" ): #8 => Amino_acids
			AAs = self.getVCFKeyIndex( values , "Amino_acids" )
			aas = AAs.split( "/" ) 
			if len( aas ) > 1:
				aas[0] = mafvariant().convertAA( aas[0] )
				aas[1] = mafvariant().convertAA( aas[1] )
			else:
				hgvsp = self.getVCFKeyIndex( values , "HGVSp" ).split( ":" )
				changep = None
				if len( hgvsp ) > 1:
					changep = re.match( "p\." , hgvsp[1] )
				if changep:
					aas = mafvariant().splitHGVSp( hgvsp[1] )
					#print( "AmilaW:aas = " + str(aas) )
					aas[0] = mafvariant().convertAA( aas[0] )
					aas[2] = mafvariant().convertAA( aas[2] )
				else:
					aas.append( None )
					needVEP = True
					preVEP.append( var )
		return aas

	def getExons( self , values ):
		exons = [None , None]
		if self.getVCFKeyIndex( values , "EXON" ): #25 => EXON
			exons = self.getVCFKeyIndex( values , "EXON" ).split( "/" )
			if len( exons ) == 1:
				exons.append(None)
		return exons

	def getIntrons( self , values ):
		introns = [None , None]
		if self.getVCFKeyIndex( values , "INTRON" ): #26 => INTRON
			introns = self.getVCFKeyIndex( values , "INTRON" ).split( "/" )
			if len( introns ) == 1:
				introns.append(None)
		return introns

	def getSIFT( self , values ):
		siftStuff = [None , None]
		if self.getVCFKeyIndex( values , "SIFT" ):
			siftStuff = self.getVCFKeyIndex( values , "SIFT" ).split( "(" ) 
			if len( siftStuff ) == 1:
				siftStuff.append( None )
			else:
				siftStuff[1] = float( siftStuff[1].rstrip( ")" ) )
		return siftStuff

	def getPolyPhen( self , values ):
		polyPhenStuff = [None , None]
		if self.getVCFKeyIndex( values , "PolyPhen" ):
			polyPhenStuff = self.getVCFKeyIndex( values , "PolyPhen" ).split( "(" ) 
			if len( polyPhenStuff ) == 1:
				polyPhenStuff.append( None )
			else:
				polyPhenStuff[1] = float( polyPhenStuff[1].rstrip( ")" ) )
		return polyPhenStuff

	def getConsequence( self , values ):
		consequence_terms = self.getVCFKeyIndex( values , "Consequence" )
		csq_terms = []
		if consequence_terms:
			csq_terms = self.getVCFKeyIndex( values , "Consequence" ).split( "&" )
		return csq_terms

	def getExAC_MAF( self , values , var ):
		#if the .vcf does not have AF
		#then check for ExAC_MAF
		if ( var.alleleFrequency is not None ):
			emaf = self.getVCFKeyIndex( values , "ExAC_MAF" )
			if emaf is not None:
				for alt in emaf.split( "&" ):
					if alt == var.alternate:
						parts = emaf.split( ":" )
						if len( parts ) > 1:
							var.alleleFrequency = parts[1]
							return True
		return False

	def getGMAF( self , values , var ):
		if ( var.alleleFrequency is not None ):
			gmaf = self.getVCFKeyIndex( values , "GMAF" )
			if gmaf is not None:
				for alt in gmaf.split( "&" ):
					if alt == var.alternate:
						parts = gmaf.split( ":" )
						if len( parts ) > 1:
							var.alleleFrequency = parts[1]
							return True
		return False

	def readMetaData( self , metadata , infos , vepInfo ):
		vepDone = False
		exacDone = False
		clinvarDone = True
		for pairs in metadata:
			if pairs == 'VEP':
				print "This .vcf has VEP annotations!"
				vepDone = True
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
		for i in infos.items():
			if i[0] == 'AF':
				print "This .vcf has AF!"
				exacDone = True
			if i[0] == 'clinvar_measureset_id':
				print( "This .vcf has ClinVar annotations!" )
				clinvarDone = True
		return [ vepDone , exacDone , clinvarDone ]

	def getCLIN_SIG( self , values , var ):
		clinical = self.getVCFKeyIndex( values , "CLIN_SIG" )
		description = ""
		if clinical:
			clinvarDone = True
			sigs = clinical.split( "&" )
			description = sigs[0]
			if description == "uncertain":
				description = chargervariant.uncertain
			elif description == "benign":
				description = chargervariant.benign
			elif description == "likely_benign":
				description = chargervariant.likelyBenign
			elif description == "likely_pathogenic":
				description = chargervariant.likelyPathogenic
			elif description == "pathogenic":
				description = chargervariant.pathogenic
			for sig in sigs:
				if sig == "uncertain":
					description = chargervariant.uncertain
					break
				if sig == "benign":
					if description == chargervariant.pathogenic or description == chargervariant.likelyPathogenic:
						description = chargervariant.uncertain
						break
					if description == chargervariant.likelyBenign:
						description = chargervariant.benign
						continue
				if sig == "likely_benign":
					if description == chargervariant.likelyPathogenic or description == chargervariant.pathogenic:
						description = chargervariant.uncertain
						break
					if description == chargervariant.benign:
						continue
					description = chargervariant.likelyBenign
				if sig == "likely_pathogenic":
					if description == chargervariant.benign or description == chargervariant.likelyBenign:
						description = chargervariant.uncertain
						break
					if description == chargervariant.pathogenic:
						continue
				if sig == "pathogenic":
					if description == chargervariant.benign or description == chargervariant.likelyBenign:
						description = chargervariant.uncertain
						break
					if description == chargervariant.likelyPathogenic:
						description = chargervariant.pathogenic
						continue
		var.clinical["description"] = description
		var.clinvarVariant.clinical["description"] = description
		var.clinical["review_status"] = ""
		var.clinvarVariant.clinical["review_status"] = ""

	def getMostSevereConsequence( self , var ):
#TODO parse "splice_region_variant&synonymous_variant" as possible consequence
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
		#print( "setting most severe consequence" )
		#print( str( len( var.vepVariant.consequences ) ) + " consequences to parse" )
		try:
			for cons in var.vepVariant.consequences:
				#cons.printVariant( "~" )
				for term in cons.terms:
					if term in severeRank:
						rank = severeRank.index( term )
					else:
						print( "CharGer Warning: Unrecognized consequence term: " + str( term ) )
						rank = 10000
					if rank < rankMostSevere:
						mostSevere = cons
						rankMostSevere = rank
						mostSevereCons = term
					elif rank == rankMostSevere:
						if cons.canonical:
							mostSevere = cons
			if mostSevere:
				#print( "setting MOST SEVERE CONSEQUENCE" )
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
			print( "CharGer::getMostSevereConsequence Warning: no consequences" + var.genomicVar() )
			pass
						
	def appendToList( self , appendTo , var , variantDict = {} ):
		if appendTo == "user":
			self.userVariants.append( var )
		elif appendTo == "pathogenic":
			pathKey = self.pathogenicKey( var )
			self.pathogenicVariants[pathKey] = var
		elif appendTo == "vepDict":
			genVar = var.vcf( )
			variantDict[genVar] = var
			return variantDict
		else:
			print( "CharGer warning: bad appendTo = " + str( appendTo ) )

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

	def readModesGeneList( self , inputFile , **kwargs ): # gene list formatted "gene", "disease", "mode of inheritance"
		print( "CharGer::readModesGeneList" )
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
				if ( charger.DOMINANT in mode_inheritance ):
					self.inheritanceGeneList[gene][disease] = charger.DOMINANT
				if ( charger.RECESSIVE in mode_inheritance ):
					self.inheritanceGeneList[gene][disease] = charger.RECESSIVE
				#print '\t'.join( [ gene , disease , mode_inheritance ] )
			inFile.close()
		except:
			print "CharGer::readModesGeneList Error: bad gene list file"
	def readPP2GeneList( self , inputFile , **kwargs ):
		print( "CharGer::readPP2GeneList" )
		try:
			inFile = self.safeOpen( inputFile , 'r' , warning = True )
			for line in inFile:
				self.PP2GeneList.add( line.strip() )
			inFile.close()
		except:
			print( "CharGer::readPP2GeneList Error: bad gene list file" )
	def readBP1GeneList( self , inputFile , **kwargs ):
		print( "CharGer::readBP1GeneList" )
		try:
			inFile = self.safeOpen( inputFile , 'r' , warning = True )
			for line in inFile:
				self.BP1GeneList.add( line.strip() )
			inFile.close()
		except:
			print( "CharGer::readBP1GeneList Error: bad gene list file" )
	def readDeNovo( self , inputFile ):
		self.readOtherMAF( inputFile , varDict = self.deNovoVariants )
	def readAssumedDeNovo( self , inputFile ):
		self.readOtherMAF( inputFile , varDict = self.assumedDeNovoVariants )
	def readCoSegregate( self , inputFile ):
		self.readOtherMAF( inputFile , varDict = self.coSegregateVariants )
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
		self.thresholdAF = kwargs.get( "thresholdAF" , 0.0005 )
		self.keepAF = kwargs.get( "keepAF" , 1 )
		t = time.time()
		self.getVEP( **kwargs )
		self.printRunTime( "VEP" , self.runTime( t ) )
		t = time.time()
		self.getClinVar( **kwargs )
		self.printRunTime( "ClinVar" , self.runTime( t ) )
		t = time.time()
		de = kwargs.get( 'exac' , "ASDF" )
		#print( de )
		self.getExAC( **kwargs )
		self.printRunTime( "exac" , self.runTime( t ) )
		self.fillMissingVariantInfo()

	def fillMissingVariantInfo( self ):
		for var in self.userVariants:
			var.fillMissingInfo( var )

	def getClinVar( self , **kwargs ):
		print( "charger::getClinVar" )
		macClinVarTSV = kwargs.get( 'macClinVarTSV' , None )
		macClinVarVCF = kwargs.get( 'macClinVarVCF' , None )
		doClinVar = kwargs.get( 'clinvar' , False )
		#print( '  '.join( [ str( doClinVar ) , str( macClinVarTSV ) , str( macClinVarVCF ) ] ) ) 
		if doClinVar:
			if macClinVarTSV:
				clinvarSet = self.getMacClinVarTSV( macClinVarTSV )
				self.userVariants = self.matchClinVar( self.userVariants , clinvarSet )
			elif macClinVarVCF:
				self.getMacClinVarVCF( macClinVarTSV )
			else:
				self.getClinVarviaREST( **kwargs )

	def getClinVarviaREST( self , **kwargs ):
		summaryBatchSize = kwargs.get( 'summaryBatchSize' , 500 )
		maxSearchBatchSize = 50
		searchBatchSize = kwargs.get( 'searchBatchSize' , maxSearchBatchSize )
		if searchBatchSize > maxSearchBatchSize:
			message = "warning: ClinVar ReST search batch size given is "
			message += "greater than max allowed ("
			message += str( maxSearchBatchSize ) + ")"
			message += ". Overriding to max search batch size."
			print( message )
			searchBatchSize = maxSearchBatchSize
		ent = entrezapi()
		i = 0
		for varsStart in range( 0 , len( self.userVariants ) , int(searchBatchSize) ):
			varsEnd = varsStart + int(searchBatchSize)
			varsSet = self.userVariants[varsStart:varsEnd]
			ent.prepQuery( varsSet )
			ent.subset = entrezapi.esearch
			ent.database = entrezapi.clinvar
			clinvarsSet = ent.doBatch( summaryBatchSize )
			varsSet = self.matchClinVar( varsSet , clinvarsSet )
			self.userVariants[varsStart:varsEnd] = varSet
			#self.clinvarVariants.update( varsBoth["clinvarVariants"] )
			#self.userVariants[varsStart:varsEnd] = self.matchClinVar( varsSet , clinvarsSet )

	def getMacClinVarTSV( self , tsvfile ):
		"""
		Function adapted to parse either version of MacArthur ClinVar file:
		- Previous version header:
		chrom	pos	ref	alt	measureset_type	measureset_id	rcv	allele_id
		symbol	hgvs_c	hgvs_p	molecular_consequence	clinical_significance
		pathogenic	benign	conflicted	review_status	gold_stars
		all_submitters	all_traits	all_pmids	inheritance_modes
		age_of_onset	prevalence	disease_mechanism	origin	xrefs

		- Current version header:
		chrom	pos	ref	alt	start	stop	strand	variation_type	variation_id
		rcv	scv	allele_id	symbol	hgvs_c	hgvs_p	molecular_consequence	clinical_significance
		clinical_significance_ordered	pathogenic	likely_pathogenic	uncertain_significance	likely_benign
		benign	review_status	review_status_ordered	last_evaluated	all_submitters	submitters_ordered
		all_traits	all_pmids	inheritance_modes	age_of_onset	prevalence	disease_mechanism	origin
		xrefs	dates_ordered	gold_stars	conflicted
		"""
		clinvarSet = dict()
		with gzip.open( tsvfile , "rb" ) as macFile:
			header = macFile.readline().strip().split("\t") # read in file header; columns are now called by name instead of index
			for line in macFile:
				fields = ( line.rstrip( ) ).split( "\t" )
				[ description , status ] = self.parseMacPathogenicity( header, fields ) # no need to specify which fields here anymore; parseMacPathogenicity now knows which specific columns to look for
				if len(header) > 27: # if yes, file is in the new format
					var = clinvarvariant( chromosome = fields[header.index("chrom")] , \
										  start = fields[header.index("pos")] , \
										  reference = fields[header.index("ref")] , \
										  alternate = fields[header.index("alt")] , \
										  uid = fields[header.index("variation_id")], \
										  gene = fields[header.index("symbol")] , \
										  clinical = { "description" : description , "review_status" : status } , \
										  trait = { fields[header.index("xrefs")] : fields[header.index("all_traits")] } )
				else: # file in the old format
					var = clinvarvariant( chromosome = fields[header.index("chrom")] , \
										  start = fields[header.index("pos")] , \
										  reference = fields[header.index("ref")] , \
										  alternate = fields[header.index("alt")] , \
										  uid = fields[header.index("measureset_id")], \
										  gene = fields[header.index("symbol")] , \
										  clinical = { "description" : description , "review_status" : status } , \
										  trait = { fields[-1] : fields[header.index("all_traits")] } )
				var.setStopFromReferenceAndAlternate( )
				var.splitHGVSc( fields[header.index("hgvs_c")] )
				var.splitHGVSp( fields[header.index("hgvs_p")] )
				#var.printVariant( "," )
				#print( var.proteogenomicVar( ) )
				#sys.exit() 
				clinvarSet[var.uid] = var
		print( "Have " + str( len( clinvarSet ) ) + " uid's from MacArthur ClinVar .tsv file: " + tsvfile )
		return clinvarSet

	@staticmethod
	def parseMacPathogenicity( header, fields ): # addded header argument, so can recognize column names from file header (defined by the modified getMacClinVarTSV function)
		named = fields[header.index("clinical_significance")]
		isPathogenic = fields[header.index("pathogenic")]
		if isPathogenic == "1;0":
			isPathogenic = 1
		elif int(isPathogenic) >= 1:
			isPathogenic = 1

		isBenign = fields[header.index("benign")]
		if isBenign == "1;0": 
			isBenign = 1
		elif int(isBenign) >= 1:
			isBenign = 1

		isConflicted = fields[header.index("conflicted")]
		if isConflicted == "1;0":
			isConflicted = 1
		elif int(isConflicted) >= 1:
			isConflicted = 1

		status = fields[header.index("review_status")]
		desc = chargervariant.uncertain

		# added snippet (next 5 lines) to solve cases where variant is categorized likely pathogenic, likely benign, 
		# or as both uncertain significance and likely pathogenic or likely benign. 
		if "likely_pathogenic" in header:
			if fields[header.index("likely_pathogenic")] >= 1:
				isPathogenic = 1
			
			if fields[header.index("likely_benign")] >= 1:
				isBenign = 1
		
		if isConflicted == 1:
			return [ desc , status ]
		
		# adjusted function to read values split by either ";" or "/"
		if isBenign == 1:
			if ";" in named: # old version
				for desc in named.split( ";" ):
					if re.match( desc.lower( ) , "ikely" ) and desc != chargervariant.benign:
						desc = chargervariant.likelyBenign
					elif re.match( desc.lower( ) , "benign" ):
						desc = chargervariant.benign
						break
			else:
				for desc in named.split( "/" ): # new version
					if re.match( desc.lower( ) , "ikely" ) and desc != chargervariant.benign:
						desc = chargervariant.likelyBenign
					elif re.match( desc.lower( ) , "benign" ):
						desc = chargervariant.benign
						break
		
		if isPathogenic == 1:
			if ";" in named:
				for desc in named.split( ";" ): # old version
					if re.match( desc.lower( ) , "ikely" ) and desc != chargervariant.pathogenic:
						desc = chargervariant.likelyPanic
					elif re.match( desc.lower( ) , "athog" ):
						desc = chargervariant.pathogenic
						break
			else:
				for desc in named.split( "/" ): # new version
					if re.match( desc.lower( ) , "ikely" ) and desc != chargervariant.pathogenic:
						desc = chargervariant.likelyPanic
					elif re.match( desc.lower( ) , "athog" ):
						desc = chargervariant.pathogenic
						break
		return [ desc , status ]

	def getMacClinVarVCF( self , vcffile ):
		print( "TODO: add macClinVarVCF read in" )
		sys.exit()

	def getExAC( self , **kwargs ):
		exacVCF = kwargs.get( 'exacVCF' , None )
		doExAC = kwargs.get( 'exac' , False )
		useHarvard = kwargs.get( 'harvard' , True )
		threshold = kwargs.get( 'thresholdAF' , 0.0005 )
		if doExAC:
			sys.stdout.write( "charger::getExAC " )
			if exacVCF:
				print( "from local file - " + exacVCF )
				exac = exacparser( )
				exacQueries = self.getUniqueGenomicVariantList( self.userVariants )
				entries = exac.searchVCF( vcf = exacVCF , queries = exacQueries )
			else:
				print( "through BioMine ReST" )
				totalVars = len( self.userVariants )
				exac = exacapi(harvard=useHarvard)
#entries by genomivVar
				exacIn = self.getUniqueGenomicVariantList( self.userVariants )
				entries = exac.getAlleleFrequencies( exacIn )
			elen = 0
			common = 0
			rare = 0
			totalVars = len( self.userVariants )
			if entries:
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
		doVEP = kwargs.get( 'vep' , False )
		preVEP = kwargs.get( 'prevep' , [] )
		doREST = kwargs.get( 'rest' , False )
		vepScript = kwargs.get( 'vepScript' , "" )
		if doVEP:
			sys.stdout.write( "charger::getVEP " )
			if not vepScript:
				print( "through BioMine ReST" )
				self.getVEPviaREST( **kwargs )
			else:
				print( "from local tool" )
				self.getVEPviaCMD( **kwargs )
		else:
			print( "charger::getVEP Warning: skipping VEP" )

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
		# pdb.set_trace()
		print "\nVEP annotated userVariants " + str( luvafter )
		print "VEP annotated " + str(aluv) + " from the original set of " + str(luv)

	def getVEPviaCMD( self , **kwargs ):
#TODO check ensemblRelease & grch for compatibility, GRCh37 -> [55,75] & GRCh38 -> [76,87]
		#print( kwargs )
		perl = kwargs.get( 'perl', "/bin/perl" ) #, defaultVEPScript )
		vepScript = kwargs.get( 'vepScript', "" ) #, defaultVEPScript )
		vepConfig = kwargs.get( 'vepConfig', "" )
		vcfFile = kwargs.get( 'vcf' , "" )
		#defaultVEPDir = "./"
		#vepDir = kwargs.get( 'vepDir' , defaultVEPDir )
		#vepCacheDir = kwargs.get( 'vepCache' , defaultVEPDir )
		defaultEnsemblRelease = str( 75 )
		defaultVEPVersion = str( 87 )
		defaultGRCh = str( 37 )
		defaultForks = str( 1 )
		##defaultVEPScript = vepDir + "/variant_effect_predictor.pl"
		##hdir = vepVersion + "_" + assembly
		#hdir = ensemblRelease + "_" + assembly
		#fa = "Homo_sapiens." + assembly + "." + ensemblRelease + \
		#	 ".dna.primary_assembly.fa.gz" 
		#defaultFastaArray = [ vepCacheDir , "homo_sapiens" , hdir , fa ] 
		#defaultFasta = '/'.join( defaultFastaArray )
		#fasta = kwargs.get( 'referenceFasta' , defaultFasta )
		#forks = kwargs.get( 'fork' , 0 )
		#print( defaultFasta )
		#print( fasta )
		defaultVEPOutput = "./"
		if vepConfig:
			with open( vepConfig, 'r') as configFH:
				configFH.next()
				for line in configFH:
					line.strip()
					match = re.match( "output_file\s+(.*)", line )
					if match:
						defaultVEPOutput = match.group(1)
						break
		outputFile = kwargs.get( 'vepOutput' , defaultVEPOutput )
		if not outputFile: # this is necessary b/c vepOutput is initialized as "None"
			outputFile = defaultVEPOutput
		#print("AmilaW:outputFile"+str(outputFile))

		if vcfFile:
			vep_command = []
			#print("AmilaW:outputFile"+str(outputFile))
			if vepConfig:
				vep_command = [ perl , vepScript , \
					"--config", vepConfig ]
				#print( "AW:inside vep config loop" )
			else:
				#defaultVEPDir = "./"
				#vepDir = kwargs.get( 'vepDir' , defaultVEPDir )
				vepCacheDir = kwargs.get( 'vepCache', "./" ) #, defaultVEPDir )
				ensemblRelease = kwargs.get( 'ensemblRelease' , defaultEnsemblRelease )
				vepVersion = kwargs.get( 'vepVersion' , defaultVEPVersion )
				grch = kwargs.get( 'grch' , defaultGRCh )
				#defaultVEPOutput = "./charger.vep.vcf"
				#defaultVEPScript = vepDir + "/variant_effect_predictor.pl"
				assembly = "GRCh" + str( grch )
				#hdir = vepVersion + "_" + assembly
				hdir = ensemblRelease + "_" + assembly
				fa = "Homo_sapiens." + assembly + "." + ensemblRelease + \
					".dna.primary_assembly.fa.gz"
				defaultFastaArray = [ vepCacheDir , "homo_sapiens" , hdir , fa ]
				defaultFasta = '/'.join( defaultFastaArray )
				fasta = kwargs.get( 'referenceFasta' , defaultFasta )
				#outputFile = kwargs.get( 'vepOutput' , defaultVEPOutput )
				#vepScript = kwargs.get( 'vepScript', "" ) #, defaultVEPScript )
				#vcfFile = kwargs.get( 'vcf' , "" )
				forks = kwargs.get( 'fork' , defaultForks )
				#vepConfig = kwargs.get( 'vepConfig', "" )
				#print( defaultFasta )
				#print( fasta )
				vep_command = [ perl , vepScript , \
					"--species" , "homo_sapiens" , \
					"--assembly" , assembly , \
					"--input_file" , vcfFile , \
					"--output_file" , outputFile , \
					"--format" , "vcf" , \
					"--fork" , forks , \
					"--fasta" , fasta , \
					"--everything" , \
					"--vcf" , \
					"--cache" , \
					"--offline" , \
					"--no_progress" , \
					"--total_length" , \
					"--no_escape" , \
					"--xref_refseq" , \
					"--force_overwrite" , \
					"--no_stats" , \
					"--dir" , vepCacheDir , \
					"--dir_cache" , vepCacheDir , \
					"--verbose" , \
					#"--quiet" , \
					#"--help" \
				]
			#print( "AW:VEP command")
			print( vep_command )
			#print("AmilaW:outputFile"+str(outputFile))
			try:
				#returnCode = subprocess.call( ' '.join( vep_command ) )
				returnCode = subprocess.call( vep_command )
				#print( str( returnCode ) )
				#print("AmilaW:outputFile")
				#print(str(outputFile))
				#print(str(os.path.getsize(outputFile)))
				if os.path.getsize( outputFile ) == 0:
					print( "CharGer ERROR: VEP did not produce output: " + outputFile )
					sys.exit( )
			except subprocess.CalledProcessError:
				print( "CharGer ERROR: VEP messed up" )
				sys.exit( )
				#pass # handle errors in the called executable
			except OSError:
				print( "CharGer ERROR: VEP not found" )
				sys.exit( )
				#pass # executable not found
			print( "local VEP done" )
			[ vepDone , preVEP , exacDone , clinvarDone , vepVariants \
			] = self.readVCF( outputFile , appendTo = "vepDict" , **kwargs )
			print( "finished reading in VEP .vcf: " + outputFile )
			print( "now matching VEP annotated variants" )
			self.matchVEP( vepVariants )

#### Helper methods for data retrieval ####
	def matchClinVar( self , userVariants , clinvarVariants ):
		print( "matchClinVar!" )
		matched = 0
		for var in userVariants:
			for uid in clinvarVariants:
				cvar = clinvarVariants[uid]
				if var.sameGenomicVariant( cvar ):
					matched += 1
					var.clinical = cvar.clinical
					cvar.fillMissingInfo( var )
					var.fillMissingInfo( cvar )
					var.clinvarVariant = cvar
					break
		print( "Matched " + str( matched ) + " ClinVar entries to user variants" )
		return userVariants

	def matchVEP( self , vepVariants ):
		# print "Charger::matchVEP"
		currentVars = str( len( self.userVariants ) )
		removedVars = 0
		for var in self.userVariants: # var is a chargervariant object
			for vepVar in vepVariants:
				if var.sameGenomicVariant( vepVariants[vepVar] ):
					# var.printVariant("|")
					charVEPVar = vepVariants[vepVar]
					# pdb.set_trace()
					# print "ACSW: charVEPVar = "
					# charVEPVar.printVariant("|")
					if not var.vepVariant:
						var.vepVariant = charVEPVar
					# charVEPVar.fillMissingInfo( var.vepVariant )
					# pdb.set_trace()
					var.fillMissingInfo( charVEPVar )
					# pdb.set_trace()
					var.copyMostSevereConsequence()
					# pdb.set_trace()
					var.vepVariant = charVEPVar
					# if var.gene == "NA" or var.gene is None:
					# 	print "ACSW:after missing gene info"
					# 	var.printVariant("|")
					# 	pdb.set_trace()
					# var.printVariant("|")
					# charVEPVar.printVariant("|")
					# pdb.set_trace()
					break
		afterFilter = len( self.userVariants )
		removedVars = int( currentVars ) - afterFilter
		print( "Removed " + str( removedVars ) + " from initial user set of " + currentVars + ", now have " + str( len( self.userVariants ) ) + " user variants left." )

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

# any level ending with C (ex. PMC1, PPC1) are CharGer-invented modules #

##### Very Strong #####
	def PVS1( self , expressionThreshold = 0.2 ):
		print "CharGer module PVS1"
		print "- truncations in genes where LOF is a known mechanism of the disease"
		print "- require the mode of inheritance to be dominant (assuming heterzygosity) and co-occurence with reduced gene expression"
		print "- run concurrently with PSC1, PMC1, PM4, PPC1, and PPC2 - "
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
		print "Recessive truncations of susceptible genes"

##### Moderate #####
	def PMC1( self ):
		print "CharGer module PMC1"
		print "Truncations of genes when no gene list provided"

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
						chromosome = fields[7]
						start = fields[8]
						stop = fields[9]
						reference = fields[10]
						alternate = fields[11]
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
		print( "CharGer module PP2: Missense variant in a gene that has low rate of benign missense and in which missense are common mechanism of disease" )
		if self.PP2GeneList: #gene
			for var in self.userVariants:
				varGene = var.gene
				varSample = var.sample
				varClass = "asdfasdf"
				if var.variantClass:
					varClass = var.variantClass
				varVEPClass = "asdfasdf"
				if var.vepVariant:
					varVEPClass = var.vepVariant.mostSevereConsequence
				if ( "missense" in varClass.lower() ) or \
					( "missense" in varVEPClass.lower() ):
					if varGene in self.PP2GeneList: # check if in gene list
						var.PP2 = True # if call is true then check expression effect
						var.addSummary( "PP2(Missense variant in gene from PP2 gene list)" )
		else: 
			print "CharGer::PP2 Error: Cannot evaluate PP2: No PP2 gene list supplied."
#		print "- "
	def PP3( self , minimumEvidence ):
		print "CharGer module PP3"
		print "- multiple lines of in silico evidence of deliterous effect"
		callSIFTdam = "damaging"
		callSIFTdel = "deleterious"
		thresholdSIFT = 0.05
		callPolyphen = "probably_damaging"
		thresholdPolyphen = 0.432
		callBlosum62 = -2
		callCompara = 2
		callImpact = "high"
		fracMaxEntScan = 0.8
		callGeneSplicer = ""
		nFound = 0
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
								nFound += 1
								break
				if var.PP3:
					var.addSummary( "PP3(Multiple (>=" + str( minimumEvidence ) \
						+ ") in silico predictions of deliterious effect=" \
						+ "|".join( case ) + ")" )
		print( "Found " + str( nFound ) + " variants with >= " + str( minimumEvidence ) + " of in silico evidence" )

	def PP4( self ):
		print "CharGer module PP4: not yet implemented"
#		print "- "
	def PP5( self ):
		print "CharGer module PP5: not yet implemented"
#		print "- "
	def PPC1( self ):
		print "CharGer module PPC1"
		print "- protein length changes due to inframe indels or nonstop variant of other, not-specificied genes - "
	def PPC2( self ):
		print "CharGer module PPC2"
		print "- protein length changes due to inframe indels or nonstop variant when no susceptibility genes given - "

### helper functions of evidence levels ###
	def runIndelModules( self ):
		maf_truncations = ["Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Splice_Site"] #,"Nonstop_Mutation"
		vep_truncations = ["transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained",\
							"frameshift_variant","start_lost"]
		lenShift = ["In_Frame_Del","In_Frame_Ins","Nonstop_Mutation"]
		vep_inframe = [	"inframe_insertion" , "inframe_deletion" , \
						"stop_lost" ]

		if not self.inheritanceGeneList:
			print( "CharGer::runIndelModules Error: Cannot evaluate PVS1 or PM4: No gene list supplied." )
		for var in self.userVariants:
			varClass = var.variantClass
			varGene = var.gene
			varVEPClass = ""
			if var.vepVariant:
				varVEPClass = var.vepVariant.mostSevereConsequence
			altPeptide = var.alternatePeptide
			if (varClass in maf_truncations) or \
					(varVEPClass in vep_truncations) or \
					altPeptide == "*" or \
					altPeptide == "fs":
				if self.inheritanceGeneList:
					if varGene in self.inheritanceGeneList: # check if in gene list
						varDisease = var.disease # no disease field in MAF; may require user input	
						if varDisease in self.inheritanceGeneList[varGene] \
						or charger.allDiseases in self.inheritanceGeneList[varGene]:
							if charger.DOMINANT in self.inheritanceGeneList[varGene][varDisease].lower() \
							or charger.DOMINANT in self.inheritanceGeneList[varGene][charger.allDiseases].lower():
								var.PVS1 = True # if call is true then check expression effect
								if self.userExpression: # consider expression data only if the user has supplied an expression matrix
									varSample = var.sample
									if self.userExpression[varSample][varGene] >= expressionThreshold:
										var.PVS1 = False
							elif charger.RECESSIVE in self.inheritanceGeneList[varGene][varDisease].lower() \
							or charger.RECESSIVE in self.inheritanceGeneList[varGene][charger.allDiseases].lower():
								var.PSC1 = True
				else: 
					# in case of no gene list, make all truncations PMC1
					var.PMC1 = True
			elif (varClass in lenShift \
			or varClass in vep_inframe):
				if self.inheritanceGeneList:
					if varGene in self.inheritanceGeneList: # check if in gene list
						varDisease = var.disease # no disease field in MAF; may require user input	
						if varDisease in self.inheritanceGeneList[varGene] \
						or charger.allDiseases in self.inheritanceGeneList[varGene]:
							if charger.DOMINANT in self.inheritanceGeneList[varGene][varDisease].lower() \
							or charger.DOMINANT in self.inheritanceGeneList[varGene][charger.allDiseases].lower():
								var.PM4 = True
							elif charger.RECESSIVE in self.inheritanceGeneList[varGene][varDisease].lower() \
							or charger.RECESSIVE in self.inheritanceGeneList[varGene][charger.allDiseases].lower():
								var.PPC1 = True
				else: 
					# in case of no gene list, make all inframes PPC1
					var.PPC2 = True

			if var.PVS1:
				var.addSummary( "PVS1(" + str( varClass ) + " in susceptible gene " + str( varGene ) + ")" )
			if var.PSC1:
				var.addSummary( "PSC1(" + str( varClass ) + " recessive in gene " + str( varGene ) + ")" )
			if var.PMC1:
				var.addSummary( "PMC1(" + str( varClass ) + " no gene list but in gene " + str( varGene ) + ")" )
			if var.PM4:
				var.addSummary( "PM4(" + str( varClass ) + " in susceptible gene " + str( varGene ) + ")" )
			if var.PPC1:
				var.addSummary( "PPC1(" + str( varClass ) + " recessive in gene " + str( varGene ) + ")" )
			if var.PPC2:
				var.addSummary( "PPC2(" + str( varClass ) + " no gene list but in gene " + str( varGene ) + ")" )

	def peptideChange( self , mod , **kwargs ):
		called = 0
		for var in self.userVariants:
			call = False
			if mod == "PS1":
				call = var.PS1
			if mod == "PM5":
				call = var.PM5
			if mod == "BSC1":
				call = var.BSC1
			if mod == "BMC1":
				call = var.BMC1
			if not call: #is already true
				CVchecked = 0
				PVchecked = 0
				if var.clinvarVariant or self.pathogenicVariants:
					if var.vepVariant:
						if var.vepVariant.consequences:
							for consequence in var.vepVariant.consequences:
								CVchecked = self.checkClinVarPC( var , mod , consequence=consequence )
								PVchecked = self.checkPathogenicVariants( var , mod , consequence=consequence )
								if CVchecked or PVchecked:
									called += 1
					else:
						CVchecked = self.checkClinVarPC( var , mod )
						PVchecked = self.checkPathogenicVariants( var , mod , var )
						if CVchecked or PVchecked:
							called += 1
			if var.PS1 and mod == "PS1":
				var.addSummary( "PS1(Peptide change is known pathogenic)" )
			if var.PM5 and mod == "PM5":
				var.addSummary( "PM5(Peptide change at the same location of a known pathogenic change)" )
			if var.BSC1 and mod == "BSC1":
				var.addSummary( "BSC1(Peptide change is known benign)" )
			if var.BMC1 and mod == "BMC1":
				var.addSummary( "BMC1(Peptide change at the same location of a known benign change)" )
		if mod == "PS1" or mod == "PM5":
			print mod + " found " + str(called) + " pathogenic variants"
		if mod == "BSC1" or mod == "BMC1":
			print mod + " found " + str(called) + " benign variants"
		
	def checkClinVarPC( self , var , mod , **kwargs ):
		called = 0
		consequence = kwargs.get( 'consequence' , var )
		if var.clinvarVariant:
			clinvarVar = var.clinvarVariant
			clin = clinvarVar.clinical
			if consequence.sameGenomicVariant( clinvarVar ):
				#if genomic change is the same, then PS1
				if clin["description"] == clinvarvariant.benign:
					if mod == "BSC1":
						#print( "also BSC1 via samGenomicVariant" )
						var.BSC1 = True # already pathogenic still suffices to be BSC1
						called = 1
				if clin["description"] == clinvarvariant.pathogenic:
					if mod == "PS1":
						#print( "also PS1 via samGenomicVariant" )
						var.PS1 = True # already pathogenic still suffices to be PS1
						called = 1
			elif consequence.sameGenomicReference( clinvarVar ):
			#if genomic change is different, but the peptide change is the same, then PS1
				if clinvarVar.alternatePeptide == consequence.alternatePeptide: #same amino acid change
					if clin["description"] == clinvarvariant.benign:
						if mod == "BSC1":
							#print( "also BSC1 via sameGenomicReference" )
							var.BSC1 = True
							called = 1
					if clin["description"] == clinvarvariant.pathogenic:
						if mod == "PS1":
							#print( "also PS1 via sameGenomicReference" )
							var.PS1 = True
							called = 1
			if consequence.samePeptideReference( clinvarVar ):
				if not consequence.samePeptideChange( clinvarVar ):
				#if peptide change is different, but the peptide reference is the same, then PM5
					if consequence.plausibleCodonFrame( clinvarVar ):
						if clin["description"] == clinvarvariant.benign:
							if mod == "BMC1":
								#print( "also PS5 via plausibleCodonFrame" )
								var.BMC1 = True # already benign still suffices to be BSC1
								called = 1
						if clin["description"] == clinvarvariant.pathogenic:
							if mod == "PM5":
								#print( "also PS5 via plausibleCodonFrame" )
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
		#Allele frequency is greater than expected for disorder
	def BS2( self ):
		print "CharGer module BS2: not yet implemented"
		# Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder, with full penetrance expected at an early age
	def BS3( self ):
		print "CharGer module BS3: not yet implemented"
		#print " - in vitro or in vivo functional studies with no damaging effect on protein function or splicing"
	def BS4( self ):
		print "CharGer module BS4: not yet implemented"
		# Lack of segregation in affected members of a family
		#Caveat: The presence of phenocopies for common phenotypes (i.e., cancer, epilepsy) can mimic lack of segregation among affected individuals. Also, families may have more than one pathogenic variant contributing to an autosomal dominant disorder, further confounding an apparent lack of segregation.
	def BSC1( self ):
		print "CharGer module BSC1"
		print "- same peptide change as a previously established benign variant"
		self.peptideChange( "BSC1" )

#### Moderate ####
	def BMC1( self ):
		print "CharGer module BMC1"
		print "- different peptide change of a benign variant at the same reference peptide"
		self.peptideChange( "BMC1" )

#### Supporting ####
	def BP1( self ):
		print( "CharGer module BP1: Missense variant in a gene for which primarily truncations cause disease" )
		if self.BP1GeneList: #gene
			for var in self.userVariants:
				varGene = var.gene
				varSample = var.sample
				varClass = "asdfasdf"
				if var.variantClass:
					varClass = var.variantClass
				varVEPClass = "asdfasdf"
				if var.vepVariant:
					varVEPClass = var.vepVariant.mostSevereConsequence
				if ( "missense" in varClass.lower() ) or \
					( "missense" in varVEPClass.lower() ):
					if varGene in self.BP1GeneList: # check if in gene list
						var.BP1 = True # if call is true then check expression effect
						var.addSummary( "BP1(Missense variant in gene from BP1 gene list)" )
		else: 
			print "CharGer::BP1 Error: Cannot evaluate BP1: No BP1 gene list supplied."
	def BP2( self ):
		print "CharGer module BP2: not yet implemented"
		#Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or observed in cis with a pathogenic variant in any inheritance pattern
	def BP3( self ):
		print "CharGer module BP3: not yet implemented"
		#In-frame deletions/insertions in a repetitive region without a known function
	def BP4( self , minimumEvidence ):
		print "CharGer module BP4"
		print " - in silico evidence of no damage"
		callSIFTdam = "damaging"
		callSIFTdel = "deleterious"
		thresholdSIFT = 0.05
		callPolyphen = "probably_damaging"
		thresholdPolyphen = 0.432
		callBlosum62 = -2
		callCompara = 2
		callImpact = "high"
		fracMaxEntScan = 0.8
		callGeneSplicer = ""
		nFound = 0
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
								nFound += 1
								break
				if var.BP4:
					var.addSummary( "BP4(Multiple (>=" + str( minimumEvidence ) \
						+ ") in silico predictions of non-deleterious effect=" \
						+ "|".join( case ) + ")" )
		print( "Found " + str( nFound ) + " variants with >= " + str( minimumEvidence ) + " with in silico evidence" )

	def BP5( self ):
		print "CharGer module BP5: not yet implemented"
		#print" - multiple lines of computational evidence suggesting no impact on gene or product"
	def BP6( self ):
		print "CharGer module BP6: not yet implemented"
		#Reputable source recently reports variant as benign, but the evidence is not available to the laboratory to perform an independent evaluation
	def BP7( self ):
		print "CharGer module BP7: not yet implemented"
		#A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site AND the nucleotide is not highly conserved

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
		print( "write " + str( len( self.userVariants ) ) + \
			   " charged user variants to " + outFile )
		annotateInput = kwargs.get( 'annotate' , False )
		skipURLTest = kwargs.get( 'skipURLTest' , True )
		if skipURLTest:
			print( "charger::writeSummary Warning: skipping pubmed link tests" )
		else:
			print( "charger::writeSummary behavior: will test pubmed links" )
		if annotateInput:
			self.annotateInputTSV( outFile , delim = delim )
			return True
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
			"CharGer_Summary" , "VCF_Details" ] )
			# "VEP_Annotations" , \
			# "VCF_Headers" , "VCF_INFO" , "CharGer_Summary"] )
		try:
			outFH.write( headLine )
			if asHTML:
				outFH.write( "</td></tr>" )
			else:
				outFH.write( "\n" )
			for var in self.userVariants:
				genVar = var.vcf()
				#if genVar in self.filtered:
				#	print( "skipping " + genVar + " -- " + var.proteogenomicVar() + " due to mutation-type or allele frequency filters." )
				#	continue
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
				self.appendStr( fields , var.alleleFrequency , emptyValue = 0 )
				#self.appendStr( fields,var.vepVariant.mostSevereConsequence) ## this line will fail you on insertions regardless of all the checks in appendStr
				mSC = var.returnMostSevereConsequence()
				self.appendStr( fields , mSC )
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
						text += var.clinvarVariant.linkPubMed( skip = skipURLTest )
						text += "\">uid="
						text += str( var.clinvarVariant.uid )
						text += "</a>"
						self.appendStr( fields, text )
					else:
						self.appendStr( fields,var.clinvarVariant.linkPubMed( skip = skipURLTest ))
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
				self.appendStr( fields , ' -- '.join( var.callSummary ) )
				self.appendStr( fields , var.vcfInfo )

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
		return True

	def annotateInputTSV( self , outFile , **kwargs ):
		annotateInput = kwargs.get( 'annotate' , False )
		if not annotateInput:
			return False
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
			i = 0
			header = inFile.readline( )
			outFile.write( header + "\t" + chargervariant.annotationsHeader( ) )
			for line in inFile:
				userVar = self.userVariants[i]
				i += 1
				line.strip()
				outFile.write( line + "\t" + userVar.annotations( ) )
			inFile.close()
			outFile.close()
		except:
			print( "TODO: fill out exception to annotateInputTSV" )
			pass

	def getVCFKeyIndex( self , values , field ):
		if field in self.vcfKeyIndex:
			if self.vcfKeyIndex[field] < len( values ):
				return values[self.vcfKeyIndex[field]]
		return None

	@staticmethod
	def appendStr( array, value , emptyValue = "NA" ):
		try:
			if value is None:
				value = emptyValue
			if value == "":
				value = emptyValue
			array.append( str( value ) )
		except:
			print "CharGer Warning: failed to append a value\n"
			array.append( emptyValue )
			pass

	@staticmethod
	def safeOpen( inputFile , rw , **kwargs ):
		errwar = kwargs.get( 'warning' , False )
		try:
			return open( inputFile , rw )
		except:
			if errwar:
				print( "CharGer Warning: could not open " + str( inputFile ) )
			else:
				print( "CharGer Error: could not open " + str( inputFile ) )
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

	class ALTEncoder( json.JSONEncoder ):
		def default( self , obj ):
			try:
				if isinstance( obj , vcf.model._Substitution ):
					return str( obj )
				if isinstance( obj , vcf.model._SV ):
					return str( obj )
				return json.JSONEncoder.default( self , obj )
			except Exception as E:
				print( E )
				sys.exit()
