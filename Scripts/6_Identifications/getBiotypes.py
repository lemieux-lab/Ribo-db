import warnings
import time
import sys
sys.path += sys.path + ["/u/ruizma/Jupyter/PipelineAligment/Scripts/PythonScripts/General/", '/u/ruizma/PIPELINE/Scripts/DataPreparation/PythonScripts/12_Compare_Canonical_Non_Canonical_db', '/u/ruizma/lib/rabadb','/u/ruizma/lib/pygeno']
warnings.filterwarnings("ignore")
import array
import math
import pickle
import itertools
import datetime
from difflib import SequenceMatcher
from pyGeno.tools.UsefulFunctions import *
import logging
import csv
import numpy as np
from scipy.stats import mannwhitneyu
from scipy.stats import kruskal
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.style
matplotlib.style.use("ggplot")
import seaborn as sns
sns.set(style="whitegrid")
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pysam
from pyGeno.Genome import * 
from pyGeno.tools.UsefulFunctions import *
import lookupOpenProtDB as openProt
import pandas as pd
import copy
import pysam
import getMCS_from_aa_seq as getMCS
from pyGeno.tools.BinarySequence import *


class biotypeMAPs:


	def __init__(self, annotation, toSaveFile, fastaFile, filepath_index, freeByesSNPs, quality ):
		self.toSaveFile = toSaveFile
		with open(annotation, 'rb') as fp:
			self.annotationTranscripts= pickle.load(fp)
		self.categories = {'Alt-TIS':0}
		self.categories_peptides = {'Alt-TIS':[]}
		self.bioTypePeptide = {}
		self.get_MCS = getMCS.getTranscriptInformation(freeByesSNPs, quality, fastaFile, filepath_index)
		print 'biotypeMAPs OK'

	def restartDictionnary(self):
		self.categories = {'Alt-TIS':0}
		self.categories_peptides = {'Alt-TIS':[]}
		self.bioTypePeptide = {}
		self.canonical_peptides = {}
		self.non_canonical_peptides = {}

	
	def getCanonicalOrNonCanonicalOrigin(self, peptides, name='MAPs'):

		only_can = 0
		only_non_can = 0
		
   		small_proteins_can = 0
		big_proteins_can = 0

		small_proteins_non_can = 0
		big_proteins_non_can = 0

		get_only_cannonique = {}
		get_only_non_cannonique = {}
		
		small_proteins_non_cannonique = {}
		big_proteins_non_cannonique = {}

		small_cannonique_proteins = {}
		big_canonique_proteins = {}

		small_novel_isoform_proteins = {}
		big_novel_isoform_proteins = {}

		summary_origin_peptide = {}
		set_peptides = set()

		start_codon_canonical = {}
		start_codon_non_canonical = {}

		tpm_source_can = {}
		tpm_source_non_can = {}

		tpm_source_can_big = {}
		tpm_source_non_can_big = {}

		tpm_source_can_small = {}
		tpm_source_non_can_small = {}

		list_prot_to_remove_can = set()
		list_prot_to_remove_non_can = set()
		cont = 0
		print 'Total Peptides ', len(peptides)
		

		for peptide, info in peptides.items():
			cannonique = info[0]
			non_cannonique = info[1]
			set_peptides.add(peptide)

			len_prot = 0
			trans_origin = '' 
				
			if non_cannonique == None:
				summary_origin_peptide[peptide] = 'Canonical'
				len_prot = cannonique[0][4]
				trans_origin = cannonique[0][0]
				codon_canonical = cannonique[0][11]
				index_prot = cannonique[0][10]
				tpm_can = cannonique[0][12]
				list_prot_to_remove_can.add(index_prot)	

				only_can += 1
				if len_prot <= 100:
					tpm_source_can_small[index_prot] = tpm_can
					small_proteins_can += 1
					small_cannonique_proteins[peptide] = [cannonique[0]]
				else:
					tpm_source_can_big[index_prot] = tpm_can
					big_proteins_can += 1
					big_canonique_proteins[peptide] = [cannonique[0]]
				
				try:
					start_codon_canonical[codon_canonical] += 1
				except:
					start_codon_canonical[codon_canonical] = 1

				tpm_source_can[index_prot] = tpm_can
				get_only_cannonique[peptide] = [cannonique[0]]

			elif cannonique == None:
				only_non_can += 1
				summary_origin_peptide[peptide] = 'Non_Canonical'
				len_prot = non_cannonique[0][4]
				codon_non_canonical = non_cannonique[0][11]
				index_prot = non_cannonique[0][10]
				tpm_non_can = non_cannonique[0][12]
				list_prot_to_remove_non_can.add(index_prot)	

				if len_prot <= 100:
					tpm_source_non_can_small[index_prot] = tpm_non_can
					small_proteins_non_can += 1
					small_proteins_non_cannonique[peptide] = [non_cannonique[0]]
				else:
					tpm_source_non_can_big[index_prot] = tpm_non_can
					big_proteins_non_can += 1
					big_proteins_non_cannonique[peptide] = [non_cannonique[0]]
				
				try:
					start_codon_non_canonical[codon_non_canonical] += 1
				except:
					start_codon_non_canonical[codon_non_canonical] = 1

				tpm_source_non_can[index_prot] = tpm_non_can
				get_only_non_cannonique[peptide] = [non_cannonique[0]]
				
			elif non_cannonique != None and cannonique != None:
				tpm_non_can = non_cannonique[0][12]
				proba_non_can = non_cannonique[0][7]
				sc_non_canonical = non_cannonique[0][3]
				kozak_non_can = non_cannonique[0][14]

				sc_canonical = cannonique[0][3]
				tpm_can = cannonique[0][12]
				proba_can = cannonique[0][7]
				codon_canonical = cannonique[0][11]
				kozak_can = cannonique[0][14]

				info_high_tpm = []
				highest_tpm = 0.0
				for info_can in cannonique:
					strand_can = info_can[2]
					trans_name_search = info_can[0]
					tis_can = info_can[3]
					info_transcript = self.annotationTranscripts[strand_can][trans_name_search]	
					tis_transcript = info_transcript['Info'][0]+':'+str(info_transcript['Info'][4][0][0])+'-'+str(info_transcript['Info'][4][0][1])+'-'+info_transcript['Info'][11]
					if tis_can == tis_transcript:
						tpm_can = info_can[12]
						if 	tpm_can > highest_tpm:
							highest_tpm = tpm_can
							info_high_tpm = info_can
					
				if highest_tpm > 0 :
					only_can += 1
					summary_origin_peptide[peptide] = 'Canonical'
					len_prot = info_high_tpm[4]
					
					index_prot = info_high_tpm[10]
					list_prot_to_remove_can.add(index_prot)	

					if len_prot <= 100:
						tpm_source_can_small[index_prot] = highest_tpm
						small_proteins_can += 1
						small_cannonique_proteins[peptide] = [info_high_tpm]
					else:
						tpm_source_can_big[index_prot] = highest_tpm
						big_proteins_can += 1
						big_canonique_proteins[peptide] = [info_high_tpm]
					
					try:
						start_codon_canonical[codon_canonical] += 1
					except:
						start_codon_canonical[codon_canonical] = 1

					tpm_source_can[index_prot] =  info_high_tpm[12]
					get_only_cannonique[peptide] = [info_high_tpm]

				else:
					only_non_can += 1
					summary_origin_peptide[peptide] = 'Non_Canonical'
					len_prot = non_cannonique[0][4]
					codon_non_canonical = non_cannonique[0][11]
					index_prot = non_cannonique[0][10]
					tpm_non_can = non_cannonique[0][12]
					
					list_prot_to_remove_non_can.add(index_prot)	

					if len_prot <= 100:
						tpm_source_non_can_small[index_prot] = tpm_non_can
						small_proteins_non_can += 1
						small_proteins_non_cannonique[peptide] = [non_cannonique[0]] 
					else:
						tpm_source_non_can_big[index_prot] = tpm_non_can
						big_proteins_non_can += 1
						big_proteins_non_cannonique[peptide] = [non_cannonique[0]]

					try:
						start_codon_non_canonical[codon_non_canonical] += 1
					except:
						start_codon_non_canonical[codon_non_canonical] = 1

					tpm_source_non_can[index_prot] = tpm_non_can
					get_only_non_cannonique[peptide] = [non_cannonique[0]]

					
			else:
				print 'Peptide with non-source protein either canonical or non canonical, Verify! ', peptide
			cont += 1

		with open(self.toSaveFile+'/summary_origin_peptide.dic', 'wb') as handle:
			pickle.dump(summary_origin_peptide, handle, protocol=pickle.HIGHEST_PROTOCOL)

		print 'Total Canonical Proteins ', only_can
		print 'small Proteins Canonical ',small_proteins_can, ' big Proteins Canonical ',big_proteins_can
		print 'Start codons Canonical', start_codon_canonical
		print 
		print 'Total Non Canonical Proteins ', only_non_can
		print 'small Proteins Non Canonical ',small_proteins_non_can, ' big Proteins Non Canonical ' , big_proteins_non_can
		print 'Start codons Non Canonical', start_codon_non_canonical
		print 
		
		total_small = small_proteins_can+small_proteins_non_can
		total_big = big_proteins_can+big_proteins_non_can
		print 'total Small Proteins ', total_small
		print 'total Big Proteins ', total_big
		total_maps = only_can + only_non_can 
		print 'Total Maps ', total_maps

		with open(self.toSaveFile+'/getOnlyNonCanonical_'+name+'.dic', 'wb') as handle:
			pickle.dump(get_only_non_cannonique, handle, protocol=pickle.HIGHEST_PROTOCOL)

		with open(self.toSaveFile+'/getOnlyCanonical_'+name+'.dic', 'wb') as handle:
			pickle.dump(get_only_cannonique, handle, protocol=pickle.HIGHEST_PROTOCOL)

		with open(self.toSaveFile+'/small_proteins_non_canonical_'+name+'.dic', 'wb') as handle:
			pickle.dump(small_proteins_non_cannonique, handle, protocol=pickle.HIGHEST_PROTOCOL)

		with open(self.toSaveFile+'/small_canonical_proteins_'+name+'.dic', 'wb') as handle:
			pickle.dump(small_cannonique_proteins, handle, protocol=pickle.HIGHEST_PROTOCOL)

		with open(self.toSaveFile+'/big_proteins_non_canonical_'+name+'.dic', 'wb') as handle:
			pickle.dump(big_proteins_non_cannonique, handle, protocol=pickle.HIGHEST_PROTOCOL)

		with open(self.toSaveFile+'/big_canonical_proteins_'+name+'.dic', 'wb') as handle:
			pickle.dump(big_canonique_proteins, handle, protocol=pickle.HIGHEST_PROTOCOL)

		return only_can, only_non_can, small_proteins_can, small_proteins_non_can, big_proteins_can, big_proteins_non_can, summary_origin_peptide, tpm_source_can, tpm_source_non_can, list_prot_to_remove_non_can, list_prot_to_remove_can, tpm_source_non_can_big, tpm_source_can_big, tpm_source_non_can_small, tpm_source_can_small, get_only_non_cannonique, get_only_cannonique 

	
	def getBioTypeAlternTis(self, peptides, name):

		toWriteBackward = ''
		toWriteForward = ''
		countPeptide = 0
		bio_to_search = 0
		open_prot_proteins = {}
		proteins_by_tis_codon = {}

		unique_proteins = {}
		small_proteins = {}
		big_proteins = {}
		peptides_in_small_prot = 0
		peptides_in_big_prot = 0
		high_number_prot = 0
		key_high_number_prot = ''
		peptides_high_number_prot = []

		keys = {}

		for peptide, info in peptides.items():
			countPeptide += 1
			transcripts_info_to_search = []
			first_proba = info[0][7]
			first_kozakScore = info[0][13]
			first_tpm = info[0][12]
			
			for ind, value in enumerate(info):
				proba_aux = value[7]
				kozakScore_aux = value[13]
				tpm = value[12]
				if (proba_aux == first_proba) and (first_kozakScore == kozakScore_aux) : 
					transcripts_info_to_search.append(value)

			final_bio = ''
			
			ordered = sorted(transcripts_info_to_search, key=lambda e:(e[0]))
			for trans in ordered:
				transcript = trans[0]
				strand = trans[2]
				places = trans[15][0][1]
				proba = trans[7]
				kozakScore = trans[14]
				start_codon = trans[3]
				len_proteine = trans[4]
				rangeTranscript = trans[6] 
				codon = trans[11]
				tpm = trans[12]
				sequence =  trans[13]
				peptide_score =  trans[16]
				indexProtein = trans[10]

				key_unique_prot = transcript+'-'+start_codon

				try:
					l = keys[key_unique_prot][0]
					l.add(indexProtein)
					p = keys[key_unique_prot][1]
					p.append(peptide)
				except KeyError:
					prot = set()
					prot.add(indexProtein)
					keys[key_unique_prot] = [prot, [peptide]]

				bio = ''
				chr = start_codon.split(':')[0]
				if 'chr' in chr:
					chr = chr.split('chr')[1]
					if chr == 'M':
						chr = 'MT'

				try:
					infoTranscript = self.annotationTranscripts[strand][transcript]
					origin = self.searchOriginPlacesPeptide(places, infoTranscript, transcript, strand, peptide)
					bio = origin[0][1][0]

					if bio == 'Alt-TIS' or bio == 'Frameshift':
						final_bio = self.getOrigin(origin, places, infoTranscript, start_codon, strand, peptide, transcript, proba)
						break
				except KeyError:
					if len(ordered) == 1:
						final_bio = 'ToSearch'
			
			if final_bio == '':
				infoTranscript = self.annotationTranscripts[strand][transcript]
				origin = self.searchOriginPlacesPeptide(places, infoTranscript, transcript, strand, peptide)
				bio = origin[0][1][0]
				final_bio = self.getOrigin(origin, places, infoTranscript, start_codon, strand, peptide, transcript, proba)

			if final_bio == 'ToSearch':
				bio_to_search += 1
			
			try:
				info_unique_prot = unique_proteins[key_unique_prot][1]
				info_unique_prot.append(peptide)
				if len(info_unique_prot) > high_number_prot:
					high_number_prot = len(info_unique_prot)
					key_high_number_prot = key_unique_prot
					peptides_high_number_prot = info_unique_prot
			except KeyError:
				unique_proteins[key_unique_prot] = [sequence, [peptide]]

			if len_proteine <= 100:
				peptides_in_small_prot += 1
				try: 
					info_small_prote = small_proteins[key_unique_prot][1]
					info_small_prote.append(peptide)
				except KeyError:
					small_proteins[key_unique_prot] = [sequence, [peptide]]
			else:
				peptides_in_big_prot += 1
				try: 
					info_big_prote = big_proteins[key_unique_prot][1]
					info_big_prote.append(peptide)
				except KeyError:
					big_proteins[key_unique_prot] = [sequence, [peptide]]

			try:
				proteins_by_tis_codon[codon] += 1
			except KeyError:
				proteins_by_tis_codon[codon] = 1

			
			if chr == 'MT':
				chr = 'M'

			if '|' in places:
				splitPlaces = places.split('|')
				for place in splitPlaces:
					if strand == '+' :
						toWriteForward += self.getStringFromPos(place, 'chr'+chr, countPeptide, peptide, transcript, rangeTranscript, start_codon, codon, len_proteine, tpm, final_bio, sequence, peptide_score, indexProtein)
					else:
						toWriteBackward += self.getStringFromPos(place, 'chr'+chr, countPeptide, peptide, transcript, rangeTranscript, start_codon, codon, len_proteine, tpm, final_bio, sequence, peptide_score, indexProtein)
			else:
				if strand == '+' :
					toWriteForward += self.getStringFromPos(places, 'chr'+chr, countPeptide, peptide, transcript, rangeTranscript, start_codon, codon, len_proteine, tpm, final_bio, sequence, peptide_score, indexProtein)
				else:
					toWriteBackward += self.getStringFromPos(places, 'chr'+chr, countPeptide, peptide, transcript, rangeTranscript, start_codon, codon, len_proteine, tpm, final_bio, sequence, peptide_score, indexProtein)
			

		print 'Categories Proteins ', self.categories
		
		fileBed = self.toSaveFile + '/toIntersectBackward_'+name+'.bed'
		fileToSave = open(fileBed, 'w')
		fileToSave.write(toWriteBackward)
		fileToSave.close()
		
		fileBed = self.toSaveFile + '/toIntersectForward_'+name+'.bed'
		fileToSave = open(fileBed, 'w')
		fileToSave.write(toWriteForward)
		fileToSave.close()

		print 'Total Peptides ', len(peptides)
		print 'Unique Proteins ', len(unique_proteins)
		print 
		print 'Total small proteins ', len(small_proteins), ' total peptides ', peptides_in_small_prot
		print 'Total big proteins ', len(big_proteins), ' total peptides ', peptides_in_big_prot
		print 'Start codons origin ', proteins_by_tis_codon

		print 'high_number_prot : ', high_number_prot 
		print 'key_high_number_prot : ', key_high_number_prot 
		print 'peptides_high_number_prot : ', peptides_high_number_prot

		dif = 0
		for key, indexes in keys.items():
			if len(indexes[0]) > 1:
				dif += 1
		print 'Len Keys ',len(keys), ' Dif ', dif

		return self.categories, unique_proteins

	def getStringFromPos(self, place, chr, bio_to_search, peptide, transcript, rangeTranscript, start_codon, codon, len_proteine, tpm, bio, sequence, peptide_score, index_proteine):
		start = place.split('-')[0]
		end = place.split('-')[1]
		string = chr+'\t'+start+'\t'+end+'\t'+str(bio_to_search)+'\t'+peptide+'\t'+transcript+'\t'+ rangeTranscript+'\t'+start_codon+'\t'+codon+'\t'+str(len_proteine)+'\t'+str(tpm)+'\t'+bio+'\t'+sequence +'\t'+str(peptide_score)+'\t'+str(index_proteine)+'\n'
		return string
	
	def getBiotypeInTranscript(self, infoTranscript, ini, fini, transcript, strand, peptide):

		toReturn = ''
		if len(infoTranscript['CDS']) > 0:
			for cds in infoTranscript['CDS']:
				if cds[0] <= ini <= cds[1] and cds[0] <= fini <= cds[1]:
					toReturn = 'In'
					break

			if toReturn == 'In':
				chr = infoTranscript['Info'][0]
				regions = infoTranscript['CDS']
				protein = self.getTranscript_and_protein( chr, regions, transcript, strand)
				for prot in protein:
					if peptide in prot:
						toReturn = 'Alt-TIS', infoTranscript['Info'][3]
						return toReturn
				toReturn = 'Frameshift', infoTranscript['Info'][3]
				return toReturn
			
		if len(infoTranscript['3UTR']) > 0:
			for utr in infoTranscript['3UTR']:
				if utr[0] <= ini <= utr[1] and utr[0] <= fini <= utr[1]:
					toReturn = '3UTR', infoTranscript['Info'][3]
					return toReturn
				elif utr[0] <= ini <= utr[1] and not (utr[0] <= fini <= utr[1]):
					toReturn = '3UTR', infoTranscript['Info'][3]
					return toReturn
					
		if len(infoTranscript['5UTR']) > 0:
			for utr in infoTranscript['5UTR']:
				if utr[0] <= ini <= utr[1] and utr[0] <= fini <= utr[1]:
					toReturn = '5UTR', infoTranscript['Info'][3]
					return toReturn
				elif utr[0] <= ini <= utr[1] and not (utr[0] <= fini <= utr[1]):
					toReturn = '5UTR', infoTranscript['Info'][3]
					return toReturn
		
		if len(infoTranscript['Introns']) > 0:
			for intron in infoTranscript['Introns']:
				if intron[0] <= ini <= intron[1] and intron[0] <= fini <= intron[1]:
					toReturn = 'Intron', infoTranscript['Info'][3]
					return toReturn
				elif intron[0] <= ini <= intron[1] and not (intron[0] <= fini <= intron[1]):
					toReturn = 'Intron', infoTranscript['Info'][3]
					return toReturn

		if len(infoTranscript['Exons']) > 0:
			for exon in infoTranscript['Exons']:
				if exon[0] <= ini <= exon[1] and exon[0] <= fini <= exon[1]:
					toReturn = infoTranscript['Info'][3], infoTranscript['Info'][3]
					return toReturn
				elif exon[0] <= ini <= exon[1] and not (exon[0] <= fini <= exon[1]):
					toReturn = infoTranscript['Info'][3], infoTranscript['Info'][3]
					return toReturn

		return 'Intergenic',


	def getTranscript_and_protein(self, chr, regions, transcript, strandTranscript):

		sequenceTranscript = ''
		ntd_sequence = self.get_MCS.getTranscriptsRegionsInformation(chr, regions, strandTranscript, '', True)
		
		if chr == 'chrM':
			proteineT = translateDNA(ntd_sequence, frame = 'f1', translTable_id='mt')
		else:
			proteineT = translateDNA(ntd_sequence, frame = 'f1', translTable_id='default')


		if '/' in proteineT:
			indexesPoly = [j for j,val in enumerate(proteineT) if val=='/']
			variants = self.get_variants(proteineT)
			return variants
		else:
			return [proteineT]

	def getInfoPosInProtein(self, proteinsList, infoTranscript, peptide, transcript):
		for prot in proteinsList:
			if peptide in prot:
				bio = ''
				info_trans = infoTranscript['Info']
				codon = info_trans[-1]
				start_codon = info_trans[4]
				if  len(start_codon) == 0:
					bio = 'Alt-TIS'
				else:
					bio = 'Canonical'
				return True, bio, infoTranscript['Info'], prot
				break
		return [False]


	def get_variants(self, protein, thresold=128):
		variants = []
		indexesPoly = [j for j,val in enumerate(protein) if val=='/']
		sequenceFirstProteine, listVariants, listPositionsVariants = self.getListVariants(protein, indexesPoly)
		c = list(itertools.product(*listVariants))
		indexes = [j for j,val in enumerate(sequenceFirstProteine) if val=='/']

		for i, pro in enumerate(c):
			copySmallSeqVariantFinal = list(sequenceFirstProteine)
			for index in range(0, len(indexes)):
				copySmallSeqVariantFinal[indexes[index]] = pro[index]

			proteine = "".join(copySmallSeqVariantFinal).split('*')[0]
			variants.append(proteine)
			if i == 128:
				break
		return variants


	def getListVariants(self, proteine, indexes):
		
		prev = indexes[0]
		listPositionsVariants = [[prev]]
		fVariant = set()
		fVariant.add(proteine[prev-1])
		fVariant.add(proteine[prev+1])
		listVariants = [fVariant]
		sequenceFirstProteine = ''
		
		for i, n in enumerate(indexes):
			if i > 0 :
				if n - 2 == indexes[i-1]:
					l = listPositionsVariants[-1]
					l.append(n)
					l2 = listVariants[-1]
					l2.add(proteine[n-1])
					l2.add(proteine[n+1])
				else:
					listPositionsVariants.append([n])
					listVariants.append(set())
					l2 = listVariants[-1]
					l2.add(proteine[n-1])
					l2.add(proteine[n+1])
		
		variant = False
		for i in range(0, len(proteine)-1):
			if i == 0 and proteine[i+1] =='/':
				sequenceFirstProteine += '/'
				variant = True
			elif i == 0 and proteine[i+1] !='/':
				sequenceFirstProteine += proteine[i]
			elif i > 0:
				if not variant and proteine[i+1] != '/':
					sequenceFirstProteine += proteine[i]
				elif variant and proteine[i+1] != '/' and proteine[i] != '/':
					variant = False
				elif not variant and proteine[i+1] == '/':
					variant = True
					sequenceFirstProteine += '/'
		if proteine[-2] != '/':
			sequenceFirstProteine += proteine[-1]

		return sequenceFirstProteine, listVariants, listPositionsVariants


	def searchOriginPlacesPeptide(self, places, infoTranscript, transcript, strand, peptide):
		origin = []
		if '|' in places:
			splitPlaces = places.split('|')
			for place in splitPlaces: 
				ini = int(place.split('-')[0])
				fini =  int(place.split('-')[1])
				origin.append((place, self.getBiotypeInTranscript(infoTranscript, ini, fini, transcript, strand, peptide)))
		else:
			ini = int(places.split('-')[0])
			fini =  int(places.split('-')[1])
			origin.append((places, self.getBiotypeInTranscript(infoTranscript, ini, fini, transcript, strand, peptide)))
		return origin

	def originPlacesAlt_TIS_Peptide(self, places):
		origin = []
		if '|' in places:
			splitPlaces = places.split('|')
			for place in splitPlaces: # chr1:28214940-28214963
				ini = int(place.split('-')[0])
				fini =  int(place.split('-')[1])
				origin.append((place, 'Alt-TIS'))
		else:
			ini = int(places.split('-')[0])
			fini =  int(places.split('-')[1])
			origin.append((places, 'Alt-TIS'))
		
		return origin

	def getOrigin(self, origin, places, infoTranscript, start_codon, strand, peptide, transcript, proba):
		if len(origin) == 1:
			type_ = origin[0][1][0]
			bio = type_
			try:
				self.categories[type_] += 1
				self.categories_peptides[type_].append((peptide, start_codon, strand, places, origin))
			except KeyError:
				self.categories[type_] = 1
				self.categories_peptides[type_] = [(peptide, start_codon, strand, places, origin)]

			try:
				infoPeptide = self.bioTypePeptide[peptide] 
				infoPeptide.append([places, transcript, strand, type_, start_codon, proba])
			except KeyError:
				self.bioTypePeptide[peptide] = [[places, transcript, strand, type_, start_codon, proba]]

		else:
			l = len(origin)
			v = 1/(l*1.0)
			nameOrigin = set()
			for o in origin:
				type_ = o[1][0]
				nameOrigin.add(type_)

			originString = ''
			if len(nameOrigin) == 1:
				originString = list(nameOrigin)[0]
				try:
					self.categories[type_] += 1
					self.categories_peptides[type_].append((peptide, start_codon, strand, places, origin))
				except KeyError:
					self.categories[type_] = 1
					self.categories_peptides[type_] = [(peptide, start_codon, strand, places, origin)]

				try:
					infoPeptide = self.bioTypePeptide[peptide] 
					infoPeptide.append([places, transcript, strand, type_, start_codon, proba])
				except KeyError:
					self.bioTypePeptide[peptide] = [[places, transcript, strand, type_, start_codon, proba]]
			else:
				for name in nameOrigin:
					originString += name
					try:
						self.categories[type_] += v
						self.categories_peptides[type_].append((peptide, start_codon, strand, places, origin))
					except KeyError:
						self.categories[type_] = v
						self.categories_peptides[type_] = [(peptide, start_codon, strand, places, origin)]
					try:
						infoPeptide = self.bioTypePeptide[peptide] 
						infoPeptide.append([places, transcript, strand, type_, start_codon, proba])
					except KeyError:
						self.bioTypePeptide[peptide] = [[places, transcript, strand, type_, start_codon, proba]]
			bio = originString
			self.bioTypePeptide[peptide] = [[places, transcript, strand, bio, start_codon, proba]]
		return bio