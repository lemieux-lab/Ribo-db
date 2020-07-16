import warnings
import time
import sys
sys.path += sys.path + ["/u/ruizma/Jupyter/PipelineAligment/Scripts/PythonScripts", "/u/ruizma/Jupyter/PipelineAligment/Scripts/PythonScripts/General/",'/u/ruizma/lib/rabadb','/u/ruizma/lib/pygeno','/u/ruizma/Jupyter/ExploringRiboSeqReads/Scripts' , '/u/ruizma/miniconda2/lib/python2.7/site-packages']
from pyGeno.Genome import * 
warnings.filterwarnings("ignore")
from pyGeno.tools.UsefulFunctions import *
from pyGeno.tools.BinarySequence import *
import array
import math
import pickle
import itertools
import datetime
from difflib import SequenceMatcher
import logging
import sys, getopt
import bisect
import numpy as np
from collections import Counter
import re
from itertools import groupby
from operator import itemgetter
from collections import OrderedDict
import operator

strongKozak = re.compile('[ATGC]GCC[AG]CC[ATGC]{3}G[ATGC]')
weakKozak = re.compile('[ATGC]{4}[AG][ATGC]{2}[ATGC]{3}G[ATGC]')



class getFromWherePeptidesDerived:

	def __init__(self, fileDBInfoCustom):

		time0 = time.time()
		
		with open (fileDBInfoCustom, 'rb') as fp:
		 	self.fileDBInfoCustom = pickle.load(fp)

		self.proteinsZero = {}
		self.peptidesByProtein = {}
		self.peptidesOrigin = {}
		self.valeurs_probas = {}
		self.start_codons = {}
		self.tpm_protein = {}

		self.getProteinVariantZero()
		
		timeFinal = time.time()
		total = (timeFinal-time0) / 60
		
		print 'getFromWherePeptidesDerived OK'


	def getProteinVariantZero(self):

		self.indexes_proteins_contained = {}

		for variant in self.fileDBInfoCustom: 	#nameProteineCandidat = [transName, sc, frame, proteineId,  tis] --> Non can tis = tis+'-'+strand
												#nameProteineCandidat = [transName, sc, frame, proteineId,  tis] --> Can
												#sc = [tis, codon, strand, proba, posIndex, kozak, positionsTrans, reference_transcript, tpm]  --> Non can
												#sc = [tis, codon, strand, proba, posIndex, kozak, positionsTransFromTis, transcript, tpm] --> Can
												#variant = [(proteine, nameProteineCandidat, contVariant, chr+':'+rangeTranscript, [], lenProteine, lenTranscript, posIndex, (posIndex, posIndex+aaLen+3))]
			proteinsContained = variant[4]
			variantCount = variant[2]
			lenProteine = variant[5]
			strand = variant[1][1][2]
			tis = variant[1][4]
			rangeTranscript = variant[3]
			indexProtein = variant[1][3]
			sequence =  variant[0]
			codon = variant[1][1][1]
			tpm = variant[1][1][8]

			reference_transcript = ''
			kozak = variant[1][1][5]
			reference_transcript = variant[1][1][7]
			
			
			if variantCount == 0:
				
				self.tpm_protein[indexProtein] = tpm
				
				transcriptProtein = variant[1][0]
				
				if 'ENST' in reference_transcript:
					transcriptProtein = reference_transcript

				try:
					self.start_codons[codon] += 1
				except:
					self.start_codons[codon] = 1

				try:
					proteins = self.proteinsZero[strand]
					proteins[indexProtein] = [transcriptProtein, lenProteine, tis, rangeTranscript, kozak, indexProtein, sequence, codon]
					
				except KeyError:
					self.proteinsZero[strand] = {indexProtein: [transcriptProtein, lenProteine, tis, rangeTranscript, kozak, indexProtein, sequence, codon]}
				
			self.indexes_proteins_contained[indexProtein] = [indexProtein]
			if len(proteinsContained) > 0:
				for variant in proteinsContained:
					self.getVariantZeroInfo(variant, indexProtein)


	def getVariantZeroInfo(self, contained, index_protein_mother):

		variantNumber = contained[2]
		indexProtein = contained[1][3]
		lenProteine = contained[5]
		transcriptProtein = contained[1][0]
		strand = contained[1][1][2]
		tis = contained[1][4]
		rangeTranscript = contained[3]
		sequence =  contained[0]
		codon = contained[1][1][1]
		tpm = contained[1][1][8]
		self.indexes_proteins_contained[index_protein_mother].append(indexProtein)
		
		
		reference_transcript = ''
		kozak = contained[1][1][5]
		reference_transcript = contained[1][1][7]

		if 'ENST' in reference_transcript:
			transcriptProtein = reference_transcript
		
		if variantNumber == 0:
			self.tpm_protein[indexProtein] = tpm
			
			try:
				self.start_codons[codon] += 1
			except:
				self.start_codons[codon] = 1

			try:
				proteins = self.proteinsZero[strand]
				proteins[indexProtein] = [transcriptProtein, lenProteine, tis, rangeTranscript, kozak, index_protein_mother, sequence, codon]
			except KeyError:
				self.proteinsZero[strand] = {indexProtein: [transcriptProtein, lenProteine, tis, rangeTranscript, kozak, index_protein_mother, sequence, codon]}
				
		if len(contained[4]) == 0:
			return 
		else:
			for variant in contained[4]:
				self.getVariantZeroInfo(variant, index_protein_mother)


	def getInfoPeptide(self, peptide, indexes_proteines, score):

		proteinesSource = []
		for indexToGetInfo in indexes_proteines:
			listStartCodonInfo = self.getStartCodonInformation(peptide, indexToGetInfo)
			proteinesSource.extend(listStartCodonInfo)
		
		# (transcript, variantNumber, strand, tis, lenProteine, index, rangeTranscript, proba, kozak, index_protein_mother, indexProtein)
		for idx, sc in enumerate(proteinesSource):
			listaPositions = []
			strand = sc[2]
			transcript = sc[0]
			rangeTranscript = sc[6]
			len_peptide = len(peptide)
			index_pos_peptide = sc[5]

			pos = self.getPositionPeptide(rangeTranscript, index_pos_peptide, len_peptide, strand, peptide)
			
			if (strand, pos) not in listaPositions:
				listaPositions.append((strand, pos))

			valueKozak = 0

			kozak = sc[8]
			try:
				m = strongKozak.match(kozak)
				if m.group():
					valueKozak = 1
			except AttributeError:
				try:
					m = weakKozak.match(kozak)
					if m.group():
						valueKozak = 0.5
				except AttributeError:
					valueKozak = 0
			proteinesSource[idx] = sc + (valueKozak, listaPositions, score,)

		proteins_final = []
		already_in = []

		for origin in proteinesSource:
			trans = origin[0]
			tis = origin[3]
			position = origin[10]
			key = trans +'-'+tis+'-'+str(position)
			if key not in already_in:
				already_in.append(key)
				proteins_final.append(origin)

		# (0:transcript, 1:variantNumber, 2:strand, 3:tis, 4:lenProteine, 5:index, 6:rangeTranscript, 7:proba, 8:kozak, 9:index_protein_mother, 10:indexProtein, 
		# 11:codon, 12:tpm, 13:sequence, 14:valueKozak, 15:listaPositions, 16:Score)

		if len(proteins_final) > 0:
			ordered = sorted(proteins_final, key=lambda e:(e[7], e[14], e[12]), reverse=True)
			name_prot = ordered[0][9]
			try:
				self.peptidesByProtein[name_prot] += 1
			except KeyError:
				self.peptidesByProtein[name_prot] = 1
			return ordered
		else:
			return proteins_final


	def getInfoPeptides(self, peptides):

		time0 = time.time()
		self.peptidesByProtein = {}
		self.peptidesOrigin = {}
		
		proteins = {}
		peptide_origin = {}
		origins_decided = {}
		max_number_peptides_by_protein = {}
		
		proteins, peptide_origin, max_number_peptides_by_protein = self.lookup_proteins(peptides)

		
		for peptide, proteins_source in peptide_origin.items():
			indexes_proteins_source = [index[0] for index in proteins_source[0]]
			peptide_score = proteins_source[1]
			peptide_aux = proteins_source[2]
			listaPositions = self.getInfoPeptide(peptide, indexes_proteins_source, peptide_score)
			if len(listaPositions) == 0 and peptide[0] == 'M':
				listaPositions = self.getInfoPeptide(peptide_aux, indexes_proteins_source, peptide_score)
			self.peptidesOrigin[peptide] = listaPositions
			

		timeFinal = time.time()
		total = (timeFinal-time0) / 60
		print 'Total runtime getInfoPeptides minutes ', total
	

	def lookup_proteins(self, peptides):
		# Return all the proteins where the peptide have been found
		proteins = {}
		peptide_origin = {}
		origins_decided = {}
		max_number_peptides_by_protein = {}
		totalPeptides = 0

		peptides_aux = {}
		if type(peptides) is list:
			for peptide in peptides:
				peptides_aux[peptide] = [peptide, '-']
		else:
			peptides_aux = peptides

		for peptide, info_peptide in peptides_aux.items():
			totalPeptides += 1
			peptide_score =  info_peptide[1]
			if peptide[0] == 'M':
				peptide_aux = peptide[1:]	
			else:
				peptide_aux = peptide

			for index, protein_info in enumerate(self.fileDBInfoCustom):
				protein = protein_info[0]
				len_prot = protein_info[5]

				if peptide in protein or peptide_aux in protein:
					inside = True
					try:
						info = proteins[index]
						info.append(peptide)
						max_number_peptides_by_protein[index] += 1
					except KeyError:
						proteins[index] = [peptide]
						max_number_peptides_by_protein[index] = 1
					try:
						info = peptide_origin[peptide][0]
						info.append((index, len_prot))
					except KeyError:
						peptide_origin[peptide] = [[(index, len_prot)], peptide_score, peptide_aux]

		print 'Total Peptides evaluated ', totalPeptides
		return proteins, peptide_origin, max_number_peptides_by_protein


	def getStartCodonInformation(self, peptideSeq, indexProt):

		completeInfoProtein_peptide = self.fileDBInfoCustom[indexProt]
		listStartCodonInfoAux = self.getDeepInfo(peptideSeq, completeInfoProtein_peptide)
		return 	listStartCodonInfoAux


	def getDeepInfo(self, peptide, lista):
		listStartCodonInfo = []
		sequence = lista[0].strip()
		
		if peptide in sequence :
												#nameProteineCandidat = [transName, sc, frame, proteineId,  tis] --> Non can tis = tis+'-'+strand
												#nameProteineCandidat = [transName, sc, frame, proteineId,  tis] --> Can
												#sc = [tis, codon, strand, proba, posIndex, kozak, positionsTrans, reference_transcript, tpm]  --> Non can
												#sc = [tis, codon, strand, proba, posIndex, kozak, positionsTransFromTis, transcript, tpm] --> Can
												# [(proteine, nameProteineCandidat, contVariant, chr+':'+rangeTranscript, [] , lenProteine, lenTranscript,  posIndex, (posIndexV, posIndexV+aaLen))]
			
			nameProteineCandidat = lista[1]
			sc = nameProteineCandidat[1]
			transcript = nameProteineCandidat[0]
			variantNumber = lista[2]
			index = lista[0].index(peptide)
			lenProteine = lista[5]
			lenTranscript = lista[6]
			strand =  sc[2]
			tis = nameProteineCandidat[4]
			proba = sc[3]
			indexProtein = nameProteineCandidat[3]
			codon = sc[1]
			tpm = sc[8]
			kozak = sc[5]
			reference_transcript = sc[7]
			
			if variantNumber == 0:
				if index == 0 and (codon != 'CTG' or codon != 'ATG'):
					if peptide[0] != 'M':
						if len(lista[4]) == 0:
							return listStartCodonInfo
						else:
							for nested in lista[4]:
								listAux = self.getDeepInfo(peptide, nested)
								if listAux:
									listStartCodonInfo.extend(listAux)
						return listStartCodonInfo	


			if 'ENST' in reference_transcript:
				transcript = reference_transcript

			try:
				posIndex =  lista[7]
				positions_variant = lista[8] 
			except IndexError:
				pass

			index_protein_mother = indexProtein

			if variantNumber != 0:
				lenProteine = self.proteinsZero[strand][indexProtein][1]
				rangeTranscript =  self.proteinsZero[strand][indexProtein][3]
				transcript =  self.proteinsZero[strand][indexProtein][0]
				index_protein_mother = self.proteinsZero[strand][indexProtein][5]
				sequence = self.proteinsZero[strand][indexProtein][6].strip()
				try:
					index = ((positions_variant[0] - posIndex)/3) + index
				except UnboundLocalError:
					pass		
			else:
				rangeTranscript = lista[3]

			if len(sequence) != lenProteine:
				print 'Different Lens'

			try:
				self.valeurs_probas[proba] += 1
			except KeyError:
				self.valeurs_probas[proba] = 1

			listStartCodonInfo.append((transcript, variantNumber, strand, tis, lenProteine, index, rangeTranscript, proba, kozak, index_protein_mother, indexProtein, codon, tpm, sequence))
			
			if len(lista[4]) == 0:
				return listStartCodonInfo
			else:
				for nested in lista[4]:
					listAux = self.getDeepInfo(peptide, nested)
					if listAux:
						listStartCodonInfo.extend(listAux)
			return listStartCodonInfo
		else:
			return listStartCodonInfo


	def getPositionPeptide(self, rangeTranscript, index_pos_peptide, len_peptide, strand, peptide):
		transcript_positions = self.decode_positions(rangeTranscript, strand)
		index_pos_peptide = index_pos_peptide * 3 
		endPos = index_pos_peptide + (len_peptide*3)
		posFinal = transcript_positions[index_pos_peptide:endPos]
		return self.find_ranges(posFinal)
		

	def decode_positions(self, positionsTrans, strand):
		postitions_list = []

		split_positions = positionsTrans.strip().split('|')
		for index, split in enumerate(split_positions):
			if index  == 0:
				split_pos = split.split('-')
				start = int(split_pos[0].split(':')[1])
				end = int(split_pos[1])
			else:
				split_pos = split.split('-')
				start = int(split_pos[0])
				end = int(split_pos[1])
			pos = range(start, end+1)
			if strand == '-':
				postitions_list = pos[::-1] + postitions_list
			else:
				postitions_list += pos
		return postitions_list


	def find_ranges(self, iterable):
		ranges = []
		ordered = sorted(iterable)
		toReturn = ''
		for k, g in groupby(enumerate(ordered), lambda (i,x):i-x):
			group = map(itemgetter(1), g)
			ranges.append((group[0], group[-1]))

		for pos in ranges:
			toReturn += str(pos[0])+'-'+str(pos[1])+'|'
			
		toReturn = toReturn[:-1]
		return toReturn

	def get_tpm_all_prot(self):
		import copy
		return copy.deepcopy(self.tpm_protein)



