import warnings
warnings.filterwarnings("ignore")
import time
import pysam
from pyGeno.Genome import * 
from pyGeno.tools.UsefulFunctions import *
import numpy as np
sys.path += sys.path + ['/u/ruizma/miniconda2/lib/python2.7/site-packages']
import pickle
from itertools import permutations 


class getTranscriptInformation:

	def __init__(self, snps, quality, fastaFile, filepath_index):
		self.snps = self.getSNPsDic(snps, quality)
		self.faFile = pysam.FastaFile(fastaFile, filepath_index)
		self.chromosomesWithSNPs = {}
		self.dico = {'AG':'R', 'CT':'Y', 'GC':'S', 'AT':'W', 'GT':'K', 'AC':'M', 'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACG':'V', 'ACTG':'N'}
		self.dico_reverse = {'A':'A','T':'T', 'C':'C', 'G':'G', 'R':'AG', 'Y':'CT', 'S':'GC', 'W':'AT', 'K':'GT', 'M':'AC', 'B':'CGT', 'D':'AGT', 'H':'ACT', 'V':'ACG', 'N':'ACTG'}
		print 'getTranscriptInformation OK'


	def getSNPsDic(self, snps_file, quality_threshold):

		time0 = time.time()
		dic = {}
		totalSnps = 0
		
		with open(snps_file) as f:
			for index, snp in enumerate(f): 
				if index > 0:
					splitLine = snp.strip().split('\t')
					chr = splitLine[0]
					if len(chr) < 3:
						chr = 'chr'+chr
					quality = float(splitLine[6])
					if quality > quality_threshold:
						totalSnps += 1
						try:
							lista = dic[chr]
							lista.append(snp)
							dic[chr] = lista
						except KeyError:
							dic[chr] = [snp]
		
		print 'Total SNPs Over the quality_threshold ', quality_threshold, ' : ', str(totalSnps)
		timeFinal = time.time()
		total = (timeFinal-time0) / 60
		print  'Total runtime getSNPsDic : '+ str(total) +' min'
		return dic


	def transformRegionsWithSNPs(self, chromosome) :
		try:
			sequence = list(self.faFile.fetch(chromosome)) 
			try:
				for snp in self.snps[chromosome]:
					splitLine = snp.strip().split('\t')
					pos = int(splitLine[2])
					letter = splitLine[4]
					letterToChange = splitLine[5]
					if len(letterToChange) == 1:
						sequence[pos] = letterToChange
					else:
						splitLetterToChange = letterToChange.split(',')
						perm = permutations(splitLetterToChange)
						for i in list(perm):
							s = self.getStringFromTuple(i)
							try:
								letterToChange = self.dico[s]
								sequence[pos] = letterToChange
								break
							except KeyError:
								pass
			except KeyError:
				return 'No_SNPs'
		except KeyError:
			return 'NoInFasta'
		return "".join(sequence)


	def getStringFromTuple(self, tuple):
		string = ''
		for t in range(len(tuple)):
			string += tuple[t]
		return string


	def getTranscriptsRegionsInformation(self, chromosomeTranscript, places, strand, protein_to_compare, equal = False):
		
		sequenceTranscript = ''				
		sequenceTranscriptPositions = []
		
		for index, place in enumerate(places):
			start = place[0]
			end = place[1]
			try:
				chromosome = self.chromosomesWithSNPs[chromosomeTranscript]
				if chromosome == 'No_SNPs' :
					sequence = self.faFile.fetch(chromosomeTranscript,start-1,end)
				elif chromosome == 'NoInFasta':
					print 'Problem with chromosome not in Fasta, verify!'
				else:
					sequence = chromosome[start-1:end]
			except KeyError :
				chromosome = self.transformRegionsWithSNPs(chromosomeTranscript)
				self.chromosomesWithSNPs[chromosomeTranscript] =  chromosome
				if chromosome == 'No_SNPs' or chromosome == 'NoInFasta':
					if chromosome == 'No_SNPs':
						sequence = self.faFile.fetch(chromosomeTranscript,start-1,end)
					elif chromosome == 'NoInFasta':
						print 'Problem with chromosome not in Fasta, verify!'
				else:
					sequence = chromosome[start-1:end]	
			
			if not equal:
				if strand == '-' :
					sequence = reverseComplement(sequence)
					sequenceTranscript = sequence + sequenceTranscript
				else:
					sequenceTranscript =  sequenceTranscript + sequence
			else:
				if strand == '-' :
					sequence = reverseComplement(sequence)
				sequenceTranscript = sequenceTranscript + sequence
			
		return sequenceTranscript


	def get_translations_sequence(self, chr, ntd_sequence, protein_to_compare, get = False):

		sequences = [ntd for ntd in self.dico_reverse[ntd_sequence[0]]]
		for ntd_1 in ntd_sequence[1:]:
			to_extend = sequences
			sequences = []
			for ntd in self.dico_reverse[ntd_1]:
				for sequence in to_extend:
					sequence += ntd
					sequences.append(sequence)
		
		sequences_that_code_for_peptide = []

		for sequence in sequences:
			if chr == 'chrM':
				proteineT = translateDNA(sequence, frame = 'f1', translTable_id='mt')
			else:
				proteineT = translateDNA(sequence, frame = 'f1', translTable_id='default')
		
			if proteineT == protein_to_compare:
				sequences_that_code_for_peptide.append(sequence)
			elif get:
				sequences_that_code_for_peptide.append(sequence)

		if len(sequences_that_code_for_peptide) == 0 and not get:
			print chr, ntd_sequence, proteineT, protein_to_compare

		return sequences_that_code_for_peptide






