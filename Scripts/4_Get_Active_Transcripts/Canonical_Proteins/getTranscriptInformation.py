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

	def __init__(self, snps, fastaFile, filepath_index):
		self.snps = snps
		self.faFile = pysam.FastaFile(fastaFile, filepath_index)
		self.chromosomesWithSNPs = {}
		self.dico = {'AG':'R', 'CT':'Y', 'GC':'S', 'AT':'W', 'GT':'K', 'AC':'M', 'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACG':'V', 'ACTG':'N'}
		print 'getTranscriptInformation OK'


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


	def getTranscriptsRegionsInformation(self, infoTranscript):
		
		nameTranscript = infoTranscript[0][0]
		chromosomeTranscript = infoTranscript[0][1]
		startTranscript = infoTranscript[0][2]
		endTranscript = infoTranscript[0][3]
		strandTranscript = infoTranscript[0][4]

		sequenceTranscript = ''				
		sequenceTranscriptPositions = []
				
		regionTranscript = chromosomeTranscript+":"+str(startTranscript)+'-'+str(endTranscript)
		
		for index, exon in enumerate(infoTranscript):
			
			if index != 0 :
				startExon = int(exon[0])
				endExon = int(exon[1])
				
				try:
					chromosome = self.chromosomesWithSNPs[chromosomeTranscript]
					if chromosome == 'No_SNPs' :
						sequence = self.faFile.fetch(chromosomeTranscript,startExon-1,endExon)
					elif chromosome == 'NoInFasta':
						print 'Problem with chromosome not in Fasta, verify!'
					else:
						sequence = chromosome[startExon-1:endExon]
				except KeyError :
					chromosome = self.transformRegionsWithSNPs(chromosomeTranscript)
					self.chromosomesWithSNPs[chromosomeTranscript] =  chromosome
					if chromosome == 'No_SNPs' or chromosome == 'NoInFasta':
						if chromosome == 'No_SNPs':
							sequence = self.faFile.fetch(chromosomeTranscript,startExon-1,endExon)
						elif chromosome == 'NoInFasta':
							print 'Problem with chromosome not in Fasta, verify!'
					else:
						sequence = chromosome[startExon-1:endExon]	
				
				if strandTranscript == '-' :
					sequence = reverseComplement(sequence)
					sequenceTranscript += sequence
					rangePositions = range(startExon, endExon + 1)[::-1]
					sequenceTranscriptPositions += rangePositions
					
				elif strandTranscript == '+':
					sequenceTranscript += sequence
					rangePositions = range(startExon, endExon + 1)
					sequenceTranscriptPositions.extend(rangePositions)

			
		return [regionTranscript, sequenceTranscript, sequenceTranscriptPositions, nameTranscript, chromosomeTranscript]

	def getTranscriptsRegionsInformation(self, range_transcript, strand):
		
		chromosomeTranscript = range_transcript.split(':')[0]
		exones = range_transcript.split(':')[1]
		
		sequenceTranscript = ''				
		
		if '|' in exones:
			exones = exones.split('|')
		else:
			exones = [exones]

		for index, exon in enumerate(exones):
			split_exon = exon.split('-')
			startExon = int(split_exon[0])
			endExon = int(split_exon[1])
			
			try:
				chromosome = self.chromosomesWithSNPs[chromosomeTranscript]
				if chromosome == 'No_SNPs' :
					sequence = self.faFile.fetch(chromosomeTranscript,startExon-1,endExon)
				elif chromosome == 'NoInFasta':
					print 'Problem with chromosome not in Fasta, verify!'
				else:
					sequence = chromosome[startExon-1:endExon]
			except KeyError :
				chromosome = self.transformRegionsWithSNPs(chromosomeTranscript)
				self.chromosomesWithSNPs[chromosomeTranscript] =  chromosome
				if chromosome == 'No_SNPs' or chromosome == 'NoInFasta':
					if chromosome == 'No_SNPs':
						sequence = self.faFile.fetch(chromosomeTranscript,startExon-1,endExon)
					elif chromosome == 'NoInFasta':
						print 'Problem with chromosome not in Fasta, verify!'
				else:
					sequence = chromosome[startExon-1:endExon]	
			
			sequenceTranscript += sequence
				
		if strand == '-' :
			sequenceTranscript = reverseComplement(sequenceTranscript)

		return sequenceTranscript
		



