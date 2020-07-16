import warnings
import time
import pysam
from pyGeno.Genome import * 
from pyGeno.tools.UsefulFunctions import *
import numpy as np
import pickle
from itertools import permutations 


class getGenomeWithSNPS:

	def __init__(self, fastaFile, filepath_index, fileSNPsFreeBayes, qualityThreshold):
		self.snps = self.getSNPsDic(fileSNPsFreeBayes, qualityThreshold)
		self.faFile = pysam.FastaFile(fastaFile, filepath_index)
		self.chromosomesWithSNPs = {}
		# IUPAC symbols
		self.dico = {'AG':'R', 'CT':'Y', 'GC':'S', 'AT':'W', 'GT':'K', 'AC':'M', 'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACG':'V', 'ACTG':'N'}
		print 'getTranscriptInformation OK'


	def getSNPsDic(self, fileSNPsFreeBayes, qualityThreshold):

		dic = {}
		totalSnps = 0
		
		with open(fileSNPsFreeBayes) as f:
			for index, snp in enumerate(f):
				if index > 0:
					splitLine = snp.strip().split('\t')
					chr = splitLine[0]
					if len(chr) < 3:
						chr = 'chr'+chr
					quality = float(splitLine[6])
					if quality > qualityThreshold:
						totalSnps += 1
						try:
							lista = dic[chr]
							lista.append(snp)
							dic[chr] = lista
						except KeyError:
							dic[chr] = [snp]
		
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
				print chromosome, ' No SNPs '
				return 'No_SNPs'
		except KeyError:
			print chromosome, ' No present in fasta '
			return 'NoInFasta'
		return "".join(sequence)


	def getStringFromTuple(self, tuple):
		string = ''
		for t in range(len(tuple)):
			string += tuple[t]
		return string


	def getRefRegion(self, chr, ini, fini):
		
		try:
			chromosome = self.chromosomesWithSNPs[chr]
			if chromosome == 'No_SNPs' :
				sequence = self.faFile.fetch(chr,ini,fini)
			elif chromosome == 'NoInFasta':
				print 'Problem with chromosome not in Fasta, to verify!'
			else:
				sequence = self.faFile.fetch(chr,ini,fini)
		except KeyError :
			chromosome = self.transformRegionsWithSNPs(chr)
			self.chromosomesWithSNPs[chr] =  chromosome
			if chromosome == 'No_SNPs' or chromosome == 'NoInFasta':
				if chromosome == 'No_SNPs':
					sequence = self.faFile.fetch(chr,ini,fini)
				elif chromosome == 'NoInFasta':
					print 'Problem with chromosome not in Fasta, to verify!'
			else:
				sequence = chromosome[ini:fini]	
			
		return sequence
		



