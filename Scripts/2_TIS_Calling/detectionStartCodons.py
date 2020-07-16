import time
import sys, getopt
from pyGeno.tools.UsefulFunctions import *
import pysam
import logging
import pickle, math
import regex as rx
import matplotlib.pyplot as plt
import re
import getGenomeSNPs as snps

listCodonsStart = [ 'ATG','CTG', 'TTG', 'GTG', 'ACG', 'ATA', 'ATT', 'ATC', 'AAG', 'AGG'] 
pattern = 'ATG|CTG|TTG|GTG|ACG|ATA|ATT|ATC|AAG|AGG'


def startCodonsByReadsFromSamFile(samFileInput, probasByLenFragment, weightsByLenght, minReads, fastaFile, filepath_index, folderToSave, fileSNPsFreeBayes, qualityThreshold):
	
	logging.info('')
	logging.info('-----------Starting detection start codons-----------')
	time_0 = time.clock()

	logging.info('probasByLenFragment path : %s', probasByLenFragment)
	with open (probasByLenFragment, 'rb') as fp:
		probasByLenFragment = pickle.load(fp)
	logging.info('probasByLenFragment : %s', str(probasByLenFragment))

	logging.info('weightsByLenght path : %s', weightsByLenght)
	with open (weightsByLenght, 'rb') as fp:
		weightsByLenght = pickle.load(fp)
	logging.info('weightsByLenght : %s', str(weightsByLenght))

	totalAlignments = 0 
	scFound = 0
	startCodons = {}
	
	toGetGenomePer = snps.getGenomeWithSNPS(fastaFile, filepath_index, fileSNPsFreeBayes, qualityThreshold)

	with open(samFileInput) as f:
		for index, line in enumerate(f):
			if '@' not in line:
				totalAlignments += 1
				splitLine = line.strip().split('\t')
				queryname = splitLine[0]
				hi = int(splitLine[12].split(':')[2]) # Hit number: aligment order assigned by STAR 
				readKey = queryname+':'+str(hi)
				cigar = splitLine[5]
				chr = splitLine[2]
				readStart = int(splitLine[3])
				strand = setStrandRead(splitLine[1])
				seq = splitLine[9]
				lenSeq = len(seq)

				realEnd, rang, operators, seqReference = getRanges(cigar, readStart, lenSeq, chr, toGetGenomePer)

				if 'D' not in operators and 'I' not in operators:
					seq = seqReference

					if strand == '-':
						seq = reverseComplement(seq)
						rang = rang[::-1]
					
					res = [(m.start(0), m.end(0), m.group()) for m in rx.finditer(pattern, seq, overlapped=True)]
					
					for element in res:
						ini = element[0]
						fini = element[1]
						codon = element[2]
						substring = rang[ini:fini]
						
						if 0 not in substring:
							if strand == '-':
								substring = substring[::-1]

							probaLenghtRead = probasByLenFragment[lenSeq]
							weightLen = weightsByLenght[lenSeq]
							weightMultimapping = 1 - (((hi-1)/(10-1))*0.99)
							weightsTotal = weightLen * weightMultimapping

							if probaLenghtRead[ini] != 0:
								probaPosition = (probaLenghtRead[ini] * weightsTotal)
								infoSC = chr+':'+str(substring[0])+'-'+str(substring[-1])
								try:
									scChr= startCodons[strand]
									try:
										chromo = scChr[chr]
										try:
											info = chromo[infoSC]
											info[3] += probaPosition
											info[4] += weightsTotal
											info[5] += 1
										except KeyError:
											scFound += 1
											chromo[infoSC] = [infoSC, substring, codon, probaPosition, weightsTotal, 1]
									except KeyError:
										scFound += 1
										scChr[chr] = {infoSC:[infoSC, substring, codon, probaPosition, weightsTotal, 1]}
								except KeyError:
									scFound += 1
									startCodons[strand] = {chr:{infoSC: [infoSC, substring, codon, probaPosition, weightsTotal, 1] }}
				
				if totalAlignments % 1000000 == 0:
					logging.info('Alignment number : %d', totalAlignments)
					logging.info('Total start codons detected : %d', scFound)

				
	logging.info('Total alignments in the bamFile : %d', totalAlignments)
	logging.info('Total start codons detected : %d', scFound)

	startCodonsFinal = []	
	for strand in startCodons.keys():
		for chr in startCodons[strand].keys():
			for sc in startCodons[strand][chr]:
				infoSC = startCodons[strand][chr][sc]
				if infoSC[5] >= minReads:
					probaPosition = safe_div(infoSC[3],infoSC[4])
					infoSC.append(probaPosition)
					startCodonsFinal.append((probaPosition, infoSC[0], infoSC[1], infoSC[2], infoSC[5], strand))
					
	sortedList = sorted(startCodonsFinal, key=lambda tup: tup[0], reverse=True ) 

	with open(folderToSave+'/SortedStartCodonsAfterApplyFilterByMinReads.list', 'wb') as fp:
		pickle.dump(sortedList, fp)	
		
	logging.info('Total start codons candidates : %d save in %s ', len(sortedList), folderToSave+'/SortedStartCodonsAfterApplyFilterByMinReads.list')

	with open(folderToSave+'/AllStartCodonsWithoutFilter.dic', 'wb') as fp:
		pickle.dump(startCodons, fp, protocol=pickle.HIGHEST_PROTOCOL)	

	logging.info('Total start without filtering by the min reads save in %s ', folderToSave+'/AllStartCodonsWithoutFilter.dic')

	timeFinal = time.clock()
	total = (timeFinal-time_0) /60
	logging.info('Total time run function getStartCodons End : %f min ', total)

		
def safe_div(x,y):
		if y == 0: return 0
		return x/(y*1.0)


def setStrandRead(strand):
	number = "{0:b}".format(int(strand))
	if len(number)>= 5:
		if '1' == number[-5]:
			return '-'
		else:
			return '+'
	else:
		return '+'
	

def getRanges(cigar, start, lenSeq, chr, toGetGenomePer):

	rang = [0]*lenSeq
	rx = re.findall('(\d+)([MISDNX=])?', cigar)
	realReadStart = start
	realEnd = start
	operators = []
	lastIndex = 0
	seqReference = ''

	for index in rx:
		operation = index[1]
		length = int(index[0])
		operators.append(operation)
		
		if ('S' in operation):
			end = length
			lastIndex += end
			realReadStart = start - length
			realEnd += end
			resSeq = toGetGenomePer.getRefRegion(chr, start-1,start+end-1)
			seqReference += resSeq

		elif ('I' in index) or ('M' in index) or ('=' in index) or ('X' in index) :
			end = length
			resSeq = toGetGenomePer.getRefRegion(chr, start-1,start+end-1)
			seqReference += resSeq

			for i in range(0,end):
				rang[lastIndex+i] = start
				start = start + 1
			
			realEnd += end
			lastIndex = lastIndex + end
			
		elif ('N' in index) or ('D' in index):
			end = start + length
			start = end
			realEnd += length	

	return realEnd-1, rang, operators, seqReference








