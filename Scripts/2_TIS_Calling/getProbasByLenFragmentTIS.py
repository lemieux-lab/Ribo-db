import time
import sys, getopt
import re
import pysam
import pickle
import random
import datetime
import logging
import getInfoFromProbas as getInfoProbas
import numpy as np


def getRanges(cigar, start, rangeRead, nameRead, strand, lenSeq):

	rang = [0]*lenSeq
	rx = re.findall('(\d+)([MISDNX=])?', cigar)
	lastIndex = 0

	for index in rx:
		if ('S' in index):
			end = int(index[0])
			lastIndex = lastIndex + end

		elif ('I' in index) or ('M' in index) or ('=' in index) or ('X' in index) :
			end = int(index[0])
			for i in range(0,end):
				rang[lastIndex+i] = start
				start = start + 1
			
			lastIndex = lastIndex + end
			intersect = rangeRead.intersection(set(rang))
			if len(intersect) == 3 :
				if strand == '-':
					orderedRange = rang[::-1]
					orderedIntersect = sorted(intersect)[::-1]
				else:
					orderedRange = rang
					orderedIntersect = sorted(intersect)
				posStart = orderedRange.index(orderedIntersect[0])
				return True, posStart, orderedRange

		elif ('N' in index) or ('D' in index):
			end = start + int(index[0])
			start = end
			
	return False, 0, []



def getProbas(all_starCodonsPositions, starCodonsPositions, bamFile, minReads, toSaveOutPuts):
	
	t0_ = time.time()

	logging.info('get Probas')
	transcripts_scs = {}
	frequenciesByLenFragment = {}
	known_scs_intersected_dic = {}
	known_scs_intersected_list = []
	splited_sc = 0 

	logging.info('Reading file starCodonsPositions')

	with open(starCodonsPositions) as f: 
		# This recovers the start codons annotated only once (A start codon genomic position can be the start
		# position of several isoforms of a gene), including the start codons spliced.
		for index, line in enumerate(f):
			splitLine = line.strip().split('\t')
			chr=splitLine[0].split(':')[0]
			start=splitLine[0].split(':')[1].split('-')[0]
			end=splitLine[0].split(':')[1].split('-')[1]
			sc = splitLine[0]
			strand = splitLine[1]
			transcript = splitLine[2]
			codon = splitLine[3]

			try:
				infoStrand = transcripts_scs[strand]
				try:
					infoTranscript = infoStrand[transcript]
					infoTranscript.append([codon, sc])
					splited_sc += 1
				except KeyError:
					infoStrand[transcript] = [[codon, sc]]
			except KeyError:
				transcripts_scs[strand] = {transcript:[[codon, sc]]}

	print 'Total Start annotated codons spliced ', splited_sc
	totalS_sc_known_with_TIS_evidence = 0
	known_intersected = []
	known_intersected_backward_keys_list = []
	known_intersected_forward_keys_list = []

	for strand, transcripts in transcripts_scs.items():
		for transcript, scs in transcripts.items():
			rangeStartCodon = []
			rows = ''
			codon = ''
			chr = scs[0][1].split(':')[0]
			names = set()
			info_to_add = ''
			key = ''

			for index, sc in enumerate(scs):
				codon += sc[0]
				if len(scs) == index + 1:
					key += sc[1]
				else:
					key += sc[1]+'|'

				if strand == '+':
					rows += pysam.view("-F0X10", bamFile, sc[1])
				elif strand == '-':
					rows += pysam.view("-f0X10", bamFile, sc[1])
				
				start=int(sc[1].split(':')[1].split('-')[0])
				end=int(sc[1].split(':')[1].split('-')[1])
				
				if start == end:
					rangeStartCodon.append(start)
				else:
					rangeStartCodon.extend(range(start,end+1))

			sorted_range_start_codon = np.sort(rangeStartCodon)
			sc_string = chr+':'+str(sorted_range_start_codon[0])+'-'+str(sorted_range_start_codon[-1])+'-'+strand+'-'+codon
			splitRows = rows.strip().split('\n')
			rangeStartCodon = set(rangeStartCodon)
			info_to_add = chr+'\t'+str(sorted_range_start_codon[0])+'\t'+str(sorted_range_start_codon[-1])+'\t'+strand+'\t'+codon

			if len(splitRows) >= minReads :
				contReads = 0
				for row in splitRows:
					if row:
						splitRow = row.split('\t')
						name = splitRow[0]
						cigar = splitRow[5]
						sequence = splitRow[9]
						start = int(splitRow[3])
						lenSeq = len(sequence)
						into, pos, rangeRead = getRanges(cigar, start, rangeStartCodon, name, strand, lenSeq)
						
						if into:
							names.add(name)
							contReads += 1
							if len(names) >= minReads:
								try:
									positionsDic = frequenciesByLenFragment[lenSeq] 
									try:
										positionsDic[pos] += 1
									except KeyError:
										positionsDic[pos] = 1
								except KeyError:
									frequenciesByLenFragment[lenSeq] = {pos:1}

				if len(names) >= minReads:
					totalS_sc_known_with_TIS_evidence += 1
					try:
						infoStrand = known_scs_intersected_dic[strand]
						infoStrand[transcript] = scs
					except KeyError:
						known_scs_intersected_dic[strand] = {transcript:scs}
					known_scs_intersected_list.append(sc_string)

					known_intersected.append(info_to_add)

					info_sc_to_add = ''
					for index, sc in enumerate(scs):
						start=int(sc[1].split(':')[1].split('-')[0])
						end=int(sc[1].split(':')[1].split('-')[1])
						info_sc_to_add += chr+'\t'+str(start)+'\t'+str(end)+'\t'+strand+'\t'+transcript+'\n'
						
					if strand == '-':
						known_intersected_backward_keys_list.append(key)
					else:
						known_intersected_forward_keys_list.append(key)

	
	logging.info('Total StartCodons retained analysis %d', totalS_sc_known_with_TIS_evidence)

	with open(toSaveOutputs+'/Known_scs_intercepted.dic', 'wb') as fp:
		pickle.dump(known_scs_intersected_dic, fp, protocol=pickle.HIGHEST_PROTOCOL)	

	logging.info('Total start codons known retained save in dic %s ', toSaveOutputs+'/Known_scs_intercepted.dic')

	with open(toSaveOutputs+'/Known_scs_intercepted.list', 'wb') as fp:
		pickle.dump(known_scs_intersected_list, fp, protocol=pickle.HIGHEST_PROTOCOL)	

	logging.info('Total start codons known retained save in list %s ', toSaveOutputs+'/Known_scs_intercepted.list')

	with open(toSaveOutputs+'/frequenciesByLenFragment.dic', 'wb') as fp:
		pickle.dump(frequenciesByLenFragment, fp, protocol=pickle.HIGHEST_PROTOCOL)

	logging.info('Dic Lens with StartCodons Retained for analysis %s', str(frequenciesByLenFragment))
	
	with open (all_starCodonsPositions, 'rb') as fp:
		all_starCodonsPositions = pickle.load(fp)

	known_intersected_backward = ''
	known_intersected_forward = ''

	for key in known_intersected_forward_keys_list:
		known_intersected_forward += all_starCodonsPositions[key]

	for key in known_intersected_backward_keys_list:
		known_intersected_backward += all_starCodonsPositions[key]

	fileToSave = open(toSaveOutputs+'/Known_SC_Intercepted_Forward.bed', 'w')
	fileToSave.write(known_intersected_forward)
	fileToSave.close()	

	fileToSave = open(toSaveOutputs+'/Known_SC_Intercepted_Backward.bed', 'w')
	fileToSave.write(known_intersected_backward)
	fileToSave.close()	

	t1 = time.time()
	total = t1-t0_
	logging.info('Total time toSeparateDataBase :  %s min', str(total/60))
	
	return frequenciesByLenFragment


def main(argv):
	
	t0 = time.time()

	bamFileInput = ''
	pathToSave = ''
	starCodonsPositions = ''
	allStarCodonsPositions = ''
	minReads = 3
	logPath =''

	try:
		opts, args = getopt.getopt(argv,"hi:o:s:a:m::l:",["bamFileInput=", "pathToSave=", "starCodonsPositions=", "allStarCodonsPositions=", 'minReads=', 'logPath='])
	except getopt.GetoptError:
		print 'main.py -i <bamFileInput> -o <pathToSave> -s <starCodonsPositions> -a <allStarCodonsPositions> -m <minReads> -l <logPath> '
		sys.exit(2)

	print opts
	for opt, arg in opts:
		if opt == '-h':
			print 'main.py -i <bamFileInput> -o <pathToSave> -s <starCodonsPositions>  -a <allStarCodonsPositions>  -m <minReads>  -l <logPath>'
			sys.exit()
		elif opt in ("-i", "--bamFileInput"):
			bamFileInput = arg
		elif opt in ("-o", "--pathToSave"):
			pathToSave = arg
		elif opt in ("-s", "--starCodonsPositions"):
			starCodonsPositions = arg
		elif opt in ("-a", "--allStarCodonsPositions"):
			allStarCodonsPositions = arg
		elif opt in ("-m", "--minReads"):
			minReads = int(arg)
		elif opt in ("-l", "--logPath"):
			logPath = arg


	now = time.strftime('%y-%m-%d_%H:%M:%S', time.localtime())
	nameLog = logPath+'gettingProbasByLenFragment.log'
	print nameLog
	logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')

	logging.info('Running Getting Info Probas By Len Fragment')

	logging.info('starCodonsPositions : %s ', starCodonsPositions)
	logging.info('bamFileInput : %s ', bamFileInput)
	logging.info('pathToSave : %s ',pathToSave)
	logging.info('Minimum reads to count : %d', minReads)

	t0_ = time.time()
	
	frequenciesByLenFragment = getProbas(allStarCodonsPositions, starCodonsPositions, bamFileInput, minReads, pathToSave)
	probas = getInfoProbas.getProbasFromDictionnaryWithFrequencies(frequenciesByLenFragment, pathToSave)
	logging.info('Probabilities \n %s, ', str(probas))
	median_mean_stdev_info = getInfoProbas.get_Median_Mean_Stdev(frequenciesByLenFragment, pathToSave)
	logging.info('Median Mean Stdev Info \n %s, ', str(median_mean_stdev_info))
	weightsByLen = getInfoProbas.getDispersion(median_mean_stdev_info, pathToSave)
	logging.info('weightsByLen \n %s, ', str(weightsByLen))

	t1 = time.time()
	total = t1-t0_
	logging.info('Total time run function getProbasByLenFragmentTis : %f min ', (total/60))


if __name__ == "__main__":
	main(sys.argv[1:])
