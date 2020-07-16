import warnings
import time
from pyGeno.Genome import * 
warnings.filterwarnings("ignore")
from pyGeno.tools.UsefulFunctions import *
from pyGeno.tools.BinarySequence import *
import array
import pickle
from getTranscriptInformation import getTranscriptInformation
import datetime
import logging
from itertools import groupby
from operator import itemgetter
import getopt
import numpy as np
import math

dico = {'M':'ATG', 'L':'CTG', 'L':'TTG', 'V':'GTG', 'T':'ACG', 'I':'ATA', 'I':'ATT', 'I':'ATC', 'K':'AAG', 'R':'AGG'}
listCodonsStart = [ 'ATG','CTG', 'TTG', 'GTG', 'ACG', 'ATA', 'ATT', 'ATC', 'AAG', 'AGG'] 
dico_iupac_ntd = {'R':['A', 'G'], 'Y': ['C', 'T'], 'S':['G', 'C'], 'W': ['A', 'T'], 'K':['G', 'T'], 'M': ['A', 'C'], 'B':['C', 'G', 'T'], 'D':['A', 'G', 'T'], 'D':['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V':['A', 'C', 'G']}

class getInfoTranscriptsIntersected:

	def __init__(self, fileSNPsFreeBayes, qualityThreshold, transcriptsFile, infoBedIntersectTranscripts, starts_annotated, dicStartCodons, strand, folderToSave, fastaFile, filepath_index):
		self.fileSNPsFreeBayes = fileSNPsFreeBayes
		self.qualityThreshold = qualityThreshold
		self.transcriptsFile = transcriptsFile
		self.infoBedIntersectTranscripts = infoBedIntersectTranscripts
		self.dicStartCodons = dicStartCodons
		self.strand = strand
		self.folderToSave = folderToSave
		self.fastaFile = fastaFile
		self.filepath_index = filepath_index
		self.starts_annotated = starts_annotated


	def getSNPsDic(self):

		time0 = time.time()
		dic = {}
		totalSnps = 0
		logging.info('Get SNPs Info::')
		logging.info('Get SNPs Info::File %s ', self.fileSNPsFreeBayes)

		with open(self.fileSNPsFreeBayes) as f:
			for index, snp in enumerate(f): 
				if index > 0:
					splitLine = snp.strip().split('\t')
					chr = splitLine[0]
					if len(chr) < 3:
						chr = 'chr'+chr
					quality = float(splitLine[6])
					if quality > self.qualityThreshold:
						totalSnps += 1
						try:
							lista = dic[chr]
							lista.append(snp)
							dic[chr] = lista
						except KeyError:
							dic[chr] = [snp]
		
		logging.info('Total SNPs Over the qualityThreshold  %d : %d ', self.qualityThreshold, totalSnps)
		timeFinal = time.time()
		total = (timeFinal-time0) / 60
		logging.info('Total runtime getSNPsDic : %f min \n', total)
		return dic, totalSnps



	def getInfoTranscripts(self):

		time0 = time.time()
		dicTranscripts = {}
		logging.info('Get Info Transcripts::')
		logging.info('Get Info Transcripts::File %s ', self.transcriptsFile)

		with open(self.transcriptsFile) as f:

			for index, transcript in enumerate(f):
				splitLine = transcript.strip().split('\t')
				infoLine = splitLine[8].split('; ')
				nameTranscript = infoLine[1].split("\"")[1]
				chromosome = splitLine[0]
				start = int(splitLine[3])
				end = int(splitLine[4])
				strand = splitLine[6]

				if splitLine[2] == 'transcript':
					dicTranscripts[nameTranscript] = [(nameTranscript, chromosome, start, end, strand)]
				elif splitLine[2] == 'exon':
					try:
						lista = dicTranscripts[nameTranscript]
						lista.append((start, end, end-start))
					except KeyError:
						logging.warning('Transcript %s not in dic verify!', nameTranscript)
		
		
		timeFinal = time.time()
		total = (timeFinal-time0) / 60
		logging.info('Total runtime getInfoTranscripts : %f min \n', total)
		return dicTranscripts



	def gettingInfoTranscriptsIntersected(self):

		time0 = time.time()
		fileNameTranscripts = '' 
		fileNameGenes = ''

		logging.info('\n\nGet Info Transcripts Intersected::')
		snpsDic, totalSnps = self.getSNPsDic()
		logging.info('Get Info Transcripts Intersected:: snpsDic OK')
		dicTranscripts = self.getInfoTranscripts()
		logging.info('Get Info Transcripts Intersected:: dicTranscripts OK')
		getTransInfo = getTranscriptInformation(snpsDic, self.fastaFile, self.filepath_index)
		logging.info('Get Info Transcripts Intersected:: getTranscriptInformation OK')

		with open(self.dicStartCodons, 'rb') as handle:
			scInfo = pickle.load(handle)

		canonical_sc = {}
		with open(self.starts_annotated) as f:
			for index, line in enumerate(f):
				splitLine = line.strip().split('\t') 
				chr = splitLine[0]
				start = splitLine[1]
				end = splitLine[2]
				strand = splitLine[3]
				tis = chr+':'+start+'-'+end+'-'+strand
				transcript = splitLine[4].strip()
				canonical_sc[tis] = transcript
		
		logging.info('Get Info Transcripts Intersected:: get SC cannoniques Information OK %d', len(canonical_sc))

		infoTranscripts = {}
		infoTranscriptsOnlyPositions = {}
		genes = {}
		
		logging.info('Summary Files : ')
		logging.info('InfoBedIntersectTranscripts %s', self.infoBedIntersectTranscripts)
		logging.info('transcriptsFile %s', self.transcriptsFile)
		logging.info('dicStartCodons %s', self.dicStartCodons)
		logging.info('qualityThreshold %d ',self.qualityThreshold)
		logging.info('folderToSave %s', self.folderToSave)
		logging.info('Strand %s ', self.strand)
		logging.info('Total SNPs Charged %d', totalSnps)
		logging.info('Total dicTranscripts %d ', len(dicTranscripts.keys()))
		notOverlap = 0

		sc_overlapped = set()
		total_info_trans_intersected = ''
		total_info_trans_not_intersected = ''

		known_sc_know_transcripts = 0
		sc_intersected = set()
		gene_positions = {}

		with open(self.infoBedIntersectTranscripts) as f:
			cont = 0
			for index, tisIntersected in enumerate(f):
				splitLine = tisIntersected.strip().split('\t')
				chrTis = splitLine[0]
				gene = splitLine[11].split(";")[0].split("\"")[1]
				tis_in_transcript = splitLine[11].split(";")[1].split("\"")[1]
				reference_transcript = splitLine[11].split(";")[2].split("\"")[1]
				tpm = float(splitLine[11].split(" TPM ")[1].split('\"')[1])
				
				if 'ENST' not in reference_transcript:
					reference_transcript = '_'

				startTis = splitLine[1]
				endTis = splitLine[2]
				tisString = chrTis+':'+str(startTis)+'-'+str(endTis)+'-'+self.strand
				sc_intersected.add(tisString)
				
				trans = ''
				if tisString in canonical_sc.keys():
					trans =  canonical_sc[tisString]
					if trans == reference_transcript:
						known_sc_know_transcripts += 1
				
				if trans != reference_transcript: 

					infoSC = scInfo[tisString] 
					sc_overlapped.add(tisString)
					codon = infoSC[3]
					pos = infoSC[1]
					strand = infoSC[5]
					range_codon = infoSC[2]
					proba = round(infoSC[0], 3)

					infoToSaveSC =  [pos, codon, strand, proba]
					sequenceTranscript_from_TIS = ''
					sequencePositions_from_TIS = ''
					codonToCompare = ''
					startTis = range_codon[0]
					middleTis = range_codon[1]
					endTis = range_codon[2]
					range_codon = set(range_codon)

					try:
						info = infoTranscripts[tis_in_transcript]
						positions = infoTranscriptsOnlyPositions[tis_in_transcript]
					except KeyError:
						infoTranscript = dicTranscripts[tis_in_transcript]
						transcript = getTransInfo.getTranscriptsRegionsInformation(infoTranscript) 
						infoTranscripts[tis_in_transcript] = [transcript[1], [], [] , []] 
						positions = transcript[2]
						infoTranscriptsOnlyPositions[tis_in_transcript] = positions
						info = infoTranscripts[tis_in_transcript]

					intersected = set(positions).intersection(range_codon)

					if len(intersected) == 3:
						if self.strand == '-':
							posIndex = positions.index(endTis)
						else:
							posIndex = positions.index(startTis)
						
						frame = posIndex % 3
						kozak = info[0][posIndex-7:posIndex+5]
						sequenceTranscript_from_TIS = info[0][posIndex:]
						codonToCompare = sequenceTranscript_from_TIS[0:3]
						positionsTrans = self.find_ranges(positions[posIndex:])
						positionsTranscript = self.find_ranges(positions)

						equal_codon = True
						if codon != codonToCompare:
							equal_codon = False
							if codonToCompare in listCodonsStart:
								equal_codon = True
							else:
								try:
									if codon[0] != codonToCompare[0]:
										r = dico_iupac_ntd[codonToCompare[0]]
										if codon[0] in r:
											equal_codon = True
									if codon[1] != codonToCompare[1]:
										r = dico_iupac_ntd[codonToCompare[1]]
										if codon[1] in r:
											equal_codon = True
									if codon[2] != codonToCompare[2]:
										r = dico_iupac_ntd[codonToCompare[2]]
										if codon[2] in r:
											equal_codon = True
								except KeyError:
									equal_codon = False

						if equal_codon:
							
							if codon != codonToCompare and codonToCompare not in listCodonsStart:
								if kozak != '':
									s = list(kozak)
									s[-3] = codon[2]
									s[-4] = codon[1]
									s[-5] = codon[0]
									kozak = "".join(s)
							
							try:
								gene_positions_list = gene_positions[gene]
								isnt_in = False
								
								for posis in gene_positions_list:
									
									if positionsTrans == posis:
										isnt_in = True
										break

								if not isnt_in:
									singleton_set = set()
									singleton_set.add(positionsTrans)
									intersect = gene_positions_list.intersection(singleton_set)
									if len(intersect) == 0 :
										gene_positions_list.add(positionsTrans)
										gene_positions[gene] = gene_positions_list
										l = info[frame+1]
										infoToSaveSC.extend([posIndex, kozak, positionsTranscript, reference_transcript, tpm])
										l.append(infoToSaveSC)
										total_info_trans_intersected += chrTis+'\t'+str(startTis)+'\t'+str(endTis)+'\t'+codon+'\t'+tis_in_transcript+'\t'+reference_transcript+'\t'+str(positions[posIndex]) +'\t'+str(positions[-1])+'\t'+self.strand+'\t'+positionsTrans+'\t'+str(proba)+'\t'+str(tpm)+'\t'+kozak+'\n'
										try:
											lGenes = genes[gene]
											if tis_in_transcript not in lGenes:
												lGenes.append(tis_in_transcript)
										except KeyError:
											genes[gene] = [tis_in_transcript]
											
							except KeyError:
								singleton_set = set()
								singleton_set.add(positionsTrans)
								gene_positions[gene] = singleton_set
								l = info[frame+1]
								infoToSaveSC.extend([posIndex, kozak, positionsTranscript, reference_transcript, tpm])
								l.append(infoToSaveSC)
								total_info_trans_intersected += chrTis+'\t'+str(startTis)+'\t'+str(endTis)+'\t'+codon+'\t'+tis_in_transcript+'\t'+reference_transcript+'\t'+str(positions[posIndex]) +'\t'+str(positions[-1])+'\t'+self.strand+'\t'+positionsTrans+'\t'+str(proba)+'\t'+str(tpm)+'\t'+kozak+'\n'
								try:
									lGenes = genes[gene]
									if tis_in_transcript not in lGenes:
										lGenes.append(tis_in_transcript)
								except KeyError:
									genes[gene] = [tis_in_transcript]
						else:
							print 'Not equal codons, codon seen ', codonToCompare, ' codon expected ', codon
							
					else:
						kozak = ''
						total_info_trans_not_intersected += tisString+'\t'+tis_in_transcript+'\t'+reference_transcript+'\t'+str(range_codon)+'\t'+str(proba)+'\t'+str(tpm)+'\t'+kozak+'\n'
						notOverlap += 1

		max_ = 0
		name_max = ''	
		for gene in gene_positions:
			if max_ < len(gene_positions[gene]):
				max_ = len(gene_positions[gene])
				name_max = gene

		logging.info('Not Overlap between SC and sequence %d', notOverlap)
		logging.info('Total Transcripts Intersected %d ', len(infoTranscripts.keys()))
		logging.info('Total Genes Intersected %d ', len(genes.keys()))
		
		fileNameTranscripts = self.folderToSave+'/totalTranscriptsIntersected'+self.strand+'.dic'
		logging.info('Saving Total Transcripts Intersected Information (dic) to %s', fileNameTranscripts)

		with open(fileNameTranscripts, 'wb') as handle:
			pickle.dump(infoTranscripts, handle, protocol=pickle.HIGHEST_PROTOCOL)

		fileTransIntercepted = self.folderToSave+'/InfoTranscriptsIntersected'+self.strand+'.gtf'
		logging.info('Saving gtf file of transcripts Intersected %s', fileTransIntercepted)

		fileToSave = open(fileTransIntercepted, 'w')
		fileToSave.write(total_info_trans_intersected)
		fileToSave.close()

		fileTransNotIntercepted = self.folderToSave+'/InfoTranscriptsNotIntersected'+self.strand+'.gtf'
		logging.info('Saving gtf file of transcripts Not Intersected %s', fileTransNotIntercepted)

		fileToSave = open(fileTransNotIntercepted, 'w')
		fileToSave.write(total_info_trans_not_intersected)
		fileToSave.close()	

		fileNameGenes = self.folderToSave+'/totalGenesIntersected'+self.strand+'.dic'
		logging.info('Saving Total Genes Intersected Information (dic) to %s', fileNameGenes)

		with open(fileNameGenes, 'wb') as handle:
			pickle.dump(genes, handle, protocol=pickle.HIGHEST_PROTOCOL)

		timeFinal = time.time()
		total = (timeFinal-time0) / 60
		logging.info('Total runtime getInfoTranscriptsIntersected : %f min \n', total)

		return fileNameTranscripts, fileNameGenes


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


	def transformation_positions_transcripts(self, infoTranscriptsOnlyPositions):
		infoPositionsTranscripts = {}

		for key, value in infoTranscriptsOnlyPositions.items():
			infoPositionsTranscripts[key] = self.find_ranges(value)

		return infoPositionsTranscripts


def main(argv):
	strand = ''
	transcriptsFile = ''
	fileSNPsFreeBayes = ''
	qualityThreshold = 20
	infoBedIntersectTranscripts = ''
	folderToSave = ''
	fastaFile = ''
	filepath_index = ''
	dicStartCodons = ''
	logPath = ''	

	try:
		opts, args = getopt.getopt(argv,"hs:t:n:q:b:o:f:i:d:k:l:",["strand=", "transcriptsFile=", "snps=", 'quality=', 'bedIntersected=','outputfile=','GenomeFasta=','GenomeIndexFasta=', 'dicStartCodons=',  'starts_annotated=', 'logPath='])
	
	except getopt.GetoptError:
		print 'main.py -s <strand>  -t <transcriptsFile> -n <snps> -q <quality> -b <bedIntersected> -o <outputfile> -f <GenomeFasta> -i <GenomeIndexFasta> -d <dicStartCodons> -k <starts_annotated> -l <logPath> '
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print 'main.py -s <strand> -t <transcriptsFile> -n <snps> -q <quality> -b <bedIntersected> -o <outputfile> -f <GenomeFasta> -i <GenomeIndexFasta> -d <dicStartCodons> -k <starts_annotated> -l <logPath> '
			sys.exit()
		elif opt in ("-s", "--strand"):
			strand = arg
		elif opt in ("-t", "--transcriptsFile"):
			transcriptsFile = arg
		elif opt in ("-n", "--snps"):
			fileSNPsFreeBayes = arg
		elif opt in ("-q", "--quality"):
			qualityThreshold = int(arg)
		elif opt in ("-b", "--bedIntersected"):
			infoBedIntersectTranscripts = arg
		elif opt in ("-o", "--outputfile"):
			folderToSave = arg
		elif opt in ("-f", "--GenomeFasta"):
			fastaFile = arg
		elif opt in ("-i", "--GenomeIndexFasta"):
			filepath_index = arg
		elif opt in ("-d", "--dicStartCodons"):
			dicStartCodons = arg
		elif opt in ("-k", "--starts_annotated"):
			starts_annotated = arg
		elif opt in ("-l", "--logPath"):
			logPath = arg
	
	now = time.strftime('%y-%m-%d_%H:%M:%S', time.localtime())

	nameLog = logPath+'GetTotalInfoTranscriptsCandidates'+strand+'.log'
	logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')
	logging.info('Date : %s', str(now))

	t0 = time.time()
	logging.info('Running Get Transcripts Sequences')
	getTranscriptsInfo = getInfoTranscriptsIntersected(fileSNPsFreeBayes, qualityThreshold, transcriptsFile, infoBedIntersectTranscripts, starts_annotated, dicStartCodons, strand, folderToSave, fastaFile, filepath_index)
	fileNameTranscripts, fileNameGenes = getTranscriptsInfo.gettingInfoTranscriptsIntersected()
	logging.info('File Transcripts %s ', fileNameTranscripts)
	logging.info('File Genes %s ', fileNameGenes)
	logging.info('Got transcript Information! \n')
	t2 = time.time()
	total = t2-t0
	logging.info('Total time to run function : %d min', total/60)

if __name__ == "__main__":
	main(sys.argv[1:])






