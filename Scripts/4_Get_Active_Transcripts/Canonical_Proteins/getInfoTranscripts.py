import warnings
import time
import sys, getopt
sys.path += sys.path + ["/u/ruizma/Jupyter/PipelineAligment/Scripts/PythonScripts", "/u/ruizma/Jupyter/PipelineAligment/Scripts/PythonScripts/General/",'/u/ruizma/lib/rabadb','/u/ruizma/lib/pygeno','/u/ruizma/Jupyter/ExploringRiboSeqReads/Scripts' , '/u/ruizma/miniconda2/lib/python2.7/site-packages', '/u/ruizma/PIPELINE/Scripts/DataPreparation/PythonScripts/11_GetTranscriptsInformation/0_Get_Proteins_Sequences/Canonical/getTranscripts/']
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
import numpy as np
import math

dico = {'M':'ATG', 'L':'CTG', 'L':'TTG', 'V':'GTG', 'T':'ACG', 'I':'ATA', 'I':'ATT', 'I':'ATC', 'K':'AAG', 'R':'AGG'}


class getInfoTranscriptsIntersected:

	def __init__(self, fileSNPsFreeBayes, tpm_dictionary, qualityThreshold, transcriptsRefFile, infoBedIntersectTranscripts, dicStartCodons, strand, folderToSave, fastaFile, filepath_index):
		self.fileSNPsFreeBayes = fileSNPsFreeBayes
		self.qualityThreshold = qualityThreshold
		self.transcriptsRefFile = transcriptsRefFile
		self.infoBedIntersectTranscripts = infoBedIntersectTranscripts
		self.strand = strand
		self.folderToSave = folderToSave
		self.fastaFile = fastaFile
		self.filepath_index = filepath_index
		self.dicStartCodons = dicStartCodons
		self.tpm_dictionary = tpm_dictionary


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
		logging.info('Get Info Transcripts::File %s ', self.transcriptsRefFile)

		with open(self.transcriptsRefFile) as f:
			for index, transcript in enumerate(f):
				splitLine = transcript.strip().split('\t')
				infoLine = splitLine[8].split('; ')
				type_ = splitLine[2]
				if type_ == 'transcript' or type_ == 'exon':
					nameTranscript = infoLine[1].split("\"")[1]
					chromosome = splitLine[0]
					start = int(splitLine[3])
					end = int(splitLine[4])
					strand = splitLine[6]
					if type_ == 'transcript':
						dicTranscripts[nameTranscript] = [(nameTranscript, chromosome, start, end, strand)]
					elif type_ == 'exon':
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

		logging.info('Get Info Transcripts Intersected::')
		snpsDic, totalSnps = self.getSNPsDic()
		logging.info('Get Info Transcripts Intersected:: snpsDic OK')
		dicTranscripts = self.getInfoTranscripts()
		logging.info('Get Info Transcripts Intersected:: dicTranscripts OK')
		getTransInfo = getTranscriptInformation(snpsDic, self.fastaFile, self.filepath_index)
		logging.info('Get Info Transcripts Intersected:: getTranscriptInformation OK')
		
		with open (self.tpm_dictionary, 'rb') as fp:
			getQuant = pickle.load(fp)
		logging.info('Get Info Quantification:: getQuantification OK')
		
		infoTranscripts = {}
		genes = {}
		
		logging.info('Summary Files : ')
		logging.info('InfoBedIntersectTranscripts %s', self.infoBedIntersectTranscripts)
		logging.info('transcriptsRefFile %s', self.transcriptsRefFile)
		logging.info('qualityThreshold %d ',self.qualityThreshold)
		logging.info('folderToSave %s', self.folderToSave)
		logging.info('Strand %s ', self.strand)
		logging.info('Total SNPs Charged %d', totalSnps)
		logging.info('Total dicTranscripts %d ', len(dicTranscripts.keys()))
		notOverlap = 0

		total_info_trans_intersected = ''

		tis_by_transcripts = {}
		with open(self.infoBedIntersectTranscripts) as f:
			for index, tisIntersected in enumerate(f):

				splitLine = tisIntersected.strip().split('\t') 
				chrTis = splitLine[0]
				startTis = int(splitLine[1])
				endTis = int(splitLine[2])+1
				range_sc = range(startTis, endTis)
				transcript_name = splitLine[4]
				gene = splitLine[14].split(";")[0].split("\"")[1]
				tis_in_transcript = splitLine[14].split(";")[1].split("\"")[1]
				
				try:
					infoTrans = tis_by_transcripts[tis_in_transcript]
					infoTrans[2].extend(range_sc)
				except KeyError:
					tis_by_transcripts[tis_in_transcript] = [chrTis, gene, range_sc]

		
		with open (self.dicStartCodons, 'rb') as fp:
			self.dicStartCodons = pickle.load(fp)
			
		trans_tpm_0 = 0
		for transcript, sc in tis_by_transcripts.items():
			chrTis = sc[0]
			gene = sc[1]
			range_sc = np.sort(sc[2])
			tisString = chrTis+':'+str(range_sc[0])+'-'+str(range_sc[2])+'-'+self.strand
			sequenceTranscript_from_TIS = ''

			if '_' in transcript:
				trans = transcript.split('_')[0]
			else:
				trans = transcript
			try:
				tpm =  getQuant[trans]
			except KeyError:
				tpm = 0
				trans_tpm_0 += 1
			
			startTis = range_sc[0]
			middleTis = range_sc[1]
			endTis = range_sc[2]
			range_sc = set(np.sort(sc[2]))
			
			try:
				sc_known = self.dicStartCodons[tisString]
				proba = round(sc_known[0], 3)
				tis = sc_known[1]
				codon = sc_known[3]
				infoToSaveSC =  [tis, codon, self.strand, proba]
			except KeyError:
				infoToSaveSC =  [tisString, 'NF', self.strand, 0]

			add_new_transcript_information = False
			
			try:
				info = infoTranscripts[transcript][0]
				positions = info[0][2]
			except KeyError:
				infoTranscript = dicTranscripts[transcript]
				transcriptInfo = getTransInfo.getTranscriptsRegionsInformation(infoTranscript) 
				add_new_transcript_information = True
				info = transcriptInfo
				positions = transcriptInfo[2]

			intersected = set(positions).intersection(range_sc)
			if len(intersected) == 3:
				if self.strand == '-':
					posIndex = positions.index(endTis)
				else:
					posIndex = positions.index(startTis)
				
				frame = posIndex % 3
				sequenceTranscript_from_TIS = info[1][posIndex:]
				codon = info[1][posIndex:posIndex+3]
				positionsTrans = self.find_ranges(info[2][posIndex:])
				kozak = info[1][posIndex-7:posIndex+5]
				infoToSaveSC.extend([posIndex, kozak, positionsTrans, transcript, tpm])

				if add_new_transcript_information:
					infoTranscripts[transcript] = [transcriptInfo, [[sequenceTranscript_from_TIS, frame, chrTis, posIndex, transcript, codon, gene, tisString, proba, infoToSaveSC, positionsTrans]]]
				else:
					infoTranscripts[transcript][1].append([sequenceTranscript_from_TIS, frame, chrTis, posIndex, transcript, codon, gene, tisString, proba, infoToSaveSC, positionsTrans])

				total_info_trans_intersected += chrTis+'\t'+str(startTis)+'\t'+str(endTis)+'\t'+codon+'\t'+transcript+'\t'+str(info[2][posIndex])+'\t'+str(info[2][-1])+'\t'+self.strand+'\t'+positionsTrans+'\t'+str(proba)+'\t'+str(tpm)+'\t'+kozak+'\n'
			else:
				notOverlap += 1


		logging.info('# not Overlap between SC and sequence %d', notOverlap)
		logging.info('Total Intersected Transcripts %d ', len(infoTranscripts.keys()))
		
		now = str(datetime.datetime.now())
		fileNameTranscripts = self.folderToSave+'/Total_Transcripts_Intersected_Canonical_'+self.strand+'.dic'
		logging.info('Saving Total Transcripts Intersected Information canonical (dic) to %s', fileNameTranscripts)

		with open(fileNameTranscripts, 'wb') as handle:
			pickle.dump(infoTranscripts, handle, protocol=pickle.HIGHEST_PROTOCOL)

		fileTransIntercepted = self.folderToSave+'/InfoTranscriptsIntersected'+self.strand+'.gtf'
		logging.info('Saving gtf file of transcripts Intersected %s', fileTransIntercepted)

		fileToSave = open(fileTransIntercepted, 'w')
		fileToSave.write(total_info_trans_intersected)
		fileToSave.close()	

		timeFinal = time.time()
		total = (timeFinal-time0) / 60
		logging.info('Total runtime getInfoTranscriptsIntersected : %f min \n', total)

		return fileNameTranscripts


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


def main(argv):
	strand = ''
	transcriptsRefFile = ''
	fileSNPsFreeBayes = ''
	qualityThreshold = 20
	infoBedIntersectTranscripts = ''
	folderToSave = ''
	fastaFile = ''
	filepath_index = ''
	dicStartCodons = ''
	logPath = ''
	tpm_dictionary = ''
	
	try:
		opts, args = getopt.getopt(argv,"hs:t:n:q:b:o:f:i:d:l:k:",["strand=", "transcriptsRefFile=", "SNPsFreeBayes=", 'QualitySNPs=', 'BedIntersectedSCTrans=', 'outputfile=', 'GenomeFasta=', 'GenomeIndexFasta=', 'dicStartCodons=', 'logPath=', 'tpm_dictionary='])
	except getopt.GetoptError:
		print 'main.py -s <strand> -t <transcriptsRefFile> -n <SNPsFreeBayes> -q <QualitySNPs> -b <BedIntersectedSCTrans> -o <outputfile> -f <GenomeFasta> -i <GenomeIndexFasta> -d <dicStartCodons> -l <logPath> -k <tpm_dictionary>'
		sys.exit(2)

	print opts
	for opt, arg in opts:
		if opt == '-h':
			print 'main.py -s <strand> -t <transcriptsRefFile> -n <SNPsFreeBayes> -q <QualitySNPs> -b <BedIntersectedSCTrans> -o <outputfile> -f <GenomeFasta> -i <GenomeIndexFasta> -d <dicStartCodons> -l <logPath> -k <tpm_dictionary>'
			sys.exit()
		elif opt in ("-s", "--strand"):
			strand = arg
		elif opt in ("-t", "--transcriptsRefFile"):
			transcriptsRefFile = arg
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
		elif opt in ("-l", "--logPath"):
			logPath = arg
		elif opt in ("-k", "--tpm_dictionary"):
			tpm_dictionary = arg

	now = datetime.datetime.now()
	nameLog = logPath+'GetTotalInfoTranscriptsCandidates'+strand+'.log'
	logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')

	logging.info('Date %s ', str(now))
	logging.info('Running Get Transcripts Sequences')
	t0 = time.time()

	getTranscriptsInfo = getInfoTranscriptsIntersected(fileSNPsFreeBayes, tpm_dictionary, qualityThreshold, transcriptsRefFile, infoBedIntersectTranscripts, dicStartCodons, strand, folderToSave, fastaFile, filepath_index)
	fileNameTranscripts = getTranscriptsInfo.gettingInfoTranscriptsIntersected()
	logging.info('Got transcript Information!')
	
	t2 = time.time()
	total = t2-t0
	logging.info('Total time run function getInfoTranscriptsIntersected : %d min', total/60)

if __name__ == "__main__":
	main(sys.argv[1:])


