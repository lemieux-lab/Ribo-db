import warnings
import time
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


class getUniqueSequences:

	def __init__(self, proteinesCandidatesForward, proteinesCandidatesBackward, toSaveOutput):
		self.proteinesCandidatesForward = proteinesCandidatesForward
		self.proteinesCandidatesBackward = proteinesCandidatesBackward
		self.toSaveOutput = toSaveOutput

	
	def getTotalCandidates(self):

		with open (self.proteinesCandidatesBackward, 'rb') as fp:
			candidatesBackward = pickle.load(fp)

		with open (self.proteinesCandidatesForward, 'rb') as fp:
			candidatesForward = pickle.load(fp)

		totalCandidates = candidatesForward + candidatesBackward
		now = str(datetime.datetime.now())
		fileNameSavingTotal = self.toSaveOutput+'/1_proteinsUniques.info'

		with open(fileNameSavingTotal, 'wb') as handle:
			pickle.dump(totalCandidates, handle, protocol=pickle.HIGHEST_PROTOCOL)

		return totalCandidates


	def writeDatabase(self, proteinesToKeep):

		time0 = time.time()
		logging.info('----------------------------------------')
		logging.info('Starting Writing Canonical DB ')
		startCodonsOrigin = {}

		if not isinstance(proteinesToKeep, list):
			with open (proteinesToKeep, 'rb') as fp:
				proteinesToKeep = pickle.load(fp)

		toWriteProt = ''
		changeInitiation = 0
		onlyATG = 0
		totalProteinsWritten = 0
		alt_prote = 0

		logging.info('Total proteins to Keep  %d', len(proteinesToKeep))
		print 'Total Proteins to Keep : ', len(proteinesToKeep)

		for index, proteine in enumerate(proteinesToKeep): 
			transcript = proteine[1][0]
			strand = proteine[1][1][2]
			seq = proteine[0]
			codon = proteine[1][1][1]
			contVariant = proteine[2]
			lon = proteine[5]
			lenProteine = len(seq)

			prots = []

			if lenProteine > 10000:
				fois =  lenProteine/9000
				ini = 0
				for i in range(0, fois+1):
					if ini != 0:
						prots.append(seq[ini-30: ini+9000])
					else:
						prots.append(seq[ini: ini+9000])
					ini += 9000

			if contVariant == 0:
				try:
					l = startCodonsOrigin[codon]
					try:
						l[lenProteine] += 1
					except KeyError:
						l[lenProteine] = 1
				except KeyError:
					startCodonsOrigin[codon] = {lenProteine:1}
				
				
				if lenProteine > 10000:
					for index2, prot in enumerate(prots) :
						if index2 == 0 and codon != 'ATG' and codon != 'CUG' :
							seq = 'M'+prot[1:]
							name = '>'+str(index)+'_'+str(index2)+'_'+strand+'_Can_ToM_splitBP\n'
							toWriteProt += name +seq +'\n'
							changeInitiation += 1
							totalProteinsWritten += 1
						else :
							name = '>'+str(index)+'_'+str(index2)+'_'+strand+'_Can_splitBP\n'
							toWriteProt += name +prot +'\n'
						if codon == 'ATG' and index == 0:
							onlyATG += 1
				else:	
					if codon != 'ATG' and codon != 'CUG' :
						seq = 'M'+seq[1:]
						name = '>'+str(index)+'_'+strand+'_Can_ToM\n'
						toWriteProt += name +seq +'\n'
						changeInitiation += 1
						totalProteinsWritten += 1
					else :
						name = '>'+str(index)+'_'+strand+'_Can_M_L\n'
						toWriteProt += name +seq +'\n'
						totalProteinsWritten += 1
						if codon == 'ATG':
							onlyATG += 1
			else:
				alt_prote += 1
				if lenProteine < 10000:
					name = '>'+str(index)+'_'+strand+'_V\n'
					toWriteProt += name +seq +'\n'
					totalProteinsWritten += 1
				else:
					for index2, prot in enumerate(prots) :
						name = '>'+str(index)+'_'+str(index2)+'_'+strand+'_Can_splitBP_Vari\n'
						toWriteProt += name +prot +'\n'
						if index2 == 0:
							totalProteinsWritten += 1

		logging.info('Total Proteins written %d ', totalProteinsWritten)
		logging.info('Total Proteins Change Initation to Methyonine %d ', changeInitiation)
		logging.info('Total Proteins with ATG start codon %d ', onlyATG)
		logging.info('Total variant proteins %d ', alt_prote)

		fileNameSavingCustom = self.toSaveOutput+'/Custom_DB_Total.fasta'
		logging.info('Saving Custom DB (fasta) to %s', fileNameSavingCustom)

		fileToSave = open(fileNameSavingCustom, 'w')
		fileToSave.write(toWriteProt)
		fileToSave.close()	

		logging.info('Saving Information of Startcodons Usage ')

		fileNameSavingStartCodonsOrigin = self.toSaveOutput+'/StartCodonsOrigin_DB.list'
		logging.info('Saving Start Codons Origin (info) to %s', fileNameSavingStartCodonsOrigin)

		with open(fileNameSavingStartCodonsOrigin, 'wb') as handle:
			pickle.dump(startCodonsOrigin, handle, protocol=pickle.HIGHEST_PROTOCOL)

		logging.info('Path to start codon usage dic %s ', fileNameSavingStartCodonsOrigin)

		timeFinal = time.time()
		total = (timeFinal-time0) / 60
		logging.info('Total runtime writeDatabase %f minutes', total)


def main(argv):
	proteinesCandidatesForward = ''
	proteinesCandidatesBackward = ''
	toSaveOutput = ''
	logPath = ''

	try:
		opts, args = getopt.getopt(argv,"hf:b:o:l:",["proteinesCandidatesForward=", "proteinesCandidatesBackward=", "toSaveOutput=", "logPath="])
	except getopt.GetoptError:
		print 'getUniqueSequences.py -f <proteinesCandidatesForward> -b <proteinesCandidatesBackward> -o <toSaveOutput> -l <logPath>'
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print 'getUniqueSequences.py -f <proteinesCandidatesForward> -b <proteinesCandidatesBackward> -o <toSaveOutput> -l <logPath>'
			sys.exit()
		elif opt in ("-f", "--proteinesCandidatesForward"):
			proteinesCandidatesForward = arg
		elif opt in ("-b", "--proteinesCandidatesBackward"):
			proteinesCandidatesBackward = arg
		elif opt in ("-o", "--toSaveOutput"):
			toSaveOutput = arg
		elif opt in ("-l", "--logPath"):
			logPath = arg

	now = time.strftime('%y-%m-%d_%H:%M:%S', time.localtime())
	nameLog = logPath+'GetCanoniqueDataBase.log'
	logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')
	logging.info('Date %s ', str(now))
	logging.info('Running Get Proteins Sequences')

	getCustomBD = getUniqueSequences(proteinesCandidatesForward, proteinesCandidatesBackward, toSaveOutput)
	proteinesToKeep = getCustomBD.getTotalCandidates()
	logging.info('Got get Info DB PP! \n')

	logging.info('-------Wiritng DB-----')
	getCustomBD.writeDatabase(proteinesToKeep)
	logging.info('Got DB! \n')

if __name__ == "__main__":
	main(sys.argv[1:])





