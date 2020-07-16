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

	def __init__(self, canonique, proteinesCandidatesForward, proteinesCandidatesBackward, toSaveOutput):
		self.canonique = canonique
		self.proteinesCandidatesForward = proteinesCandidatesForward
		self.proteinesCandidatesBackward = proteinesCandidatesBackward
		self.toSaveOutput = toSaveOutput
	

	def addingProteinsFrom_Cannonique(self):

		time0 = time.time()
		startCodonsOrigin = {}
		namesPP = [] 
		proteinesPP = []
		
		with open(self.canonique) as f:
			for index, line in enumerate(f):
				if '>' not in line:
					proteinesPP.append(line)
				else:
					namesPP.append(line)

		totalCandidates = self.getTotalCandidates(False) 
		
		logging.info('Summary Files : ')
		logging.info('File DB canonical db %s \n',self.canonique)
		
		logging.info('Total canonical proteins %d ',len(proteinesPP))
		logging.info('Total non canonical proteins candidates  %d  ',len(totalCandidates))

		notContained = []
		notContainedNames = []
		containedProteins = []
		contained_ = 0

		for index, proteine in enumerate(proteinesPP):
			seq = proteine.strip()
			contained = self.removeRepeatedVariants(seq, totalCandidates)
			if not contained:
				notContained.append(seq)
				notContainedNames.append(namesPP[index].strip())
			else:
				contained_ += 1
				containedProteins.append(index)
			if index % 5000 == 0:
				logging.info('Protein Contained : %d',contained_) 
				logging.info('Protein not Contained : %d', len(notContained)) 
				logging.info('Index : %d',index)

		return fileProteinsToADD, fileNamesProteinsToADD


	def removeRepeatedVariants(self, proteine, proteinesDB):
		contained = False
		for proteineDB in proteinesDB:
			if proteine == proteineDB[0]:
				contained = True
				return contained
		return contained

		
	def getTotalCandidates(self, save):

		with open (self.proteinesCandidatesBackward, 'rb') as fp:
			candidatesBackward = pickle.load(fp)

		with open (self.proteinesCandidatesForward, 'rb') as fp:
			candidatesForward = pickle.load(fp)

		logging.info('Summary Files : ')
		logging.info('Proteins candidatesBackward %s ',self.proteinesCandidatesBackward)
		logging.info('Proteins candidatesForward %s ',self.proteinesCandidatesForward)
		logging.info('Total proteins Backward %d  ',len(candidatesBackward))
		logging.info('Total proteins Forward %d  ',len(candidatesForward))
		
		totalCandidates = candidatesForward + candidatesBackward
		logging.info('Total proteins Candidates  %d  ',len(totalCandidates))
		if save:		
			fileNameSavingTotal = self.toSaveOutput+'/1_proteinsUniques.info'

			with open(fileNameSavingTotal, 'wb') as handle:
				pickle.dump(totalCandidates, handle, protocol=pickle.HIGHEST_PROTOCOL)

		return totalCandidates


	def writeDatabase(self, proteinesToKeep, notContained, notContainedNames, add):

		time0 = time.time()
		logging.info('----------------------------------------')
		logging.info('Starting Writing Custom DB writeDatabase')
		startCodonsOrigin = {}

		if add:
			with open (notContained, 'rb') as fp:
				notContained = pickle.load(fp)
			
			with open (notContainedNames, 'rb') as fp:
				notContainedNames = pickle.load(fp)

		toWriteProt = ''
		changeInitiation = 0
		onlyATG = 0
		totalProteinsWritten = 0

		logging.info('Total proteins to Keep  %d', len(proteinesToKeep))
		
		index_prot = 0
		for index, proteine in enumerate(proteinesToKeep):
			transcript = proteine[1][0]
			strand = proteine[1][1][2]
			seq = proteine[0]
			codon = proteine[1][1][1]
			contVariant = proteine[2]
			lenProteine = len(seq)
			lon = proteine[5]
			protein_id = proteine[1][3]

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
						if index2 == 0 and codon != 'ATG' and codon != 'CTG' :
							seq = 'M'+prot[1:]
							name = '>'+str(index_prot)+'_'+str(protein_id)+'_'+str(index2)+'_'+strand+'_ToM_splitBP\n'
							toWriteProt += name +seq +'\n'
							changeInitiation += 1
							totalProteinsWritten += 1
						else :
							name = '>'+str(index_prot)+'_'+str(protein_id)+'_'+str(index2)+'_'+strand+'_splitBP\n' 
							toWriteProt += name +prot +'\n'
						if codon == 'ATG' and index == 0:
							onlyATG += 1
				else:	
					if codon != 'ATG' and codon != 'CTG' : 
						seq = 'M'+seq[1:]
						name = '>'+str(index_prot)+'_'+str(protein_id)+'_'+strand+'_ToM\n'
						toWriteProt += name +seq +'\n'
						changeInitiation += 1
						totalProteinsWritten += 1
					else :
						name = '>'+str(index_prot)+'_'+str(protein_id)+'_'+strand+'_M_L\n' 
						toWriteProt += name +seq +'\n'
						totalProteinsWritten += 1
						if codon == 'ATG':
							onlyATG += 1
			else:
				if lenProteine < 10000:
					name =  '>'+str(index_prot)+'_'+str(protein_id)+'_'+strand+'_V\n' 
					toWriteProt += name +seq +'\n'
					totalProteinsWritten += 1
				else:
					for index2, prot in enumerate(prots) :
						name = '>'+str(index_prot)+'_'+str(protein_id)+'_'+str(index2)+'_'+strand+'_splitBP_Vari\n'
						toWriteProt += name +prot +'\n'
						if index2 == 0:
							totalProteinsWritten += 1

			index_prot += 1

		if add:
			for index, protein in enumerate(notContained):
				name = notContainedNames[index]+'\n'
				toWriteProt += name +protein +'\n'

		logging.info('Total Proteins written %d ', totalProteinsWritten)
		logging.info('Total Proteins Change Initiation to Methyonine %d ', changeInitiation)
		logging.info('Total Proteins with ATG start codon %d ', onlyATG)

		if add:
			logging.info('Total Proteins added from canonique %d ', len(notContained))
		
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



def str_to_bool(s):
	if s == 'True':
		return True
	elif s == 'False':
		return False
	else:
		raise ValueError

def main(argv):
	getNonCanonical = False
	cannonique = ''
	proteinesCandidatesForward = ''
	proteinesCandidatesBackward = ''
	toSaveOutput = ''
	logPath = ''
	proteinsFromCann = '' 
	namesProteinsFromCann = ''
	concatenation = False
	maxLength = 10000

	try:
		opts, args = getopt.getopt(argv,"hc:f:b:o:p:n:a:u:l:m:",["cannonique=","proteinesCandidatesForward=","proteinesCandidatesBackward=", "toSaveOutput=", "proteinsFromCann=", "namesProteinsFromCann=", "getNonCanonical=", "union=", "logPath=", "maxLength="])
	except getopt.GetoptError:
		print 'getUniqueSequences -c <cannonique> -f <proteinesCandidatesForward> -b <proteinesCandidatesBackward> -o <toSaveOutput> -p <proteinsFromCann> -n <namesProteinsFromCann> -a <getNonCanonical> -u <union> -l <logPath> -m <maxLength>'
		sys.exit(2)

	print opts
	for opt, arg in opts:
		if opt == '-h':
			print 'getUniqueSequences -c <cannonique> -f <proteinesCandidatesForward> -b <proteinesCandidatesBackward> -o <toSaveOutput> -p <proteinsFromCann> -n <namesProteinsFromCann> -a <getNonCanonical> -u <union> -l <logPath> -m <maxLength>'
			sys.exit()
		elif opt in ("-c", "--cannonique"):
			cannonique = arg
		elif opt in ("-f", "--proteinesCandidatesForward"):
			proteinesCandidatesForward = arg
		elif opt in ("-b", "--proteinesCandidatesBackward"):
			proteinesCandidatesBackward = arg
		elif opt in ("-o", "--toSaveOutput"):
			toSaveOutput = arg
		elif opt in ("-l", "--logPath"):
			logPath = arg
		elif opt in ("-p", "--proteinsFromCann"):
			proteinsFromCann = arg
		elif opt in ("-n", "--namesProteinsFromCann"):
			namesProteinsFromCann = arg
		elif opt in ("-a", "--getNonCanonical"):
			getNonCanonical = str_to_bool(arg)
		elif opt in ("-u", "--Concatenation"):
			concatenation = str_to_bool(arg)
		elif opt in ("-m", "--maxLength"):
			maxLength = int(arg)


	now = time.strftime('%y-%m-%d_%H:%M:%S', time.localtime())
	if getNonCanonical:
		nameLog = logPath+'GetNonCanonicalDB.log'
		logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')
		logging.info('Date %s',str(now))
		t0_ = time.time()
		getCustomBD = getUniqueSequences('', proteinesCandidatesForward, proteinesCandidatesBackward, toSaveOutput)
		
		logging.info('Running Get Proteins Sequences')
		logging.info('---------------------Writing DB---------------------------')
		proteinesToKeep = getCustomBD.getTotalCandidates(True)
		getCustomBD.writeDatabase(proteinesToKeep, '', '', False)
		logging.info('Got DB! \n')


	if concatenation:
		nameLog = logPath+'GetConcatenationDB_Canonical_And_NonCanonical.log'
		logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')
		logging.info('Date %s',str(now))
		
		t0_ = time.time()
	
		getCustomBD = getUniqueSequences(cannonique, proteinesCandidatesForward, proteinesCandidatesBackward, toSaveOutput)
	
		logging.info('Running Get Sequences In cannonique to Add to Custom')
		proteinsFromCann, namesProteinsFromCann = getCustomBD.addingProteinsFrom_Cannonique()
		logging.info('Got Info DB PP! \n')

		logging.info('Running Get Proteins Sequences')
		logging.info('---------------------Writing DB---------------------------')
		proteinesToKeep = getCustomBD.getTotalCandidates(False)
		getCustomBD.writeDatabase(proteinesToKeep, proteinsFromCann, namesProteinsFromCann, True)
		logging.info('Got DB! \n')

	t2 = time.time()
	total = t2-t0_
	logging.info('Total time run function getProteines : %f min', total/60)

if __name__ == "__main__":
	main(sys.argv[1:])





