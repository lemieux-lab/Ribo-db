
import sys, getopt
import time
import pickle
import datetime
import logging
import detectionStartCodons as detection
import auxFunctionsDetectionSC as aux


def runStartCodonsCandidats(folder, probasByLenFragment, weightsByLenght, minReads, fastaFile, filepath_index, samFile, freeBayesFile, qualityThreshold):

	logging.info('Running mode : getting Start Codons ')
	
	time0 = time.time()
	
	detection.startCodonsByReadsFromSamFile(samFile, probasByLenFragment, weightsByLenght, minReads, fastaFile, filepath_index, folder, freeBayesFile, qualityThreshold)
	logging.info('done ! ')
	
	timeFinal = time.time()
	total = (timeFinal-time0) / 60
	logging.info('Total time getting all start codons candidats mode : %f min ', total)


def runResults(folder, retainedStartsKnown):

	time0 = time.time()
	logging.info('Running mode : getting results')
	aux.getResultsStartCodonsCandidats(folder+'/SortedStartCodonsFilteredByMinReads.list', retainedStartsKnown, folder)

	timeFinal = time.time()
	total = (timeFinal-time0) / 60
	logging.info('Total time to run Only Results mode : %f min', total) 


def str_to_bool(s):
	if s == 'True':
		return True
	elif s == 'False':
		return False
	else:
		raise ValueError


def main(argv):
	
	t0 = time.time()

	getResults = False
	folderToSave = ''
	probasByLenFragment = ''
	weightsByLenght = ''
	minReads = 3
	fastaFile = ''
	filepath_index = ''
	samFile = ''
	logPath = ''
	retainedStartsApprentissage = ''
	retainedStartsValidation = ''
	freeBayesFile = ''
	qualityThreshold = 20
	retainedStartsKnown = ''

	try:
		opts, args = getopt.getopt(argv,"hr:o:p:w:m:g:i:s:k:a:f:q:l:",["getResults=", "folderToSave=", "probasByLenFragment=", "weightsByLenght=", 'minReads=', 'fastaFile=', 'filepath_index=', 'samFile=', 'retainedStartsKnown=', 'freeBayesFile=', 'qualityThreshold=', 'logPath='])
	except getopt.GetoptError:
		print 'main.py -r <getResults> -o <folderToSave> -p <probasByLenFragment> -w <weightsByLenght> -m <minReads> -g <fastaFile> -i <filepath_index> -s <bamFile>  -k <retainedStartsKnown>   -f <freeBayesFile>  -q <qualityThreshold>  -l <logPath> '
		sys.exit(2)

	print opts
	for opt, arg in opts:
		if opt == '-h':
			print 'main.py -r <getResults> -o <folderToSave> -p <probasByLenFragment> -w <weightsByLenght> -m <minReads> -g <fastaFile> -i <filepath_index> -s <bamFile>  -k <retainedStartsKnown>   -f <freeBayesFile>  -q <qualityThreshold> -l <logPath> '
			sys.exit()
		elif opt in ("-r", "--getResults"):
			getResults = str_to_bool(arg)
		elif opt in ("-o", "--folderToSave"):
			folderToSave = arg
		elif opt in ("-p", "--probasByLenFragment"):
			probasByLenFragment = arg
		elif opt in ("-w", "--weightsByLenght"):
			weightsByLenght = arg
		elif opt in ("-m", "--minReads"):
			minReads = int(arg)
		elif opt in ("-g", "--fastaFile"):
			fastaFile = arg
		elif opt in ("-i", "--filepath_index"):
			filepath_index = arg
		elif opt in ("-s", "--samFile"):
			samFile = arg
		elif opt in ("-k", "--retainedStartsKnown"):
			retainedStartsKnown = arg
		elif opt in ("-l", "--logPath"):
			logPath = arg
		elif opt in ("-f", "--freeBayesFile"):
			freeBayesFile = arg
		elif opt in ("-q", "--qualityThreshold"):
			qualityThreshold = int(arg)
		
		
	now = time.strftime('%y-%m-%d_%H:%M:%S', time.localtime())
	if getResults:
		nameLog = logPath+'StartCodons_search_Detection.log'
	else:
		nameLog = logPath+'StartCodons_search_GetResults.log'
	
	logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')
	logging.info('Date : %s', str(now))

	logging.info('Folder to Log : %s', nameLog)
	t0_ = time.time()
	if getResults:
		logging.info('Running mode : getting Start Codons %s ', str(getResults))
		logging.info('Parametters :')
		logging.info('folderToSave : %s', folderToSave)
		logging.info('probasByLenFragment path : %s ', probasByLenFragment)
		logging.info('weightsByLenght : %s ', str(weightsByLenght))
		logging.info('minReads : %d ', minReads)
		logging.info('samFile path : %s',samFile)
		logging.info('freebayes : %s',freeBayesFile)
		logging.info('qualityThreshold : %d',qualityThreshold)
		runStartCodonsCandidats(folderToSave, probasByLenFragment, weightsByLenght, minReads, fastaFile, filepath_index, samFile, freeBayesFile, qualityThreshold)
		
	else:
		t0_ = time.time()
		logging.info('Running mode : getting results')
		logging.info('folderToSave path : %s ', folderToSave)
		logging.info('retainedStartsKnown path : %s', retainedStartsKnown)
		runResults(folderToSave, retainedStartsKnown)
		t1 = time.time()
		total = t1-t0_
		logging.info('Total time getting Only Results : %f min ', (total/60))
		
	t1 = time.time()
	total = t1-t0
	logging.info('Total time whole Searching Start Codons : %f min ', (total/60))

if __name__ == "__main__":
	main(sys.argv[1:])



