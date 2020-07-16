import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import random
import math
import pickle
import time as t
import logging
import getInfoToPrintROC as getInfo


def getResultsStartCodonsCandidats(startCodonsCandidats, knownStartsCodon, fileToSaveRes):

	time0 = t.time()

	logging.info('=========== Resultats ==========')

	with open (knownStartsCodon, 'rb') as fp:
		retainedStartsKnown = pickle.load(fp)

	status, totalListe, isHere, len_retainedSC, coverage, startCodonsKnown, indexStartCodonsKnown = getInfo.getCourbeRoc(startCodonsCandidats, retainedStartsKnown)
	saveFileCourbe = fileToSaveRes+'/ROC_courbe.png'

	logging.info('Total start known codons retained %d ', len_retainedSC)
	logging.info('First Proba %s  Last Proba %s' , str(totalListe[0]),  str(totalListe[-1]))
	logging.info('True Positives %d False Positives %d Total %d', sum(status),  (len(status) - sum(status)), len(status))

	indexDisMin, xThreshold, yThreshold, auc = printCourbeROC(status, saveFileCourbe)
	logging.info('ROC_courbe AUC %f ', auc)
	if indexDisMin == 0:
		indexDisMin = len(totalListe) -1

	logging.info('')
	logging.info('Total SC found in genome %d', len(totalListe))
	logging.info('Total SC saved after getting threshold %d ', len(totalListe[:indexDisMin]))
	logging.info('Threshold %s ', str(totalListe[indexDisMin]))
	logging.info('Coordinate  X : %f Y : %f',  xThreshold, yThreshold)

	logging.info('Validation : ')
	getInfo.validation(startCodonsCandidats, retainedStartsKnown, fileToSaveRes)

	startCodonsCandidats = totalListe[:indexDisMin]
	totalList_with_SC_Reads_OverMinimum = {}

	for item in startCodonsCandidats:
		proba = item[0]
		countReads = item[4]
		strand = item[5]
		codon = item[3]
		infoStart = item[1]+'-'+strand
		totalList_with_SC_Reads_OverMinimum[infoStart] = item

	with open(fileToSaveRes+'/StartCodonsAfterThreshold.dic', 'wb') as handle:
		pickle.dump(totalList_with_SC_Reads_OverMinimum, handle, protocol=pickle.HIGHEST_PROTOCOL)

	with open(fileToSaveRes+'/StartCodonsAfterThreshold.list', 'wb') as fp:
		pickle.dump(startCodonsCandidats, fp)

	
	logging.info('Total SC found in genome %d ', len(totalListe))
	
	getInfo.getStartCodonsSavedByThreshold(startCodonsCandidats, fileToSaveRes, False)

	timeFinal = t.time()
	total = (timeFinal-time0) / 60
	logging.info('Total time running Resultats : %f min ', total)


def printCourbeROC(status, saveFileCourbe):

	logging.info('printCourbeROC')
	x = [0] 
	y = [0]
	
	try:
		mT = 1/(sum(status)*1.0) # total number of start codons known, true positives, increment courbe roc Y
	except ZeroDivisionError:
		mT = 0

	try:
		mF = 1/((len(status) - sum(status))*1.0) # total number of start codons unknown, False positives, increment courbe roc X
	except ZeroDivisionError:
		mF = 0

	distMin = 100000
	indexDisMin = 0
	xThreshold = 0
	yThreshold = 0
	
	for index, state in enumerate(status) :
		if state == 1:
			valueY = y[-1] + mT
			valueX = x[-1]

			y.append(valueY)
			x.append(valueX)
		else:
			valueX = x[-1] + mF
			valueY = y[-1]
			x.append(valueX)
			y.append(valueY)

		dist = math.sqrt(math.pow(0-valueX,2)+math.pow(1-valueY,2))
		if dist < distMin:
			indexDisMin = index
			distMin = dist
			xThreshold = valueX
			yThreshold = valueY
	
	auc = np.trapz(y,x)
	plt.figure()
	lw = 2
	plt.plot(x, y, color='darkorange', lw=lw, label='ROC curve (area = %0.5f)' % auc)
	plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic - ROC')
	plt.legend(loc="lower right")
	plt.savefig(saveFileCourbe)
	logging.info('Done!')
	return indexDisMin, xThreshold, yThreshold, auc
