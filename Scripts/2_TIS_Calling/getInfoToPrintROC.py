import warnings
warnings.filterwarnings("ignore")
import time
import array
import math
import numpy as np
import re
import sys
import pickle
import datetime
import pandas as pd
import seaborn as sns
import logging


def validation(startCodonsCandidates, retainedStartsKnown, folder):
	
	
	with open (startCodonsCandidates, 'rb') as fp:
		startCodonsCandidates = pickle.load(fp)
	
	logging.info('Total startCodonsCandidates %d  ',len(startCodonsCandidates))
	
	dicStartCodonsCandidates =  getDicFromListSC(startCodonsCandidates)

	logging.info('Total codons detected %d  Total start codons known retained to validation %d ', len(dicStartCodonsCandidates), len(retainedStartsKnown))
	isHere = 0
	
	notFoundSC = ""
	foundSC = []
	found_sc_known_with_proba = {}
	for item in list(retainedStartsKnown):
		
		try:
			present = dicStartCodonsCandidates[item]
			foundSC.append(present)
			found_sc_known_with_proba[item[:-4]] = present
			isHere += 1
		except KeyError:
			notFoundSC += item+'\n'
	
	logging.info('Codons start known found %d Codons start known not found %d ', isHere,  (len(retainedStartsKnown) - isHere))
	with open(folder+'/Canonical_SC.dic', 'wb') as handle:
		pickle.dump(found_sc_known_with_proba, handle, protocol=pickle.HIGHEST_PROTOCOL)


def getDicFromListSC(itemlist):

	listeItems = {}
	
	for index, item in enumerate(itemlist):
		proba = item[0]
		countReads = item[4]
		strand = item[5]
		codon = item[3]
		infoStart = item[1]+'-'+strand+'-'+codon
		listeItems[infoStart] = item

	return listeItems


def getCourbeRoc(startCodonsCandidates, retainedStartsKnown):
	
	isHere = 0
	with open (startCodonsCandidates, 'rb') as fp:
		startCodonsCandidates = pickle.load(fp)

	coverage = []

	retainedSC = {}
	for item in retainedStartsKnown:
		retainedSC[item] = item
	
	print 'Retained SC known', len(retainedSC)

	status = []
	startCodonsKnown = []
	indexStartCodonsKnown = []
	totalList_with_SC_Reads_OverMinimum = []
	cont = 0
	
	for index, item in enumerate(startCodonsCandidates):
		proba = item[0]
		countReads = item[4]
		strand = item[5]
		cont += 1
		coverage.append(countReads)
		codon = item[3]
		infoStart = item[1]+'-'+strand+'-'+codon
	
		try:
			present = retainedSC[infoStart] 
			status.append(1)
			isHere += 1
			startCodonsKnown.append(item[5])
			indexStartCodonsKnown.append(cont)
		except  KeyError:
			status.append(0)

		totalList_with_SC_Reads_OverMinimum.append(item)
	
	logging.info('Total start codons : %d  Status : %d ', len(totalList_with_SC_Reads_OverMinimum), len(status))
	return status, totalList_with_SC_Reads_OverMinimum, isHere, len(retainedSC), coverage, startCodonsKnown, indexStartCodonsKnown


def createDicForStartCodons():
	dic = {}
	listCodonsStart = [ 'ATG','CTG', 'TTG', 'GTG', 'ACG', 'ATA', 'ATT', 'ATC', 'AAG', 'AGG']
	for sc in listCodonsStart:
		dic[sc] = []
	return dic

def getStartCodonsSavedByThreshold(startCodonsCandidats, folderToSave, file):

	if file:
		with open (startCodonsCandidats, 'rb') as fp:
			startCodonsRetained = pickle.load(fp)
	else:
		startCodonsRetained = startCodonsCandidats

	dic = {}
	probas = []
	codonsOrder = []
	dicCodons = createDicForStartCodons()

	for idx, item in enumerate(startCodonsRetained):
		proba = item[0]
		codon = item[3]
		probas.append(proba)
		codonsOrder.append(codon)
		dicCodons[codon].append(proba)

		try:
			dic[codon] += 1
			
		except KeyError:
			dic[codon] = 1
	
	getBedFileStartCodonsCandidats(startCodonsCandidats, folderToSave, file)
	

def getBedFileStartCodonsCandidats(startCodonsCandidatsFile, folderToSave, file):

	if file:
		with open (startCodonsCandidatsFile, 'rb') as fp:
			itemlist1 = pickle.load(fp)
	else:
		itemlist1 = startCodonsCandidatsFile

	fileForward = open(folderToSave+'/BedFileForward.bed', 'w') 
	fileBackward = open(folderToSave+'/BedFileBackward.bed', 'w') 

	toWriteForward = ''
	toWriteBackward = ''

	
	for item in itemlist1:
		proba = item[0]
		countReads = item[4]
		strand = item[5]
		codon = item[3]
		chr = item[1].split(':')[0]
		ini = item[1].split(':')[1].split('-')[0]
		fini =  item[1].split(':')[1].split('-')[1]
		
		if strand == '+':
			toWriteForward += chr+"\t"+ini+"\t"+fini+'\n'
		else:
			toWriteBackward += chr+"\t"+ini+"\t"+fini+'\n'

	fileForward.write(toWriteForward)
	fileForward.close()

	fileBackward.write(toWriteBackward)
	fileBackward.close()



