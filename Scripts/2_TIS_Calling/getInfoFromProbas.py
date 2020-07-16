import time
import numpy as np
import pickle
import logging


def getProbasFromDictionnaryWithFrequencies(dicFrequencies, toSaveOutputs):
	logging.info('getProbasFromDictionnaryWithFrequencies')
	dicWithProbas = {}
	for i in dicFrequencies.keys():
		dicWithProbas[i] = {}
		sumValue = sum(dicFrequencies[i].values())
		for j in dicFrequencies[i].keys():
			dicWithProbas[i][j] = dicFrequencies[i][j]/(sumValue*1.0)
	toSaveFile = toSaveOutputs+'/probasByLenFragmentFromFrequencies.dic'
	with open(toSaveFile, 'wb') as fp:
		pickle.dump(dicWithProbas, fp, protocol=pickle.HIGHEST_PROTOCOL)
	logging.info('saveDic in: %s', toSaveFile)
	return dicWithProbas


def get_Median_Mean_Stdev(dicFrequencies, toSaveOutputs):
	
	logging.info('getMoyenne_Stdev')
	mean_median_std_info = {}
	keys = sorted(dicFrequencies.keys())

	for k in keys:
		dicAux = dicFrequencies[k]
		aux_freq = []
		for key in dicFrequencies[k].keys():
			aux_freq.extend([key] * dicAux[key])
		info = sum(dicAux.values()), np.mean(aux_freq), np.median(aux_freq), np.std(aux_freq)
		mean_median_std_info[k] = info
	
	logging.info('L_read\t#_reads_L \t Mean \t Median\t STD')

	keys = sorted(mean_median_std_info.keys())
	for key in keys:
		logging.info('%s \t %s ', str(key), str(mean_median_std_info[key]))

	toSaveFile = toSaveOutputs+'/infoStdev.dic'
	with open(toSaveFile, 'wb') as fp:
		pickle.dump(mean_median_std_info, fp, protocol=pickle.HIGHEST_PROTOCOL)
	logging.info('saveDic in: %s', toSaveFile)
	return mean_median_std_info


def getDispersion(mean_median_std_info, toSaveOutputs):

	# Heuristic H_1(l): assign a normalized weight to each read length (26-34 nt), 
	# computed through the standard deviation of the read positions acting as start codons
	logging.info('Get Dispersion')
	dic_weight = {}
	stdev = [i[3] for i in mean_median_std_info.values()]
	maxStdv = max(stdev)
	minStdv = min(stdev)
	keys = sorted(mean_median_std_info.keys())

	for key in keys:
	    valeur = (1 - (((mean_median_std_info[key][3] - minStdv)/(maxStdv-minStdv))*0.99))
	    dic_weight[key] = valeur

	logging.info(' Dic : %s ', str(dic_weight))

	toSaveFile = toSaveOutputs+'/weightsByLenght.dic'
	with open(toSaveFile, 'wb') as fp:
		pickle.dump(dic_weight, fp, protocol=pickle.HIGHEST_PROTOCOL)
	logging.info('saveDic weightsByLenght in: %s', toSaveFile)
	return dic_weight








