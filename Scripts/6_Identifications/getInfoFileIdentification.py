import warnings
import time
import sys
warnings.filterwarnings("ignore")
import math
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
matplotlib.style.use("ggplot")
import seaborn as sns
sns.set(style="whitegrid")
import pysam
import pandas as pd
	

class getInfoFileIdentification:

	def __init__(self, fileIdentification, nameDB):
		self.nameDB = nameDB
		self.scans, self.minScore, self.maxScore, self.totalPeptides = self.getListFromIdentifiedPeptidesFile(fileIdentification)
		self.drawHistogramme()
		print 'getInfoFileIndentification OK'


	def getListFromIdentifiedPeptidesFile(self, fileIdentification):

		print 'Running getListFromIdentifiedPeptidesFile ',self.nameDB
		scans = {}
		totalPeptides = 0
		minScore = 10000
		maxScore = 0
		peptides_by_lenght = {8:0, 9:0, 10:0, 11:0}
		capaDf = pd.read_csv(fileIdentification, header=1, skiprows=1)

		for index, row in capaDf.iterrows():
			file = row['File']
			scan = int(row['Scan number'])
			peptide = row['Peptide sequence']
			score = float(row['Score'])
			massError = float(row['Mass error (ppm)'])
			url = row['MS/MS url']

			try:
				proteineAccesion = row['Protein accessions / positions']
			except KeyError:
				proteineAccesion = 'NAN'

			if score < minScore:
				minScore = score 
			if score > maxScore:
				maxScore = score 

			try:
				info = scans[peptide]
				 
			except KeyError:
				scans[peptide] = [peptide, score, massError, scan, proteineAccesion, file, url]
				peptides_by_lenght[len(peptide)] += 1
		
		print 'Minimun Score ', minScore
		print 'Maximun Score ', maxScore
		print 'Total peptides ', len(scans)
		print 'Peptides by lenght ' , peptides_by_lenght
		return scans, minScore, maxScore, len(scans)


	def drawHistogramme(self):
		name = 'Distribution Score '+self.nameDB
		x = []
		for key, value in self.scans.items():
   			x.append(value[1])
		mu = np.mean(x)
		sigma = np.std(x)
		
		plt.figure(figsize=(8,5))
		n, bins, patches = plt.hist(x, bins='auto', facecolor='blue', alpha=0.5)
		
		title = '{0} $\mu= {1:.3f}$, $\sigma= {2:.3f}$ \n'.format(name, mu, sigma)
		plt.title(title)
		plt.xlabel('Score')
		plt.ylabel('Counts')
		plt.subplots_adjust(left=0.15)
		plt.show()