import warnings
import time
import sys
sys.path += sys.path + ["/u/ruizma/Jupyter/PipelineAligment/Scripts/PythonScripts", "/u/ruizma/Jupyter/PipelineAligment/Scripts/PythonScripts/General/",'/u/ruizma/lib/rabadb','/u/ruizma/lib/pygeno','/u/ruizma/Jupyter/ExploringRiboSeqReads/Scripts' , '/u/ruizma/miniconda2/lib/python2.7/site-packages']
from pyGeno.Genome import * 
warnings.filterwarnings("ignore")
from pyGeno.tools.UsefulFunctions import *
from pyGeno.tools.BinarySequence import *
import array
import math
import pickle
import itertools
import datetime
import logging
import getopt
from itertools import groupby
from operator import itemgetter
import getopt

class getProteinsSequences:

	def __init__(self, infoTranscripts1, infoTranscripts2, genesInfo1, genesInfo2, strand, folderToSave, lenghtSidesSequence, minLengProteine):
		self.infoTranscripts1 = infoTranscripts1
		self.genesInfo1 = genesInfo1
		self.infoTranscripts2 = infoTranscripts2
		self.genesInfo2 = genesInfo2
		self.folderToSave = folderToSave
		self.strand = strand
		self.lenghtSidesSequence = lenghtSidesSequence
		self.minLengProteine = minLengProteine
		self.total_info_proteins = 'Protein_id'+'\t'+'Protein_location'+'\t'+'TIS'+'\t'+'Codon'+'\t'+'Motif'+'\t'+'TPM_transcript'+'\t'+'Score_TIS'+'\t'+'Reference_Transcript'+'\t'+'StringTie_transcript'+'\t'+'Exons_range_protein'+'\t'+'Len_transcript'+'\t'+'Len_protein'+'\t'+'Variant_number'+'\t'+'Seq_proteine'+'\n'
		print 'getProteins OK'


	def getProteins(self):

		time0 = time.time()
		logging.info('Getting Proteins::')
		logging.info('Getting Proteins::lenghtSidesSequence %d', self.lenghtSidesSequence)
		logging.info('Getting Proteins::minLengProteine %d ', self.minLengProteine)
		logging.info('Folder to get infoTranscript1 %s ', self.infoTranscripts1)
		logging.info('Folder to get genes1 %s ', self.genesInfo1)
		logging.info('Folder to get infoTranscript2 %s ', self.infoTranscripts2)
		logging.info('Folder to get genes2 %s ', self.genesInfo2)

		with open (self.infoTranscripts1, 'rb') as fp:
			transcripts1 = pickle.load(fp)

		with open (self.genesInfo1, 'rb') as fp:
			genes1 = pickle.load(fp)

		with open (self.infoTranscripts2, 'rb') as fp:
			transcripts2 = pickle.load(fp)

		with open (self.genesInfo2, 'rb') as fp:
			genes2 = pickle.load(fp)

		logging.info('Summary Files : ')
		logging.info('infoTranscripts %s', self.infoTranscripts1)
		logging.info('genesInfo %s', self.genesInfo1)
		logging.info('infoTranscripts %s', self.infoTranscripts2)
		logging.info('genesInfo %s', self.genesInfo2)
		logging.info('\nInfo Transcripts and Genes collected ...')

		print 'Info Transcripts and Genes collected ...'
		self.startCodonsOrigin = {}

		logging.info('Total de Transcripts 1 Intersected %d', len(transcripts1.keys()))
		logging.info('Total de Genes 1 Intersected %d', len(genes1.keys()))

		logging.info('Total de Transcripts 2 Intersected %d', len(transcripts2.keys()))
		logging.info('Total de Genes 2 Intersected %d', len(genes2.keys()))

		self.proteinesToKeep = []
		gCont = 0
		totalGenes = len(genes1.keys())
		print '\nProcessing each gene of ', totalGenes ,' genes 1'
		logging.info('Processing each gene of %d genes', totalGenes)

		self.proteineId = 0
		geneNumber = 0

		self.getProteinsToKeep(genes1, transcripts1, '1')
		timeFinal = time.time()
		total = (timeFinal-time0) / 60
		logging.info('Total runtime getProteins first round%f min ', total)
		logging.info('Total proteines candidates to Evaluate first round %d', len(self.proteinesToKeep))
		logging.info('Start Codons Origin first round %s', str(self.startCodonsOrigin))

		self.getProteinsToKeep(genes2, transcripts2, '2')
				
		timeFinal = time.time()
		total = (timeFinal-time0) / 60

		logging.info('Total runtime getProteins final round %f min ', total)
		logging.info('Total proteines candidates to Evaluate final round %d', len(self.proteinesToKeep))
		logging.info('Start Codons Origin final round %s', str(self.startCodonsOrigin))

		print 'Total runtime getProteins : '+ str(total) +' min'
		print 'Total proteines candidates to Evaluate '+str(len(self.proteinesToKeep))+'\n'
		
		fileName = self.folderToSave+'/ProteinsCandidates_'+self.strand+'.dic'
		logging.info('Saving Total Proteins Candidates (dic) to %s', fileName)

		with open(fileName, 'wb') as handle:
			pickle.dump(self.proteinesToKeep, handle, protocol=pickle.HIGHEST_PROTOCOL)

		fileName = self.folderToSave+'/StartCodonsOrigin_'+self.strand+'.dic'
		logging.info('Saving Total Start Codons Proteins Candidates (dic) to %s', fileName)

		with open(fileName, 'wb') as handle:
			pickle.dump(self.startCodonsOrigin, handle, protocol=pickle.HIGHEST_PROTOCOL)
		
		fileProtsIntercepted = self.folderToSave+'/InfoProteinsKeept'+self.strand+'.gtf'
		logging.info('Saving gtf file of Proteins Intersected %s', fileProtsIntercepted)

		fileToSave = open(fileProtsIntercepted, 'w')
		fileToSave.write(self.total_info_proteins)
		fileToSave.close()


	def getProteinsToKeep(self, genes, transcripts, type_):

		time0 = time.time()
		geneNumber = 0
		totalGenes = len(genes.keys())
		print '\nProcessing each gene of ', totalGenes ,' genes '
		logging.info('Processing each gene of %d genes', totalGenes)

		positions_transcript = {}

		for keyG, transcriptsInGene in genes.items():
			geneNumber += 1
			logging.info('Total Transcripts in Gene %s : %d ', keyG, len(transcriptsInGene))
			print 'Total Transcripts in Gene ', keyG, ' : ', len(transcriptsInGene)
			print 'Gene ', geneNumber, '/', totalGenes, ' Proteins to Keep ', len(self.proteinesToKeep)
			logging.info('Gene %d/%d Proteins to Keep %d ', geneNumber, totalGenes, len(self.proteinesToKeep))
			proteinesInGene = []
			proteinesInAllTranscripts = []
			
			for indexTransName, transName in enumerate(transcriptsInGene):

				transcriptInfo = transcripts[transName]
				transName = transName+'_'+type_
				transcript = transcriptInfo[0]
				proteinesInAllFrames = []

				for fr in range(1,4):
					frame = transcriptInfo[fr]
					proteinesInFrameTranscript = []
					for indexSc , sc in enumerate(frame): 
					#['chrX:107605024-107605026', 'CTG', '+', proba, posIndex, kozak, positionsTrans, reference_transcript, tpm] 
					#[pos, codon, strand,proba, posIndex, kozak, positionsTrans, reference_transcript, tpm] 
						tis = sc[0]
						chrTis = sc[0].split(':')[0]
						sc_pos = sc[0]
						codon = sc[1]
						posIndex = sc[4]
						strand = sc[2]
						positionsTranscript = sc[6]
						
						sequenceTranscript_from_TIS = transcript[posIndex:]
						
						try:
							positions_list = positions_transcript[transName]
						except KeyError:
							positions_list = self.decode_positions(positionsTranscript)
							positions_transcript[transName] = positions_list

						positions_list = positions_list[posIndex:]
						nameProteineCandidat = [transName, sc, fr, self.proteineId,  tis+'-'+strand]
						variants, lenProteine = self.getTranslationAndVariants(chrTis, sequenceTranscript_from_TIS, nameProteineCandidat, codon, posIndex, positions_list)
						
						tmp_lss = self.lenghtSidesSequence
						if len(variants) > 2000:
							self.lenghtSidesSequence = tmp_lss - 5
							variants, lenProteine = self.getTranslationAndVariants(chrTis, sequenceTranscript_from_TIS, nameProteineCandidat, codon, posIndex, positions_list)
							self.lenghtSidesSequence = tmp_lss

						if len(variants) > 0:
							try:
								l = self.startCodonsOrigin[sc_pos]
								l[1].append(lenProteine)
							except KeyError:
								self.startCodonsOrigin[sc_pos] = [codon, [lenProteine]]
						
						self.proteineId += 1
						if len(proteinesInFrameTranscript) == 0:
							proteinesInFrameTranscript.extend(variants)
						else:
							for variant in variants:
								proteinesInFrameTranscript = self.removeRepeatedVariants(variant, proteinesInFrameTranscript, 'Same Frame '+str(fr))
					
					if len(proteinesInAllFrames) == 0:
						proteinesInAllFrames.extend(proteinesInFrameTranscript)
					else:
						for proteineInFrame in proteinesInFrameTranscript:
							proteinesInAllFrames = self.removeRepeatedVariants(proteineInFrame, proteinesInAllFrames, 'All Frames ')
				
				if len(proteinesInAllTranscripts) == 0:
					proteinesInAllTranscripts.extend(proteinesInAllFrames)
				else:
					for proteine in proteinesInAllFrames:
						proteinesInAllTranscripts = self.removeRepeatedVariants(proteine, proteinesInAllTranscripts, 'Same Transcript '+transName)
			
			if len(proteinesInGene) == 0:
				proteinesInGene.extend(proteinesInAllTranscripts)
			else:
				for proteine in proteinesInAllTranscripts:
					proteinesInGene = self.removeRepeatedVariants(proteine, proteinesInGene, 'Same Gene')
			
			if len(self.proteinesToKeep) == 0:
				self.proteinesToKeep.extend(proteinesInGene)
			else:
				for proteine in proteinesInGene:
					self.proteinesToKeep = self.removeRepeatedVariants(proteine, self.proteinesToKeep, 'All genes' )
			
			if geneNumber % 5000 == 0:
				timeFinal = time.time()
				total = (timeFinal-time0) / 60
				logging.info('Genes Observed : %d out of %d', geneNumber, totalGenes)
				print 'Genes Observed : ', geneNumber, ' out of ', totalGenes
				logging.info('Proteines keeped : %d', len(self.proteinesToKeep))
				print 'Proteines keeped : ', len(self.proteinesToKeep)
				logging.info('Time %f', total)
				print 'Time : ', total
	

	def getTranslationAndVariants(self, chr, sequenceTranscript_from_TIS, nameProteineCandidat, codon, posIndex, postitions_list):
		
		sequenceTranscript_from_TIS = sequenceTranscript_from_TIS.split('N')[0]
		if chr == 'chrM':
			proteineT = translateDNA(sequenceTranscript_from_TIS, frame = 'f1', translTable_id='mt')
		else:
			proteineT = translateDNA(sequenceTranscript_from_TIS, frame = 'f1', translTable_id='default')

		indexesCodonStop = [i for i,val in enumerate(proteineT) if val=='*']
		proteine = proteineT
		
		if indexesCodonStop:
			for indCS in indexesCodonStop:
				if len(indexesCodonStop) == 1 and indCS == 0 :
					proteine = ''
				elif indCS == 0 and proteine[1] != '/':
					proteine = proteine[:indCS]
					break
				elif indCS+1 < len(proteineT) and indCS-1 >= 0:
					if proteineT[indCS-1] != '/' and proteineT[indCS+1] != '/':
						proteine = proteine[:indCS]
						break
				elif indCS == len(proteine) - 1 and proteineT[indCS-1] != '/':
					proteine = proteine[:indCS]
					break
		if len(proteine) > 0:
			variants, lenProteine = self.getTranslations(proteine, nameProteineCandidat, sequenceTranscript_from_TIS, posIndex, postitions_list)
			return variants, lenProteine
		else:
			return [], 0


	def getTranslations(self, proteine, nameProteineCandidat, sequenceTranscript_from_TIS, posIndex, postitions_list):

		indexes = [i for i,val in enumerate(proteine) if val=='/']
		contVariant = 0

		# nameProteineCandidat = [transName, sc, fr, self.proteineId,  tis+'-'+strand]
		# [pos, codon, strand, proba, posIndex, kozak, positionsTrans, tpm] 

		chr = nameProteineCandidat[4].split(':')[0]
		tis = nameProteineCandidat[4]
		transcript =  nameProteineCandidat[0]
		codon = nameProteineCandidat[1][1]
		proba = nameProteineCandidat[1][3]
		kozak = nameProteineCandidat[1][5]
		tpm = str(nameProteineCandidat[1][7])

		if len(indexes) == 0:
			lenProteine = len(proteine)
			if lenProteine >= self.minLengProteine : 
				aaLen = lenProteine * 3
				lenTranscript = len(postitions_list[:aaLen+3])
				rangeTranscript = self.find_ranges(postitions_list[:aaLen+3])
				self.total_info_proteins += chr+'\t'+str(postitions_list[:aaLen+3][0])+'\t'+str(postitions_list[:aaLen+3][-1])+'\t'+tis+'\t'+codon+'\t'+kozak+'\t'+tpm+'\t'+str(proba)+'\t'+transcript+'\t'+rangeTranscript+'\t'+str(lenTranscript)+'\t'+str(lenProteine)+'\n'
				return [(proteine, nameProteineCandidat, contVariant, chr+':'+rangeTranscript, [] , lenProteine, lenTranscript, posIndex, (posIndex, posIndex+aaLen+3))], lenProteine
			else:
				return [], len(proteine)

		sequenceFirstProteine, listVariants, listPositionsVariants  = self.getListVariants(proteine, indexes)
		
		proteinesVariants = []
		lenProteine = len(sequenceFirstProteine)

		if lenProteine >= self.minLengProteine:
			f = 0
			stringFirstProteine = ''
			for p in sequenceFirstProteine:
				if p == '/':
					variant = list(listVariants[f])
					if variant[0] != '*':
						stringFirstProteine += variant[0]
					else:
						stringFirstProteine += variant[1]
					f += 1
				else:
					stringFirstProteine += p

			aaLen = len(stringFirstProteine) * 3
			if lenProteine != len(stringFirstProteine):
				logging.warning('Len sequences "first Proteins" not equal, verify!')
				logging.warning('Protein : %s', proteine)
				logging.warning('stringFirstProteine : %s', stringFirstProteine)
				logging.warning('sequenceFirstProteine : %s', sequenceFirstProteine)
				return
			if '*' in stringFirstProteine:
				logging.warning('Stop Codon in stringFirstProteine Sequence, verify! ')
				logging.warning('stringFirstProteine %s', stringFirstProteine)
				logging.warning('sequenceFirstProteine %s', sequenceFirstProteine)
				logging.warning('proteine %s', proteine)
				logging.warning('listVariants %s', str(listVariants))
				return

			lenTranscript = len(postitions_list[:aaLen+3])
			rangeTranscript = self.find_ranges(postitions_list[:aaLen+3])
			self.total_info_proteins += chr+'\t'+str(postitions_list[:aaLen+3][0])+'\t'+str(postitions_list[:aaLen+3][-1])+'\t'+tis+'\t'+codon+'\t'+kozak+'\t'+tpm+'\t'+str(proba)+'\t'+transcript+'\t'+rangeTranscript+'\t'+str(lenTranscript)+'\t'+str(lenProteine)+'\n'
		
			proteinesVariants.append((stringFirstProteine, nameProteineCandidat, contVariant, chr+':'+rangeTranscript, [], lenProteine, lenTranscript, posIndex, (posIndex, posIndex+aaLen+3)))
			contVariant += 1
			indexesPoly = [j for j,val in enumerate(sequenceFirstProteine) if val=='/']
			control = [False] * len(indexesPoly)
			
			for ind, value in enumerate(indexesPoly):
				if not control[ind]:

					variantsSmallSeq = []
					posPoly = 0
					addMoreTail = 0
					if value > self.lenghtSidesSequence :
						start = value - self.lenghtSidesSequence
						posPoly = self.lenghtSidesSequence
					else:
						start =  0	
						posPoly = value
						addMoreTail = self.lenghtSidesSequence - posPoly

					if len(sequenceFirstProteine)-value > self.lenghtSidesSequence + addMoreTail:
						if addMoreTail > 0 :
							end =  value + self.lenghtSidesSequence + addMoreTail
						else:
							end =  value + self.lenghtSidesSequence
					else:
						end = len(sequenceFirstProteine)
				 	
				 	smallSeqVariant = sequenceFirstProteine[start:end]
				 	
					while(True):
						if '/' in smallSeqVariant[-self.lenghtSidesSequence+2:]:
				 			if (end + self.lenghtSidesSequence) < len(sequenceFirstProteine):
					 			smallSeqVariant = sequenceFirstProteine[start:(end + self.lenghtSidesSequence)]
					 			end = end + self.lenghtSidesSequence
					 		else:
					 			smallSeqVariant = sequenceFirstProteine[start:]
					 			end = len(sequenceFirstProteine)
					 			break
					 	elif (end + self.lenghtSidesSequence) <= len(sequenceFirstProteine) and '/' in sequenceFirstProteine[end:(end + self.lenghtSidesSequence)]:
					 		smallSeqVariant = sequenceFirstProteine[start:(end + self.lenghtSidesSequence)]
					 		end = end + self.lenghtSidesSequence
					 	elif (end + self.lenghtSidesSequence) > len(sequenceFirstProteine) and '/' in sequenceFirstProteine[end:]:
					 		smallSeqVariant = sequenceFirstProteine[start:]
					 		end = len(sequenceFirstProteine)
					 		break
					 	else:
					 		break
					
					occu = smallSeqVariant.count("/")
					for i in range(ind, ind+occu):
						control[i] = True
						variantsSmallSeq.append(listVariants[i])

					posIndexV = posIndex + (start * 3)

					c = list(itertools.product(*variantsSmallSeq))
					posFinalPoly = [j for j,val in enumerate(smallSeqVariant) if val=='/']
					variants = []
					
					for i, pro in enumerate(c):
						copySmallSeqVariantFinal = list(smallSeqVariant)
						for index in range(0, len(posFinalPoly)):
							copySmallSeqVariantFinal[posFinalPoly[index]] = pro[index]

						proteine = "".join(copySmallSeqVariantFinal).split('*')[0]
						lenVariant = len(proteine)
						aaLen = lenVariant * 3
						
						if lenVariant >= self.minLengProteine:
							if proteine not in stringFirstProteine:
								if '*' in copySmallSeqVariantFinal:
									lenProt = ((posIndexV+aaLen) - posIndex)/3 
									lenTranscript = len(postitions_list[:lenProt*3])
									rangeTranscript = self.find_ranges(postitions_list[:(lenProt*3)+3])
									self.total_info_proteins += chr+'\t'+str(postitions_list[:(lenProt*3)+3][0])+'\t'+str(postitions_list[:(lenProt*3)+3][-1])+'\t'+tis+'\t'+codon+'\t'+kozak+'\t'+tpm+'\t'+str(proba)+'\t'+transcript+'\t'+rangeTranscript+'\t'+str(lenTranscript)+'\t'+str(lenProteine)+'\n'
								else:
									lenProt = lenProteine
								
								variants.append((proteine, nameProteineCandidat, contVariant, chr+':'+rangeTranscript, [], lenProt, lenTranscript, posIndex, (posIndexV, posIndexV+aaLen)))
								contVariant += 1
					
					proteinesVariants.extend(variants)
		
		return proteinesVariants, lenProteine


	def removeRepeatedVariants(self, proteine, variants, description):
		nameProteineCandidat = proteine[1]
		contained = False
		v = False
		for index, vari in enumerate(variants):
			if proteine[0] in vari[0]:
				variants[index][4].append(proteine)
				contained = True
				break
			elif vari[0] in proteine[0]:
				v = True
				proteine[4].append(vari)
				variants.remove(vari)
		if not contained:
			variants.append(proteine)
			
		return variants


	def getListVariants(self, proteine, indexes):
		
		prev = indexes[0]
		listPositionsVariants = [[prev]]
		fVariant = set()
		fVariant.add(proteine[prev-1])
		fVariant.add(proteine[prev+1])
		listVariants = [fVariant]
		sequenceFirstProteine = ''
		
		for i, n in enumerate(indexes):
			if i > 0 :
				if n - 2 == indexes[i-1]:
					l = listPositionsVariants[-1]
					l.append(n)
					l2 = listVariants[-1]
					l2.add(proteine[n-1])
					l2.add(proteine[n+1])
				else:
					listPositionsVariants.append([n])
					listVariants.append(set())
					l2 = listVariants[-1]
					l2.add(proteine[n-1])
					l2.add(proteine[n+1])
		
		variant = False
		for i in range(0, len(proteine)-1):
			if i == 0 and proteine[i+1] =='/':
				sequenceFirstProteine += '/'
				variant = True
			elif i == 0 and proteine[i+1] !='/':
				sequenceFirstProteine += proteine[i]
			elif i > 0:
				if not variant and proteine[i+1] != '/':
					sequenceFirstProteine += proteine[i]
				elif variant and proteine[i+1] != '/' and proteine[i] != '/':
					variant = False
				elif not variant and proteine[i+1] == '/':
					variant = True
					sequenceFirstProteine += '/'
		if proteine[-2] != '/':
			sequenceFirstProteine += proteine[-1]

		return sequenceFirstProteine, listVariants, listPositionsVariants

	
	def decode_positions(self, positionsTrans):
		postitions_list = []

		split_positions = positionsTrans.strip().split('|')
		for split in split_positions:
			split_pos = split.split('-')
			start = int(split_pos[0])
			end = int(split_pos[1])
			pos = range(start, end+1)
			if self.strand == '-':
				postitions_list = pos[::-1] + postitions_list
			else:
				postitions_list += pos
		return postitions_list


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
	folderToSave = ''
	logPath = ''
	fileNameTranscripts1 = ''
	fileNameGenes1 = ''
	fileNameTranscripts2 = ''
	fileNameGenes2 = ''

	try:
		opts, args = getopt.getopt(argv,"hs:t:g:x:y:f:l:",["strand=", "transcripts1=", "genes1=", "transcripts2=", "genes2=", 'folderToSave=', 'logPath='])
	except getopt.GetoptError:
		print 'main.py -s <strand> -t <transcripts1> -g <genes1> -x <transcripts2> -y <genes2> -f <folderToSave> -l <logPath>'
		sys.exit(2)

	print opts
	for opt, arg in opts:
		if opt == '-h':
			print 'main.py -s <strand> -t <transcripts1> -g <genes1> -x <transcripts2> -y <genes2> -f <folderToSave> -l <logPath>'
			sys.exit()
		elif opt in ("-s", "--strand"):
			strand = arg
		elif opt in ("-t", "--transcripts1"):
			fileNameTranscripts1 = arg
		elif opt in ("-x", "--transcripts2"):
			fileNameTranscripts2 = arg
		elif opt in ("-g", "--genes1"):
			fileNameGenes1 = arg
		elif opt in ("-y", "--genes2"):
			fileNameGenes2 = arg
		elif opt in ("-f", "--folderToSave"):
			folderToSave = arg
		elif opt in ("-l", "--logPath"):
			logPath = arg
	
	now = time.strftime('%y-%m-%d_%H:%M:%S', time.localtime())

	nameLog = logPath+'GetTotalInfoTranscriptsCandidates'+strand+'.log'
	print nameLog

	logging.info('Running Get Proteins Sequences')
	t0 = time.time()
	
	logging.info('Running Get Proteins Sequences')
	getProteinsSequencesInfo = getProteinsSequences(fileNameTranscripts1, fileNameTranscripts2, fileNameGenes1, fileNameGenes2, strand, folderToSave, 20, 8)
	getProteinsSequencesInfo.getProteins()
	logging.info('Got Proteins Sequences Information! \n')

	t2 = time.time()
	total = t2-t0
	logging.info('Total time run function getStartCodonsSavedByThreshold End : %d min', total/60)

if __name__ == "__main__":
	main(sys.argv[1:])










