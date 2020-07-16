import warnings
import time
import sys,getopt
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
from itertools import groupby
from operator import itemgetter


class getProteinsSequences:

	def __init__(self, infoTranscripts, strand, folderToSave, lenghtSidesSequence, minLengProteine):
		self.infoTranscripts = infoTranscripts
		self.folderToSave = folderToSave
		self.strand = strand
		self.lenghtSidesSequence = lenghtSidesSequence
		self.minLengProteine = minLengProteine
		self.toWriteTranscripts = ''
		self.total_info_proteins = ''
		self.total_info_proteins = 'Protein_id'+'\t'+'Protein_location'+'\t'+'TIS'+'\t'+'Codon'+'\t'+'Motif'+'\t'+'TPM_transcript'+'\t'+'Score_TIS'+'\t'+'Reference_Transcript'+'\t'+'StringTie_transcript'+'\t'+'Exons_range_protein'+'\t'+'Len_transcript'+'\t'+'Len_protein'+'\t'+'Variant_number'+'\t'+'Seq_proteine'+'\t'+'Seq_RNA\n'
		

	def getProteins(self):

		time0 = time.time()
		logging.info('Getting Proteins::')
		logging.info('Getting Proteins::lenghtSidesSequence %d', self.lenghtSidesSequence)
		logging.info('Getting Proteins::minLengProteine %d ', self.minLengProteine)
		logging.info('Folder to get infoTranscript %s ', self.infoTranscripts)
		
		with open (self.infoTranscripts, 'rb') as fp:
			transcripts = pickle.load(fp)

		logging.info('Summary Files : ')
		logging.info('infoTranscripts %s', self.infoTranscripts)
		logging.info('Info Transcripts collected ...')

		logging.info('Total Transcripts Intersected %d', len(transcripts.keys()))
		
		trans = 0
		proteineId = 0
		
		proteinesToKeep = []
		knownStartCodonsRetained = []
		transcriptsCanoniquesRetained = {}
		
		for key, into in transcripts.items():
			for info in into[1]:
				transcript = info[0]
				frame = info[1]
				chrTis = info[2]
				posIndex = info[3]
				transcriptName = info[4]
				codon = info[5]
				gene = info[6]
				tis = info[7]
				proba = info[8]
				sc = info[9]
				positionsTrans = info[10]
				
				transcriptsCanoniquesRetained[transcriptName] = [tis, frame]
				
				logging.info('-----------------Transcript------------- %s', transcriptName)
				
				nameProteineCandidat = [transcriptName, sc, frame, proteineId, tis]
				postitions_list = self.decode_positions(positionsTrans)

				variants, lenProteine = self.getTranslationAndVariants(chrTis, transcript, nameProteineCandidat, posIndex, postitions_list)
				if len(variants) > 2000:
					self.lenghtSidesSequence = 15
					variants, lenProteine = self.getTranslationAndVariants(chrTis, transcript, nameProteineCandidat, posIndex, postitions_list)
					self.lenghtSidesSequence = 20

				if len(variants) > 0:
					knownStartCodonsRetained.append((tis, self.strand, transcriptName))

				proteineId += 1
				
				if len(proteinesToKeep) == 0:
					proteinesToKeep = variants
				else:
					for proteine in variants:
						proteinesToKeep = self.removeRepeatedVariants(proteine, proteinesToKeep, 'Canonical' )
		
		timeFinal = time.time()
		total = (timeFinal-time0) / 60

		logging.info('Total runtime getProteins canonical %f min ', total)
		logging.info('Total proteins candidates to evaluate %d', len(proteinesToKeep))
		
		now = time.strftime('%y-%m-%d_%H:%M:%S', time.localtime())
		fileName = self.folderToSave+'/Proteins_Candidates_Canonical_'+self.strand+'.dic'
		logging.info('Saving Total Proteins Candidates Canonical (dic) to %s', fileName)

		with open(fileName, 'wb') as handle:
			pickle.dump(proteinesToKeep, handle, protocol=pickle.HIGHEST_PROTOCOL)

		fileName = self.folderToSave+'/Start_Codons_Retained_'+self.strand+'.dic'
		logging.info('Saving Start Codons Retained %s', fileName)

		with open(fileName, 'wb') as handle:
			pickle.dump(knownStartCodonsRetained, handle, protocol=pickle.HIGHEST_PROTOCOL)

		fileName = self.folderToSave+'/Transcripts_Canonical_Retained_'+self.strand+'.dic'
		logging.info('Transcripts Canonical Retained %s', fileName)

		with open(fileName, 'wb') as handle:
			pickle.dump(transcriptsCanoniquesRetained, handle, protocol=pickle.HIGHEST_PROTOCOL)

		fileProtsIntercepted = self.folderToSave+'/Info_Proteins_Kept'+self.strand+'.gtf'
		logging.info('Saving gtf file of Proteins Intersected %s', fileProtsIntercepted)

		fileToSave = open(fileProtsIntercepted, 'w')
		fileToSave.write(self.total_info_proteins)
		fileToSave.close()	

		
	def getTranslationAndVariants(self, chr, sequenceTranscript_from_TIS, nameProteineCandidat, posIndex, postitions_list):
		
		if chr == 'chrM':
			proteineT = translateDNA(sequenceTranscript_from_TIS, frame = 'f1', translTable_id='mt')
		else:
			proteineT = translateDNA(sequenceTranscript_from_TIS, frame = 'f1', translTable_id='default')

		indexesCodonStop = [i for i,val in enumerate(proteineT) if val=='*']
		proteine = proteineT
		
		if indexesCodonStop:
			for indCS in indexesCodonStop:
				if indCS == 0 and proteine[1] != '/':
					proteine = proteine[:indCS]
					break
				elif indCS+1 < len(proteineT) and indCS-1 >= 0:
					if proteineT[indCS-1] != '/' and proteineT[indCS+1] != '/':
						proteine = proteine[:indCS]
						break
				elif indCS == len(proteine) - 1 and proteineT[indCS-1] != '/':
					proteine = proteine[:indCS]
					break
		
		indexes = [i for i,val in enumerate(proteine) if val=='/']
		if len(proteine) - (len(indexes)*2) >= self.minLengProteine :
			variants, lenProteine = self.getTranslations(proteine, nameProteineCandidat, sequenceTranscript_from_TIS, posIndex, postitions_list, indexes)
			return variants, lenProteine
		else:
			return [], 0


	def getTranslations(self, proteine, nameProteineCandidat, sequenceTranscript_from_TIS, posIndex, postitions_list, indexes):

		contVariant = 0
		chr = nameProteineCandidat[4].split(':')[0]
		tis = nameProteineCandidat[4]
		transcript =  nameProteineCandidat[0]
		codon = nameProteineCandidat[1][1]
		proba = nameProteineCandidat[1][3]
		protein_id = nameProteineCandidat[3]
		kozak = nameProteineCandidat[1][5]
		tpm = str(nameProteineCandidat[1][8])
		

		if len(indexes) == 0:
			lenProteine = len(proteine)
			if lenProteine >= self.minLengProteine : 
				aaLen = lenProteine * 3
				lenTranscript = len(postitions_list[:aaLen+3])
				rangeTranscript = self.find_ranges(postitions_list[:aaLen+3])
				transcript_sequence = sequenceTranscript_from_TIS[:aaLen+3]

				self.total_info_proteins += str(protein_id)+'\t'+chr+':'+str(postitions_list[:aaLen+3][0])+'-'+str(postitions_list[:aaLen+3][-1])+'\t'+tis+'\t'+codon+'\t'+kozak+'\t'+tpm+'\t'+str(proba)+'\t'+transcript+ '\t'+transcript+'\t'+chr+':'+rangeTranscript+'\t'+str(lenTranscript)+'\t'+str(lenProteine)+'\t'+str(contVariant)+'\t'+proteine+'\t'+transcript_sequence+'\n'
				return [(proteine, nameProteineCandidat, contVariant, chr+':'+rangeTranscript, [], lenProteine, lenTranscript, posIndex, (posIndex, posIndex+aaLen+3))], lenProteine
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
			transcript_sequence = sequenceTranscript_from_TIS[:aaLen+3]

			self.total_info_proteins += str(protein_id)+'\t'+chr+':'+str(postitions_list[:aaLen+3][0])+'-'+str(postitions_list[:aaLen+3][-1])+'\t'+tis+'\t'+codon+'\t'+kozak+'\t'+tpm+'\t'+str(proba)+'\t'+transcript+ '\t'+transcript+'\t'+chr+':'+rangeTranscript+'\t'+str(lenTranscript)+'\t'+str(lenProteine)+'\t'+str(contVariant)+'\t'+stringFirstProteine+'\t'+transcript_sequence+'\n'
				
			proteinesVariants.append((stringFirstProteine, nameProteineCandidat, contVariant, chr+':'+rangeTranscript, [] , lenProteine, lenTranscript, posIndex, (posIndex, posIndex+aaLen+3)))
			
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
					variants_aux = []
					for i, pro in enumerate(c):
						copySmallSeqVariantFinal = list(smallSeqVariant)
						for index in range(0, len(posFinalPoly)):
							copySmallSeqVariantFinal[posFinalPoly[index]] = pro[index]

						proteine = "".join(copySmallSeqVariantFinal).split('*')[0]
						lenVariant = len(proteine)
						aaLen = lenVariant * 3

						if lenVariant >= self.minLengProteine:
							if proteine not in stringFirstProteine:
								if proteine not in variants_aux:
									if '*' in copySmallSeqVariantFinal:
										lenProt = ((posIndexV+aaLen) - posIndex)/3 
										lenTranscript = len(postitions_list[:lenProt*3])
										rangeTranscript = self.find_ranges(postitions_list[:(lenProt*3)+3])
										self.total_info_proteins += str(protein_id)+'\t'+chr+':'+str(postitions_list[:(lenProt*3)+3][0])+'-'+str(postitions_list[:(lenProt*3)+3][-1])+'\t'+tis+'\t'+codon+'\t'+kozak+'\t'+tpm+'\t'+str(proba)+'\t'+transcript+ '\t'+transcript+'\t'+chr+':'+rangeTranscript+'\t'+str(lenTranscript)+'\t'+str(lenProteine)+'\t'+str(contVariant)+'\t'+proteine+'\t'+''+'\n'
									else:
										lenProt = lenProteine
										aaLen_aux = len(stringFirstProteine) * 3
										self.total_info_proteins += str(protein_id)+'\t'+chr+':'+str(postitions_list[:aaLen_aux+3][0])+'-'+str(postitions_list[:aaLen_aux+3][-1])+'\t'+tis+'\t'+codon+'\t'+kozak+'\t'+tpm+'\t'+str(proba)+'\t'+transcript+ '\t'+transcript+'\t'+chr+':'+rangeTranscript+'\t'+str(lenTranscript)+'\t'+str(lenProteine)+'\t'+str(contVariant)+'\t'+proteine+'\t'+''+'\n'

									variants.append((proteine, nameProteineCandidat, contVariant, chr+':'+rangeTranscript, [], lenProt, lenTranscript, posIndex, (posIndexV, posIndexV+aaLen)))
									contVariant += 1
									variants_aux.append(proteine)

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
	fileNameTranscripts = ''
	fileSNPsFreeBayes = ''
	qualityThreshold = 20
	infoBedIntersectTranscripts = ''
	folderToSave = ''
	fastaFile = ''
	filepath_index = ''
	dicStartCodons = ''
	logPath = ''
	
	try:
		opts, args = getopt.getopt(argv,"hs:t:o:l:",["strand=", "fileNameTranscripts=", 'outputfile=','logPath='])
	except getopt.GetoptError:
		print 'main.py -s <strand> -t <fileNameTranscripts> -o <outputfile> -l <logPath> '
		sys.exit(2)

	print opts
	for opt, arg in opts:
		if opt == '-h':
			print 'main.py -s <strand> -t <fileNameTranscripts> -o <outputfile> -l <logPath> '
			sys.exit()
		elif opt in ("-s", "--strand"):
			strand = arg
		elif opt in ("-t", "--fileNameTranscripts"):
			fileNameTranscripts = arg
		elif opt in ("-o", "--outputfile"):
			folderToSave = arg
		elif opt in ("-l", "--logPath"):
			logPath = arg

	now = datetime.datetime.now()
	nameLog = logPath+'GetTotalInfoProteinsCandidates'+strand+'.log'
	logging.basicConfig(filename=nameLog, filemode='w', level=logging.INFO, format='%(asctime)s %(message)s')

	logging.info('Date %s ', str(now))
	logging.info('Running Get Proteins Sequences')
	t0 = time.time()
	getProteinsSequencesInfo = getProteinsSequences(fileNameTranscripts, strand, folderToSave, 20, 8)
	getProteinsSequencesInfo.getProteins()
	logging.info('Got Proteins Sequences Information!')
	t2 = time.time()
	total = t2-t0
	logging.info('Total time run function : %d min', total/60)


if __name__ == "__main__":
	main(sys.argv[1:])

