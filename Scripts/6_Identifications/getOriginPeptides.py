import warnings
import time
import sys
sys.path += sys.path + ["../../Scripts/4_Get_Active_Transcripts/Canonical_Proteins/"]
warnings.filterwarnings("ignore")
import array
import math
import pickle
import itertools
import datetime
import logging
import csv
import numpy as np
import pandas as pd
import getInfoIdentifications as getInfo
from collections import Counter
from collections import OrderedDict
from operator import itemgetter
import getopt
import pysam
from pyGeno.Genome import * 
from pyGeno.tools.UsefulFunctions import *
import getMCS_from_aa_seq as getMCS
from Bio import pairwise2
import textdistance
from difflib import SequenceMatcher
from getTranscriptInformation import getTranscriptInformation

class Origin_peptides:

	def __init__(self, peptides_bed_file, peptides_bed_file_intersected, infoTranscripts, strand, genome, genome_index, name, fileSNPsFreeBayes, qualityThreshold):
		self.peptides_bed_file = peptides_bed_file
		self.peptides_bed_file_intersected = peptides_bed_file_intersected
		self.infoTranscripts = infoTranscripts
		self.strand = strand
		self.genome = genome
		self.genome_index = genome_index
		self.faFile = pysam.FastaFile(genome, genome_index)
		self.name = name
		self.fileSNPsFreeBayes = fileSNPsFreeBayes
		self.qualityThreshold = qualityThreshold
		self.references_transcript_info = {}


	def get_origin_peptides(self):

		total_peptides = {}
		peptides = set()

		with open(self.peptides_bed_file) as f:
			for index, line in enumerate(f):
				splitLine = line.strip().split('\t')
				chr = splitLine[0]
				start = splitLine[1]
				end = splitLine[2]
				serie = int(splitLine[3])
				peptide = splitLine[4]
				peptides.add(peptide)

				len_proteine = int(splitLine[9])
				transcript = splitLine[5]
				rangeTranscript = splitLine[6] 
				codon = splitLine[8]
				tis_string = splitLine[7]
				bio = splitLine[11]
				tpm = float(splitLine[10])
				sequence =  splitLine[12]
				try:
					score = float(splitLine[13]) 
				except ValueError:
					score = 0
				index_protein = int(splitLine[14])

				key_peptide_position = chr+':'+start+'-'+end
				info_peptide = [ transcript, rangeTranscript, tis_string, codon, len_proteine, tpm, bio, sequence, score, index_protein] 
				try:
					positions_peptides = total_peptides[peptide][1]
					dic_trans_serie = {}
					positions_peptides[key_peptide_position] = dic_trans_serie
				except KeyError:
					dic_trans_serie = {}
					total_peptides[peptide] = [info_peptide, {key_peptide_position: dic_trans_serie}]

		with open(self.peptides_bed_file_intersected) as f:
			for index, line in enumerate(f):
				splitLine = line.strip().split('\t')
				chr = splitLine[0]
				start = splitLine[1]
				end = splitLine[2]
				peptide = splitLine[4]
				transcript = ''
				transcript_support_level = ''
				gene_type = ''

				try:
					transcript = splitLine[23].split(' transcript_id ')[1].split('\"')[1]
					transcript_support_level = splitLine[23].split(' transcript_support_level ')[1].split('\"')[1]
					gene_type = splitLine[23].split(' gene_type ')[1].split('\"')[1]
				except IndexError:
					pass

				key_peptide_position = chr+':'+start+'-'+end
				
				try:
					positions_peptides = total_peptides[peptide]
					dic_trans_serie = positions_peptides[1][key_peptide_position]
					dic_trans_serie[transcript] = [transcript_support_level, gene_type]
					
				except KeyError:
					pass
				
		print 'Total peptides in ',self.strand,' strand =', len(total_peptides)

		origin_by_type, origin_by_peptide, origin_of_each_peptide, final_origin_alt_tis_peptides = self.get_origin_information_for_peptide(total_peptides)
		return origin_by_type, origin_by_peptide, origin_of_each_peptide, final_origin_alt_tis_peptides



	def get_origin_information_for_peptide(self, total_peptides):

		for peptide, info_intersection in total_peptides.items():

			info_peptide = info_intersection[0]
			information_series = info_intersection[1]
			if len(information_series) == 1:
				total_peptides = self.get_biotype_from_transcripts(information_series.keys(), information_series.values()[0], total_peptides, peptide)
			else:
				keys = information_series.keys()
				transcripts_keys = []
				merged_dic = []

				for key in keys:
					transcripts_keys.append(information_series[key].keys())
					merged_dic += information_series[key].items()

				merged_dic = dict(merged_dic)

				set_intersection_keys = set(transcripts_keys[0])

				for keys_ in transcripts_keys:
					set_intersection_keys = set_intersection_keys.intersection(keys_)

				dic_keeped_trans = {}
				for key in information_series:
					dic_keeped_trans[key] = {}

				for key in set_intersection_keys:
					for key_ in dic_keeped_trans.keys():
						dic_keeped_trans[key_][key] = merged_dic[key]

				total_peptides[peptide][1] = dic_keeped_trans
				
				total_peptides = self.get_biotype_from_transcripts(dic_keeped_trans.keys(), dic_keeped_trans.values()[0], total_peptides, peptide)


		origin_by_type, origin_by_peptide, origin_of_each_peptide, final_origin_alt_tis_peptides = self.get_final_origin_each_peptide(total_peptides)
		return origin_by_type, origin_by_peptide, origin_of_each_peptide, final_origin_alt_tis_peptides


	def get_biotype_from_transcripts(self, info_transcripts_keys, info_transcripts_values, total_peptides, peptide):

		set_to_set = []

		for transcript, information_trans in info_transcripts_values.items():
			
			transcript_support_level = information_trans[0]
			gene_type = information_trans[1]
			
			if (transcript_support_level == 'NA' and gene_type == 'protein_coding') or transcript_support_level != 'NA':

				for position in info_transcripts_keys:
					set_to_set_aux = []
					chr = position.split(':')[0]
					start = int(position.split(':')[1].split('-')[0])
					end = int(position.split(':')[1].split('-')[1])
					info_transcript = self.infoTranscripts[self.strand][transcript]

					if len(info_transcript['CDS']) == 0:
						keys = ['Exons', '3UTR', '5UTR', 'Introns']
					else:
						keys = ['3UTR', '5UTR', 'Introns', 'CDS']

					in_cds = False

					for key in keys:
						if key == 'CDS':
							if len(info_transcript['CDS']) > 0:
								for cds in info_transcript['CDS']:
									start_elem = int(cds[0])
									end_elem = int(cds[1])
									if start_elem <= start <= end_elem and start_elem <= end <= end_elem:
										in_cds = True
										break
							if in_cds:
								chr = info_transcript['Info'][0]
								regions = info_transcript['CDS']
								protein = self.getTranscript_and_protein( chr, regions, transcript, self.strand)
								
								if peptide in protein:
									try:
										peptide_info_dic = self.references_transcript_info['Alt-TIS']
										try:
											peptide_info_dic[peptide].append((transcript, info_transcript, protein))
										except KeyError:
											peptide_info_dic[peptide] = [(transcript, info_transcript, protein)]
									except KeyError:
										self.references_transcript_info['Alt-TIS'] = {peptide: [(transcript, info_transcript, protein)]}

									if 'Alt-TIS' not in set_to_set_aux:
										set_to_set_aux.append('Alt-TIS')
								else:
									into = self.get_presence_peptide_in_protein(protein, peptide)
									if into:
										try:
											peptide_info_dic = self.references_transcript_info['Alt-TIS']
											try:
												peptide_info_dic[peptide].append((transcript, info_transcript, protein))
											except KeyError:
												peptide_info_dic[peptide] = [(transcript, info_transcript, protein)]
										except KeyError:
											self.references_transcript_info['Alt-TIS'] = {peptide: [(transcript, info_transcript, protein)]}
											
										if 'Alt-TIS' not in set_to_set_aux:
											set_to_set_aux.append('Alt-TIS')
									elif 'Frameshift' not in set_to_set_aux:
										set_to_set_aux.append('Frameshift')
										try:
											peptide_info_dic = self.references_transcript_info['Frameshift']
											try:
												peptide_info_dic[peptide].append((transcript, info_transcript, protein))
											except KeyError:
												peptide_info_dic[peptide] = [(transcript, info_transcript, protein)]
										except KeyError:
											self.references_transcript_info['Frameshift'] = {peptide: [(transcript, info_transcript, protein)]}
						else:

							for elem in info_transcript[key]:
								start_elem = int(elem[0])
								end_elem = int(elem[1])
								if start_elem <= start <= end_elem and start_elem <= end <= end_elem:
									try:
										peptide_info_dic = self.references_transcript_info[key]
										try:
											peptide_info_dic[peptide].append((transcript, info_transcript))
										except KeyError:
											peptide_info_dic[peptide] = [(transcript, info_transcript)]
									except KeyError:
										self.references_transcript_info[key] = {peptide: [(transcript, info_transcript)]}
									
									if key not in set_to_set_aux:
										set_to_set_aux.append(key)
									break

					if len(set_to_set_aux) == 0:
						keys = ['3UTR', '5UTR', 'Introns']
						try:
							for key in keys:
								for elem in info_transcript[key]:
									start_elem = int(elem[0])
									end_elem = int(elem[1])
									
									if self.strand == '-' and key == '5UTR':
										if end >= start_elem and end <= end_elem:
											try:
												peptide_info_dic = self.references_transcript_info[key]
												try:
													peptide_info_dic[peptide].append((transcript, info_transcript))
												except KeyError:
													peptide_info_dic[peptide] = [(transcript, info_transcript)]
											except KeyError:
												self.references_transcript_info[key] = {peptide: [(transcript, info_transcript)]}

											if key not in set_to_set_aux:
												set_to_set_aux.append(key)
											break

									elif self.strand == '+' and key == '5UTR':
										if start >= start_elem and start <= end_elem:
											try:
												peptide_info_dic = self.references_transcript_info[key]
												try:
													peptide_info_dic[peptide].append((transcript, info_transcript))
												except KeyError:
													peptide_info_dic[peptide] = [(transcript, info_transcript)]
											except KeyError:
												self.references_transcript_info[key] = {peptide: [(transcript, info_transcript)]}
											if key not in set_to_set_aux:
												set_to_set_aux.append(key)
											break

									if self.strand == '-' and key == '3UTR':
										if start >= start_elem and start <= end_elem:
											try:
												peptide_info_dic = self.references_transcript_info[key]
												try:
													peptide_info_dic[peptide].append((transcript, info_transcript))
												except KeyError:
													peptide_info_dic[peptide] = [(transcript, info_transcript)]
											except KeyError:
												self.references_transcript_info[key] = {peptide: [(transcript, info_transcript)]}
											if key not in set_to_set_aux:
												set_to_set_aux.append(key)
											break

									elif self.strand == '+' and key == '3UTR':
										if end >= start_elem and end <= end_elem:
											try:
												peptide_info_dic = self.references_transcript_info[key]
												try:
													peptide_info_dic[peptide].append((transcript, info_transcript))
												except KeyError:
													peptide_info_dic[peptide] = [(transcript, info_transcript)]
											except KeyError:
												self.references_transcript_info[key] = {peptide: [(transcript, info_transcript)]}
											if key not in set_to_set_aux:
												set_to_set_aux.append(key)
											break

									if key == 'Introns':
										if (start >= start_elem and start <= end_elem) or (end >= start_elem and end <= end_elem):
											try:
												peptide_info_dic = self.references_transcript_info[key]
												try:
													peptide_info_dic[peptide].append((transcript, info_transcript))
												except KeyError:
													peptide_info_dic[peptide] = [(transcript, info_transcript)]
											except KeyError:
												self.references_transcript_info[key] = {peptide: [(transcript, info_transcript)]}
											if key not in set_to_set_aux:
												set_to_set_aux.append(key)
											break
						except KeyError:
							pass

				if len(set_to_set_aux) > 1:
					origin = ''
					for indx, ori in enumerate(set_to_set_aux):
						if indx == 0:
							origin += ori
						else:
							origin += "|" + ori
					set_to_set.append(set_to_set_aux[0])
				else:
					set_to_set.extend(set_to_set_aux)
		
		total_peptides[peptide].append(set_to_set)
		return total_peptides

	def get_presence_peptide_in_protein(self, protein, peptide):
		if len(protein) > len(peptide):
			len_peptide = len(peptide)
			for ini in range(0, len(protein)-len_peptide):
				matchs = 0
				sub_seq = protein[ini:ini+len_peptide]
				for i,c in enumerate(sub_seq):
					if peptide[i] == c:
						matchs += 1
				if matchs == len_peptide - 1:
					return True
			return False
		else:
			False

	def get_final_origin_each_peptide(self, total_peptides):

		origin_by_type = {'3UTR':0, '5UTR':0, 'Exons':0, 'Alt-TIS': 0, 'Frameshift': 0, 'Intergenic':0, 'Introns':0}
		origin_by_peptide = {}
		origin_of_each_peptide = {}
		final_origin_alt_tis_peptides = {} 
		

		for peptide, value in total_peptides.items():
			cont = 0
			info_peptide = value[0]
			bio = info_peptide[6]
			positions =  value[1].keys()
			set_to_set = value[2]
			count = []
			max_key = ''

			positions_concat = ''
			for index, position in enumerate(positions):
				if index == 0:
					positions_concat += position
				else:
					positions_concat += '|'+position
			
			if len(set_to_set) > 0:
				count = Counter(set_to_set)
				d = OrderedDict(sorted(count.items(), key=itemgetter(1), reverse=True))
				max_count = d[d.keys()[0]]
				max_key = d.keys()[0]

				if max_key == 'Introns' :
					info = count['CDS']
					if info > 0:
						max_key = 'CDS'
					else:
						info = count['Exons']
						if info > 0:
							max_key = 'Exons'

				if bio != max_key:
					info = count[bio]
					if info > 0:
						max_key = bio

				if bio == '3UTR' and max_key != '3UTR'  :
					info = count['3UTR']
					if info > 0:
						max_key = '3UTR'

				if bio == '5UTR' and max_key != '5UTR'  :
					info = count['5UTR']
					if info > 0:
						max_key = '5UTR'

				if bio == 'Alt-TIS' and max_key != 'Alt-TIS'  :
					info = count['Alt-TIS']
					if info > 0:
						max_key = 'Alt-TIS'

				if bio not in  origin_by_type.keys() and bio != 'ToSearch':
					info = count['Exons']
					if info > 0:
						max_key = 'Exons'

				if bio == 'ToSearch'  :
					info = count['Alt-TIS']
					if info > 0:
						max_key = 'Alt-TIS'

			if max_key == '':
				max_key = 'Intergenic'
			
			origin_by_type[max_key] += 1

			try:
				info = origin_by_peptide[max_key] 
				info.append((peptide, info_peptide, positions_concat, self.strand, count))
			except KeyError:
				origin_by_peptide[max_key]  = [(peptide, info_peptide, positions_concat, self.strand, count)]

			transcript = info_peptide[0]
			sequence = info_peptide[7]

			if max_key != 'Intergenic':
				info_pept = self.references_transcript_info[max_key][peptide]
				final_origin_alt_tis_peptides[peptide] = info_pept
				reference_transcript = ''

				if 'ENS' not in transcript:
					if max_key == 'Alt-TIS':
						names_transcript = ''
						for i, prot_inf in enumerate(info_pept):
							names_transcript += prot_inf[0]+'|'
							if sequence in prot_inf[2]:
								reference_transcript += prot_inf[0]+'|'
						if len(reference_transcript) > 0:
							reference_transcript = reference_transcript[:-1]
						else:
							reference_transcript = names_transcript[:-1]
					else:
						names_transcript = ''
						for prot_inf in info_pept:
							names_transcript += prot_inf[0]+'|'
						reference_transcript = names_transcript[:-1]
				else:
					reference_transcript = transcript
			else:
				reference_transcript = transcript
			info_peptide.append(reference_transcript)
			origin_of_each_peptide[peptide] = (peptide, info_peptide, positions_concat, self.strand, count, max_key)
		return origin_by_type, origin_by_peptide, origin_of_each_peptide, final_origin_alt_tis_peptides


	def getTranscript_and_protein(self, chr, regions, transcript, strandTranscript):

		sequenceTranscript = ''
		for cds in regions:
			startExon = cds[0]
			endExon = cds[1]
			sequence = self.faFile.fetch(chr,startExon-1,endExon)

			if strandTranscript == '-' :
				sequence = reverseComplement(sequence)
				sequenceTranscript = sequenceTranscript + sequence
				
			elif strandTranscript == '+':
				sequenceTranscript =  sequenceTranscript + sequence

		if chr == 'chrM':
			proteineT = translateDNA(sequenceTranscript, frame = 'f1', translTable_id='mt')
		else:
			proteineT = translateDNA(sequenceTranscript, frame = 'f1', translTable_id='default')

		if '/' in proteineT:
			print proteineT
		return proteineT

	def getSNPsDic(self):

		time0 = time.time()
		dic = {}
		totalSnps = 0
		logging.info('Get SNPs Info::')
		logging.info('Get SNPs Info::File %s ', self.fileSNPsFreeBayes)

		with open(self.fileSNPsFreeBayes) as f:
			for index, snp in enumerate(f): # 1	1	19216	19217	G	G,C	0.141377	freebayes
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
		
		print 'Total SNPs Over the qualityThreshold ',self.qualityThreshold, ' : ', str(totalSnps)
		logging.info('Total SNPs Over the qualityThreshold  %d : %d ', self.qualityThreshold, totalSnps)
		timeFinal = time.time()
		total = (timeFinal-time0) / 60
		print  'Total runtime getSNPsDic : '+ str(total) +' min'
		logging.info('Total runtime getSNPsDic : %f min \n', total)
		return dic, totalSnps

	def writeCSV(self, peptidesToWrite, pathToSave, name, reference_transcript_peptides_backward, reference_transcript_peptides_forward ):
		import csv
		snpsDic, totalSnps = self.getSNPsDic()
		getTransInfo = getTranscriptInformation(snpsDic, self.genome, self.genome_index)

		sequences_final = {}
		csvData = [['Origin', 'Peptide', 'Position_peptide', 'Strand', 'Transcript', 'Reference Transcript', 'Range_transcript', 'Tis_Transcript', 'SC', 'Len_Protein', 'MS_Score', 'TPM','Bio', 'Type Transcript', 'Index_Protein','Count', 'Sequence', 'RNA_Sequence']]
		
		csvData_info_prot = [['Protein','Origin', 'Peptides', 'Position_peptides', 'Strand', 'Transcript', 'Reference Transcript', 'Range_transcript', 'Tis_Transcript', 'SC', 'Len_Protein', 'MS_Score','TPM', 'Type Transcript', 'Index_Protein', 'Sequence', 'RNA_Sequence']]
		
		big_proteins = 0
		info_proteins = {}
		peptides_by_protein = {}
		small_proteins = 0
		origin_by_protein = {'3UTR':0, '5UTR':0, 'Exons':0, 'Alt-TIS': 0, 'Frameshift': 0, 'Intergenic':0, 'Introns':0}
		len_proteins = []
		len_proteins_by_type = {}
		total_proteins = 0
		tpm_proteins = {}
		to_write_alt_tis = ''
		cont_found_ref_trans_alt_tis = 0

		for origin, info in peptidesToWrite.items():

			for peptide_info in info:
				peptide = peptide_info[0]

				info = peptide_info[1]
				transcript = info[0]
				rangeTranscript = info[1]
				tis_string = info[2]
				sc = info[3]
				len_proteine = info[4]
				tpm = float(info[5])
				bio = info[6]
				sequence =  info[7]
				key_prot = transcript+'-'+tis_string
				index_protein = info[9]

				if len(sequence) != int(len_proteine):
					print 'Differente Lens probably because variant ?', len(sequence), int(len_proteine)

				try:
					info_prot = sequences_final[sequence]
					info_prot[0].append(info)
					info_prot[1].append(key_prot)
					info_prot[2].append(peptide)
				except KeyError:
					sequences_final[sequence] = [[info],[key_prot], [peptide]]

				
				peptide_score = info[8]
				
				strand = peptide_info[3]
				positions_peptides = peptide_info[2]
				count = peptide_info[4]

				try:
					info_transcript = self.infoTranscripts[strand][transcript]
					type_transcript = info_transcript['Info'][3]
				except KeyError:
					type_transcript = 'Novel Transcript'
				
				transcript_sequence = getTransInfo.getTranscriptsRegionsInformation(rangeTranscript, strand)

				reference_transcript = ''
				if 'ENS' not in transcript:
					if origin != 'Intergenic':
						dic_aux = reference_transcript_peptides_forward
						if strand == '-':
							dic_aux = reference_transcript_peptides_backward
						to_print = dic_aux[peptide]

						if origin == 'Alt-TIS':
							names_transcript = ''
							for i, prot_inf in enumerate(to_print):
								names_transcript += prot_inf[0]+'|'
								if sequence in prot_inf[2]:
									reference_transcript += prot_inf[0]+'|'
									cont_found_ref_trans_alt_tis += 1
								if i == 0:
									to_write_alt_tis += '\n'+peptide+' '+type_transcript+' '+strand+'\t'+sequence+'\n'
								to_write_alt_tis += prot_inf[0]+'\t'+prot_inf[2]+'\t'+str(sequence in prot_inf[2])+'\n'
							if len(reference_transcript) > 0:
								reference_transcript = reference_transcript[:-1]
							else:
								reference_transcript = names_transcript[:-1]
						else:
							names_transcript = ''
							for prot_inf in to_print:
								names_transcript += prot_inf[0]+'|'
							reference_transcript = names_transcript[:-1]
					else:
						reference_transcript = transcript
				else:
					reference_transcript = transcript

				to_add = [origin, peptide, positions_peptides, strand, transcript, reference_transcript, rangeTranscript, tis_string, sc, len_proteine, peptide_score, tpm, bio, type_transcript, index_protein, str(count), sequence, transcript_sequence]
				
				csvData.append(to_add)
				
				try:
					info_prot = info_proteins[key_prot]
					info_prot[1].append(origin)
					info_prot[2].append(peptide)
					info_prot[3].append(positions_peptides)
					info_prot[11].append(peptide_score)
					info_prot[13].append(bio)
					info_prot[16].append(str(count))
				except KeyError:
					total_proteins += 1
					info_proteins[key_prot] = [key_prot, [origin], [peptide], [positions_peptides], strand, transcript, reference_transcript, rangeTranscript, tis_string, sc, len_proteine, [peptide_score], tpm, [bio], type_transcript, index_protein, [str(count)], sequence, transcript_sequence ]
					try:
						tpm_add = tpm_proteins[origin]
						tpm_add.append(tpm)
					except KeyError:
						tpm_proteins[origin] = [tpm]

				peptides_by_protein[peptide] = [key_prot, origin, peptide, positions_peptides, strand, transcript, reference_transcript, rangeTranscript, tis_string, sc, len_proteine, peptide_score, tpm, [bio], type_transcript, index_protein, str(count), sequence, transcript_sequence ]

		print 'Total Proteins ', total_proteins
		print 'Total set Proteins ', len(sequences_final)
		
		with open(pathToSave+'peptides_origin_'+name+'.csv', 'w') as csvFile:
			writer = csv.writer(csvFile)
			writer.writerows(csvData)

		csvFile.close()

		with open(pathToSave+'/peptides_by_protein_'+name+'.dic', 'wb') as handle:
			pickle.dump(peptides_by_protein, handle, protocol=pickle.HIGHEST_PROTOCOL)

		small_protein = ''
		big_protein = ''
		supported_peptides_by_protein = {}
		peptides_small_proteins = 0
		peptides_big_proteins = 0
		
		file_to_save = open(pathToSave+'Get_info_novel_isoform_'+name+'.info', 'w')
		file_to_save.write(to_write_alt_tis)
		file_to_save.close()

		for protein, info_protein in info_proteins.items():
			n_peptides = len(info_protein[2])
			if n_peptides > 1:
				supported_peptides_by_protein[protein] = n_peptides
			origin_count = Counter(info_protein[1])
			d = OrderedDict(sorted(origin_count.items(), key=itemgetter(1), reverse=True))
			max_count = d[d.keys()[0]]
			max_key = d.keys()[0]
			len_prot = int(info_protein[10])

			if len(origin_count) > 1:
				print origin_count, protein, max_key
			
			len_proteins.append(len_prot)
			try:
				len_proteins_by_type[max_key].append(len_prot)
			except KeyError:
				len_proteins_by_type[max_key] = [len_prot]
			
			scores_aux = [str(score) for score in info_protein[11]]					
			sequence_protein = info_protein[17]
			codon_tis = info_protein[9]

			max_key_aux = max_key
			if 'Canonical' in self.name:
				max_key_aux = 'Canonical'
			if 'Cryptic' in self.name:
				if 'Alt-TIS' in max_key:
					max_key_aux = 'Novel Isoform'
				elif 'Exons' in max_key:
					max_key_aux = 'Annotated Non-coding Transcript'
				if codon_tis != 'ATG' and codon_tis != 'CTG':
					sequence_protein = 'M'+ sequence_protein[1:]

			to_add = [info_protein[0], max_key_aux, ",".join(info_protein[2]), ",".join(info_protein[3]), info_protein[4], info_protein[5], info_protein[6], info_protein[7], info_protein[8], codon_tis, info_protein[10], ",".join(scores_aux), info_protein[12], info_protein[14], info_protein[15], sequence_protein, info_protein[18]]
			csvData_info_prot.append(to_add)

			if len_prot <= 100:
				small_protein += '>'+info_protein[0]+'|'+",".join(info_protein[2])+'\n'+info_protein[17]+'\n'
				small_proteins += 1
				peptides_small_proteins += n_peptides
		
			else:
				big_protein += '>'+info_protein[0]+'|'+",".join(info_protein[2])+'\n'+info_protein[17]+'\n'
				big_proteins += 1
				peptides_big_proteins += n_peptides

			try: 
				origin_by_protein[max_key] += 1
			except KeyError:
				origin_by_protein[max_key] = 1

		with open(pathToSave+'protein_origin_'+name+'.csv', 'w') as csvFile:
			writer = csv.writer(csvFile)
			writer.writerows(csvData_info_prot)

		file_to_save = open(pathToSave+'Small_proteins_source_peptide_'+name+'.fastq', 'w')
		file_to_save.write(small_protein)
		file_to_save.close()

		file_to_save = open(pathToSave+'Big_proteins_source_peptide_'+name+'.fastq', 'w')
		file_to_save.write(big_protein)
		file_to_save.close()


		csvFile.close()

		with open(pathToSave+'/info_proteins_'+name+'.dic', 'wb') as handle:
			pickle.dump(info_proteins, handle, protocol=pickle.HIGHEST_PROTOCOL)

		print 
		print '===========Info Proteins ==========='
		print name+'_Total_small_proteins =', small_proteins
		print name+'_Total_big_proteins =', big_proteins

		print name+'_Total_peptides_small_proteins =', peptides_small_proteins
		print name+'_Total_peptides_big_proteins =', peptides_big_proteins


		print
		print 'Total Proteins with one or more support peptides ', len(supported_peptides_by_protein)
		print 'Min Number of peptides ', min(supported_peptides_by_protein.values()), ' Max number of peptides ', max(supported_peptides_by_protein.values())
		print
		print 'origin_by_protein = ', origin_by_protein
		total_exons = sum(origin_by_protein.values()) - origin_by_protein['Intergenic'] - origin_by_protein['Introns']
		sizes = [total_exons, origin_by_protein['Intergenic'], origin_by_protein['Introns']] 
		print 'sizes = ',sizes
		print
		return len_proteins, len_proteins_by_type, origin_by_protein, tpm_proteins
		

def main(argv):
	
	time0 = time.time()

	saveFile = ''

	peptides_bed_file_backward = '' 
	peptides_bed_file_intersected_backward = '' 
	annotations_file = ''
	

	peptides_bed_file_forward = '' 
	peptides_bed_file_intersected_forward = ''
	genome = ''
	genome_index = ''

	name = ''
	freeByesSNPs = ''
	quality = ''

	try:
		opts, args = getopt.getopt(argv,"hs:b:i:a:f:d:t:n:g:x:l:q:",["saveFile=", 'peptides_bed_file_backward=', 'peptides_bed_file_intersected_backward=', 'annotations_file=' , 'peptides_bed_file_forward=', 'peptides_bed_file_intersected_forward=', 'name=', 'genome=', 'genome_index=', 'freeByesSNPs=', 'quality='])
	except getopt.GetoptError:
		print 'main.py -s <saveFile> -b <peptides_bed_file_backward> -i <peptides_bed_file_intersected_backward> -a <annotations_file> -f <peptides_bed_file_forward> -d <peptides_bed_file_intersected_forward> -n <name> -g <genome> -x <genome_index> -l <freeByesSNPs> -q <quality>'
		print argv
		sys.exit(2)

	#print opts
	for opt, arg in opts:
		if opt == '-h':
			print 'main.py -s <saveFile> -b <peptides_bed_file_backward> -i <peptides_bed_file_intersected_backward> -a <annotations_file> -f <peptides_bed_file_forward> -d <peptides_bed_file_intersected_forward> -n <name> -g <genome> -x <genome_index> -l <freeByesSNPs> -q <quality>'
			sys.exit()
		elif opt in ("-s", "--saveFile"):
			saveFile = arg
		elif opt in ("-b", "--peptides_bed_file_backward"):
			peptides_bed_file_backward = arg
		elif opt in ("-i", "--peptides_bed_file_intersected_backward"):
			peptides_bed_file_intersected_backward = arg
		elif opt in ("-a", "--annotations_file"):
			annotations_file = arg
		elif opt in ("-f", "--peptides_bed_file_forward"):
			peptides_bed_file_forward = arg
		elif opt in ("-d", "--peptides_bed_file_intersected_forward"):
			peptides_bed_file_intersected_forward = arg
		elif opt in ("-n", "--name"):
			name = arg
		elif opt in ("-g", "--genome"):
			genome = arg
		elif opt in ("-x", "--genome_index"):
			genome_index = arg
		elif opt in ("-l", "--freeByesSNPs"):
			freeByesSNPs = arg
		elif opt in ("-q", "--quality"):
			quality = arg


	now = time.strftime('%y-%m-%d_%H:%M:%S', time.localtime())

	with open(annotations_file, 'rb') as fp:
		infoTranscripts= pickle.load(fp)
	
	origin_peptide_backward = Origin_peptides(peptides_bed_file_backward, peptides_bed_file_intersected_backward, infoTranscripts, '-', genome, genome_index, name, freeByesSNPs, quality)
	origin_by_type_backward, origin_by_peptide_backward, origin_of_each_peptide_backward, final_origin_alt_tis_peptides_backward = origin_peptide_backward.get_origin_peptides()

	
	origin_peptide_forward = Origin_peptides(peptides_bed_file_forward, peptides_bed_file_intersected_forward, infoTranscripts, '+', genome, genome_index, name, freeByesSNPs, quality)
	origin_by_type_forward, origin_by_peptide_forward, origin_of_each_peptide_forward, final_origin_alt_tis_peptides_forward = origin_peptide_forward.get_origin_peptides()

	
	to_calcule_pie = {} 
	for key, value in origin_by_peptide_backward.items():
		try:
			origin_by_peptide_forward[key] += value
		except KeyError:
			origin_by_peptide_forward[key] = value
	
	for key, value in origin_by_type_backward.items():
		try:
			origin_by_type_forward[key] += value
			if key != 'Intergenic' and key != 'Introns':
				to_calcule_pie[key] = origin_by_type_forward[key]
		except KeyError:
			origin_by_type_forward[key] = value
			if key != 'Intergenic' and key != 'Introns':
				to_calcule_pie[key] = origin_by_type_forward[key]

	print 'Total Peptides ', sum(origin_by_type_forward.values())
	print 
	print 'Biotype Protein source peptides  = ', origin_by_type_forward
	sizes = [sum(to_calcule_pie.values()), origin_by_type_forward['Intergenic'], origin_by_type_forward['Introns']] 
	print 'sizes = ',sizes

	len_proteins, len_proteins_by_type, origin_by_type, tpm_proteins = origin_peptide_forward.writeCSV(origin_by_peptide_forward, saveFile, name, final_origin_alt_tis_peptides_backward, final_origin_alt_tis_peptides_forward)
	t1 = time.time()
	total = t1-time0
	print 'Total time running function ', total/60.0, ' min'

	with open(saveFile+'/origin_by_peptide_'+name+'.dic', 'wb') as handle:
		pickle.dump(origin_by_peptide_forward, handle, protocol=pickle.HIGHEST_PROTOCOL)

	merged_dic = dict(origin_of_each_peptide_forward.items() + origin_of_each_peptide_backward.items())
	print 'Total origin_of_each_peptide ', len(merged_dic)
	with open(saveFile+'/origin_of_each_peptide_'+name+'.dic', 'wb') as handle:
		pickle.dump(merged_dic, handle, protocol=pickle.HIGHEST_PROTOCOL)

	return origin_by_type, origin_by_peptide_forward, merged_dic, len_proteins, len_proteins_by_type, tpm_proteins

if __name__ == "__main__":
	main(sys.argv[1:])

