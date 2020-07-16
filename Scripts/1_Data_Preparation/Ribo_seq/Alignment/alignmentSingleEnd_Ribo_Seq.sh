#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: Alignment Ribo-seq reads to the genome using STAR aligner
'''

genomeDirectory=$1  				# Path to the index of the Genome
outputFilterMatchInt=$2				# Number of minimum matched bases 
inputFilesR1_1=$3					# Paths to files that contain input reads
inputFilesR1_2=$4					# Paths to files that contain input reads
outputFile=$5						# Path where output will be store
anchorMM=$6							# Anchor MM
nameTask=$7							# To identify the task 
saveOutputQsub=$8                   # To save output Qsub
logPath=$9							# Where to save the log of this Pipelines step

module add torque
module add star

echo '#################### Alignment to Genome ######################'

echo ' > Path to Input R1_1: '$inputFilesR1_1
echo ' > Path to Input R1_2 : '$inputFilesR1_2
echo ' > Path to Save : '$outputFile

echo "STAR --runThreadN 10 --genomeDir $genomeDirectory --outSAMtype BAM SortedByCoordinate --alignEndsType EndToEnd --quantMode TranscriptomeSAM GeneCounts --seedSearchStartLmax 15 --winAnchorMultimapNmax $anchorMM --outFilterMatchNmin $outputFilterMatchInt --outFilterMismatchNoverLmax 0.05 --outReadsUnmapped Fastx --readFilesIn $inputFilesR1_1,$inputFilesR1_2 --outSAMattributes NH HI MD --outFileNamePrefix $outputFile" | qsub -V -l nodes=1:ppn=32,mem=40gb,walltime=48:00:00 -j oe -N $nameTask -d $saveOutputQsub
echo 'Alignment for '$outputFile' Done!'


nameFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'_Alignment_Genome.log'

exec 3<> $fileToLog
    echo "$nameFileLog" >&3
    echo "" >&3
    echo '#################### Alignment to Genome ######################' >&3
    echo ' > Path to Input R1_1: '$inputFilesR1_1 >&3
	echo ' > Path to Input R1_2 : '$inputFilesR1_2 >&3
    echo ' > Path to Save : '$outputFile >&3
    echo "" >&3
    echo "RUN STAR --runThreadN 10 --genomeDir $genomeDirectory --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --seedSearchStartLmax 15 --winAnchorMultimapNmax $anchorMM --outFilterMatchNmin $outputFilterMatchInt --outFilterMismatchNoverLmax 0.05 --outReadsUnmapped Fastx --readFilesIn $inputFilesR1_1,$inputFilesR1_2 --outSAMattributes NH HI MD --outFileNamePrefix $outputFile" >&3

# Close fd 3
exec 3>&-

# STAR Options

# genomeDir								: Directory where is found the Index for the genome of Interest
# outFilterMismatchNoverLmax		 	: this is the max number of mismatches per read relative to read length:  for 1 read x 23b, max number of mismatches is 0.05*23 = 1 mismatch for a read	#This is what we typically use for ENCODE small RNA-seq. This means that for reads <20b no mismatches are allowed, 20-39b: 1 mismatch, 40-59b 2 mismatches and so on..
# outSAMtype BAM SortedByCoordinate		: output sorted by coordinate Aligned.sortedByCoord.out.bam file, which means sorted by chromosomic order
# quantMode TranscriptomeSAM			: option STAR will outputs alignments translated into transcript coordinates in the Aligned.toTranscriptome.out.bam 
#										: file (in addition to alignments in genomic coordinates in Aligned.*.sam/bam files). 
# readFilesIn							: paths to files that contain input read1 (and, if needed, read2)
# outReadsUnmapped Fastx				: output of unmapped reads (besides SAM). output in separate fasta/fastq files, Unmapped.out.mate1/2
# outFilterMatchNmin					: int:  alignment will be output only if the number of matched bases is higher than this value
# outFileNamePrefix 					: Path to save the outputs
