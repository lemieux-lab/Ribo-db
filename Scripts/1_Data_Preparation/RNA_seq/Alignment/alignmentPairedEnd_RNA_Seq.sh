#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: Alignment paired end reads of RNA-seq
'''

genomeDirectory=$1  				# Path to the index of the Genome
outputFilterMatchInt=$2				# Number of minimum matched bases 

inputFilesR1_1=$3					# Paths to files that contain input reads R1_1
inputFilesR1_2=$4					# Paths to files that contain input reads R1_2

inputFilesR2_1=$5                   # Paths to files that contain input reads R2_1
inputFilesR2_2=$6                   # Paths to files that contain input reads R2_2


outputFile=7					   # Path where output will be store
nameTask=8						   # To identify the task 
saveOutputQsub=9                   # To save output Qsub
logPath=${10}					   # Where to save the log of this Pipelines step

module add torque
module add star

echo '#################### Alignment to Genome ######################'

echo ' > Path to Input R1_1 : '$inputFilesR1_1
echo ' > Path to Input R1_2 : '$inputFilesR1_2

echo ' > Path to Input R2_1 : '$inputFilesR2_1
echo ' > Path to Input R2_2 : '$inputFilesR2_2

echo ' > Path to Save : '$outputFile

echo "STAR --runThreadN 10 --genomeDir $genomeDirectory --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMatchNmin $outputFilterMatchInt --outFilterMismatchNoverLmax 0.05 --outReadsUnmapped Fastx --readFilesIn $inputFilesR1_1,$inputFilesR1_2 $inputFilesR2_1,$inputFilesR2_2 --outSAMattributes All --outFileNamePrefix $outputFile" | qsub -V -l nodes=1:ppn=32,mem=40gb,walltime=48:00:00 -j oe -N $nameTask -d $saveOutputQsub
echo 'Alignment for '$outputFile' Done!'


nameFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'_Alignment_Genome.log'
echo $fileToLog
exec 3<> $fileToLog
    echo "$nameFileLog" >&3
    echo "" >&3
    echo '#################### Alignment to Genome ######################' >&3
    echo ' > Path to Input R1_1 : '$inputFilesR1_1 >&3
    echo ' > Path to Input R1_2 : '$inputFilesR1_2 >&3
    
    echo ' > Path to Input R2_1 : '$inputFilesR2_1 >&3
    echo ' > Path to Input R2_2 : '$inputFilesR2_2 >&3
    
    echo ' > Path to Save : '$outputFile >&3
    echo "" >&3
    echo "RUN STAR --runThreadN 10 --genomeDir $genomeDirectory --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMatchNmin $outputFilterMatchInt --outFilterMismatchNoverLmax 0.05 --outReadsUnmapped Fastx --readFilesIn $inputFilesR1_1,$inputFilesR1_2 $inputFilesR2_1,$inputFilesR2_2 --outSAMattributes All --outFileNamePrefix $outputFile" >&3

# Close fd 3
exec 3>&-

# STAR Options

# genomeDir								: Directory where is found the Index for the genome of Interest
# outFilterMismatchNoverLmax		 	: this is the max number of mismatches per read relative to read length:  for 1 read x 23b, max number of mismatches is 0.05*23 = 1 mismatch for a read	
# outSAMtype BAM SortedByCoordinate		: output sorted by coordinate Aligned.sortedByCoord.out.bam file, which means sorted by chromosomic order
# quantMode TranscriptomeSAM			: option STAR will outputs alignments translated into transcript coordinates in the Aligned.toTranscriptome.out.bam 
#										: file (in addition to alignments in genomic coordinates in Aligned.*.sam/bam files). 
# readFilesIn							: paths to files that contain input read1 (and, if needed, read2)
# outReadsUnmapped Fastx				: output of unmapped reads (besides SAM). output in separate fasta/fastq files, Unmapped.out.mate1/2
# outFilterMatchNmin					: int:  alignment will be output only if the number of matched bases is higher than this value
# outFileNamePrefix 					: Path to save the outputs


