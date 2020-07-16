#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: Trimmomatic script call
'''


r1=$1 
r2=$2 
outputFile=$3

echo 'File Input R1 ' $1 					#gets the path where the input data can be found for R1
echo 'File Input R2 ' $2 					#gets the path where the input data can be found for R2
echo 'File OutPut Trimmed ' $3				#gets the path where the output data will be saved

java -Xms4g -Xmx4g -jar .../trimmomatic-0.35.jar PE \
             -threads 10 -phred33 \
             $r1 \
             $r2 \
             $outputFile/R1_trimmed.fastq \
             $outputFile/R1_trimmed.unpaired.fastq \
             $outputFile/R2_trimmed.fastq \
             $outputFile/R2_trimmed.unpaired.fastq \
             ILLUMINACLIP:.../adapters/TruSeq_and_nextera_adapters.fa:2:30:10:8:true \
             LEADING:20 TRAILING:20 MINLEN:20


# Trimmomatic Options

#-phred33 : Specify that reads having a quality score below 33 must be discarded

#$r1 \  : Specify input files
#$r2 \ 	: Specify input files

#$outputFile/R1_trimmed.fastq \				: Specify output files 
#$outputFile/R1_trimmed.unpaired.fastq \	: Specify output files 
#$outputFile/R2_trimmed.fastq \				: Specify output files 
#$outputFile/R2_trimmed.unpaired.fastq \	: Specify output files 


#ILLUMINACLIP:.../adapters/TruSeq3-PE-2.fa:2:30:10:8:true \ : Specify the file where the adpaters are listed.
#ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
#    fastaWithAdaptersEtc: specifies the path to a fasta file containing all the adapters, PCR sequences etc. The naming of the various sequences within this file determines how they are used. See below.
#    seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
#    palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
#    simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.


#LEADING:20 	: LEADING: Cut bases off the start of a read, if below a threshold quality
#TRAILING:20 	: TRAILING: Cut bases off the end of a read, if below a threshold quality
#MINLEN:20		: MINLEN: Drop the read if it is below a specified length
