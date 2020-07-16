#!/bin/bash

inputFiles=$1						# Path to files that contain the BAM File
outputFile=$2						# Path where output assemblies will be store
library=$3							# Type of library between fr et rf
gtfGuide=$4							# Guidance to the aseemblie transcript, annotated GTF
numberMinReads=$5					# Minimun number of reads to the assemblie
minimunLenght=$6					# Sets the minimum length allowed for the predicted transcripts. Default: 200
minimunGap=$7						# Minimum locus gap separation value. Reads that are mapped closer than this distance are merged together in the same processing bundle
nameTask=$8							# To identify the task 
saveOutputQsub=$9                   # To save output Qsub
logPath=${10}						# Where to save the log of this Pipelines step
gene_abund=${11}					# Path to save gene abudance
cov_refs=${12}						# Path to save coverage

module add torque
module load stringtie/1.3.6
echo 
echo '#################### Transcript Assemblies  ######################'

echo ' > Path to Input : '$inputFiles
echo ' > Path to Save : '$outputFile

if [ $library -eq 1 ]
then # Ribosome Profiling
	echo "stringtie $inputFiles -G $gtfGuide -o $outputFile -p 32 --fr -m $minimunLenght -g $minimunGap -c $numberMinReads -C $cov_refs -A $gene_abund" | qsub -V -l nodes=1:ppn=32,mem=100gb,walltime=48:00:00 -j oe -N $nameTask -d $saveOutputQsub
else # RNA
	echo "stringtie $inputFiles -G $gtfGuide -o $outputFile -p 32 --rf -m $minimunLenght -g $minimunGap -c $numberMinReads -C $cov_refs -A $gene_abund" | qsub -V -l nodes=1:ppn=32,mem=100gb,walltime=48:00:00 -j oe -N $nameTask -d $saveOutputQsub
fi
echo 'StringTie for '$inputFiles' Done!'


dateFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'.log'
echo $fileToLog

exec 3<> $fileToLog
    echo "$dateFileLog" >&3
    echo "" >&3
    echo '#################### Transcript Assemblies ######################' >&3
    echo ' > Path to Input : '$inputFiles >&3
    echo ' > Path to Save : '$outputFile >&3
    echo "" >&3
    if [ $library -eq 1 ]
	then # Ribosome Profiling
		echo "RUN StringTie with this options : stringtie $inputFiles -G $gtfGuide -o $outputFile -p 32 --fr -m $minimunLenght -g $minimunGap -c $numberMinReads -C $cov_refs -A $gene_abund"  >&3
	else # RNA
		echo "RUN StringTie with this options : stringtie $inputFiles -G $gtfGuide -o $outputFile -p 32 --rf -m $minimunLenght -g $minimunGap -c $numberMinReads -C $cov_refs -A $gene_abund"  >&3
	fi

exec 3>&-

