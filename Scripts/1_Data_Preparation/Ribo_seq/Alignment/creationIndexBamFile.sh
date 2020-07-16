#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: Generation Index BAM file .bai
'''

inputFile=$1						# Paths to files that contain the BAM File
outputFile=$2						# Path where output index will be store
nameTask=$3							# To identify the task 
saveOutputQsub=$4                   # To save output Qsub
logPath=$5							# Where to save the log of this Pipelines step

module add torque
module add samtools

echo '#################### Index for BAM file ######################'

echo ' > Path to Input : '$inputFile
echo ' > Path to Save : '$outputFile

echo "samtools index $inputFile $outputFile" | qsub -V -l nodes=1:ppn=32,mem=40gb,walltime=48:00:00 -j oe -N $nameTask -d $saveOutputQsub
echo 'Index for '$inputFile' Done!'


nameFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'.log'

exec 3<> $fileToLog
    echo "$nameFileLog" >&3
    echo "" >&3
    echo '#################### Index_BAM_File ######################' >&3
    echo ' > Path to Input : '$inputFile >&3
    echo ' > Path to Save : '$outputFile >&3
    echo "" >&3
    echo "RUN Samtools Index with this options : samtools index $inputFile $outputFile" >&3

# Close fd 3
exec 3>&-
