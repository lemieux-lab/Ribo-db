#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: Deduplication BAM File
'''

entry=$1						# Bam File Input
outPutBamFile=$2				# Bam File Output
outPutFolder=$3					# Bam File Folder
nameTask=$4						# To identify the task 
saveOutputQsub=$5				# Where to save the log of this Pipelines step
logPath=$6						# Where to save the log of this Pipelines step

module add torque

echo '#################### Deduplication ######################'

echo ' > Path to Input  : '$entry
echo ' > Path to Output Bam File : '$outPutBamFile

echo "umi_tools dedup -I $entry --multimapping-detection-method=NH --output-stats=$outPutFolder.DD.stats -S $outPutBamFile"  | qsub -V -l nodes=1:ppn=32,mem=40gb,walltime=48:00:00 -j oe -N $nameTask -d $saveOutputQsub 


nameFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'.log'
echo $fileToLog

exec 3<> $fileToLog
    echo "$nameFileLog" >&3
    echo "" >&3
    echo '#################### Deduplication ######################' >&3

	echo ' > Path to Input  : '$entry >&3
	echo ' > Path to Output Bam File : '$outPutBamFile >&3

    echo "" >&3
    echo "RUN umi_tools dedup -I $entry --multimapping-detection-method=NH --output-stats=$outPutFolder.DD.stats -S $outPutBamFile" >&3

# Close fd 3
exec 3>&-
