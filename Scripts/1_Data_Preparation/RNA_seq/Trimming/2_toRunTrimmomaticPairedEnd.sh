#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: Call the Trimmomatic paired end script
'''


inputFileR1=$1  				# Path to the R1 reads
inputFileR2=$2  				# Path to the R2 reads
outputFile=$3					# Path to sqve Outputs
nameTask=$4						# To identify the task 
saveOutputQsub=$5           	# To save output Qsub
logPath=$6						# Where to save the log of this Pipelines step

module add torque

echo '#################### Trimmomatic Paired END ######################'

echo ' > Path to Input R1 : '$inputFileR1
echo ' > Path to Input R2 : '$inputFileR2
echo ' > Path to Save : '$outputFile

echo "sh ../../../Scripts/1_Data_Preparation/RNA_seq/Trimming/3_trimmer_TruSeqAdapters_PairedEnd.sh $inputFileR1 $inputFileR2 $outputFile" | qsub -V -l nodes=1:ppn=10,mem=50gb,walltime=48:00:00 -j oe -N $nameTask -d $saveOutputQsub
echo 'Trimming : output '$outputFile' Done!'


nameFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'_TrimmingReads.log'
echo $fileToLog
exec 3<> $fileToLog
    echo "$nameFileLog" >&3
    echo "" >&3
	echo '#################### Trimmomatic Paired END ######################' >&3
    echo ' > Path to Input R1 : '$inputFileR1 >&3
	echo ' > Path to Input R2 : '$inputFileR2 >&3
    echo ' > Path to Save : '$outputFile >&3
    echo "" >&3
    echo "RUN sh ../../../Scripts/1_Data_Preparation/RNA_seq/Trimming/3_trimmer_TruSeqAdapters_PairedEnd.sh $inputFileR1 $inputFileR2 $outputFile" >&3

# Close fd 3
exec 3>&-

