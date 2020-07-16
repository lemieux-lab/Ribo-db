#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: To call script that filter Ribo-seq reads according to their length 26-34 ntds
'''

input=$1  							# Path to Reads to clip
output=$2							# Path to save clipped reads
nameTask=$3							# To identify the task 
saveOutputQsub=$4                   # To save output Qsub
logPath=$5							# Where to save the log of this Pipelines step

module add torque

echo 'Run Filtering for '$input', Save Output: '$output', Task Name Qsub: '$nameTask ', Save Qsub output: '$saveOutputQsub

#[ -e $output ] && rm $output

echo "sh /u/ruizma/PIPELINE/Scripts/DataPreparation/BashScripts/2_Filter_by_Length/1_toFilter.sh $input $output" | qsub -V -l nodes=1:ppn=32,mem=40gb,walltime=48:00:00 -j oe -N $nameTask -d $saveOutputQsub
echo 'Launched! '
echo


# This save a log file in the path logPath, with the name fo the date + Alignment_Genome.log
nameFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'.log'

exec 3<> $fileToLog
    echo "$nameFileLog" >&3
    echo "" >&3
    echo '#################### Alignment to Genome ######################' >&3
    echo ' > Path to Input : '$input >&3
    echo ' > Path to Save : '$output >&3
    echo ' > Name Task Qsub : '$nameTask >&3
    echo ' > Path to save file output Qsub : '$saveOutputQsub >&3
    echo "" >&3
    echo "RUN sh /u/ruizma/PIPELINE/Scripts/DataPreparation/BashScripts/2_Filter_by_Length/1_toFilter.sh $input $output " >&3

# Close fd 3
exec 3>&-

#sh /u/ruizma/PIPELINE/Scripts/DataPreparation/BashScripts/2_Filter_by_Length/0_toRunFilter.sh 