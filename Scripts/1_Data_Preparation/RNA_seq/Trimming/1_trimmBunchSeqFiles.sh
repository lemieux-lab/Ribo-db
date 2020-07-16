#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: This scripts calls to the Trimmomatic paired end script
'''

masterDirectory=$1              # Path to find the samples (fastq) 
pathToSaveOutput=$2             # Path to save output after trimming
saveOutputQsub=$3               # Path to save qsub output
logPath=$4                      # Path to save logs

# This scripts finds the path of the paired end reads in the directory masterDirectory
cd $masterDirectory
echo $masterDirectory

for d in */ ; do
    echo
    echo "$d"
    path=$masterDirectory/$d
    cd $path
    outPut=$pathToSaveOutput/$d
    rm -rf $outPut
    mkdir $outPut
    r1=''
    r2=''
    IFS='/' read -ra ADDR <<< "$d"
    nameTask='Trimmomatic'${ADDR[0]}
    for entry in "$path"*
	do
        if [[ $entry == *"R1"* ]] && [[ $entry != *"gz"* ]]; then
            r1=$entry
        elif [[ $entry == *"R2"* ]] && [[ $entry != *"gz"* ]]; then
            r2=$entry
        fi
	done
    sh ../../../Scripts/1_Data_Preparation/RNA_seq/Trimming/2_toRunTrimmomaticPairedEnd.sh $r1 $r2 $outPut $nameTask $saveOutputQsub $logPath
done
echo
echo 'Trimmming all files : OK!'
