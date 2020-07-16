#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: UMI WhiteList
'''

entry=$1						# Fasta File Entry
outPutWhiteList=$2				# White List Output
nameTask=$3						# To identify the task 
saveOutputQsub=$4				# Where to save the log of this Pipelines step
logPath=$5						# Where to save the log of this Pipelines step

module add torque

echo '#################### UMIS WhiteList ######################'

echo ' > Path to Input  : '$entry
echo ' > Path to Output WhiteList : '$outPutWhiteList

bcPattern='"(?P<umi_1>.{2}).{15,40}(?P<umi_2>.{5})(?P<cell_1>.{5})"'

cd $masterDirectory

echo "umi_tools whitelist -I $entry --extract-method=regex --bc-pattern=$bcPattern --set-cell-number=1 > $outPutWhiteList"  | qsub -V -l nodes=1:ppn=32,mem=40gb,walltime=48:00:00 -j oe -N $nameTask -d $saveOutputQsub 


nameFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'.log'
echo $fileToLog

exec 3<> $fileToLog
    echo "$nameFileLog" >&3
    echo "" >&3
    echo '#################### UMIS WhiteList ######################' >&3
    echo ' > Path to Input  : '$entry >&3
	echo ' > Path to Output WhiteList : '$outPutWhiteList >&3

    echo "" >&3
    echo "RUN umi_tools whitelist -I $entry --extract-method=regex --bc-pattern=$bcPattern --set-cell-number=1 > $outPutWhiteList" >&3


exec 3>&-
