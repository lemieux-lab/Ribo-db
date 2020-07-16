#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: UMI extract pattern
'''

entry=$1						# Fasta File Entry
outPut=$2						# Path Output
whitelist=$3					# Path to WhiteList
nameTask=$4						# To identify the task 
saveOutputQsub=$5				# Where to save the log of this Pipelines step
logPath=$6						# Where to save the log of this Pipelines step

module add torque

echo '#################### UMIS- Extraction ######################'

echo ' > Path to Input  : '$entry
echo ' > Path to WhiteList : '$whitelist
echo ' > Path to Output : '$outPut

bcPattern='"(?P<umi_1>.{2}).{20,40}(?P<umi_2>.{5})(?P<cell_1>.{5})"'

echo "umi_tools extract -I $entry -S $outPut --extract-method=regex --whitelist=$whitelist --filter-cell-barcode --error-correct-cell --bc-pattern=$bcPattern"  | qsub -V -l nodes=1:ppn=32,mem=40gb,walltime=48:00:00 -j oe -N $nameTask -d $saveOutputQsub 

 
nameFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'.log'
echo $fileToLog

exec 3<> $fileToLog
    echo "$nameFileLog" >&3
    echo "" >&3
    echo '#################### UMIS- Extraction ######################' >&3
    echo ' > Path to Input  : '$entry >&3
	echo ' > Path to WhiteList : '$whitelist >&3
	echo ' > Path to Output : '$outPut 

    echo "" >&3
    echo "RUN umi_tools extract -I $entry -S $outPut --extract-method=regex --whitelist=$whitelist --filter-cell-barcode --error-correct-cell --bc-pattern=$bcPattern" >&3

# Close fd 3
exec 3>&-
