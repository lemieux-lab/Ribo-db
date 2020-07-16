#!/bin/bash

input=$1						# Bed file Start Codons Candidates
input2=$2						# GTF transcripts Assemblies after merge
output=$3						# Path to the output intersec bed file
nameTask=$4						# To identify the task 
saveOutputQsub=$5				# To save output Qsub
logPath=$6						# Where to save the log of this Pipelines step


module add bedtools
module add torque

echo $saveOutputQsub
echo $output

echo "bedtools intersect -a $input -b $input2  -wao -f 1 | grep -w transcript > $output" | qsub -V -l nodes=1:ppn=8,mem=12gb,walltime=5:00:00 -j oe -N $nameTask -d $saveOutputQsub


nameFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'.log'
echo $fileToLog

exec 3<> $fileToLog
    echo "$nameFileLog" >&3
    echo "" >&3
    echo '#################### Bedtools Intersection SC-Transcripts  ######################' >&3
    echo ' > Date : '$nameFileLog >&3
    echo ' > Path to Input A : '$input >&3
	echo ' > Path to Input B : '$input2 >&3
	echo ' > Path to Save A intersect B : '$output >&3
    echo "" >&3
    echo "RUN Bedtools with this options : bedtools intersect -a $input -b $input2  -wao -f 1 | grep -w exon > $output" >&3

exec 3>&-

