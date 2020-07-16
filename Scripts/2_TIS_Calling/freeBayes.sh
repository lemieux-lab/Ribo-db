#!/bin/bash

listFile=$1						# List Bam Files in which search the SNPs
genomePathFai=$2				# Fai file for the genome
genomePath=$3					# Path to the reference genome
outputFile=$4					# Output Files
nameTask=$5						# To identify the task 
saveOutputQsub=$6				# To save output Qsub
logPath=$7						# Where to save the log of this Pipelines step

module add torque
module load freebayes/1.2.0

echo "freebayes-parallel ../../../Ribo_db/Data_Input_Scripts/ref.fa.1mbp.regions 32 -f $genomePath -L $listFile > $outputFile/output.var.5X.vcf" | qsub  -V -l nodes=1:ppn=32,mem=50gb,walltime=150:00:00 -j oe -N $nameTask -d $saveOutputQsub
echo 'Launched! '

nameFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'.log'
echo $fileToLog

exec 3<> $fileToLog
    echo "$nameFileLog" >&3
    echo "" >&3
    echo '#################### SNPs Calling : FreeBayes  ######################' >&3
    echo ' > Date ', $nameFileLog >&3
    echo ' > Path to Input : '$listFile >&3
    echo ' > Path to Save : '$outputFile >&3
    echo "" >&3
    echo "RUN FreeBayes with this options : freebayes-parallel ../../../Ribo_db/Data_Input_Scripts/ref.fa.1mbp.regions 32 -f $genomePath -L $listFile > $outputFile/output.var.5X.vcf" >&3

# Close fd 3
exec 3>&-

