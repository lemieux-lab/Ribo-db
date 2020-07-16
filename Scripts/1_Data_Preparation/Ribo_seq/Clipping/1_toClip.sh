#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: fastx_clipper tool applied to Ribo-seq reads
'''

input=$1  							# Path to Reads to clip
output=$2							# Path to save clipped reads
lN=$3								# Discard sequences shorter than N nucleotides
seq=$4								# Sequence to clip
nameTask=$5							# To identify the task 
saveOutputQsub=$6                   # To save output Qsub
logPath=$7							# Where to save the log of this Pipelines step

module add torque
module add fastx-toolkit/0.0.14

echo 'Run Clipper for '$input', Save Output: '$output', Minimun length: '$lN', Sequence to Clip: '$seq', Task Name Qsub: '$nameTask ', Save Qsub output: '$saveOutputQsub

[ -e $output ] && rm $output

echo "fastx_clipper -Q33 -a $seq -l $lN -v -i $input -o $output" | qsub -V -l nodes=1:ppn=32,mem=50gb,walltime=48:00:00 -j oe -N $nameTask -d $saveOutputQsub
echo 'Launched! '
echo


nameFileLog=$( date -u | tr " " .)
fileToLog=$logPath'/'$nameTask'.log'

exec 3<> $fileToLog
    echo "$nameFileLog" >&3
    echo "" >&3
    echo '#################### Clipper reads ######################' >&3
    echo ' > Path to Input : '$input >&3
    echo ' > Path to Save : '$output >&3
    echo ' > Name Task Qsub : '$nameTask >&3
    echo ' > Path to save file output Qsub : '$saveOutputQsub >&3
    echo "" >&3
    echo "RUN fastx_clipper -Q33 -a $seq -l $lN -v -i $input -o $output " >&3

# Close fd 3
exec 3>&-

# Options fastx_clipper

# Q33			: Keep only reads quality 33
# -a $seq 		: Sequence to be clipped
# -l $lN 		: Minimum lenght that a sequence can have in order to be retained. If a lenght sequence is lower than a threshold the sequence is discarded
# -i $input 	: Input file. Path to the sequences file
# -o $output 	: Output file. Path in which the clipped sequences will be stored

# See more documentation for qsub at:
# https://www.nas.nasa.gov/hecc/support/kb/commonly-used-qsub-options-in-pbs-scripts-or-in-the-qsub-command-line_175.html

# See more documentqtion for fastx_clipper at :
# http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_clipper_usage
