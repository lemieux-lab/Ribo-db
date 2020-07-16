#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: To call bash script that will clip the reads according to data found in the manifest.ini file
'''


module add torque
manifestFile=$1                     # Input manifest file, where paths to every file of replicate is defiend along 
                                    # with the sequence to clip and the path to save the clipped reads
logPath=$2  

# Common variables
outputQsub=''                       # To save output Qsub
nameTask=''							# To identify the task 

##############################################################################
#                           Clipping DATA                                    #
##############################################################################

#Init variables for clipping
input=''  							# Path to Reads to clip
output=''							# Path to save clipped reads
lN=''								# Discard sequences shorter than N nucleotides
seq=''								# Sequence to clip

lastCommand=''
while IFS='' read -r line || [[ -n "$line" ]]; do
    if [[ $line == *"outputQsub"* ]]; then
        IFS='=' read -ra ADDR <<< "$line"
        outputQsub=${ADDR[1]}
        echo 'To Save Ouput QSub '$outputQsub
        echo
    fi
    if [[ $line == *"Command"* ]]; then
        IFS='::' read -ra ADDR <<< "$line"
        lastCommand=${ADDR[2]}
        echo '>> Command :: '$lastCommand
    fi
    ## Command to Clip
    if [[ $lastCommand == "InitToClip" ]]; then
        IFS='=' read -ra ADDR <<< "$line"
        variable=${ADDR[0]}
        value=${ADDR[1]}
        if [[ $variable == "input" ]]; then
            input=$value
        elif [[ $variable == "output" ]]; then
            output=$value
        elif [[ $variable == "lN" ]]; then
            lN=$value
        elif [[ $variable == "seq" ]]; then
            seq=$value
        elif [[ $variable == "nameTask" ]]; then
            nameTask=$value
        fi
        if [[ $line == "EndToClip" ]]; then
            if [[ $input != '' ]] && [[ $output != '' ]] && [[ $lN != '' ]] && [[ $seq != '' ]] && [[ $nameTask != '' ]]; then
                sh ../../../Scripts/1_Data_Preparation/Ribo_seq/Clipping/1_toClip.sh $input $output $lN $seq $nameTask $outputQsub $logPath
                input=''
                output=''
                lN=''
                seq=''
                nameTask=''
            else
                echo 'Any of demanding parameters is not fullfilled'
                echo $outputQsub $input $output $lN $seq $nameTask
            fi
        fi
    fi
done < $manifestFile
