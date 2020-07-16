#!/bin/bash

'''
# Author: Maria Virginia Ruiz
# Objective: Filtering Ribo-seq reads according to their length 26-34 ntds
'''

input=$1  							# Path to Reads to filter
output=$2							# Path to save Filtered reads
logPath=$3							# Where to save the log of this Pipelines step

echo 'Filter reads for '$input', Save Output: '$output

awk '!last { last = $0; next } length($0)<=34 && length($0)>=26 { print last; print } { last = "" }' $input > $output


# Explanation Options
# This script filter out the reads that are larger than 34 ntds and smaller than 26.