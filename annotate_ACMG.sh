#!/bin/bash
#$ -cwd

#source activate py36

INPUT=${1}   # input vcf file 
OUTPUT=${2}  # output file

time python ACMGv14.py -i $INPUT -o $OUTPUT
