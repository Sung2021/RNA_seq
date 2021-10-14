#!/bin/bash

path=~/Desktop/Sung_work/fastq/19092_50_01
mkdir fastqc_result
cd $path/

for input in $(ls *.fastq.gz) ;do fastqc $path/${input} ; done
