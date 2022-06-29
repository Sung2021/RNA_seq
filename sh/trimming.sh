#!/bin/bash

path=.
mkdir trimmed
for input in ko1 ko2 ko3 wt1 wt2 wt3 ; do java -jar ~/Desktop/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 \
$path/${input}.R1.fastq.gz $path/${input}.R2.fastq.gz \
$path/trimmed/${input}.R1.paired.fastq.gz $path/trimmed/${input}.R1.unpaired.fastq.gz \
$path/trimmed/${input}.R2.paired.fastq.gz $path/trimmed/${input}.R2.unpaired.fastq.gz CROP:110; done
