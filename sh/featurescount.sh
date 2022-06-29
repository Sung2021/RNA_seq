#!/bin/bash

path=.
mkdir $path/feature
for i in wt1 wt2 wt3 ko1 ko2 ko3;
do featureCounts -a ~/Desktop/ref/mm10_refgenes.gtf \
-o $path/feature/${i}_featurecounts.txt \
$path/star_out/${i}.Aligned.sortedByCoord.out.bam 2> $path/feature/${i}.log 
cut -f 1,7,8,9,10,11,12 $path/feature/${i}_featurecounts.txt > $path/feature/${i}.feat.mat.txt ; done
