#!/bin/bash

path=~/Desktop/Sung_work/fastq/19092_50_01
mkdir fastqc_result
cd $path/

for input in $(ls *.fastq.gz) ;do fastqc $path/${input} ; done


### trimming 
## gzip files only
cd $path/
mkdir trimmed
for input in ko1 ko2 ko3 wt1 wt2 wt3 ; do java -jar ~/Desktop/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 $path/${input}.R1.fastq.gz $path/${input}.R2.fastq.gz $path/trimmed/${input}.R1.paired.fastq.gz $path/trimmed/${input}.R1.unpaired.fastq.gz $path/trimmed/${input}.R2.paired.fastq.gz $path/trimmed/${input}.R2.unpaired.fastq.gz CROP:110; done


### STAR alignment for RNA-seq data
path=~/Desktop/Sung_work/fastq/19092_50_01

cd $path/star_out
ulimit -n 1000
for i in wt1 wt2 wt3 ko1 ko2 ko3;
do
 sampleid="${i}"
mkdir $path/${sampleid}
gzip -d $path/trimmed/${sampleid}.R1.paired.fastq.gz
gzip -d $path/trimmed/${sampleid}.R2.paired.fastq.gz
~/Desktop/software/STAR-2.7.9a/bin/MacOSX_x86_64/STAR --runThreadN 8 \
--readFilesIn $path/trimmed/${sampleid}.R1.paired.fastq $path/trimmed/${sampleid}.R2.paired.fastq   \
--genomeDir ~/Desktop/ref/mm10/star_index \
--outFileNamePrefix $path/star_out/${sampleid} \
--outSAMtype BAM SortedByCoordinate ;
done 


cd $path/star_out

### sorting star output
mkdir $path/sorted
cd $path/star_out
for i in wt1 wt2 wt3 ko1 ko2 ko3;
do
 sampleid="${i}"
 samtools sort -@ 10 $path/star_out/${sampleid} -o $path/sorted/${sampleid}.sorted.bam
 samtools index -@ 10 $path/sorted/${sampleid}.sorted.bam 
done

### featureCounts
mkdir $path/features
cd $path/sorted
for i in wt1 wt2 wt3 ko1 ko2 ko3;
do
 sampleid="${i}"
featureCounts -a ~/Desktop/ref/mm10_refgenes.gtf  -o $path/features/${sampleid}_featurecounts.txt $path/sorted/${sampleid}.sorted.bam 2> $path/features/${sampleid}.log
done


### cleaning up the featurecounts matrix
cd $path/features
for i in wt1 wt2 wt3 ko1 ko2 ko3;
do
 sampleid="${i}"
cut -f 1,7,8,9,10,11,12 $path/features/${sampleid}_featurecounts.txt > $path/features/${sampleid}_featurecounts.matrix.txt 
done
