#!/bin/bash

path=~/Desktop/Sung_work/fastq/19092_50_01_07_18
mkdir $path/fastqc_result
cd $path/

for input in $(ls *.fastq.gz) ;do fastqc $path/${input} ; done

### trimming 
crop=110
cd $path/
mkdir trimmed
for input in d5.wt1 d5.wt2 d5.wt3 d5.mt1 d5.mt2 d5.mt3 d9.wt1 d9.wt2 d9.wt3 d9.mt1 d9.mt2 d9.mt3 ; do java -jar ~/Desktop/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 $path/Tfh.${input}.R1.fastq.gz $path/Tfh.${input}.R2.fastq.gz $path/trimmed/${input}.R1.paired.fastq.gz $path/trimmed/${input}.R1.unpaired.fastq.gz $path/trimmed/${input}.R2.paired.fastq.gz $path/trimmed/${input}.R2.unpaired.fastq.gz CROP:${crop}; done


### STAR alignment for RNA-seq data
path=~/Desktop/Sung_work/fastq/19092_50_01_07_18

ulimit -n 1000
for input in d5.wt1 d5.wt2 d5.wt3 d5.mt1 d5.mt2 d5.mt3 d9.wt1 d9.wt2 d9.wt3 d9.mt1 d9.mt2 d9.mt3 ;
do 
gzip -d $path/trimmed/${input}.R1.paired.fastq.gz
gzip -d $path/trimmed/${input}.R2.paired.fastq.gz

~/Desktop/software/STAR-2.7.9a/bin/MacOSX_x86_64/STAR --runThreadN 8 \
--readFilesIn $path/trimmed/${input}.R1.paired.fastq $path/trimmed/${input}.R2.paired.fastq   \
--genomeDir ~/Desktop/ref/mm10/star_index \
--outFileNamePrefix $path/star_out/${input} \
--outSAMtype BAM SortedByCoordinate ;
done 

cd $path/star_out

### sorting star output
mkdir $path/sorted
cd $path/star_out
for input in d5.wt1 d5.wt2 d5.wt3 d5.mt1 d5.mt2 d5.mt3 d9.wt1 d9.wt2 d9.wt3 d9.mt1 d9.mt2 d9.mt3 ;
do
 samtools sort -@ 10 $path/star_out/${input}Aligned.sortedByCoord.out.bam -o $path/sorted/${input}.sorted.bam
 samtools index -@ 10 $path/sorted/${input}.sorted.bam ;
done


### featureCounts
mkdir $path/features
cd $path/sorted
for input in d5.wt1 d5.wt2 d5.wt3 d5.mt1 d5.mt2 d5.mt3 d9.wt1 d9.wt2 d9.wt3 d9.mt1 d9.mt2 d9.mt3 ;
do
featureCounts -a ~/Desktop/ref/mm10_refgenes.gtf  -o $path/features/${input}_featurecounts.txt $path/sorted/${input}.sorted.bam 2> $path/features/${input}.log ;
done


### cleaning up the featurecounts matrix
cd $path/features
for input in *txt;
do
cut -f 1,7,8,9,10,11,12 $path/features/${input} > $path/features/${input}_matrix.txt ;
done

