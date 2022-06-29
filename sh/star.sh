#!/biin/bash

path=.
mkdir star_out
### STAR alignment for RNA-seq data
ulimit -n 1000
for i in wt1 wt2 wt3 ko1 ko2 ko3;
do
gzip -d $path/trimmed/${i}.R1.paired.fastq.gz
gzip -d $path/trimmed/${i}.R2.paired.fastq.gz
~/Desktop/software/STAR-2.7.9a/bin/MacOSX_x86_64/STAR --runThreadN 8 \
--readFilesIn $path/trimmed/${i}.R1.paired.fastq $path/trimmed/${i}.R2.paired.fastq   \
--genomeDir ~/Desktop/ref/mm10/star_index \
--outFileNamePrefix $path/star_out/${i}. \
--outSAMtype BAM SortedByCoordinate ;
done 
