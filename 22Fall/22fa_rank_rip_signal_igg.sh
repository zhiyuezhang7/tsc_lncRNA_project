#!/bin/bash

#SBATCH -J igg_%j
#SBATCH -p general
#SBATCH --mem=50G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH -o igg_%j.out
#SBATCH -e igg_%j.err
#SBATCH --mail-type ALL 
#SBATCH --mail-user zhiyue@live.unc.edu

#### Goal: do everything required of igg

## Example command:
    # sbatch 22fa_rank_rip_signal_igg.sh

## load modules
module load star/2.5.4b
module load samtools

# input files:
    # /proj/calabrlb/users/Zhiyue/22_sp/rip/*igg* .

### 1. alignment: generate SAM, BAM, WIG
    
## put rip files and do following work in directory named igg
mkdir igg
cd igg
cp /proj/calabrlb/users/Zhiyue/22_sp/rip/*igg* .
gunzip *.fastq.gz
    # concatenate all fastq file names, separated by comma
        # remove spaces and newlines
rip_filenames=$(ls -m *.fastq | sed -r 's/\s+//g' | tr -d '\n')
echo $rip_filenames

# calculate read length
fastq_lin_num_igg=$(wc -l *igg*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
(( reads_igg = fastq_lin_num_igg / 4 ))
echo "Reads of igg is:"
echo $reads_igg

## use STAR aligner to align reads to mm10 mouse genome.
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix igg_ --readFilesIn $rip_filenames

## Filter and output to sam file
samtools view -h -Sq 30 igg_Aligned.out.sam > igg_Aligned_filteredsq30.out.sam

## sam to wig
# split sam by strand
samtools view -h -F 0x10 igg_Aligned_filteredsq30.out.sam > igg_Aligned_neg.out.sam
samtools view -h -f 0x10 igg_Aligned_filteredsq30.out.sam > igg_Aligned_pos.out.sam

# Generate wiggle tracks
    # Positive strand is blue; negative strand is red
    # All colors usable: pink, babyblue, grey, black, red, green, blue, yellow, lightgreen, purple, navyblue maroon, brown, orange, or magenta
    # move bigsam_to_wig_mm10.pl to directory
cp /proj/calabrlb/users/Zhiyue/22_sp/code/bigsam_to_wig_mm10.pl .
perl bigsam_to_wig_mm10.pl igg_Aligned_neg.out.sam igg_neg red
perl bigsam_to_wig_mm10.pl igg_Aligned_pos.out.sam igg_pos blue

## convert SAM to BAM
samtools view -S -b igg_Aligned_neg.out.sam > igg_neg.bam
samtools view -S -b igg_Aligned_pos.out.sam > igg_pos.bam

# index and sort BAM
samtools sort igg_neg.bam -o igg_sorted_neg.bam
samtools sort igg_pos.bam -o igg_sorted_pos.bam
samtools index igg_sorted_neg.bam
samtools index igg_sorted_pos.bam

# # save sam, wig, bam files
# cp igg_Aligned_filteredsq30.out.sam /proj/calabrlb/users/Zhiyue/22_sp/rip_sam
# cp igg_neg.wig /proj/calabrlb/users/Zhiyue/22_sp/rip_wig
# cp igg_pos.wig /proj/calabrlb/users/Zhiyue/22_sp/rip_wig

# cp igg_sorted_neg.bam /proj/calabrlb/users/Zhiyue/22_sp/rip_bam
# cp igg_sorted_pos.bam /proj/calabrlb/users/Zhiyue/22_sp/rip_bam
# cp igg_sorted_neg.bam.bai /proj/calabrlb/users/Zhiyue/22_sp/rip_bam
# cp igg_sorted_pos.bam.bai /proj/calabrlb/users/Zhiyue/22_sp/rip_bam