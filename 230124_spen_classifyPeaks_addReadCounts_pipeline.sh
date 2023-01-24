#!/bin/bash

#SBATCH -J peak_%j
#SBATCH -p general
#SBATCH --mem=50G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH -o peak_%j.out
#SBATCH -e peak_%j.err
#SBATCH --mail-type ALL 
#SBATCH --mail-user zhiyue@live.unc.edu

## input files:
    # /proj/calabrlb/users/Zhiyue/23_sp/0112_spen_fastq/
        # JBT20221213_delSpen1_SPEN_RIP_S18_R1_001.fastq.gz
        # JBT20221213_delSpen2_SPEN_RIP_S20_R1_001.fastq.gz
        # JBT20221213_WT1_SPEN_RIP_S17_R1_001.fastq.gz
        # JBT20221213_WT2_SPEN_RIP_S19_R1_001.fastq.gz

## load modules
module load star/2.5.4b
module load samtools
module load macs
module load subread
module load bedtools

#### Part 1: from rip, STAR genome alignment and MACS peak-calling 

### 1A. STAR alignment
## All WT
# move rip files here
mkdir rip_fastq
cp /proj/calabrlb/users/Zhiyue/23_sp/0112_spen_fastq/JBT20221213_WT1_SPEN_RIP_S17_R1_001.fastq.gz rip_fastq
cp /proj/calabrlb/users/Zhiyue/23_sp/0112_spen_fastq/JBT20221213_WT2_SPEN_RIP_S19_R1_001.fastq.gz rip_fastq
gunzip rip_fastq/*.fastq.gz

# use STAR aligner to align reads to mm10 mouse genome.
rip_filenames=$(ls -m rip_fastq/*.fastq | sed -r 's/\s+//g' | tr -d '\n')
echo $rip_filenames
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix spen_ --readFilesIn $rip_filenames

# Filter and output to sam file
samtools view -h -Sq 30 spen_Aligned.out.sam > spen_Aligned_filteredsq30.out.sam

## WT1
wt1_filename="rip_fastq/JBT20221213_WT1_SPEN_RIP_S17_R1_001.fastq"

# use STAR aligner to align reads to mm10 mouse genome.
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix spen_wt1_ --readFilesIn $wt1_filename

# Filter and output to sam file
samtools view -h -Sq 30 spen_wt1_Aligned.out.sam > spen_wt1_Aligned_filteredsq30.out.sam

## WT2
wt2_filename="rip_fastq/JBT20221213_WT2_SPEN_RIP_S19_R1_001.fastq"

# use STAR aligner to align reads to mm10 mouse genome.
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix spen_wt2_ --readFilesIn $wt2_filename

# Filter and output to sam file
samtools view -h -Sq 30 spen_wt2_Aligned.out.sam > spen_wt2_Aligned_filteredsq30.out.sam


## control
# move control files here
mkdir ctrl_fastq
cp /proj/calabrlb/users/Zhiyue/23_sp/0112_spen_fastq/JBT20221213_delSpen1_SPEN_RIP_S18_R1_001.fastq.gz ctrl_fastq
cp /proj/calabrlb/users/Zhiyue/23_sp/0112_spen_fastq/JBT20221213_delSpen2_SPEN_RIP_S20_R1_001.fastq.gz ctrl_fastq
gunzip ctrl_fastq/*.fastq.gz

# use STAR aligner to align reads to mm10 mouse genome.
ctrl_filenames=$(ls -m ctrl_fastq/*.fastq | sed -r 's/\s+//g' | tr -d '\n')
echo $ctrl_filenames
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix ctrl_ --readFilesIn $ctrl_filenames

# Filter and output to sam file
samtools view -h -Sq 30 ctrl_Aligned.out.sam > ctrl_Aligned_filteredsq30.out.sam


### 1B. MACS peak-calling

# split sam by strand
samtools view -h -F 0x10 spen_Aligned_filteredsq30.out.sam > spen_Aligned_neg.out.sam
samtools view -h -f 0x10 spen_Aligned_filteredsq30.out.sam > spen_Aligned_pos.out.sam

# load randomization code
cp /proj/calabrlb/users/Zhiyue/22_sp/code/macs_strand_rand_sam.pl .

# randomize strand
perl macs_strand_rand_sam.pl spen_Aligned_pos.out.sam spen_Aligned_pos.out.rand.sam
perl macs_strand_rand_sam.pl spen_Aligned_neg.out.sam spen_Aligned_neg.out.rand.sam

# macs
    # merge size = 304 = 4 * 76
macs2 callpeak -t spen_Aligned_pos.out.rand.sam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir spen_pos_peaks -n spen_pos
macs2 callpeak -t spen_Aligned_neg.out.rand.sam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir spen_neg_peaks -n spen_neg 
    # add trackline for UCSC
cat spen_neg_peaks/spen_neg_peaks.broadPeak | sed "1 i\track type=broadPeak name="spen_peaks_neg" description="spen_peaks_neg" nextItemButton=on" > spen_neg_peaks_ucsc.broadPeak
cat spen_pos_peaks/spen_pos_peaks.broadPeak | sed "1 i\track type=broadPeak name="spen_peaks_pos" description="spen_peaks_pos" nextItemButton=on" > spen_pos_peaks_ucsc.broadPeak



#### Part 2: featureCounts count number of reads in each peak

## convert broadPeak to SAF file in preparation for featureCounts
cat spen_neg_peaks_ucsc.broadPeak | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"-"}' > spen_saf.txt
cat spen_pos_peaks_ucsc.broadPeak | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"+"}' >> spen_saf.txt
cat spen_saf.txt | awk -v OFS="\t" '{print NR,$0}' | sed '1 i\GeneID\tChr\tStart\tEnd\tStrand' > spen_peaks.saf

## featureCounts with SAF - count reads in each peak for All_WT, WT1, WT2, and ctrl
featureCounts -s 2 -F SAF -a spen_peaks.saf -o spen_fc spen_Aligned_filteredsq30.out.sam
featureCounts -s 2 -F SAF -a spen_peaks.saf -o spen_wt1_fc spen_wt1_Aligned_filteredsq30.out.sam
featureCounts -s 2 -F SAF -a spen_peaks.saf -o spen_wt2_fc spen_wt2_Aligned_filteredsq30.out.sam
featureCounts -s 2 -F SAF -a spen_peaks.saf -o spen_ctrl_fc ctrl_Aligned_filteredsq30.out.sam

## make summary file
# extract raw counts
cut -f7 spen_wt1_fc > spen_wt1_counts.txt
cut -f7 spen_wt2_fc > spen_wt2_counts.txt
cut -f7 spen_ctrl_fc > spen_ctrl_counts.txt

# concatenate files
paste spen_fc spen_wt1_counts.txt spen_wt2_counts.txt spen_ctrl_counts.txt | sed 1d > spen_fc.txt



#### Part 3: find true peaks enriched over ctrl (rpm > 2 * ctrl_rpm)

## calculate rpm = (count * 1000,000 / reads)
    # counts: $7 (All WT) $8 (WT1) $9 (WT2) $10 (ctrl) of spen_fc.txt

# var: reads = number of lines in fastq files / 4
fastq_lin_num=$(wc -l rip_fastq/* | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
(( reads = fastq_lin_num / 4 ))

fastq_lin_num_wt1=$(wc -l rip_fastq/JBT20221213_WT1_SPEN_RIP_S17_R1_001.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
(( reads_wt1 = fastq_lin_num_wt1 / 4 ))

fastq_lin_num_wt2=$(wc -l rip_fastq/JBT20221213_WT2_SPEN_RIP_S19_R1_001.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
(( reads_wt2 = fastq_lin_num_wt2 / 4 ))

fastq_lin_num_ctrl=$(wc -l ctrl_fastq/* | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
(( reads_ctrl = fastq_lin_num_ctrl / 4 ))

# print rpm to $11-$14 columns
cat spen_fc.txt | sed 1d | awk -v OFS="\t" -v reads=$reads -v reads_wt1=$reads_wt1 -v reads_wt2=$reads_wt2 -v reads_ctrl=$reads_ctrl '{print $0, (($7*1000000)/reads), (($8*1000000)/reads_wt1), (($9*1000000)/reads_wt2), (($10*1000000)/reads_ctrl)}' | sed "1 i\GeneID\tChr\tStart\tEnd\tStrand\tLength\tspen_counts\tspen_wt1_counts\tspen_wt2_counts\tctrl_counts\tspen_rpm\tspen_wt1_rpm\tspen_wt2_rpm\tctrl_rpm" > spen_fc_rpm.txt
    # counts: $7 (All WT) $8 (WT1) $9 (WT2) $10 (ctrl) of spen_fc_rpm.txt
    # rpm: $11 (All WT) $12 (WT1) $13 (WT2) $14 (ctrl) of spen_fc_rpm.txt

# filter peaks with rpm 2 folds over ctrl
cat spen_fc_rpm.txt | awk '{if(NR == 1) print $0; else if($11 > 2 * $14) print $0}' > spen_fc_rpm_2ctrl.txt

# convert enricheded peaks to BED format
cat spen_fc_rpm_2ctrl.txt | sed 1d | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > spen_peaks_2ctrl.bed

# add trackline for UCSC
cat spen_peaks_2ctrl.bed | sed "1 i\track name=spen_peaks_2ctrl description=spen_peaks_2ctrl colorByStrand='<0,0,255 255,0,0>'" > spen_peaks_2ctrl_ucsc.bed

# Summary of Part 3: got spen_fc_rpm_2ctrl.txt (with counts and rpm) and spen_peaks_2ctrl.bed



#### Part 4: classify peaks
    # classify each peak into protein-coding/lncRNA/other, and exonic/intronic/UTR5/UTR3

### A. Classify peaks as protein-coding (C), lncRNA (NC), or other 
        # (classify as protein-coding if overlapping with both)

## Find lncRNA and protein-coding intervals
    # for increased efficiency, you can only perform this part once

# var: gtf input file with utr annotation
gtf='/proj/calabrlb/users/Zhiyue/22_sp/gtf/gencode.vM25.basic.annotation.utr.gtf'

# extract gtf_nc.bed and gtf_c.bed
    # var: nc_pattern and c_pattern
nc_pattern='transcript_type "bidirectional_promoter_lncRNA"\|transcript_type "macro_lncRNA"\|transcript_type "antisense"\|transcript_type "3prime_overlapping_ncRNA"\|transcript_type "lincRNA"\|transcript_type "processed_transcript"\|transcript_type "sense_intronic"\|transcript_type "sense_overlapping"'
c_pattern='transcript_type "protein_coding"'

cat $gtf | grep "$nc_pattern" | awk -v OFS="\t" '{print $1,$4-1,$5,".","0",$7}' > gtf_nc.bed
cat $gtf | grep "$c_pattern" | awk -v OFS="\t" '{print $1,$4-1,$5,".","0",$7}' > gtf_c.bed

# sort and merge NC/C intervals
sortBed -i gtf_nc.bed > gtf_nc_sorted.bed
bedtools merge -i gtf_nc_sorted.bed -s -c 6 -o distinct | awk -v OFS="\t" '{print $1,$2,$3,".","0",$4}'  > gtf_nc_merged.bed
sortBed -i gtf_c.bed > gtf_c_sorted.bed
bedtools merge -i gtf_c_sorted.bed -s -c 6 -o distinct | awk -v OFS="\t" '{print $1,$2,$3,".","0",$4}'  > gtf_c_merged.bed

## Classify peaks as NC/C

# report peak intersection with NC/C intervals
bedtools intersect -a spen_peaks_2ctrl.bed -b gtf_c_merged.bed -s -wao > spen_peaks_2ctrl_c_intersect.txt
bedtools intersect -a spen_peaks_2ctrl.bed -b gtf_nc_merged.bed -s -wao > spen_peaks_2ctrl_nc_intersect.txt

# merge peak intersection with NC/C intervals
    # ex. peak 50 may have 2 lines in the intersection file
    # This ensures intersect_len.txt has only 1 line for each peak
cat spen_peaks_2ctrl_c_intersect.txt | awk '{arr[$4]+=$13}END{for (i in arr) print i, arr[i]}' | sort -V | cut -d ' ' -f2  > spen_peaks_2ctrl_c_intersect_len.txt
cat spen_peaks_2ctrl_nc_intersect.txt | awk '{arr[$4]+=$13}END{for (i in arr) print i, arr[i]}' | sort -V | cut -d ' ' -f2  > spen_peaks_2ctrl_nc_intersect_len.txt

# label peaks as NC/C
    # (classify as protein-coding if overlapping with both)
paste spen_peaks_2ctrl_c_intersect_len.txt spen_peaks_2ctrl_nc_intersect_len.txt | awk '{if($1 > 0) print "protein_coding"; else if($2 > 0) print "lncRNA"; else print "other"}' > spen_peaks_2ctrl_codingClass.txt



### B. Classify peaks as exonic / intronic / UTR5 / UTR3 
    # (label as exonic if >50% exon intersection)

## Find lncRNA and protein-coding intervals
    # for increased efficiency, you can only perform this part once

# make gtf_exon.bed, gtf_utr5.bed, gtf_utr3.bed
cat $gtf | awk '{if ($3 ~ /exon/) print $0}' | awk -v OFS="\t" '{print $1,$4-1,$5,".","0",$7}' > gtf_exon.bed
cat $gtf | awk '{if ($3 ~ /UTR5/) print $0}' | awk -v OFS="\t" '{print $1,$4-1,$5,".","0",$7}' > gtf_utr5.bed
cat $gtf | awk '{if ($3 ~ /UTR3/) print $0}' | awk -v OFS="\t" '{print $1,$4-1,$5,".","0",$7}' > gtf_utr3.bed

## Classify peaks as exonic / intronic / UTR5 / UTR3 

# report peak intersection with 3 groups
sortBed -i gtf_exon.bed > gtf_exon_sorted.bed
bedtools merge -i gtf_exon_sorted.bed -s -c 6 -o distinct | awk -v OFS="\t" '{print $1,$2,$3,".","0",$4}'  > gtf_exon_merged.bed
bedtools intersect -a spen_peaks_2ctrl.bed -b gtf_exon_merged.bed -s -wao > spen_peaks_2ctrl_exon_intersect.txt

sortBed -i gtf_utr5.bed > gtf_utr5_sorted.bed
bedtools merge -i gtf_utr5_sorted.bed -s -c 6 -o distinct | awk -v OFS="\t" '{print $1,$2,$3,".","0",$4}'  > gtf_utr5_merged.bed
bedtools intersect -a spen_peaks_2ctrl.bed -b gtf_utr5_merged.bed -s -wao > spen_peaks_2ctrl_utr5_intersect.txt

sortBed -i gtf_utr3.bed > gtf_utr3_sorted.bed
bedtools merge -i gtf_utr3_sorted.bed -s -c 6 -o distinct | awk -v OFS="\t" '{print $1,$2,$3,".","0",$4}'  > gtf_utr3_merged.bed
bedtools intersect -a spen_peaks_2ctrl.bed -b gtf_utr3_merged.bed -s -wao > spen_peaks_2ctrl_utr3_intersect.txt

# merge exonic intersection within each peak
cat spen_peaks_2ctrl_exon_intersect.txt | awk '{arr[$4]+=$13}END{for (i in arr) print i, arr[i]}' | sort -V | cut -d ' ' -f2  > spen_peaks_2ctrl_exon_intersect_len.txt
cat spen_peaks_2ctrl_utr5_intersect.txt | awk '{arr[$4]+=$13}END{for (i in arr) print i, arr[i]}' | sort -V | cut -d ' ' -f2  > spen_peaks_2ctrl_utr5_intersect_len.txt
cat spen_peaks_2ctrl_utr3_intersect.txt | awk '{arr[$4]+=$13}END{for (i in arr) print i, arr[i]}' | sort -V | cut -d ' ' -f2  > spen_peaks_2ctrl_utr3_intersect_len.txt

# label peaks as exonic/intronic/utr5/utr3
    # (label as exonic if >50% exon intersection)
    # Note: assuming no header
paste spen_peaks_2ctrl.bed spen_peaks_2ctrl_exon_intersect_len.txt spen_peaks_2ctrl_utr5_intersect_len.txt spen_peaks_2ctrl_utr3_intersect_len.txt | awk '{if ($8 > 0) print "UTR5"; else if ($9 > 0) print "UTR3"; else if( $7/($3-$2) > 0.5) print "exonic"; else print "intronic"}' > spen_peaks_2ctrl_exonClass.txt

### Summary:
# add header to each label
cat spen_peaks_2ctrl_codingClass.txt | sed "1 i\coding_class" > spen_peaks_2ctrl_codingClass_wHeader.txt
cat spen_peaks_2ctrl_exonClass.txt | sed "1 i\exonic_class" > spen_peaks_2ctrl_exonClass_wHeader.txt

# concatenate files
paste spen_fc_rpm_2ctrl.txt spen_peaks_2ctrl_codingClass_wHeader.txt spen_peaks_2ctrl_exonClass_wHeader.txt > spen_peaks_2ctrl_classified_wIndCounts.txt
cat spen_peaks_2ctrl_classified_wIndCounts.txt | tr '\t' ',' > spen_peaks_2ctrl_classified_wIndCounts.csv