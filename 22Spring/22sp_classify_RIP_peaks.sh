#!/bin/bash

#SBATCH -J peakCls_%j
#SBATCH -p general
#SBATCH --mem=50G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH -o peakCls_%j.out
#SBATCH -e peakCls_%j.err
#SBATCH --mail-type ALL 
#SBATCH --mail-user zhiyue@live.unc.edu

## argument: RBP_name (all lower cases!!!)
    # example command: (TO EDIT)
    # sbatch name.sh alyref

## input files:
    # /proj/calabrlb/users/Zhiyue/22_sp/rip/*$1* = all rip files needed
    # /proj/calabrlb/users/Zhiyue/22_sp/gtf/gencode.vM25.basic.annotation.utr.gtf

## load modules
module load star/2.5.4b
module load samtools
module load macs
module load subread
module load bedtools


#### Part 1: from rip, STAR genome alignment and MACS peak-calling 
    # (important note: run this part with keyword igg first before running the entire pipeline!)

## put rip files and do following work in directory named $1
mkdir $1
cd $1
cp /proj/calabrlb/users/Zhiyue/22_sp/rip/*$1* .
gunzip *.fastq.gz
    # concatenate all fastq file names, separated by comma
        # remove spaces and newlines
rip_filenames=$(ls -m *.fastq | sed -r 's/\s+//g' | tr -d '\n')
echo $rip_filenames

### 1A. STAR alignment (merging all input files)
## use STAR aligner to align reads to mm10 mouse genome.
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix $1_ --readFilesIn $rip_filenames

## Filter and output to sam file
samtools view -h -Sq 30 $1_Aligned.out.sam > $1_Aligned_filteredsq30.out.sam

## sam to wig
# split sam by strand
samtools view -h -F 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_neg.out.sam
samtools view -h -f 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_pos.out.sam

# Generate wiggle tracks
    # Positive strand is blue; negative strand is red
    # All colors usable: pink, babyblue, grey, black, red, green, blue, yellow, lightgreen, purple, navyblue maroon, brown, orange, or magenta
    # move bigsam_to_wig_mm10.pl to directory
cp /proj/calabrlb/users/Zhiyue/22_sp/code/bigsam_to_wig_mm10.pl .
perl bigsam_to_wig_mm10.pl $1_Aligned_neg.out.sam $1_neg red
perl bigsam_to_wig_mm10.pl $1_Aligned_pos.out.sam $1_pos blue


### 1B. MACS

# load randomization code
cp /proj/calabrlb/users/Zhiyue/22_sp/code/macs_strand_rand_sam.pl .

# randomize strand
perl macs_strand_rand_sam.pl $1_Aligned_pos.out.sam $1_Aligned_pos.out.rand.sam
perl macs_strand_rand_sam.pl $1_Aligned_neg.out.sam $1_Aligned_neg.out.rand.sam

# macs
    # merge size = 304 = 4 * 76
            # NOTE: check if trackline is correct
macs2 callpeak -t $1_Aligned_pos.out.rand.sam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir $1_pos_peaks -n $1_pos_peaks
macs2 callpeak -t $1_Aligned_neg.out.rand.sam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir $1_neg_peaks -n $1_neg_peaks 
    # add trackline for UCSC
cat $1_neg_peaks/$1_neg_peaks.broadPeak | sed "1 i\track type=broadPeak name="$1_peaks_neg" description="$1_peaks_neg" nextItemButton=on" > $1_neg_peaks_ucsc.broadPeak
cat $1_pos_peaks/$1_pos_peaks.broadPeak | sed "1 i\track type=broadPeak name="$1_peaks_pos" description="$1_peaks_pos" nextItemButton=on" > $1_pos_peaks_ucsc.broadPeak



#### Part 2: featureCounts count number of reads in each peak

## convert broadPeak to SAF file in preparation for featureCounts
cat $1_neg_peaks_ucsc.broadPeak | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"-"}' > $1_saf.txt
cat $1_pos_peaks_ucsc.broadPeak | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"+"}' >> $1_saf.txt
cat $1_saf.txt | awk -v OFS="\t" '{print NR,$0}' | sed '1 i\GeneID\tChr\tStart\tEnd\tStrand' > $1_peaks.saf

## featureCounts with SAF (both rip and igg)
featureCounts -s 2 -F SAF -a $1_peaks.saf -o $1_fc $1_Aligned_filteredsq30.out.sam
featureCounts -s 2 -F SAF -a $1_peaks.saf -o $1_igg_fc /proj/calabrlb/users/Zhiyue/22_sp/rip_sam/igg_Aligned_filteredsq30.out.sam
    # to edit: location of igg file

# make summary file
cut -f7 $1_igg_fc > $1_igg_counts.txt
paste $1_fc $1_igg_counts.txt | sed 1d > $1_fc.txt
cat $1_fc.txt | tr "\t" ","  > $1_fc.csv



#### Part 3: find true peaks enriched over igg (rpm > 2 * igg_rpm)

## calculate rpm = (count * 1000,000 / reads)
    # count = $7 and $8 of fc table
   
# var: reads = number of lines in fastq files / 4
fastq_lin_num=$(wc -l *$1*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
(( reads = fastq_lin_num / 4 ))

# reads_igg = 124344953
    # only calculate once:
        # cp /proj/calabrlb/users/Zhiyue/22_sp/rip/*igg* .
        # gunzip *igg*.fastq.gz
        # fastq_lin_num_igg=$(wc -l *igg*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
        # (( reads_igg = fastq_lin_num_igg / 4 ))
        # echo $reads_igg

# print rpm to $9 and $10 columns
cat $1_fc.txt | sed 1d | awk -v OFS="\t" -v reads=$reads -v reads_igg=124344953 '{print $0, (($7*1000000)/reads), (($8*1000000)/reads_igg)}' | sed "1 i\GeneID\tChr\tStart\tEnd\tStrand\tLength\t$1_counts\tigg_counts\t$1_rpm\tigg_rpm" > $1_fc_rpm.txt
cat $1_fc_rpm.txt | tr "\t" ","  > $1_fc_rpm.csv

# filter peaks with rpm 2 folds over igg
cat $1_fc_rpm.txt | awk '{if(NR == 1) print $0; else if($9 > 2 * $10) print $0}' > $1_fc_rpm_2igg.txt

# convert enricheded peaks to BED format
cat $1_fc_rpm_2igg.txt | sed 1d | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > $1_peaks_2igg.bed

# add trackline for UCSC
cat $1_peaks_2igg.bed | sed "1 i\track name=$1_peaks_2igg description=$1_peaks_2igg colorByStrand='<0,0,255 255,0,0>'" > $1_peaks_2igg_ucsc.bed



#### Part 4: classify peaks
    # classify each peak into protein-coding/lncRNA/other, and exonic/intronic/UTR5/UTR3

### A. Classify peaks as protein-coding (C), lncRNA (NC), or other 
        # (classify as protein-coding if overlapping with both)

# for increased efficiency, only perform this part once

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

# report peak intersection with NC/C intervals
bedtools intersect -a $1_peaks_2igg.bed -b gtf_c_merged.bed -s -wao > $1_peaks_2igg_c_intersect.txt
bedtools intersect -a $1_peaks_2igg.bed -b gtf_nc_merged.bed -s -wao > $1_peaks_2igg_nc_intersect.txt

# merge peak intersection with NC/C intervals
    # ex. peak 50 may have 2 lines in the intersection file
    # This ensures intersect_len.txt has only 1 line for each peak
cat $1_peaks_2igg_c_intersect.txt | awk '{arr[$4]+=$13}END{for (i in arr) print i, arr[i]}' | sort -V | cut -d ' ' -f2  > $1_peaks_2igg_c_intersect_len.txt
cat $1_peaks_2igg_nc_intersect.txt | awk '{arr[$4]+=$13}END{for (i in arr) print i, arr[i]}' | sort -V | cut -d ' ' -f2  > $1_peaks_2igg_nc_intersect_len.txt

# label peaks as NC/C
    # (classify as protein-coding if overlapping with both)
paste $1_peaks_2igg_c_intersect_len.txt $1_peaks_2igg_nc_intersect_len.txt | awk '{if($1 > 0) print "protein_coding"; else if($2 > 0) print "lncRNA"; else print "other"}' > $1_peaks_2igg_codingClass.txt



### B. Classify peaks as exonic / intronic / UTR5 / UTR3 
    # (label as exonic if >50% exon intersection)

# for increased efficiency, only perform this part once
# only need

# make gtf_exon.bed, gtf_utr5.bed, gtf_utr3.bed
cat $gtf | awk '{if ($3 ~ /exon/) print $0}' | awk -v OFS="\t" '{print $1,$4-1,$5,".","0",$7}' > gtf_exon.bed
cat $gtf | awk '{if ($3 ~ /UTR5/) print $0}' | awk -v OFS="\t" '{print $1,$4-1,$5,".","0",$7}' > gtf_utr5.bed
cat $gtf | awk '{if ($3 ~ /UTR3/) print $0}' | awk -v OFS="\t" '{print $1,$4-1,$5,".","0",$7}' > gtf_utr3.bed

# report peak intersection with 3 groups
sortBed -i gtf_exon.bed > gtf_exon_sorted.bed
bedtools merge -i gtf_exon_sorted.bed -s -c 6 -o distinct | awk -v OFS="\t" '{print $1,$2,$3,".","0",$4}'  > gtf_exon_merged.bed
bedtools intersect -a $1_peaks_2igg.bed -b gtf_exon_merged.bed -s -wao > $1_peaks_2igg_exon_intersect.txt

sortBed -i gtf_utr5.bed > gtf_utr5_sorted.bed
bedtools merge -i gtf_utr5_sorted.bed -s -c 6 -o distinct | awk -v OFS="\t" '{print $1,$2,$3,".","0",$4}'  > gtf_utr5_merged.bed
bedtools intersect -a $1_peaks_2igg.bed -b gtf_utr5_merged.bed -s -wao > $1_peaks_2igg_utr5_intersect.txt

sortBed -i gtf_utr3.bed > gtf_utr3_sorted.bed
bedtools merge -i gtf_utr3_sorted.bed -s -c 6 -o distinct | awk -v OFS="\t" '{print $1,$2,$3,".","0",$4}'  > gtf_utr3_merged.bed
bedtools intersect -a $1_peaks_2igg.bed -b gtf_utr3_merged.bed -s -wao > $1_peaks_2igg_utr3_intersect.txt


# merge exonic intersection within each peak
cat $1_peaks_2igg_exon_intersect.txt | awk '{arr[$4]+=$13}END{for (i in arr) print i, arr[i]}' | sort -V | cut -d ' ' -f2  > $1_peaks_2igg_exon_intersect_len.txt
cat $1_peaks_2igg_utr5_intersect.txt | awk '{arr[$4]+=$13}END{for (i in arr) print i, arr[i]}' | sort -V | cut -d ' ' -f2  > $1_peaks_2igg_utr5_intersect_len.txt
cat $1_peaks_2igg_utr3_intersect.txt | awk '{arr[$4]+=$13}END{for (i in arr) print i, arr[i]}' | sort -V | cut -d ' ' -f2  > $1_peaks_2igg_utr3_intersect_len.txt

# label peaks as exonic/intronic/utr5/utr3
    # (label as exonic if >50% exon intersection)
    # Note: assuming no header
paste $1_peaks_2igg.bed $1_peaks_2igg_exon_intersect_len.txt $1_peaks_2igg_utr5_intersect_len.txt $1_peaks_2igg_utr3_intersect_len.txt | awk '{if ($8 > 0) print "UTR5"; else if ($9 > 0) print "UTR3"; else if( $7/($3-$2) > 0.5) print "exonic"; else print "intronic"}' > $1_peaks_2igg_exonClass.txt

### Summary:
paste $1_peaks_2igg.bed $1_peaks_2igg_codingClass.txt $1_peaks_2igg_exonClass.txt > $1_peaks_2igg_classified.txt
