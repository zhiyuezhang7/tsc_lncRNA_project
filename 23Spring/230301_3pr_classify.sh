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

### Goals: find true peaks of rip enriched over 2*igg;
    # count reads within each peak, including reads from all rips, igg, and individual rips (if rip is more than one);
    # classify peaks based on exonic/intronic/utr5/utr3, and coding/lncRNA/other

## input files needed: 
    # rips in $2
        # If you supply >1 rip files, they will be indexed in alphabetical order (the order of "ls" command).
        # Ex. For alyref.fastq and ciz1.fastq, alyref will be rip1 and ciz1 will be rip2.
    # igg alignment:
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_sam/igg_Aligned_filteredsq30.out.sam
    # randomization code:
        # /proj/calabrlb/users/Zhiyue/22_sp/code/macs_strand_rand_sam.pl
    # gencode annotation (with utf info):
        # /proj/calabrlb/users/Zhiyue/22_sp/gtf/gencode.vM25.basic.annotation.utr.gtf

## Example command:
# sbatch 230301_3pr_classify.sh ezh2 "/proj/calabrlb/users/Zhiyue/23_sp/0221_3pr_rip/*ezh2*"
    # $1 = rbp name (all lower cases) e.g., alyref
    # $2 = all rip files; you can put them in someDirectory and use "someDirectory/*"

## load modules
module load star/2.5.4b
module load samtools
module load macs
module load subread
module load bedtools

#### Part 1: from rip, STAR genome alignment and MACS peak-calling 

## do following work in directory named $1_classify
mkdir $1_classify
cd $1_classify

### 1A. STAR alignment

## STAR alignment of all rips
# move rip files here
mkdir $1
cp ${2} ${1}
gunzip $1/*.fastq.gz

# use STAR aligner to align reads to mm10 mouse genome.
rip_filenames=$(ls -m $1/*.fastq | sed -r 's/\s+//g' | tr -d '\n')
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix $1_ --readFilesIn $rip_filenames

# Filter and output to sam file
samtools view -h -Sq 30 $1_Aligned.out.sam > $1_Aligned_filteredsq30.out.sam

## STAR alignment of individual rips
# decide if there are more than 1 rip
rip_file_num=$(ls $1 | wc -l)
if [ $rip_file_num -gt 1 ]
then

    IFS=',' read -r -a rip_arr <<< "$rip_filenames" # make array with rip file names

    # align each rip file
    for i in ${!rip_arr[@]}; do # loop on indices of array
        rip_name="${rip_arr[$i]}" # for each file name
        echo The rip${i} file is $rip_name. >> $1_rip_info.txt # NOTE: check in $1_rip_info.txt to see which rip gets assigned to which number!
        j=$(($i + 1)) # index+1 so that it's 1 based for the file name

        # align to mm10 genome
        STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix $1_rip${j}_ --readFilesIn $rip_name
        
        # filter to sam file
        samtools view -h -Sq 30 $1_rip${j}_Aligned.out.sam > $1_rip${j}_Aligned_filteredsq30.out.sam
    done

fi


### 1B. MACS peak-calling

# split sam by strand
samtools view -h -F 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_neg.out.sam
samtools view -h -f 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_pos.out.sam

# load randomization code
cp /proj/calabrlb/users/Zhiyue/22_sp/code/macs_strand_rand_sam.pl .

# randomize strand
perl macs_strand_rand_sam.pl $1_Aligned_pos.out.sam $1_Aligned_pos.out.rand.sam
perl macs_strand_rand_sam.pl $1_Aligned_neg.out.sam $1_Aligned_neg.out.rand.sam

# macs
    # merge size = 304 = 4 * 76
macs2 callpeak -t $1_Aligned_pos.out.rand.sam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir $1_pos_peaks -n $1_pos
macs2 callpeak -t $1_Aligned_neg.out.rand.sam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir $1_neg_peaks -n $1_neg 
    # add trackline for UCSC
cat $1_neg_peaks/$1_neg_peaks.broadPeak | sed "1 i\track type=broadPeak name="$1_peaks_neg" description="$1_peaks_neg" nextItemButton=on" > $1_neg_peaks_ucsc.broadPeak
cat $1_pos_peaks/$1_pos_peaks.broadPeak | sed "1 i\track type=broadPeak name="$1_peaks_pos" description="$1_peaks_pos" nextItemButton=on" > $1_pos_peaks_ucsc.broadPeak



#### Part 2: featureCounts count number of reads in each peak

## convert broadPeak to SAF file in preparation for featureCounts
cat $1_neg_peaks_ucsc.broadPeak | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"-"}' > $1_saf.txt
cat $1_pos_peaks_ucsc.broadPeak | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"+"}' >> $1_saf.txt
cat $1_saf.txt | awk -v OFS="\t" '{print NR,$0}' | sed '1 i\GeneID\tChr\tStart\tEnd\tStrand' > $1_peaks.saf

## featureCounts with SAF
# count reads in each peak for all_rips and igg
featureCounts -s 2 -F SAF -a $1_peaks.saf -o $1_fc $1_Aligned_filteredsq30.out.sam
featureCounts -s 2 -F SAF -a $1_peaks.saf -o $1_igg_fc /proj/calabrlb/users/Zhiyue/22_sp/rip_sam/igg_Aligned_filteredsq30.out.sam

# extract raw reads
cut -f7 $1_igg_fc > $1_igg_counts.txt # TO-DO: 2 LINE HEADER

# count reads in each peak for individual rips
if [ $rip_file_num -gt 1 ]
then

    for i in ${!rip_arr[@]}; do # loop on indices
        j=$(($i + 1)) # index+1 so that it's 1 based
        featureCounts -s 2 -F SAF -a $1_peaks.saf -o $1_rip${j}_fc $1_rip${j}_Aligned_filteredsq30.out.sam
        
        # extract raw reads
        cut -f7 $1_rip${j}_fc > $1_rip${j}_counts.txt
        
        # prep for concatenate count files
        if [ $j == 1 ]
        then
        rip_count_filename="$1_rip1_counts.txt"
        else
        rip_count_filename="${rip_count_filename} $1_rip${j}_counts.txt"
        fi
    done

    # make summary file: all_rip, igg, individual rips
    paste $1_fc $1_igg_counts.txt $rip_count_filename | sed 1d > $1_fc.txt

else

    # make summary file: rip, igg (when only one rip file)
    paste $1_fc $1_igg_counts.txt | sed 1d > $1_fc.txt

fi



#### Part 3: find true peaks enriched over ctrl (rpm > 2 * ctrl_rpm)

## calculate rpm = (count * 1000,000 / reads)

# rpm conversion for all rips:
# var: reads = number of lines in fastq files / 4
fastq_lin_num=$(wc -l $1/*$1*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
(( reads = fastq_lin_num / 4 ))
cat $1_fc.txt | sed 1d | awk -v OFS="\t" -v reads=$reads '{print $0, (($7*1000000)/reads)}' > $1_fc_allRpm.txt

# TO-DO: OKAY UP UNTIL THIS POINT. PROBLEM IN THIS LINE? 2 header only removed one!!!

# rpm conversion for igg:
    # reads_igg 
        # only calculate once:
            # cp /proj/calabrlb/users/Zhiyue/22_sp/rip/*igg* .
            # gunzip *igg*.fastq.gz
            # fastq_lin_num_igg=$(wc -l *igg*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
            # (( reads_igg = fastq_lin_num_igg / 4 ))
            # echo $reads_igg
cat $1_igg_counts.txt  | sed 1,2d | awk -v reads_igg=124344953 '{print (($1*1000000)/reads_igg)}' > $1_igg_rpm.txt

# rpm conversion for individual rips:

if [ $rip_file_num -gt 1 ]
then

    for i in ${!rip_arr[@]}; do # loop on indices
        j=$(($i + 1)) # index+1 so that it's 1 based
        # count reads
        fastq_lin_num=$(wc -l ${rip_arr[$i]} | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
        (( reads = fastq_lin_num / 4 ))
        # rpm conversion
        cat $1_rip${j}_counts.txt | sed 1,2d | awk -v reads=$reads '{print (($1*1000000)/reads)}' > $1_rip${j}_rpm.txt

        # prep for concatenate rpm files and file header
        if [ $j == 1 ]
        then
        rip_count_header="$1_rip1_counts"
        rip_rpm_header="$1_rip1_rpm"
        rip_rpm_filename="$1_rip1_rpm.txt"
        else
        rip_count_header="${rip_count_header}\t$1_rip${j}_counts"
        rip_rpm_header="${rip_rpm_header}\t$1_rip${j}_rpm"
        rip_rpm_filename="${rip_rpm_filename} $1_rip${j}_rpm.txt"
        fi
    done

    # make summary of rpm: all rips, igg, individual rips
    paste $1_fc_allRpm.txt $1_igg_rpm.txt $rip_rpm_filename | sed "1 i\GeneID\tChr\tStart\tEnd\tStrand\tLength\t$1_counts\tigg_counts\t${rip_count_header}\t$1_rpm\tigg_rpm\t${rip_rpm_header}" > $1_fc_rpm.txt

else

    # make summary of rpm: rip, igg
    paste $1_fc_allRpm.txt $1_igg_rpm.txt | sed "1 i\GeneID\tChr\tStart\tEnd\tStrand\tLength\t$1_counts\tigg_counts\t$1_rpm\tigg_rpm" > $1_fc_rpm.txt

fi

## filter peaks with rpm 2 folds over igg

if [ $rip_file_num -gt 1 ]
then
    ((rip_all_col = rip_file_num + 9)) # calculate column number of rip_rpm
    ((igg_col = rip_all_col + 1)) # calculate column number of igg_rpm
    cat $1_fc_rpm.txt | awk -v c1=$rip_all_col -v c2=$igg_col  '{if(NR == 1) print $0; else if($c1 > 2 * $c2) print $0}' > $1_fc_rpm_2igg.txt
else
    cat $1_fc_rpm.txt | awk '{if(NR == 1) print $0; else if($9 > 2 * $10) print $0}' > $1_fc_rpm_2igg.txt
fi

# Note: all following commands are about peaks and do not take into account how many rip files you have.

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
sed -i '1 i\coding_class' $1_peaks_2igg_codingClass.txt
sed -i '1 i\exon_class' $1_peaks_2igg_exonClass.txt
paste $1_fc_rpm_2igg.txt $1_peaks_2igg_codingClass.txt $1_peaks_2igg_exonClass.txt > $1_peaks_2igg_classified.txt
cat $1_peaks_2igg_classified.txt | tr '\t' ',' > $1_peaks_2igg_classified.csv