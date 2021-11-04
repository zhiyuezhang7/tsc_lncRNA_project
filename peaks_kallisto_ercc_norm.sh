#!/bin/bash

#SBATCH -J pken_%j
#SBATCH -p general
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=12:00:00
#SBATCH -o pken_%j.out
#SBATCH -e pken_%j.err
#SBATCH --mail-type ALL 
#SBATCH --mail-user zhiyue@live.unc.edu

## BEFORE GETTING STARTED
    # Change path to directories appropriately

## Files needed (in current directory)
    # 0. scripts
        # peaks_kallisto_ercc_norm.sh
        # ercc_norm.py
    # 1. 3 RIP-seq output files (with specific antibody, IgG control, RNA input) 
        # ex: 
            # tsc_B5_hnrnpk_rip_S10_R1_001.fastq
            # tsc_B5_igg_mouse-rip_S2_R1_001.fastq.gz
            # tsc_B5_input-rip_S1_R1_001.fastq.gz
    # 2. peaks table 
        # ex: hnrnpk_pirhana_rpkm_ercc_normalization.txt
    # 3. Kallisto transcriptome index
        # ex: gencode.vM25.basic.annotation.complete.ERCC.idx
            # To build it, run the commands below individually (takes approximately 40 min):
                # module load kallisto
                # sbatch --mem 250g -N 1 -n 16 -t 3:00:00 -J kallisto_idx -o kallisto_idx.out -e kallisto_idx.err --wrap='kallisto index -i gencode.vM25.basic.annotation.complete.ERCC.idx gencode.vM25.basic.annotation.complete.ERCC.fa'
    
## Command format: sbatch peaks_kallisto_ercc_norm.sh $1 $2 $3 $4 $5
    # ex: sbatch peaks_kallisto_ercc_norm.sh tsc_B5_hnrnpk_rip_S10_R1_001 tsc_B5_igg_mouse-rip_S2_R1_001 tsc_B5_input-rip_S1_R1_001 hnrnpk_pirhana_rpkm_ercc_normalization hnrnpk

    # $1 = rootname of RIP-seq output with specific antibody 
        # ex: tsc_B5_hnrnpk_rip_S10_R1_001
    # $2 = rootname of RIP-seq output with IgG control 
        # ex: tsc_B5_igg_mouse-rip_S2_R1_001
    # $3 = rootname of RIP-seq output with input RNA 
        # ex: tsc_B5_input-rip_S1_R1_001
    # $4 = rootname of peak table
        # ex: hnrnpk_pirhana_rpkm_ercc_normalization
    # $5 = name of the target RBP
        # ex: hnrnpk

## loading all the tool modules for processing pipeline
    # must use star/2.5.4b (version 20101) due to the STAR indexed files for mm9 being generated with version 20101.
module load star/2.5.4b
module load samtools
module load subread
module load kallisto

## unzip fastq files if needed.
gunzip $1.fastq.gz
gunzip $2.fastq.gz
gunzip $3.fastq.gz

### A. Align peaks to the transcriptome

### A1. Align FASTQ to genome with STAR
## use STAR aligner to align reads to mm10 mouse genome.
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix $1_ --readFilesIn $1.fastq
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix $2_ --readFilesIn $2.fastq
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix $3_ --readFilesIn $3.fastq

## filter with MAPQ>30 and output to SAM file
samtools view -h -Sq 30 $1_Aligned.out.sam > $1_Aligned_filteredsq30.out.sam
samtools view -h -Sq 30 $2_Aligned.out.sam > $2_Aligned_filteredsq30.out.sam
samtools view -h -Sq 30 $3_Aligned.out.sam > $3_Aligned_filteredsq30.out.sam

### A2. Filter alignment_to_genome.SAM with peaks.BED
# only output alignments overlapping the peak BED file
# must split strands first for strand-specific filtering

## convert peaks table to peaks.bed
cat $4.txt | sed 1d | awk '{print $2"\t"$3-1"\t"$4"\t"$1"\t""0""\t"$5}' | sed 's/"//g' > $5_peaks.bed

## split BED by strand
cat $5_peaks.bed | awk '{if ($6=="+") print $0}' > $5_peaks_pos.bed
cat $5_peaks.bed | awk '{if ($6=="-") print $0}' > $5_peaks_neg.bed

## split SAM by strand
    # -h: keep header
samtools view -h -F 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_neg.out.sam
samtools view -h -f 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_pos.out.sam

samtools view -h -F 0x10 $2_Aligned_filteredsq30.out.sam > $2_Aligned_neg.out.sam
samtools view -h -f 0x10 $2_Aligned_filteredsq30.out.sam > $2_Aligned_pos.out.sam

samtools view -h -F 0x10 $3_Aligned_filteredsq30.out.sam > $3_Aligned_neg.out.sam
samtools view -h -f 0x10 $3_Aligned_filteredsq30.out.sam > $3_Aligned_pos.out.sam

## convert SAM to BAM
    # -S: input is SAM (by default it expects BAM)
    # -b: output is BAM (by default it produces BAM) 
samtools view -S -b $1_Aligned_neg.out.sam > $1_neg.bam
samtools view -S -b $1_Aligned_pos.out.sam > $1_pos.bam

samtools view -S -b $2_Aligned_neg.out.sam > $2_neg.bam
samtools view -S -b $2_Aligned_pos.out.sam > $2_pos.bam

samtools view -S -b $3_Aligned_neg.out.sam > $3_neg.bam
samtools view -S -b $3_Aligned_pos.out.sam > $3_pos.bam
 
## sort and index BAM
samtools sort $1_neg.bam -o $1_sorted_neg.bam
samtools sort $1_pos.bam -o $1_sorted_pos.bam
samtools index $1_sorted_neg.bam
samtools index $1_sorted_pos.bam

samtools sort $2_neg.bam -o $2_sorted_neg.bam
samtools sort $2_pos.bam -o $2_sorted_pos.bam
samtools index $2_sorted_neg.bam
samtools index $2_sorted_pos.bam

samtools sort $3_neg.bam -o $3_sorted_neg.bam
samtools sort $3_pos.bam -o $3_sorted_pos.bam
samtools index $3_sorted_neg.bam
samtools index $3_sorted_pos.bam

## filtering: extract SAM entries overlapping BED
    # -L: only output alignments overlapping the input BED FILE.
    # -h: keep header
samtools view $1_sorted_neg.bam -h -L $5_peaks_neg.bed > $1_peaks_neg.out.sam
samtools view $1_sorted_pos.bam -h -L $5_peaks_pos.bed > $1_peaks_pos.out.sam

samtools view $2_sorted_neg.bam -h -L $5_peaks_neg.bed > $2_peaks_neg.out.sam
samtools view $2_sorted_pos.bam -h -L $5_peaks_pos.bed > $2_peaks_pos.out.sam
 
samtools view $3_sorted_neg.bam -h -L $5_peaks_neg.bed > $3_peaks_neg.out.sam
samtools view $3_sorted_pos.bam -h -L $5_peaks_pos.bed > $3_peaks_pos.out.sam

### A3. Align filtered peaks.SAM to transcriptome with Kallisto

## convert peaks.sam back to FASTQ
cat $1_peaks_neg.out.sam | sed '/^@/d' | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $1_peaks.fastq
cat $1_peaks_pos.out.sam | sed '/^@/d' | awk '{print "@"$1"\n"$10"\n+\n"$11}' >> $1_peaks.fastq

cat $2_peaks_neg.out.sam | sed '/^@/d' | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $2_peaks.fastq
cat $2_peaks_pos.out.sam | sed '/^@/d' | awk '{print "@"$1"\n"$10"\n+\n"$11}' >> $2_peaks.fastq
 
cat $3_peaks_neg.out.sam | sed '/^@/d' | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $3_peaks.fastq
cat $3_peaks_pos.out.sam | sed '/^@/d' | awk '{print "@"$1"\n"$10"\n+\n"$11}' >> $3_peaks.fastq

## align both peaks (filtered) and original treatment file (unfiltered) to gencode.vM25.basic.annotation.complete.ERCC.fa
    # -i: location of the transcriptome index
    # -o: kallisto output folder for each dataset
    # --single -l 200 -s 50: single-end data
    # -l and -s are required for SE reads: estimated average and standard deviation of fragment length
kallisto quant -i gencode.vM25.basic.annotation.complete.ERCC.idx -o $1_filtered --single -l 200 -s 50 $1_peaks.fastq
kallisto quant -i gencode.vM25.basic.annotation.complete.ERCC.idx -o $1_unfiltered --single -l 200 -s 50 $1.fastq

kallisto quant -i gencode.vM25.basic.annotation.complete.ERCC.idx -o $2_filtered --single -l 200 -s 50 $2_peaks.fastq
kallisto quant -i gencode.vM25.basic.annotation.complete.ERCC.idx -o $2_unfiltered --single -l 200 -s 50 $2.fastq

kallisto quant -i gencode.vM25.basic.annotation.complete.ERCC.idx -o $3_filtered --single -l 200 -s 50 $3_peaks.fastq
kallisto quant -i gencode.vM25.basic.annotation.complete.ERCC.idx -o $3_unfiltered --single -l 200 -s 50 $3.fastq

### B. Normalize kallisto results

### B1. Make table summarizing estimated counts of each transcript from each group
mkdir table
cat $1_filtered/abundance.tsv | cut -f1,2,3 | sed 1d > table/id_and_lengths.txt

cat $1_filtered/abundance.tsv | cut -f4 | sed 1d > table/$1_filtered_est_cts.txt
cat $1_unfiltered/abundance.tsv | cut -f4 | sed 1d > table/$1_unfiltered_est_cts.txt
cat $2_filtered/abundance.tsv | cut -f4 | sed 1d > table/$2_filtered_est_cts.txt
cat $2_unfiltered/abundance.tsv | cut -f4 | sed 1d > table/$2_unfiltered_est_cts.txt
cat $3_filtered/abundance.tsv | cut -f4 | sed 1d > table/$3_filtered_est_cts.txt
cat $3_unfiltered/abundance.tsv | cut -f4 | sed 1d > table/$3_unfiltered_est_cts.txt

paste table/id_and_lengths.txt table/$3_unfiltered_est_cts.txt table/$3_filtered_est_cts.txt table/$1_unfiltered_est_cts.txt table/$1_filtered_est_cts.txt table/$2_unfiltered_est_cts.txt table/$2_filtered_est_cts.txt > est_cts_all.txt
sed -i  '1i gene_ID\tlength\teff_length\tinput\tinput_peaks\trip\trip_peaks\tigg\tigg_peaks' est_cts_all.txt
cat est_cts_all.txt | tr '\t' ',' > est_cts_all.csv

### B2. ERCC normalization

# extract ERCC est_cts
cat est_cts_all.txt | grep "ERCC" > est_cts_ercc.txt

# convert ERCC est_cts to FPKM
    # FPKM = Estimated_counts / (length / 1000) / (read_counts / 1,000,000)
    # = Estimated_counts * 10^9 / length / read_counts

# calculate read counts of each sample
printf '%s\t' $3 > read_counts.txt
var=$(cat $3.fastq | wc -l) && expr $var / 4 >>  read_counts.txt
printf '%s\t' $1 >> read_counts.txt
var=$(cat $1.fastq | wc -l) && expr $var / 4 >>  read_counts.txt
printf '%s\t' $2 >> read_counts.txt
var=$(cat $2.fastq | wc -l) && expr $var / 4 >>  read_counts.txt

# Normalization
    # calculate ERCC and all FPKM
    # find upper quartile for non-zero ERCC FPKM in each unfiltered sample
    # multiply all filtered rip and igg FPKM by ERCC scaling factor
        # all_rip_ft(FPKM) * ( ercc_input_unft(Q3) / ercc_rip_unft(Q3) )
        # all_igg_ft(FPKM) * ( ercc_input_unft(Q3) / ercc_igg_unft(Q3) )
    # Subtract background noise for each gene
        # all_rip_ft(scaled_FPKM) - all_igg_ft(scaled_FPKM)
    # normalize to input
        # (all_rip_ft(scaled_FPKM) - all_igg_ft(scaled_FPKM)) / all_input_unft(FPKM)
python ercc_norm.py $5

    # output (all peak_filtered except for unfiltered input)
        # fpkm_input.txt
        # fpkm_rip.txt
        # fpkm_igg.txt
        # ercc_rip.txt
        # ercc_igg.txt
        # rip_less_igg.txt
        # rip_over_input.txt

# Now we have the ERCC-normalized FPKM counts for the transcriptome alignment of RIP-seq!!!

## Make summarizing table with gene_ID, length, eff_length + the above 7 columns
cat est_cts_all.txt | cut -f1,2,3 > transcript_info.txt
paste transcript_info.txt fpkm_input.txt fpkm_rip.txt fpkm_igg.txt ercc_rip.txt ercc_igg.txt rip_less_igg.txt rip_over_input.txt > $5_kallisto_fpkm_ercc_normalization.txt
cat $5_kallisto_fpkm_ercc_normalization.txt | tr '\t' ',' > $5_kallisto_fpkm_ercc_normalization.csv

## Final output: proteinName_kallisto_fpkm_ercc_normalization.csv