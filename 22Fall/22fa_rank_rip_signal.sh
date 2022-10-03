#!/bin/bash

#SBATCH -J ripRank_%j
#SBATCH -p general
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=12:00:00
#SBATCH -o ripRank_%j.out
#SBATCH -e ripRank_%j.err
#SBATCH --mail-type ALL 
#SBATCH --mail-user zhiyue@live.unc.edu

## argument: RBP_name (all lower cases!)
    # example command:
        # sbatch 22fa_rank_rip_signal.sh $1
        # $1 = rbp name (all lower cases) e.g., alyref

## input files:
    # rip files:
        # /proj/calabrlb/users/Zhiyue/22_sp/rip/*$1*
    # bam alignment of igg (sorted and indexed):
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_neg.bam
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_pos.bam
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_neg.bam.bai
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_pos.bam.bai
    # kallisto reference index:
        # /proj/calabrlb/users/Zhiyue/22_sp/ref_genome/gencode.vM25.basic.annotation.complete.ERCC.idx

## load modules
module load star/2.5.4b
module load samtools
module load macs
module load subread
module load bedtools
module load kallisto

### 1. Align RIP-seq data to the genome; generate wiggle tracks

## put rip files and do following work in directory named $1
mkdir $1
cd $1
cp /proj/calabrlb/users/Zhiyue/22_sp/rip/*$1* .
gunzip *.fastq.gz
    # concatenate all fastq file names, separated by comma; remove spaces and newlines
rip_filenames=$(ls -m *.fastq | sed -r 's/\s+//g' | tr -d '\n')
echo $rip_filenames

## use STAR aligner to align reads to mm10 mouse genome.
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix $1_ --readFilesIn $rip_filenames
# Filter and output to sam file
samtools view -h -Sq 30 $1_Aligned.out.sam > $1_Aligned_filteredsq30.out.sam

## make wiggle tracks
# split sam by strand
samtools view -h -F 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_neg.out.sam
samtools view -h -f 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_pos.out.sam

# Generate wiggle tracks
    # Positive strand is blue; negative strand is red
    # All colors usable: pink, babyblue, grey, black, red, green, blue, yellow, lightgreen, purple, navyblue maroon, brown, orange, or magenta
cp /proj/calabrlb/users/Zhiyue/22_sp/code/bigsam_to_wig_mm10.pl .
perl bigsam_to_wig_mm10.pl $1_Aligned_neg.out.sam $1_neg red
perl bigsam_to_wig_mm10.pl $1_Aligned_pos.out.sam $1_pos blue

# save sam and wig files
# cp $1_Aligned_filteredsq30.out.sam /proj/calabrlb/users/Zhiyue/22_sp/rip_sam
# cp $1_neg.wig /proj/calabrlb/users/Zhiyue/22_sp/rip_wig
# cp $1_pos.wig /proj/calabrlb/users/Zhiyue/22_sp/rip_wig


### 2. Call peaks of protein-binding

# randomize strand
cp /proj/calabrlb/users/Zhiyue/22_sp/code/macs_strand_rand_sam.pl .
perl macs_strand_rand_sam.pl $1_Aligned_pos.out.sam $1_Aligned_pos.out.rand.sam
perl macs_strand_rand_sam.pl $1_Aligned_neg.out.sam $1_Aligned_neg.out.rand.sam

# macs
    # merge size = 304 = 4 * 76
macs2 callpeak -t $1_Aligned_pos.out.rand.sam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir $1_pos_peaks -n $1_pos
macs2 callpeak -t $1_Aligned_neg.out.rand.sam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir $1_neg_peaks -n $1_neg 
    # add trackline for UCSC
cat $1_neg_peaks/$1_neg_peaks.broadPeak | sed "1 i\track type=broadPeak name="$1_peaks_neg" description="$1_peaks_neg" nextItemButton=on" > $1_neg_peaks_ucsc.broadPeak
cat $1_pos_peaks/$1_pos_peaks.broadPeak | sed "1 i\track type=broadPeak name="$1_peaks_pos" description="$1_peaks_pos" nextItemButton=on" > $1_pos_peaks_ucsc.broadPeak

# save .broadPeak files
# cp $1_neg_peaks_ucsc.broadPeak /proj/calabrlb/users/Zhiyue/22_sp/rip_peaks_ucsc/
# cp $1_pos_peaks_ucsc.broadPeak /proj/calabrlb/users/Zhiyue/22_sp/rip_peaks_ucsc/


### 3. Find true peaks enriched over IgG

## convert broadPeak to SAF file in preparation for featureCounts
cat $1_neg_peaks_ucsc.broadPeak | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"-"}' > $1_saf.txt
cat $1_pos_peaks_ucsc.broadPeak | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"+"}' >> $1_saf.txt
cat $1_saf.txt | awk -v OFS="\t" '{print NR,$0}' | sed '1 i\GeneID\tChr\tStart\tEnd\tStrand' > $1_peaks.saf

## featureCounts with SAF (both rip and igg): count reads assigned to each peak
featureCounts -s 2 -F SAF -a $1_peaks.saf -o $1_fc $1_Aligned_filteredsq30.out.sam
featureCounts -s 2 -F SAF -a $1_peaks.saf -o $1_igg_fc /proj/calabrlb/users/Zhiyue/22_sp/rip_sam/igg_Aligned_filteredsq30.out.sam

# make summary file
cut -f7 $1_igg_fc > $1_igg_counts.txt
paste $1_fc $1_igg_counts.txt | sed 1d > $1_fc.txt
cat $1_fc.txt | tr "\t" ","  > $1_fc.csv

# save summary file
# cp $1_fc.txt /proj/calabrlb/users/Zhiyue/22_sp/rip_peaks_fc/
# cp $1_fc.csv /proj/calabrlb/users/Zhiyue/22_sp/rip_peaks_fc/

# normalize read count to rpm
    # rpm = (count * 1000,000 / reads)
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

# save files
# cp $1_fc_rpm.txt /proj/calabrlb/users/Zhiyue/22_sp/rip_peaks_fc_rpm
# cp $1_fc_rpm.csv /proj/calabrlb/users/Zhiyue/22_sp/rip_peaks_fc_rpm
# cp $1_peaks_2igg.bed /proj/calabrlb/users/Zhiyue/22_sp/rip_peaks_2igg
# cp $1_peaks_2igg_ucsc.bed /proj/calabrlb/users/Zhiyue/22_sp/rip_peaks_2igg_ucsc


### 4. Find RIP-seq reads at true peaks

# split peak strands for strand-specific filtering
cat $1_peaks_2igg.bed | awk '{if ($6=="+") print $0}' > $1_peaks_2igg_pos.bed
cat $1_peaks_2igg.bed | awk '{if ($6=="-") print $0}' > $1_peaks_2igg_neg.bed

# split sam by strand
samtools view -h -F 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_neg.out.sam
samtools view -h -f 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_pos.out.sam

## convert SAM to BAM
    # -S: take in SAM (by default it expects BAM)
    # -b: output is BAM (by default it produces BAM) 
samtools view -S -b $1_Aligned_neg.out.sam > $1_neg.bam
samtools view -S -b $1_Aligned_pos.out.sam > $1_pos.bam

## sort and index BAM
samtools sort $1_neg.bam -o $1_sorted_neg.bam
samtools sort $1_pos.bam -o $1_sorted_pos.bam
samtools index $1_sorted_neg.bam
samtools index $1_sorted_pos.bam

# save files to proj: bam
# cp $1_sorted_neg.bam /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/
# cp $1_sorted_pos.bam /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/
# cp $1_sorted_neg.bam.bai /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/
# cp $1_sorted_pos.bam.bai /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/

## filtering: extract SAM entries overlapping peaks.bed from both rip and igg
    # -L: only output alignments overlapping the BED FILE.
    # -h: keep header
samtools view $1_sorted_neg.bam -h -L $1_peaks_2igg_neg.bed > $1_peaks_neg.out.sam
samtools view $1_sorted_pos.bam -h -L $1_peaks_2igg_pos.bed > $1_peaks_pos.out.sam

samtools view igg_sorted_neg.bam -h -L $1_peaks_2igg_neg.bed > igg_$1_peaks_neg.out.sam
samtools view igg_sorted_pos.bam -h -L $1_peaks_2igg_pos.bed > igg_$1_peaks_pos.out.sam

## To convert peaks.sam and igg_peaks.sam back to FASTQ, must do bam processing first
# convert SAM to BAM
samtools view -S -b $1_peaks_neg.out.sam > $1_peaks_neg.bam
samtools view -S -b $1_peaks_pos.out.sam > $1_peaks_pos.bam

samtools view -S -b igg_$1_peaks_neg.out.sam > igg_$1_peaks_neg.bam
samtools view -S -b igg_$1_peaks_pos.out.sam > igg_$1_peaks_pos.bam

# convert BAM to FASTQ
samtools fastq $1_peaks_neg.bam > $1_peaks.fastq
samtools fastq $1_peaks_pos.bam >> $1_peaks.fastq

samtools fastq igg_$1_peaks_neg.bam > igg_$1_peaks.fastq
samtools fastq igg_$1_peaks_pos.bam >> igg_$1_peaks.fastq

### 5. Align RIP-seq reads at true peaks to the transcriptome

# var: kallisto index (input file)
    # derived from gencode.vM25.basic.annotation.complete.ERCC.fa:
    # sbatch --mem 250g -N 1 -n 16 -t 3:00:00 -J kallisto_idx -o kallisto_idx.out -e kallisto_idx.err --wrap='kallisto index -i gencode.vM25.basic.annotation.complete.ERCC.idx gencode.vM25.basic.annotation.complete.ERCC.fa'
kallisto_idx='/proj/calabrlb/users/Zhiyue/22_sp/ref_genome/gencode.vM25.basic.annotation.complete.ERCC.idx'

## kallisto alignment
    # -i: location of the transcriptome index
    # -o: kallisto output folder for each dataset
    # --single -l 200 -s 50: single-end data
    # -l and -s are required for SE reads: estimated average and standard deviation of fragment length
    # rf-stranded: reverse stranded
kallisto quant -i $kallisto_idx -o $1_kallisto --single -l 200 -s 50 --rf-stranded $1_peaks.fastq
kallisto quant -i $kallisto_idx -o igg_$1_kallisto --single -l 200 -s 50 --rf-stranded igg_$1_peaks.fastq

## Make table summarizing estimated counts
cat $1_kallisto/abundance.tsv | cut -f1,2,3,4 | sed 1d > $1_kallisto_id_est_cts.txt
cat igg_$1_kallisto/abundance.tsv | cut -f4 | sed 1d > igg_$1_kallisto_est_cts.txt
paste $1_kallisto_id_est_cts.txt igg_$1_kallisto_est_cts.txt > $1_kallisto_igg.txt

## Normalization: convert estimated counts to rpm
# calculate rip reads
fastq_lin_num=$(wc -l rip_fastq/*$1*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
(( reads = fastq_lin_num / 4 ))

# col6 = rip rpm; col7 = igg rpm; col8 = rip_rpm - igg_rpm = 6 - 7
cat $1_kallisto_igg.txt | awk -v OFS="\t" -v reads=$reads -v reads_igg=124344953 '{print $0, (($4*1000000)/reads), (($5*1000000)/reads_igg)}' | awk -v OFS="\t" '{print $0, $6-$7}' > $1_kallisto_igg_rpm.txt
sed -i "1i gene_ID\tlength\teff_length\t$1_est_cts\tigg_est_cts\t$1_rpm\tigg_rpm\t$1_rpm_over_igg" $1_kallisto_igg_rpm.txt
cat $1_kallisto_igg_rpm.txt | tr '\t' ',' > $1_kallisto_igg_rpm.csv

# # save files to proj: kallisto tables
# cp $1_kallisto_igg_rpm.txt /proj/calabrlb/users/Zhiyue/22_sp/kallisto/
# cp $1_kallisto_igg_rpm.csv /proj/calabrlb/users/Zhiyue/22_sp/kallisto/