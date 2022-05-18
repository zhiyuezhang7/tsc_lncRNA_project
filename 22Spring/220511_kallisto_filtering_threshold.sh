#!/bin/bash

#SBATCH -J thre_%j
#SBATCH -p general
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=3:00:00
#SBATCH -o thre_%j.out
#SBATCH -e thre_%j.err
#SBATCH --mail-type ALL 
#SBATCH --mail-user zhiyue@live.unc.edu

module load kallisto

### Goal: Choose threshold for expressed AND chromatin enriched transcripts

# Files needed:
    # total RNAseq files (x4)
    # fractionation files (chromatin, cytoplasm)

cp /proj/calabrlb/users/Zhiyue/21_02_25/total_RNAseq/*fastq* .
cp /proj/calabrlb/users/Zhiyue/21_fall/fractionation/*chromatin* .
cp /proj/calabrlb/users/Zhiyue/21_fall/fractionation/*cytoplasm* .

# 1. run kallisto

kallisto_idx='/proj/calabrlb/users/Zhiyue/22_sp/ref_genome/gencode.vM25.basic.annotation.complete.ERCC.idx'

kallisto quant -i $kallisto_idx -o total_rna_kallisto --single -l 200 -s 50 --rf-stranded d165_krabrtta_rna_S10_R1_001.fastq.gz d165_vprtta_rna_S8_R1_001.fastq.gz d91_krab_rtta_rna_S12_R1_001.fastq.gz d91_vp_rtta_rna_S10_R1_001.fastq.gz
kallisto quant -i $kallisto_idx -o chr_kallisto --rf-stranded tsc_b16_chromatin_rna_S2_R1_001.fastq.gz  tsc_b16_chromatin_rna_S2_R2_001.fastq.gz
kallisto quant -i $kallisto_idx -o cyto_kallisto --rf-stranded tsc_b16_cytoplasm_rna_S4_R1_001.fastq.gz  tsc_b16_cytoplasm_rna_S4_R2_001.fastq.gz

# 2. make table of results (gene info and tpm)
cut_col_num="1,2,3,5,10,15"
paste total_rna_kallisto/abundance.tsv chr_kallisto/abundance.tsv cyto_kallisto/abundance.tsv | cut -f"$cut_col_num" | sed 1d > total_chr_cyto_tpm.txt

## Calculate for each transcript:
    # %chromatin = (chrom_tpm)/(chrom_tpm+cyto_tpm)
cat total_chr_cyto_tpm.txt | awk -v OFS="\t" '{if ($5+$6 > 0) print $0, $5/($5+$6); else print $0, "-1"}' > total_chr_cyto_tpm_perc.txt

cat total_chr_cyto_tpm_perc.txt | cut -f4 > total_tpm.txt
cat total_chr_cyto_tpm_perc.txt | cut -f7 > chr_perc.txt

sed -i  '1i gene_ID\tlength\teff_length\ttotal_tpm\tchr_tpm\tcyto_tpm\tchr%' total_chr_cyto_tpm_perc.txt

# save kallisto results of total rna and fractionation data:
# cp total_rna_kallisto/abundance.tsv /proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/threshold/total_rna_kallisto.txt
# cp chr_kallisto/abundance.tsv /proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/threshold/chr_kallisto.txt
# cp cyto_kallisto/abundance.tsv /proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/threshold/cyto_kallisto.txt
# cp total_chr_cyto_tpm_perc.txt /proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/threshold/

# 3. find threshold with python histogram
    # 220511_kallisto_filtering_threshold_histogram.ipynb
# 4. filtering with thresholds
