#!/bin/bash

#SBATCH -J fil_%j
#SBATCH -p general
#SBATCH --time=0:15:00
#SBATCH -o fil_%j.out
#SBATCH -e fil_%j.err
#SBATCH --mail-type ALL 
#SBATCH --mail-user zhiyue@live.unc.edu

# Goal: make table with kallisto results from all 27 RBPs; add kallisto results from totalRNA and fractionation (chromatin and cytoplasm) data
    # then filter with total_RNA_TPM > 0.25; chromatin% > 0.75

# extract both rip (rpm) and rip_less_igg (rpm) from individual kallisto tables
    # just putting these 2 columns side by side
    # we have 27 files; we want the 6th and 8th column from each file
cut_col_num=$(seq 8 8 216 | awk -v OFS="" '{ print $1-2","$1","}' | tr -d '\n' | sed 's/.$//')
cut_col_num="${cut_col_num}"
paste /proj/calabrlb/users/Zhiyue/22_sp/kallisto/*_kallisto_igg_rpm.txt | cut -f"$cut_col_num" > all27rbp_kallisto_igg_rpm_0.txt
paste /proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/threshold/total_chr_cyto_tpm_perc.txt all27rbp_kallisto_igg_rpm_0.txt > all27rbp_kallisto_igg_rpm.txt

cat all27rbp_kallisto_igg_rpm.txt | tr '\t' ',' > all27rbp_kallisto_igg_rpm.csv

# # save to proj
# cp all27rbp_kallisto_igg_rpm.* /proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/

# filter kallisto matrix with selected rna satisfying thresholds
cp /proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/threshold/total_chr_cyto_tpm_perc.txt bg_list.txt
cp /proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/all27rbp_kallisto_igg_rpm.txt our_list.txt
    # copy header to output
head -1 our_list.txt > our_list_filtered.txt
sed -i 1d our_list.txt

# 0.25 = expression level; 0.75 = chr% level
cat bg_list.txt | sed 1d | awk -v OFS="\t" '{if ($4 > 0.25 && $7 > 0.75) print 1; else print 0}' > rna_filter.txt
# our list has 3+2*27 = 57 col; paste the T/F of rna select before 1st col
paste rna_filter.txt our_list.txt | awk -v OFS="\t" '{if ($1 == 1) print $0}' | cut -f2- >> our_list_filtered.txt

# # save to proj
# cp our_list_filtered.txt /proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/all27rbp_kallisto_igg_rpm_filtered.txt
# cat our_list_filtered.txt | tr '\t' ',' > /proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/all27rbp_kallisto_igg_rpm_filtered.csv
