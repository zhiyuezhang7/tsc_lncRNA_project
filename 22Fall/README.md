# Identification of Repressive RNAs with Xist-Like Functions in the Mouse Transcriptome: Materials and Methods
Zhiyue Zhang's Honors Thesis project in the Calabrese Lab at UNC School of Medicine


## Quantifying RNA-Protein Interactions
### RNA Immunoprecipitation Data
Mickey Murvin - another member of the Calabrese Lab - performed RNA Immunoprecipitation (RIP) on mouse trophoblast stem cells (TSCs)  with antibodies of 27 RNA-binding proteins important for Xist function: `Aly/Ref,G9a,HnrnpC,HNRNPK,HnrnpM,HnrnpU,Jarid2,LBR,MAtr3,Nudt21,PABPN1,PTBP1,RBM15,Ring1b,RYBP,SAFB,SPEN,SRSF1,SUPT16H,SUZ12,U2AF35,XRN2,Tia1,Ciz1,U2AF65,Nxf1,SFPQ`. `IgG` antibodies were used as control.

I stored the RIP-seq data in `/proj/calabrlb/users/Zhiyue/22_sp/rip/` on Longleaf:
```
tsc_b12_alyref_rip_S6_R1_001.fastq.gz
tsc_b12_g9a_rip_S6_R1_001.fastq.gz
tsc_b12_hnrnpc_rip_S4_R1_001.fastq.gz
tsc_b12_hnrnpk_rip_S7_R1_001.fastq.gz
tsc_b12_hnrnpm_rip_S5_R1_001.fastq.gz
tsc_b12_hnrnpu_rip_S3_R1_001.fastq.gz
tsc_b12_igg_goat_rip_S2_R1_001.fastq.gz
tsc_b12_igg_mouse_rip_S1_R1_001.fastq.gz
tsc_b12_igg_rb_rip_S2_R1_001.fastq.gz
tsc_b12_jarid2_rip_S5_R1_001.fastq.gz
tsc_b12_lbr_rip_S7_R1_001.fastq.gz
tsc_b12_matrin3_rip_S3_R1_001.fastq.gz
tsc_b12_pabpn1_rip_S8_R1_001.fastq.gz
tsc_b12_ptbp1_rip_S9_R1_001.fastq.gz
tsc_b12_rbm15_rip_S8_R1_001.fastq.gz
tsc_b12_rybp_rip_S12_R1_001.fastq.gz
tsc_b12_srsf1_rip_S4_R1_001.fastq.gz
tsc_b15_nudt21_b1_rip_S8_R1_001.fastq.gz
tsc_b15_supt16h_b1_rip_S7_R1_001.fastq.gz
tsc_b15_u2af35_b2_rip_S13_R1_001.fastq.gz
tsc_b15_u2af65_b2_rip_S16_R1_001.fastq.gz
tsc_b15_xrn2_b2_rip_S14_R1_001.fastq.gz
tsc_B5_alyref-rip_S8_R1_001.fastq.gz
tsc_B5_hnrnpk-rip_S10_R1_001.fastq.gz
tsc_B5_igg_mouse-rip_S2_R1_001.fastq.gz
tsc_B5_matrin3-rip_S14_R1_001.fastq.gz
tsc_B5_rybp-rip_S15_R1_001.fastq.gz
tsc_B5_spen-rip_S13_R1_001.fastq.gz
tsc_b8_igg_goat_rip_S3_R1_001.fastq.gz
tsc_b8_igg_rb_rip_S2_R1_001.fastq.gz
tsc_b8_jarid2_rip_S10_R1_001.fastq.gz
tsc_b8_ring1b_rip_S12_R1_001.fastq.gz
tsc_b8_suz12_rip_S16_R1_001.fastq.gz
tsc_c3_g9a_rip_S2_R1_001.fastq.gz
tsc_c3_hnrnpc_rip_S10_R1_001.fastq.gz
tsc_c3_hnrnpm_rip_S9_R1_001.fastq.gz
tsc_c3_hnrnpu_rip_S8_R1_001.fastq.gz
tsc_c3_lbr_rip_S4_R1_001.fastq.gz
tsc_c3_nudt21_rip_S14_R1_001.fastq.gz
tsc_c3_pabpn1_rip_S3_R1_001.fastq.gz
tsc_c3_ptbp1_rip_S5_R1_001.fastq.gz
tsc_c3_rbm15_rip_S7_R1_001.fastq.gz
tsc_c3_srsf1_rip_S6_R1_001.fastq.gz
tsc_c3_supt16h_rip_S13_R1_001.fastq.gz
tsc_c3_tia1_1_rip_S15_R1_001.fastq.gz
tsc_c3_tia1_2_rip_S16_R1_001.fastq.gz
tsc_c3_u2af35_rip_S12_R1_001.fastq.gz
tsc_c3_xrn2_rip_S11_R1_001.fastq.gz
tsc_c9_ciz1_2_rip_S18_R1_001.fastq.gz
tsc_c9_hnrnpc_rip_S15_R1_001.fastq.gz
tsc_c9_hnrnpm_rip_S14_R1_001.fastq.gz
tsc_c9_hnrnpu_rip_S16_R1_001.fastq.gz
tsc_c9_spen_rip_S21_R1_001.fastq.gz
tsc_c9_srsf1_rip_S19_R1_001.fastq.gz
tsc_c9_u2af35_rip_S20_R1_001.fastq.gz
tsc_nxf1_S8_R1_001.fastq.gz
tsc_ring1b_S7_R1_001.fastq.gz
tsc_safb1-rip_mm_S19_R1_001.fastq
tsc_sfpq_S10_R1_001.fastq.gz
tsc_spen_novus_S15_R1_001.fastq.gz
tsc_suz12_S6_R1_001.fastq.gz
```

### Computing Protein-Binding Signals
I performed the following steps for RIP-seq data with each of the 27 antibodies:
1. Align RIP-seq data to the mm10 genome
2. Call peaks of protein-binding
3. Find true peaks enriched over IgG
4. Find RIP-seq reads at true peaks
5. Align RIP-seq reads at true peaks to the transcriptome (reference used: `gencode.vM25.basic.annotation.complete.ERCC.fa`) 

I first ran `22fa_rank_rip_signal_igg.sh` to prepare IgG data.

Then I ran `22fa_rank_rip_signal.sh` for each of the 27 proteins:
```
sbatch 22fa_rank_rip_signal.sh pabpn1
sbatch 22fa_rank_rip_signal.sh hnrnpk
sbatch 22fa_rank_rip_signal.sh alyref
sbatch 22fa_rank_rip_signal.sh g9a
sbatch 22fa_rank_rip_signal.sh hnrnpc
sbatch 22fa_rank_rip_signal.sh hnrnpm
sbatch 22fa_rank_rip_signal.sh hnrnpu
sbatch 22fa_rank_rip_signal.sh jarid2
sbatch 22fa_rank_rip_signal.sh lbr
sbatch 22fa_rank_rip_signal.sh matrin3
sbatch 22fa_rank_rip_signal.sh nudt21
sbatch 22fa_rank_rip_signal.sh ptbp1
sbatch 22fa_rank_rip_signal.sh rbm15
sbatch 22fa_rank_rip_signal.sh ring1b
sbatch 22fa_rank_rip_signal.sh rybp
sbatch 22fa_rank_rip_signal.sh safb
sbatch 22fa_rank_rip_signal.sh spen
sbatch 22fa_rank_rip_signal.sh srsf1
sbatch 22fa_rank_rip_signal.sh supt16h
sbatch 22fa_rank_rip_signal.sh suz12
sbatch 22fa_rank_rip_signal.sh u2af35
sbatch 22fa_rank_rip_signal.sh xrn2
sbatch 22fa_rank_rip_signal.sh tia1
sbatch 22fa_rank_rip_signal.sh ciz1
sbatch 22fa_rank_rip_signal.sh u2af65
sbatch 22fa_rank_rip_signal.sh nxf1
sbatch 22fa_rank_rip_signal.sh sfpq
```

To show the importance of using IgG as a control, the wiggle track of RIP-seq data and the true peaks called can be viewed at [this UCSC session](https://genome.ucsc.edu/s/zhiyue/zhiyue_220511_rip_peaks).

## Extracting Expressed, Chromatin-Enriched RNAs
### Total RNA-seq and Fractionation Data
Previous analysis was on the whole transcriptome. To select only RNAs that are highly expressed and chromatin enriched in TSCs, I used total RNA-seq and fractionation data to determine thresholds of expression and chromatin enrichment.
I used the total RNA-seq data of TSCs published by the lab in [Molecular Cell](https://pubmed.ncbi.nlm.nih.gov/31256989/) (Schertzer et al., 2019), as stored in `/proj/calabrlb/users/Zhiyue/21_02_25/total_RNAseq/`:
```
d165_krabrtta_rna_S10_R1_001.fastq.gz
d165_vprtta_rna_S8_R1_001.fastq.gz
d91_krab_rtta_rna_S12_R1_001.fastq.gz
d91_vp_rtta_rna_S10_R1_001.fastq.gz
```
I used the fractionation data generated by Mickey Murvin, as stored in `/proj/calabrlb/users/Zhiyue/21_fall/fractionation/`:
```
tsc_b16_chromatin_rna_S2_R1_001.fastq.gz
tsc_b16_chromatin_rna_S2_R2_001.fastq.gz
tsc_b16_cytoplasm_rna_S4_R1_001.fastq.gz
tsc_b16_cytoplasm_rna_S4_R2_001.fastq.gz
```

### Filtering Expressed, Chromatin-Enriched RNAs 
I ran `22fa_filtering_threshold.sh` to compute the total, cytoplasmic, and chromatinic expression of transcripts.
Then I ran `22fa_filtering_threshold.ipynb` to plot histogram of total expression and chromatin enrichment: %chromatin = chromatin_TPM /(chromatin_TPM + cytoplasm_TPM). I determined the thresholds based on the inflection points.
I ran `22fa_filtering.sh` to filter the protein-binding profile data to only select expressed, chromatin-enriched RNAs. The resulting data of selected RNAs is `all27rbp_kallisto_igg_rpm_filtered.csv`.

## Quantifying K-mer Similarities between Selected RNAs and Xist, Airn, Kcnq1ot1
I extracted the sequence data of the RNAs selected in previous steps:
```
# linearize all.fa
FA="/proj/calabrlb/users/Zhiyue/22_sp/ref_genome/gencode.vM25.basic.annotation.complete.ERCC.fa"
cat $FA | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' > gencode_all_linear.fa

# get gene_name_13788.txt
kallisto_dat="/proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/all27rbp_kallisto_igg_rpm_filtered.txt"
cat $kallisto_dat | cut -f1 | sed '1d' > gene_name_13788.txt

# extract 13788.fa
cat gencode_all_linear.fa | grep -f gene_name_13788.txt > RNA_13788_linear.fa
cat RNA_13788_linear.fa | tr "\t" "\n" > RNA_13788.fa
```
In my local terminal, I ran [SEEKR](https://github.com/CalabreseLab/seekr) to compute the k-mer similarities between these RNAs and Xist, Airn, Kcnq1ot1. Input files `gencode.vM25.basic.annotation.complete.ERCC.fa` and `mm10_XKA_mc.fa` were needed:
```
seekr_norm_vectors gencode.vM25.basic.annotation.complete.ERCC.fa  -k 6  -mv mean_6mers.npy -sv std_6mers.npy
seekr_kmer_counts RNA_13788.fa -o RNA_13788_6mers.csv -k 6 -mv mean_6mers.npy -sv std_6mers.npy
seekr_kmer_counts mm10_XKA_mc.fa -o xka_6mers.csv -k 6  -l pre -mv mean_6mers.npy -sv std_6mers.npy
seekr_pearson RNA_13788_6mers.csv xka_6mers.csv -o RNA_13788_vs_xka.csv
```
Overall, I have selected 13788 RNAs with their protein-binding information stored in `all27rbp_kallisto_igg_rpm_filtered.csv` and k-mer information stored in `RNA_13788_vs_xka.csv`.

## Computing Correlation
I ran `22fa_correlation.Rmd` to evaluate correlation between k-mer and protein-binding similarities.

I added annotation of whether each transcript is spliced or unspliced, lncRNA or other types of RNAs. To add information of lncRNA, I found the names of lncRNA transcripts on Unix:
```
# get lncRNA_gene_name.txt
GTF="/proj/calabrlb/users/Zhiyue/21_02_25/gencode.vM25.basic.annotation.gtf"
grep 'transcript_type "bidirectional_promoter_lncRNA"\|transcript_type "macro_lncRNA"\|transcript_type "antisense"\|transcript_type "3prime_overlapping_ncRNA"\|transcript_type "lincRNA"\|transcript_type "processed_transcript"\|transcript_type "sense_intronic"\|transcript_type "sense_overlapping"' $GTF > gtf_lncRNA_extract.txt
cat gtf_lncRNA_extract.txt | cut -f9 | cut -f4 -d';' | grep -o '".*"' | sed 's/"//g' | sort |  uniq -d > lncRNA_gene_name.txt

# get kallisto_lncRNA_gene_name.txt
kallisto_dat="/proj/calabrlb/users/Zhiyue/22_sp/kallisto_filtering/all27rbp_kallisto_igg_rpm_filtered.txt"
cat $kallisto_dat | cut -f1 | sed '1d' > gene_name_13788.txt
cat gene_name_13788.txt | grep -f lncRNA_gene_name.txt > gene_name_13788_lncRNA.txt
```
