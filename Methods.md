# Methods and Materials

## Quantification of Gene Expression

To determine what lncRNAs count as "expressed" in TSCs, I used the total RNA-seq data of TSCs published by the lab in [Molecular Cell](https://pubmed.ncbi.nlm.nih.gov/31256989/) (Schertzer et al., 2019).
I used the following files:
```
d165_krabrtta_rna_S10_R1_001.fastq.gz
d165_vprtta_rna_S8_R1_001.fastq.gz
d91_krab_rtta_rna_S12_R1_001.fastq.gz
d91_vp_rtta_rna_S10_R1_001.fastq.gz
```
To align reads to the mouse genome, I used [RSEM and bowtie](https://pubmed.ncbi.nlm.nih.gov/21816040/).
The reference mouse genome used in alignment below was derived from [GENCODE](https://www.gencodegenes.org/).
```
gencode.vM25.basic.annotation.complete.ERCC.fa
```
I ran RSEM commands in UNIX. (Unless otherwise specified, all commands listed in this file are UNIX commands.)
```
module load rsem
module load bowtie/1.2.2
sbatch --mem 50g --mail-type ALL --mail-user zhiyue@live.unc.edu -N 1 -n 16 -o bb_ref.out -e bb_ref.err --wrap='bowtie-build --threads 16 gencode.vM25.basic.annotation.complete.ERCC.fa gencode.vM25_ercc_rsem'
sbatch --mem 50g --mail-type ALL --mail-user zhiyue@live.unc.edu -N 1 -n 8 -o rsem_ref.out -e rsem_ref.err --wrap='rsem-prepare-reference gencode.vM25.basic.annotation.complete.ERCC.fa gencode.vM25_ercc_rsem'
sbatch -N 1 -n 16 --mail-type ALL --mail-user zhiyue@live.unc.edu --mem 50g -t 24:00:00 -J rsemgc -o rsemgc.out -e rsemgc.err --wrap='rsem-calculate-expression -p 16 --forward-prob 0 d91_krab_rtta_rna_S12_R1_001.fastq,d165_krabrtta_rna_S10_R1_001.fastq,d91_vp_rtta_rna_S10_R1_001.fastq,d165_vprtta_rna_S8_R1_001.fastq gencode.vM25_ercc_rsem tsc_dcas9_gcvM25_rs'
squeue -u zhiyue
```
