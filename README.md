# TSC lncRNA project: Methods and Materials
Undergraduate research project in the Calabrese Lab, Spring 2021

## Quantification of gene expression
To determine what lncRNAs count as "expressed" in mouse TSCs, I used the total RNA-seq data of TSCs published by the lab in [Molecular Cell](https://pubmed.ncbi.nlm.nih.gov/31256989/) (Schertzer et al., 2019).
I used the following files:
```
d165_krabrtta_rna_S10_R1_001.fastq.gz
d165_vprtta_rna_S8_R1_001.fastq.gz
d91_krab_rtta_rna_S12_R1_001.fastq.gz
d91_vp_rtta_rna_S10_R1_001.fastq.gz
```
### Alignment to whole genome
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
### Extraction of lncRNA
I extracted only lncRNAs from the isoforms.results file produced by RSEM. The transcript types of lncRNAs were specified in [the Ensembl genome browser](https://www.gencodegenes.org/pages/biotypes.html), excluding “non-coding”. The gene names of all RNAs were derived from the vM25 Basic gene annotation (CHR) gtf file provided by [GENCODE](https://www.gencodegenes.org/mouse/).
```
GTF=$(ls *.gtf)
echo $GTF 
```
The output above should be the name of the gtf file used.
```
grep 'transcript_type "bidirectional_promoter_lncRNA"\|transcript_type "macro_lncRNA"\|transcript_type "antisense"\|transcript_type "3prime_overlapping_ncRNA"\|transcript_type "lincRNA"\|transcript_type "processed_transcript"\|transcript_type "sense_intronic"\|transcript_type "sense_overlapping"\|' $GTF > gtf_lncRNA_extract.txt
cat gtf_lncRNA_extract.txt | cut -f5 -d';' | sort | uniq -c
```
The command above should only print lncRNA types. 
```
cat gtf_lncRNA_extract.txt | cut -f9 | cut -f4 -d';' | grep -o '".*"' | sed 's/"//g' | sort |  uniq -d > lncRNA_gene_name.txt
cat tsc_dcas9_gcvM25_rs.isoforms.results | grep -f lncRNA_gene_name.txt > tsc_rsem_lncRNA_all.txt
cat tsc_dcas9_gcvM25_rs.isoforms.results | cut -f6 | grep -v 'TPM' > TPM_All.txt
cat tsc_rsem_lncRNA_all.txt | cut -f6 > TPM_w.txt
```
### Determination of "expressed" TPM value
I made histograms of the expression levels (TPM values) of all genes and only lncRNA genes in TSCs, with Matplotlib and NumPy in Python.
To download the graphs, comment *plt.show()* and decomment the last two lines of code.
For all genes:
```
import matplotlib.pyplot as plt
import numpy as np
import math
 
file = open("TPM_All.txt")
d = file.read().splitlines() # list
e = [math.log(float(z)+0.001,2) for z in d]
 
plt.figure(figsize=(100,30))
n, bins, patches = plt.hist(x=e, bins="auto", color='#4B9CD3', alpha=0.5, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xticks(fontsize=80)
plt.yticks(fontsize=80)
plt.xlabel('log2(TPM)',fontsize=100)
plt.ylabel('Frequency',fontsize=100)
plt.title('Expression Levels of All Genes in Mouse TSCs\n',fontsize=100,fontweight='bold')
maxfreq = n.max()
plt.ylim(ymax=10000)
plt.xlim(xmax=10)
 
plt.axvline(x=-2, ymin=0, ymax=1, linewidth=4, color='r')
plt.text(x=-2, y=-400, s='−2', fontsize=80, color='r', horizontalalignment='center')
 
plt.show()
#plt.savefig('Histogram_TPM_all_TSC.png')
#plt.savefig('Histogram_TPM_all_TSC.pdf')
```
For lncRNA genes:
```
import matplotlib.pyplot as plt
import numpy as np
import math
 
file = open("TPM_w.txt")
d = file.read().splitlines() # list
e = [math.log(float(z)+0.001,2) for z in d]
 
plt.figure(figsize=(100,20))
n, bins, patches = plt.hist(x=e, bins="auto", color='#4B9CD3', alpha=0.5, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xticks(fontsize=60)
plt.yticks(fontsize=60)
plt.xlabel('TPM value',fontsize=80)
plt.ylabel('Frequency',fontsize=80)
plt.title('Expression Levels of lncRNA Genes in Mouse TSCs',fontsize=80)
maxfreq = n.max()
plt.ylim(ymax=10000)
plt.axvline(x=-2, ymin=0, ymax=1, linewidth=4, color='r')
plt.text(x=-2, y=-500, s='−2', fontsize=60, color='r', horizontalalignment='center')
 
plt.show()
#plt.savefig('Histogram_lncRNA.png')
#plt.savefig('Histogram_lncRNA.pdf')
```
The rough inflection point of TPM_All.txt determined the TPM benchmark for “expressed” lncRNAs to be 0.25. I extracted the RSEM results of all expressed lncRNAs.
```
cat tsc_rsem_lncRNA_all.txt | awk '$6 >= 0.25' > tsc_rsem_lncRNA_expressed.txt
```

## Sequence Comparison with SEEKR
