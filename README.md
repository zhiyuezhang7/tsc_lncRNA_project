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
### Alignment to the whole genome
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
### Extraction of lncRNA expression levels
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

## Sequence comparison with SEEKR
The second step is to run [SEEKR](https://pubmed.ncbi.nlm.nih.gov/30224646/) on the sequences of expressed lncRNAs versus the sequences of Xist, Airn, and Kcnq1ot1 (derived from [the UCSC genome browser](https://genome.ucsc.edu/)), and partition expressed lncRNAs into groups based on their *k*-mer profiles.

### Extraction of lncRNA sequences
I first generated two files with the sequences of lncRNAs and expressed lncRNAs, respectively. The fa file used containing all mouse sequences was still gencode.vM25.basic.annotation.complete.ERCC.fa, with its name changed to gencode_all.fa below.
```
cat gencode_all.fa | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' > gencode_all_linear.fa
cat gencode_all_linear.fa | grep -f lncRNA_gene_name.txt > gencode_lncRNAs_linear.fa
cat gencode_lncRNAs_linear.fa | tr "\t" "\n" > gencode_lncRNAs.fa
cat tsc_rsem_lncRNA_expressed.txt | cut -f1 | uniq > tsc_rsem_lncRNA_expressed_info.txt
cat gencode_lncRNAs_linear.fa | grep -f tsc_rsem_lncRNA_expressed_info.txt > gencode_lncRNAs_expressed_linear.fa
cat gencode_lncRNAs_expressed_linear.fa | tr "\t" "\n" > gencode_expressed_lncRNAs.fa
```
### SEEKR
After downloading gencode_lncRNAs.fa and gencode_expressed_lncRNAs.fa, I installed SEEKR from my local Anaconda, and ran SEEKR on my local terminal.
```
seekr_norm_vectors gencode_lncRNAs.fa  -k 6  -l pre -mv mean_6mers.npy -sv std_6mers.npy
seekr_kmer_counts gencode_expressed_lncRNAs.fa -o expressed_6mers.csv -k 6  -l pre -mv mean_6mers.npy -sv std_6mers.npy
seekr_kmer_counts mm10_XKA_mc.fa.txt -o xka_6mers.csv -k 6  -l pre -mv mean_6mers.npy -sv std_6mers.npy
seekr_pearson expressed_6mers.csv xka_6mers.csv -o expressed_v_xka.csv
```

## Analysis of SEEKR results
The result file of SEEKR was expressed_v_xka.csv.
### Sequential similarity to *Xist*, *Airn*, and *Kcnq1ot1*
In UNIX, I extracted expressed lncRNAs' Pearson's R values to Xist, Airn, and Kcnq1ot1, respectively (referred to as XAK below).
```
cat expressed_v_xka.csv |  cut -f2 -d','  | sed 1d > R_X.txt
cat expressed_v_xka.csv |  cut -f3 -d','  | sed 1d > R_A.txt
cat expressed_v_xka.csv |  cut -f4 -d','  | sed 1d > R_K.txt
```
In excel, I calculated the mean, standard deviation, and 99.9th percentile (z-score>=3) for each list. 
I made histograms for each list with python.

For expressed vs *Xist*:
```
import matplotlib.pyplot as plt
import numpy as np
 
file = open("R_X.txt")
d = file.read().splitlines() # list
e = [round(float(z), 4) for z in d]
 
plt.figure(figsize=(100,30))
 
n, bins, patches = plt.hist(x=e, bins="auto", color='#4B9CD3', alpha=0.5, rwidth=0.85)
 
plt.grid(axis='y', alpha=0.75)
 
plt.xticks(fontsize=80)
plt.yticks(fontsize=80)
plt.xlabel("Pearson's R-Values",fontsize=100)
plt.ylabel('Frequency',fontsize=100)
plt.title('Pearson Correlations between k-mer Profiles of        and Expressed lncRNAs\n',fontsize=100,fontweight='bold')
plt.text(x=0.314, y=325, s='Xist', style='italic', fontsize=100,fontweight='bold')
 
maxfreq = n.max()
 
plt.ylim(ymax=300)
plt.xlim(xmax=0.8)
plt.xlim(xmin=-0.42)
 
plt.axvline(x=0.2777, ymin=0, ymax=1, linewidth=4, color='r')
plt.text(x=0.2777, y=-12, s='0.2777', fontsize=80, color='r', horizontalalignment='center')
 
plt.text(x=0.4, y=175, s='Mean = 0.0038', style='italic', fontsize=100)
plt.text(x=0.4, y=150, s='Standard Deviation = 0.0913', style='italic', fontsize=100)
 
plt.show()
#plt.savefig('Histogram_R_Expressed_vs_Xist.png')
#plt.savefig('Histogram_R_Expressed_vs_Xist.pdf')
```
For expressed vs *Airn*:
```
import matplotlib.pyplot as plt
import numpy as np
 
file = open("R_A.txt")
d = file.read().splitlines() # list
e = [round(float(z), 4) for z in d]
 
plt.figure(figsize=(100,30))
 
n, bins, patches = plt.hist(x=e, bins="auto", color='#4B9CD3', alpha=0.5, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xticks(fontsize=80)
plt.yticks(fontsize=80)
plt.xlabel("Pearson's R-Values",fontsize=100)
plt.ylabel('Frequency',fontsize=100)
plt.title('Pearson Correlations between k-mer Profiles of        and Expressed lncRNAs\n',fontsize=100,fontweight='bold')
plt.text(x=0.314, y=325, s='Airn', style='italic', fontsize=100,fontweight='bold')
 
plt.ylim(ymax=300)
plt.xlim(xmax=0.8)
plt.xlim(xmin=-0.42)
 
plt.axvline(x=0.2573, ymin=0, ymax=1, linewidth=4, color='r')
plt.text(x=0.2573, y=-12, s='0.2573', fontsize=80, color='r', horizontalalignment='center')
 
plt.text(x=0.4, y=175, s='Mean = -0.0036', style='italic', fontsize=100)
plt.text(x=0.4, y=150, s='Standard Deviation = 0.0870', style='italic', fontsize=100)
 
plt.show()
#plt.savefig('Histogram_R_Expressed_vs_Airn.png')
#plt.savefig('Histogram_R_Expressed_vs_Airn.pdf')
```
For expressed vs *Kcnq1ot1*:
```
import matplotlib.pyplot as plt
import numpy as np
 
file = open("R_K.txt")
d = file.read().splitlines() # list
e = [round(float(z), 4) for z in d]
 
plt.figure(figsize=(100,30))
 
n, bins, patches = plt.hist(x=e, bins="auto", color='#4B9CD3', alpha=0.5, rwidth=0.85)
 
plt.grid(axis='y', alpha=0.75)
plt.xticks(fontsize=80)
plt.yticks(fontsize=80)
plt.xlabel("Pearson's R-Values",fontsize=100)
plt.ylabel('Frequency',fontsize=100)
plt.title('Pearson Correlations between k-mer Profiles of                 and Expressed lncRNAs\n',fontsize=100,fontweight='bold')
plt.text(x=0.28, y=325, s='Kcnq1ot1', style='italic', fontsize=100,fontweight='bold')
 
plt.ylim(ymax=300)
plt.xlim(xmax=0.8)
plt.xlim(xmin=-0.42)
 
plt.axvline(x=0.4644, ymin=0, ymax=1, linewidth=4, color='r')
plt.text(x=0.4644, y=-12, s='0.4644', fontsize=80, color='r', horizontalalignment='center')
 
plt.text(x=0.4, y=175, s='Mean = -0.0094', style='italic', fontsize=100)
plt.text(x=0.4, y=150, s='Standard Deviation = 0.1579', style='italic', fontsize=100)
 
plt.show()
#plt.savefig('Histogram_R_Expressed_vs_Kcnq1ot1.png')
#plt.savefig('Histogram_R_Expressed_vs_Kcnq1ot1.pdf')
```
### Partition into groups
I partitioned these 3 lists of expressed lncRNAs into 4 groups, respectively - the 99.9th percentile group, and 3 equal-sized tersiles.

Note on file names:

R_X_1.txt -> 99.9th percentile group for *Xist*

R_X_ter1.txt -> first tersile group for *Xist* (In this case, each tersile has 924 or 925 lncRNAs.)
```
cat expressed_v_xka.csv | tr ',' '\t' |sed 1d > expressed_v_xka.txt
cat expressed_v_xka.txt | cut -f1,2 | awk '{if($2>=0.277711657){print $0}}' >R_X_1.txt
cat expressed_v_xka.txt | cut -f1,3 | awk '{if($2>=0.257312462){print $0}}' >R_A_1.txt
cat expressed_v_xka.txt | cut -f1,4 | awk '{if($2>=0.464400786){print $0}}' >R_K_1.txt

cat expressed_v_xka.txt | cut -f1,2 | sort -r -g -k 2 | sed -n '32,955p' > R_X_ter1.txt
cat expressed_v_xka.txt | cut -f1,2 | sort -r -g -k 2 |sed -n '956,1879 p' > R_X_ter2.txt
cat expressed_v_xka.txt | cut -f1,2 | sort -r -g -k 2 |sed -n '1880,2803 p' > R_X_ter3.txt
cat expressed_v_xka.txt | cut -f1,3 | sort -r -g -k 2 |sed -n '31,955 p' > R_A_ter1.txt
cat expressed_v_xka.txt | cut -f1,3 | sort -r -g -k 2 |sed -n '956,1879 p' > R_A_ter2.txt
cat expressed_v_xka.txt | cut -f1,3 | sort -r -g -k 2 |sed -n '1880,2803 p' > R_A_ter3.txt
cat expressed_v_xka.txt | cut -f1,4 | sort -r -g -k 2 |sed -n '31,955 p' > R_K_ter1.txt
cat expressed_v_xka.txt | cut -f1,4 | sort -r -g -k 2 |sed -n '956,1879 p' > R_K_ter2.txt
cat expressed_v_xka.txt | cut -f1,4 | sort -r -g -k 2 |sed -n '1880,2803 p' > R_K_ter3.txt
```
For each group, I then generated the number of spliced versus unspliced RNAs, GC-content, mean RNA length, and made box and whisker plots of gene expression.
### Number of Unspliced RNAs
The results would be printed to the terminal.
```
cat R_X_1.txt | grep '.*unspliced.*' | wc -l
cat R_X_2.txt | grep '.*unspliced.*' | wc -l
cat R_A_1.txt | grep '.*unspliced.*' | wc -l
cat R_A_2.txt | grep '.*unspliced.*' | wc -l
cat R_K_1.txt | grep '.*unspliced.*' | wc -l
cat R_K_2.txt | grep '.*unspliced.*' | wc -l
cat R_X_ter1.txt | grep '.*unspliced.*' | wc -l
cat R_X_ter2.txt | grep '.*unspliced.*' | wc -l
cat R_X_ter3.txt | grep '.*unspliced.*' | wc -l
cat R_A_ter1.txt | grep '.*unspliced.*' | wc -l
cat R_A_ter2.txt | grep '.*unspliced.*' | wc -l
cat R_A_ter3.txt | grep '.*unspliced.*' | wc -l
cat R_K_ter1.txt | grep '.*unspliced.*' | wc -l
cat R_K_ter2.txt | grep '.*unspliced.*' | wc -l
cat R_K_ter3.txt | grep '.*unspliced.*' | wc -l
```
### GC-content
First derive the sequences for each group.
```
cat R_X_1.txt | cut -f1 > R_X_1names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_X_1names.txt > R_X_1seq.txt
wc -l R_X_1seq.txt

cat R_A_1.txt | cut -f1 > R_A_1names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_A_1names.txt > R_A_1seq.txt
wc -l R_A_1seq.txt

cat R_K_1.txt | cut -f1 > R_K_1names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_K_1names.txt > R_K_1seq.txt
wc -l R_K_1seq.txt
```
The GC-content for each group would be printed to the terminal.
```
cat R_X_1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }' 
cat R_A_1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }'  
cat R_K_1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }' 

cat R_X_ter1.txt | cut -f1 > R_X_ter1names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_X_ter1names.txt > R_X_ter1seq.txt
cat R_X_ter1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }'
 
cat R_X_ter2.txt | cut -f1 > R_X_ter2names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_X_ter2names.txt > R_X_ter2seq.txt
cat R_X_ter2seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }'
 
cat R_X_ter3.txt | cut -f1 > R_X_ter3names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_X_ter3names.txt > R_X_ter3seq.txt
cat R_X_ter3seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }'
 
cat R_A_ter1.txt | cut -f1 > R_A_ter1names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_A_ter1names.txt > R_A_ter1seq.txt
cat R_A_ter1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }'
 
cat R_A_ter2.txt | cut -f1 > R_A_ter2names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_A_ter2names.txt > R_A_ter2seq.txt
cat R_A_ter2seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }'
 
cat R_A_ter3.txt | cut -f1 > R_A_ter3names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_A_ter3names.txt > R_A_ter3seq.txt
cat R_A_ter3seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }'
 
cat R_K_ter1.txt | cut -f1 > R_K_ter1names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_K_ter1names.txt > R_K_ter1seq.txt
cat R_K_ter1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }'
 
cat R_K_ter2.txt | cut -f1 > R_K_ter2names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_K_ter2names.txt > R_K_ter2seq.txt
cat R_K_ter2seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }'
 
cat R_K_ter3.txt | cut -f1 > R_K_ter3names.txt
cat gencode_lncRNAs_expressed_linear.fa | grep -f R_K_ter3names.txt > R_K_ter3seq.txt
cat R_K_ter3seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f%%\n", (gc*100)/(gc+at) }'
```

### Mean RNA length
The results would be printed to the terminal.
```
cat R_X_1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'
cat R_A_1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'
cat R_K_1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'

cat R_X_ter1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'
 
cat R_X_ter2seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'
 
cat R_X_ter3seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'
 
cat R_A_ter1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'
 
cat R_A_ter2seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'
 
cat R_A_ter3seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'
 
cat R_K_ter1seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'
 
cat R_K_ter2seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'
 
cat R_K_ter3seq.txt | cut -f2 | tr '/t' '/n' | awk '!/^>/{line++; char+=gsub(/[aAtTnNcCgG]/,"");} END{ printf "%s%s %s%.2f\n", "Number of sequences: ", line, " Mean RNA length: ", (char/line) }'
```
### Expression levels
First, I extracted the TPM values of each group.
```
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_X_1names.txt | cut -f6 > R_X_1TPM.txt
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_A_1names.txt | cut -f6 > R_A_1TPM.txt
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_K_1names.txt | cut -f6 > R_K_1TPM.txt
 
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_X_ter1names.txt | cut -f6 > R_X_ter1TPM.txt
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_X_ter2names.txt | cut -f6 > R_X_ter2TPM.txt
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_X_ter3names.txt | cut -f6 > R_X_ter3TPM.txt
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_A_ter1names.txt | cut -f6 > R_A_ter1TPM.txt
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_A_ter2names.txt | cut -f6 > R_A_ter2TPM.txt
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_A_ter3names.txt | cut -f6 > R_A_ter3TPM.txt
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_K_ter1names.txt | cut -f6 > R_K_ter1TPM.txt
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_K_ter2names.txt | cut -f6 > R_K_ter2TPM.txt
cat tsc_rsem_lncRNA_expressed.txt | sed -e 's/^/>/' | grep -f R_K_ter3names.txt | cut -f6 > R_K_ter3TPM.txt
```
Then I made boxplots of each list of TPM values with python.
```
import matplotlib.pyplot as plt
import numpy as np
import math
 
file = open("R_X_1TPM.txt")
d0 = file.read().splitlines() # list
x0 = [math.log(float(z)+0.001,2) for z in d0]
 
file = open("R_X_ter1TPM.txt")
d1 = file.read().splitlines() # list
x1 = [math.log(float(z)+0.001,2) for z in d1]
 
file = open("R_X_ter2TPM.txt")
d2 = file.read().splitlines() # list
x2 = [math.log(float(z)+0.001,2) for z in d2]
 
file = open("R_X_ter3TPM.txt")
d3 = file.read().splitlines() # list
x3 = [math.log(float(z)+0.001,2) for z in d3]
 
file = open("R_A_1TPM.txt")
d4 = file.read().splitlines() # list
a0 = [math.log(float(z)+0.001,2) for z in d4]
 
file = open("R_A_ter1TPM.txt")
d5 = file.read().splitlines() # list
a1 = [math.log(float(z)+0.001,2) for z in d5]
 
file = open("R_A_ter2TPM.txt")
d6 = file.read().splitlines() # list
a2 = [math.log(float(z)+0.001,2) for z in d6]
 
file = open("R_A_ter3TPM.txt")
d7 = file.read().splitlines() # list
a3 = [math.log(float(z)+0.001,2) for z in d7]
 
file = open("R_K_1TPM.txt")
d8 = file.read().splitlines() # list
k0 = [math.log(float(z)+0.001,2) for z in d8]
 
file = open("R_K_ter1TPM.txt")
d9 = file.read().splitlines() # list
k1 = [math.log(float(z)+0.001,2) for z in d9]
 
file = open("R_K_ter2TPM.txt")
d10 = file.read().splitlines() # list
k2 = [math.log(float(z)+0.001,2) for z in d10]
 
file = open("R_K_ter3TPM.txt")
d11 = file.read().splitlines() # list
k3 = [math.log(float(z)+0.001,2) for z in d11]
 
x = [x0, x1, x2, x3, a0, a1, a2, a3, k0, k1, k2, k3]
 
fig = plt.figure(1, figsize =(100, 30))
 
# Create an axes instance
ax = fig.add_subplot(111)
ax.set_xticklabels(['X_99.9', 'X_ter1', 'X_ter2', 'X_ter3', 'A_99.9', 'A_ter1', 'A_ter2', 'A_ter3', 'K_99.9', 'K_ter1', 'K_ter2', 'K_ter3'],fontsize=80)
plt.yticks(fontsize=80)
ax.set_ylabel('log2(TPM)', fontsize=100)
#ax.set_xlabel('Grouping based on Sequential Similarity with Xist, Airn, and Kcnq1ot1', fontsize=13)
 
# Create the boxplot
bp = ax.boxplot(x,showfliers=False,patch_artist=True,boxprops=dict(linewidth=4,facecolor="white", color="black"),
                medianprops=dict(linewidth=4,color="red"),capprops=dict(linewidth=4),
                whiskerprops=dict(linewidth=4),flierprops=dict(linewidth=4))
#plt.boxplot(data[:,:3], positions=[1,2,3], notch=True, patch_artist=True,boxprops=dict(facecolor=c, color=c),capprops=dict(color=c),whiskerprops=dict(color=c),flierprops=dict(color=c, markeredgecolor=c),medianprops=dict(color=c))
 
plt.axvline(x=4.5, linewidth=4, color='grey')
plt.axvline(x=8.5, linewidth=4, color='grey')
 
plt.title('Expression Levels of Expressed lncRNA Genes in Mouse TSCs\n',fontsize=100,fontweight='bold')
plt.show()
#plt.savefig('Boxplot_TPM_Expressed_lncRNA.png')
#plt.savefig('Boxplot_TPM_Expressed_lncRNA.pdf')
```

## Summary of findings
Lastly, I made a masterfile, which is the SEEKR results file plus the following information for each expressed lncRNA: spliced or unspliced, GC-contnet, RNA length, TPM value.
```
cat expressed_v_xka.txt | cut -f1,2,3,4 | awk '{if ($1 ~ /unspliced/) {print $1, $2, $3, $4, "unspliced";} else {print $1, $2, $3, $4, "spliced"}}' > expressed_masterfile.txt
cat gencode_lncRNAs_expressed_linear.fa | sort | cut -f2 | awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,""); printf "%.2f%%\n", (gc*100)/(gc+at)}' > gc_masterfile.txt
cat gencode_lncRNAs_expressed_linear.fa | sort | cut -f2 | awk '!/^>/{char+=gsub(/[aAtTnNcCgG]/,""); printf "%d\n", char}' > length_masterfile.txt
cat tsc_rsem_lncRNA_expressed.txt | sort | cut -f6 > TPM_masterfile.txt

sort -o expressed_masterfile.txt expressed_masterfile.txt | less
paste expressed_masterfile.txt gc_masterfile.txt length_masterfile.txt TPM_masterfile.txt | tr ' ' '\t' > masterfile.txt
sed  -i '1i gene_info\t>Xist\t>Airn\t>Kcnq1ot1\tspliced/unspliced\tGC-content\tgene_length\tTPM'  masterfile.txt
```
