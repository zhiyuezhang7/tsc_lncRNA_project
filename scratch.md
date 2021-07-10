# Classify RIP signals
From a RIP file, identify if each peak overlaps with lncRNA genes, protein-coding genes (If both, will label as protein-coding), or other DNA regions; calculate the proportion of each class.
I will label lncRNA genes as NC genes, and protein-coding genes as C genes.
The RIP file used here is named `hnrnpk_pirhana_rpkm_ercc_normalization.txt`.

## 1. Convert file format (Python)
Convert mouse NC genes, mouse C genes, and the RIP file into .bed format in preparation of using [bedtools](https://bedtools.readthedocs.io/en/latest/).
The NC genes and C genes are derived from [GENCODE](https://www.gencodegenes.org/)'s basic .gtf file: 'gencode.vM25.basic.annotation.complete.ERCC.fa`.

### 1a. Extract NC genes and C genes from .gtf file
```
import re
file = open("gencode.vM25.basic.annotation.gtf", "r")
lines = file.readlines()
file.close()

nc_pattern = ['transcript_type "3prime_overlapping_ncRNA"', 'transcript_type "antisense"', 
 'transcript_type "bidirectional_promoter_lncRNA"', 'transcript_type "lincRNA"', 
 'transcript_type "macro_lncRNA"', 'transcript_type "processed_transcript"', 
 'transcript_type "sense_intronic"', 'transcript_type "sense_overlapping"']

c_pattern = 'transcript_type "protein_coding"'

nc_gtf = open("nc.gtf", "w")
c_gtf = open("c.gtf", "w")

for line in lines:
    if re.search(r"(?=("+'|'.join(nc_pattern)+r"))", line):
        nc_gtf.write(line)
    else:
        if re.search(c_pattern, line):
            c_gtf.write(line)
    
nc_gtf.close()
c_gtf.close()
```

### 1b. Convert .gtf to .bed
For NC:
```
nc_gtf = open("nc.gtf", "r")
nc_lines = nc_gtf.readlines()
nc_gtf.close()

nc_bed = open("nc.bed", "w")
for line in nc_lines:
    arr = line.split("\t")
    nc_bed.write(arr[0] + '\t' + str(int(arr[3])-1) + '\t' + arr[4] + '\t' + '.' + '\t' 
    + '.' + '\t' + arr[6] + '\n')
nc_bed.close()
```

