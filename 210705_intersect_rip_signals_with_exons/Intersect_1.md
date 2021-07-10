### Unix

Files needed:
GENCODE's basic .gtf file: 'gencode.vM25.basic.annotation.gtf'
RIP file: 'hnrnpk_pirhana_rpkm_ercc_normalization.txt'

Goals: Report the percentage of exonic intersection within each peak, and the genomic coordinates of such intersection (if any).

1. Format conversion: convert mouse exons, and the RIP file into .bed format in preparation of using bedtools. 
```
cat gencode.vM25.basic.annotation.gtf | awk '{if ($3 ~ /exon/) {print $0}}' > exon.gtf
cat exon.gtf | awk  '{print $1,$4-1,$5,".",".",$7}' | tr ' ' '\t'  > exon.bed
cat hnrnpk_pirhana_rpkm_ercc_normalization.txt | sed 1d | awk  '{print $2"\t"$3-1"\t"$4"\t"$1"\t"".""\t"$5}' | sed 's/"//g' > rip.bed
```
2. Merge exon intervals and report intersection with rip data.
```
module load bedtools
sortBed -i exon.bed > exon_sorted.bed
bedtools merge -i exon_sorted.bed -s -c 6 -o distinct | awk  '{print $1,$2,$3,".",".",$4}' | tr ' ' '\t'  > exon_merged.bed
bedtools intersect -a rip.bed -b exon_merged.bed -s -wao > intersect_exon.txt
```
