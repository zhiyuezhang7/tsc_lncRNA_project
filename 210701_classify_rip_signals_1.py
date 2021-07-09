# Python

# Files needed:
# GENCODE's basic .gtf file: 'gencode.vM25.basic.annotation.gtf'
# RIP file: 'hnrnpk_pirhana_rpkm_ercc_normalization.txt'

# Goals:
# From a RIP file, identify if each peak overlaps with lncRNA genes, protein-coding genes (If both, will label as protein-coding), or other DNA regions
# Calculate the proportion of each class.

# I will label lncRNA genes as NC genes, and protein-coding genes as C genes. 

# 1. Format conversion: convert mouse NC genes, mouse C genes, and the RIP file into .bed format in preparation of using bedtools. 

# 1a. extract NC genes and C genes from .gtf file
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

# 1b. convert .gtf to .bed format
# for NC
nc_gtf = open("nc.gtf", "r")
nc_lines = nc_gtf.readlines()
nc_gtf.close()

nc_bed = open("nc.bed", "w")
for line in nc_lines:
    arr = line.split("\t")
    nc_bed.write(arr[0] + '\t' + str(int(arr[3])-1) + '\t' + arr[4] + '\t' + '.' + '\t' 
    + '.' + '\t' + arr[6] + '\n')
nc_bed.close()

# for C
c_gtf = open("c.gtf", "r")
c_lines = c_gtf.readlines()
c_gtf.close()

c_bed = open("c.bed", "w")
for line in c_lines:
    arr = line.split("\t")
    c_bed.write(arr[0] + '\t' + str(int(arr[3])-1) + '\t' + arr[4] + '\t' + '.' + '\t' 
    + '.' + '\t' + arr[6] + '\n')
c_bed.close()

# 1c. convert rip data to .bed format
rip = open("hnrnpk_pirhana_rpkm_ercc_normalization.txt", "r")
rip_lines = rip.readlines()
del rip_lines[0]
rip.close()

rip_bed = open("rip.bed", "w")
for line in rip_lines:
    arr = line.split("\t")
    rip_bed.write(arr[1].strip('"') + "\t" + str(int(arr[2])-1) + "\t" + arr[3] + "\t" + arr[0].strip('"') + "\t" + "." + "\t" + arr[4].strip('"') + "\n")
rip_bed.close()
