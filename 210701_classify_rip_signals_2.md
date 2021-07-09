### Unix

2. Use bedtools to merge gene intervals, and intersect RIP data with such intervals
2a. merge NC and C intervals
```
module load bedtools
sortBed -i nc.bed > nc_sorted.bed
bedtools merge -i nc_sorted.bed -s -c 6 -o distinct | awk  '{print $1,$2,$3,".",".",$4}' | tr ' ' '\t'  > nc_merged.bed
	-s (strandness) -c 6 -o distinct (report strand)

sortBed -i c.bed > c_sorted.bed
bedtools merge -i c_sorted.bed -s -c 6 -o distinct | awk  '{print $1,$2,$3,".",".",$4}' | tr ' ' '\t'  > c_merged.bed
```

2b. report intersection with RIP data
```
bedtools intersect -a rip.bed -b c_merged.bed -s -wao > c_intersect.txt
bedtools intersect -a rip.bed -b nc_merged.bed -s -wao > nc_intersect.txt
```
