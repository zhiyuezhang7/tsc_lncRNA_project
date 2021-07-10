### Unix

4. Count number of exonic coordinates for each peak.
```
cat intersect_exon.txt | cut -f4,13 | awk '{a[$1]+=$2}END{for (i in a) print i, a[i]}' | sort -V | cut -d ' ' -f2  > intersect_merged.txt
sed -i '1i exonic_length' intersect_merged.txt
```
5. Make a new file reporting intersection (and classification, if wanted). The result would be in `hnrnpk_pirhana_rpkm_ercc_normalization_exonic_intersection.csv`.
```
paste hnrnpk_pirhana_rpkm_ercc_normalization.txt intersection_points.txt intersect_merged.txt | tr '\t' ',' > hnrnpk_pirhana_rpkm_ercc_normalization_exonic_intersection.csv
```
Note: if you have run the [code](https://github.com/zhiyuezhang7/tsc_lncRNA_project/tree/main/210701_classify_rip_signals) for classifying rip signals (into protein-coding/lncRNA/other), you can report classification too by running the following commands.

The result would be in `hnrnpk_pirhana_rpkm_ercc_normalization_summary.csv`.
```
sed -i '1i classification' peak_labels.txt
paste hnrnpk_pirhana_rpkm_ercc_normalization.txt peak_labels.txt intersection_points.txt intersect_merged.txt | tr '\t' ',' > hnrnpk_pirhana_rpkm_ercc_normalization_summary.csv
```
### Excel
6. Within Excel, calculate the exonic percentage of each peak and the (normalized) number of exonic signals.

a. Start a new column labeled `exonic_percentage`. Apply exonic_percentage = exonic_length / length (`=S2/F2`) to whole column.

b. Start another column labeled `exonic_signals`. Apply exonic_signals = exonic_percentage * ercchnrnpk (`=T2*O2`) to whole column. ("ercchnrnpk" column represents the normalized RIP signal.)

### Done!
