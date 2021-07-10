# Python

# 3. Classify all the peaks into 3 classes (C, NC, other) and calculate percentage of each

# A list of labels for each peak (based on the order in the RIP file) will be printed to "peak_labels.txt".
# The percentage of each class will be printed to the terminal.

c_number = 0;
nc_number = 0;
all_number = 0;

# 3a. Label peaks that overlaps with C
c_intersect = open("c_intersect.txt", "r")
c_lines = c_intersect.readlines()
c_intersect.close()

labels = []

for line in c_lines:
    arr = line.split("\t")
    if arr[12] == '0\n':
        labels.append("")
    else:
        labels.append("protein-coding")
        c_number += 1

# 3b. label peaks that overlaps with NC (but not C)
nc_intersect = open("nc_intersect.txt", "r")
nc_lines = nc_intersect.readlines()
nc_intersect.close()

for line in nc_lines:
    arr = line.split("\t")
    index = int(arr[3].strip("Peak_")) - 1
    if arr[12] != '0\n':
        if labels[index] == '':
            labels[index] = "lncRNA"
            nc_number += 1

# 3c. produce .txt of labels
peak_labels = open("peak_labels.txt", "w")

del labels[len(labels) - 1] # delete the last newline ('\n') which has no meaning

for label in labels:
    peak_labels.write(label + "\n")

peak_labels.close()

# 3d. calculate percentage of each class
all_number = len(labels)
c_percentage = float(c_number) / float(all_number)
c_percentage = round(c_percentage, 4) * 100
nc_percentage = float(nc_number) / float(all_number)
nc_percentage = round(nc_percentage, 4) * 100
print('Within ' + str(all_number) + ' total peaks, there are ' + str(c_number) + ' peaks that overlap with protein-coding genes, and ' + str(nc_number) + ' peaks that overlap only with lncRNA genes')
print('Percentage of protein-coding peaks: ' + str(c_percentage) + "%")
print('Percentage of lncRNA peaks: ' + str(nc_percentage) + "%")
print('Percentage of other peaks: ' + str(100.0 - nc_percentage - c_percentage) + "%" )
