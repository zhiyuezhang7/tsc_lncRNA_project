import numpy as np

## 1. convert estimated counts to FPKM for all_est_cts and ERCC_est_cts

# calculate read counts for 6 groups
read_counts_file = open("read_counts.txt", "r")
read_counts = read_counts_file.readlines()
read_counts_file.close()

def est_cts_to_fpkm(est_cts_file_name):
    cts_file = open(est_cts_file_name, "r")
    cts = cts_file.readlines()
    cts_file.close()

    input_FPKM_unft = []
    input_FPKM_ft = []
    rip_FPKM_unft = []
    rip_FPKM_ft = []
    igg_FPKM_unft = []
    igg_FPKM_ft = []
    all_group_FPKM = [input_FPKM_unft, input_FPKM_ft, rip_FPKM_unft, rip_FPKM_ft, igg_FPKM_unft, igg_FPKM_ft]

    i = 0
    for sample in read_counts:
        sample_arr = sample.split("\t")
        print(sample_arr[0])
        reads = int(sample_arr[1])
        est_cts_col_unft = 3 + 2 * i
        est_cts_col_ft = 3 + 2 * i + 1

        for line in cts:
            arr = line.split("\t")
            if arr[0] != "gene_ID":
                # FPKM = Estimated_counts / (eff_length / 1000) / (read_counts / 1,000,000)
                FPKM_unft = float(arr[est_cts_col_unft]) * 1000 / float(arr[2])
                FPKM_unft = FPKM_unft * 1000000 / reads
                all_group_FPKM[2 * i].append(FPKM_unft)

                FPKM_ft = float(arr[est_cts_col_ft]) * 1000 / float(arr[2])
                FPKM_ft = FPKM_ft * 1000000 / reads
                all_group_FPKM[2 * i + 1].append(FPKM_ft)
        i+=1

    return all_group_FPKM

# ERCC FPKM
ERCC_FPKM = est_cts_to_fpkm("est_cts_ercc.txt")
    # delete filtered groups
del ERCC_FPKM[5] 
del ERCC_FPKM[3]
del ERCC_FPKM[1]

# All FPKM
all_FPKM = est_cts_to_fpkm("est_cts_all.txt")

## 2. scale rip_ft(FPKM) and igg_ft(FPKM) to ERCC

# find upper quartile of non-zero FPKMs in ERCC unft groups
    # note: delete printing statements after checking

def find_upper_quartile(FPKM_groups, output_header):
    q3_arr = []
    q3_out = open("fpkm_q3_" + output_header + ".txt", "w")
    for FPKMarr in FPKM_groups: # for each array in the 2D array
        non_zero_FPKM = []
        for FPKM in FPKMarr:
            if FPKM > 0 :
                non_zero_FPKM.append(FPKM)
        q3 = np.percentile(non_zero_FPKM, 75)
        q3_arr.append(q3)
        q3_out.write(str(q3) + "\t")
    q3_out.close() # output upper quartiles in each columns
    return q3_arr

# q3 of ERCC FPKM
ERCC_FPKM_q3 = find_upper_quartile(ERCC_FPKM, "ercc")

# use ERCC to generate scaling factor
scaling_factor_rip = ERCC_FPKM_q3[0] / ERCC_FPKM_q3[1] # input_unft(Q3) / rip_unft(Q3)
scaling_factor_igg = ERCC_FPKM_q3[0] / ERCC_FPKM_q3[2]  # input_unft(Q3) / igg_unft(Q3)

# multiply all filtered rip and igg FPKM by ERCC scaling factor
    # all_rip_ft(FPKM) * ( ercc_input_unft(Q3) / ercc_rip_unft(Q3) )
    # all_igg_ft(FPKM) * ( ercc_input_unft(Q3) / ercc_igg_unft(Q3) )

# 3. Subtract background noise (igg signal) for each transcript
    # all_rip_ft(FPKM) - all_igg_ft(FPKM)

# 4. Output (all ft except for input)
    # fpkm_input.txt
    # fpkm_rip.txt
    # fpkm_igg.txt
    # ercc_rip.txt
    # ercc_igg.txt
    # rip_less_igg.txt
    # rip_over_input.txt

output0 = open("fpkm_input.txt", "w")
output1 = open("fpkm_rip.txt", "w")
output2 = open("fpkm_igg.txt", "w")
output3 = open("ercc_rip.txt", "w")
output4 = open("ercc_igg.txt", "w")
output5 = open("rip_less_igg.txt", "w")
output6 = open("rip_over_input.txt", "w")

output0.write("fpkm_input_unfiltered\n")
output1.write("fpkm_rip\n")
output2.write("fpkm_igg\n")
output3.write("ercc_rip\n")
output4.write("ercc_igg\n")
output5.write("rip_less_igg\n")
output6.write("rip_over_input\n")

all_FPKM_input_unft = all_FPKM[0]
all_FPKM_rip_ft = all_FPKM[3]
all_FPKM_igg_ft = all_FPKM[5]
all_FPKM_rip_ft_scaled = []
all_FPKM_igg_ft_scaled = []

for z in all_FPKM_input_unft:
    output0.write(str(z) + "\n")

for x in all_FPKM_rip_ft:
    output1.write(str(x) + "\n")

    x_scaled = x * scaling_factor_rip
    all_FPKM_rip_ft_scaled.append(x_scaled)

    output3.write(str(x_scaled) + "\n")

for y in all_FPKM_igg_ft:
    output2.write(str(y) + "\n")

    y_scaled = y * scaling_factor_igg
    all_FPKM_igg_ft_scaled.append(y_scaled)

    output4.write(str(y_scaled) + "\n")

i = 0
while i < len(all_FPKM_rip_ft):
    FPKM_igg_norm = all_FPKM_rip_ft_scaled[i] - all_FPKM_igg_ft_scaled[i]
    output5.write(str(FPKM_igg_norm) + "\n")

    if all_FPKM_input_unft[i] > 0:
        FPKM_input_norm = FPKM_igg_norm / all_FPKM_input_unft[i]
    else:
        FPKM_input_norm = 0
    
    output6.write(str(FPKM_input_norm) + "\n")

    i+=1

output0.close()
output1.close()
output2.close()
output3.close()
output4.close()
output5.close()
output6.close()
