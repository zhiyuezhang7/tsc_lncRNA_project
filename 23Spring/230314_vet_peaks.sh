#!/bin/bash

#SBATCH -J vet_%j
#SBATCH -p general
#SBATCH --mem=50G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=0:30:00
#SBATCH -o vet_%j.out
#SBATCH -e vet_%j.err
#SBATCH --mail-type ALL 
#SBATCH --mail-user zhiyue@live.unc.edu

## IMPORTANT: In the same directory where you run 230301_3pr_classify.sh, where all the "$1_classify" folders are

## Example command: sbatch 230314_vet_peaks.sh $1
    # $1 = protein name

## Input needed :
    # "$1_classify" folder
    # python correlation script: /proj/calabrlb/users/Zhiyue/23_sp/code/230313_corr.py


module load samtools

cd $1_classify

## print header
echo $1 > $1_vet_peaks.txt # output file


### 1. how many reads per replicate

# get individual rip filenames
rip_filenames=$(ls -m $1/*.fastq | sed -r 's/\s+//g' | tr -d '\n')
# get replicate number
rip_file_num=$(ls $1 | wc -l)

## count reads per replicate
IFS=',' read -r -a rip_arr <<< "$rip_filenames" # make array with rip file names

for i in ${!rip_arr[@]}; do # loop on indices of array

    # count reads
    fastq_lin_num=$(wc -l ${rip_arr[$i]} | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
    (( reads = fastq_lin_num / 4 ))

    j=$(($i + 1)) # index+1 so that it's 1 based

    if [ ${j} -eq 1 ] # for the format of output
    then
        if [ ${j} -eq ${rip_file_num} ]
        then
            # first also last; print with \n
            echo $reads >> $1_vet_peaks.txt
        else
            # first not last; print without \n
            echo -n $reads >> $1_vet_peaks.txt
        fi
    else
        if [ ${j} -eq ${rip_file_num} ]
        then
            # last; print with \n
            echo ';'$reads >> $1_vet_peaks.txt
        else
            # middle; print without \n
            echo -n ';'$reads >> $1_vet_peaks.txt
        fi
    fi

done



### 2. how many reads aligned per replicate (MAPQ > 30)

if [ ${rip_file_num} -gt 1 ]
then

    # there is more than 1 replicate
    for i in $(seq 1 $rip_file_num); do # for each replicate ripi
    
        # count reads aligned (MAPQ > 30)
            # -c: count alignments
            # -F 260: exclude unmapped and secondary alignments
                # A read may map ambiguously to multiple locations, e.g. due to repeats. The primary alignment is chosen to have the highest alignment score.
        reads=$(samtools view -c -F 260 $1_rip${i}_Aligned_filteredsq30.out.sam)

        if [ ${i} -eq 1 ] # for the format of output
        then
            # first; print without \n
            echo -n $reads >> $1_vet_peaks.txt
        else
            if [ ${i} -eq ${rip_file_num} ]
            then
                # last; print with \n
                echo ';'$reads >> $1_vet_peaks.txt
            else
                # middle; print without \n
                echo -n ';'$reads >> $1_vet_peaks.txt
            fi
        fi

    done

else

    # only one replicate
    reads=$(samtools view -c -F 260 $1_Aligned_filteredsq30.out.sam)
    echo $reads >> $1_vet_peaks.txt

fi



### 3. how many total peaks
wc -l < $1_saf.txt >> $1_vet_peaks.txt



### 4. how many >2 igg peaks (in >= 2 replicates)

if [ ${rip_file_num} -gt 1 ]
then

    # there is more than 1 replicate
    # add count = 0 to last column indicating #(replicates >2igg)
    cat $1_fc_rpm.txt | sed 1d | awk -v OFS="\t" '{print $0, 0}' > $1_ripi_2igg.txt

    # When rip_file_num > 1, 
            # column number of $1_fc_rpm.txt = 6 (peak info) + 2 + rip_file_num (counts) + 2 + rip_file_num (rpm)
    ((igg_col = 10 + rip_file_num)) # column number of igg_rpm

    for i in $(seq 1 $rip_file_num) # for each replicate ripi
    do
        ((rip_col = 10 + rip_file_num + i)) # column number of ripi_rpm
        # (ripi_rpm > 2*igg_rpm) -> if TRUE, increment count in the last column
            # awk -i inplace: edit file in place; don't use cat, must put file name at the end!
        awk -i inplace -v OFS="\t" -v c1=$rip_col -v c2=$igg_col  '{if($c1 > 2 * $c2) {$NF=$NF+1;print} else {print $0}}' $1_ripi_2igg.txt
    done
    
    # find peaks with >2igg in >= replicates
        # the last column of $1_ripi_2igg.txt indicates count of replicates >2igg for a certain peak
    cat $1_ripi_2igg.txt | awk -v OFS="\t" '{if($NF > 1) {print $0}}' > $1_fc_rpm_2igg_ind.txt
    # count peaks with >2igg in >= replicates
    cat $1_fc_rpm_2igg_ind.txt | wc -l >> $1_vet_peaks.txt

else

    # there is only one replicate, print NA
    echo NA >> $1_vet_peaks.txt

fi



### 5. how many >4 igg peaks (in >= 2 replicates)

if [ ${rip_file_num} -gt 1 ]
then

    # there is more than 1 replicate
    # add count = 0 to last column indicating #(replicates >2igg)
    cat $1_fc_rpm.txt | sed 1d | awk -v OFS="\t" '{print $0, 0}' > $1_ripi_4igg.txt

    # When rip_file_num > 1, 
            # column number of $1_fc_rpm.txt = 6 (peak info) + 2 + rip_file_num (counts) + 2 + rip_file_num (rpm)
    ((igg_col = 10 + rip_file_num)) # column number of igg_rpm

    for i in $(seq 1 $rip_file_num) # for each replicate ripi
    do
        ((rip_col = 10 + rip_file_num + i)) # column number of ripi_rpm
        # (ripi_rpm > 4*igg_rpm) -> if TRUE, increment count in the last column
            # awk -i inplace: edit file in place; don't use cat, must put file name at the end!
        awk -i inplace -v OFS="\t" -v c1=$rip_col -v c2=$igg_col  '{if($c1 > 4 * $c2) {$NF=$NF+1;print} else {print $0}}' $1_ripi_4igg.txt
    done
    
    # find peaks with >4igg in >= replicates
        # the last column of $1_ripi_4igg.txt indicates number of replicates >4igg
    cat $1_ripi_4igg.txt | awk -v OFS="\t" '{if($NF > 1) {print $0}}' > $1_fc_rpm_4igg_ind.txt
    # count peaks with >4igg in >= replicates
    cat $1_fc_rpm_4igg_ind.txt | wc -l >> $1_vet_peaks.txt

else

    # there is only one replicate, print NA
    echo NA >> $1_vet_peaks.txt

fi



### 6. correlation between reads from individual replicates in >2 igg peaks

## In $1_fc_rpm_2igg_ind.txt, get corr($1_rip1_rpm, $1_rip2_rpm) for every pair
    # Output order: r1,2; r1,3; ... r1,n; r2,3; r2,4;... r2,n; ... r(n-1),n

if [ $rip_file_num -gt 1 ]
then
    # there is more than 1 replicate
    
    for i in $(seq 1 $(($rip_file_num - 1))) # for each ripi with i < number(rip)
    do
        ((ripi_col = 10 + rip_file_num + i)) # column number of ripi_rpm

        for j in $(seq $(($i + 1)) $rip_file_num ) # for each ripj with i < j <= number(rip)
        do

            ((ripj_col = 10 + rip_file_num + j)) # column number of ripj_rpm

            ## correlation between ripi_rpm and ripj_rpm
            # get 2 lists
            cat $1_fc_rpm_2igg_ind.txt | cut -f${ripi_col} > $1_ripi_rpm.txt
            cat $1_fc_rpm_2igg_ind.txt | cut -f${ripj_col} > $1_ripj_rpm.txt
            # corr in python:
            pearson_r=$(python /proj/calabrlb/users/Zhiyue/23_sp/code/230313_corr.py $1_ripi_rpm.txt $1_ripj_rpm.txt)

            ## print i,j correlation to output
            if [ ${i} -eq 1 -a ${j} -eq 2 ] # for the format of output
            then
                if [ $rip_file_num -eq 2 ]
                then
                    # first also last; print with \n
                    echo $pearson_r >> $1_vet_peaks.txt
                else
                    # first not last; print without \n
                    echo -n $pearson_r >> $1_vet_peaks.txt
                fi
            else
                rip_file_num_minus1=$(($rip_file_num - 1))
                if [ ${i} -eq ${rip_file_num_minus1} -a ${j} -eq ${rip_file_num} ]
                then
                    # last; print with \n
                    echo ';'$pearson_r >> $1_vet_peaks.txt
                else
                    # middle; print without \n
                    echo -n ';'$pearson_r >> $1_vet_peaks.txt
                fi
            fi

        done
    done

else

    # there is only one replicate, print NA
    echo NA >> $1_vet_peaks.txt

fi



### 7. correlation between reads from individual replicates in >4 igg peaks

## In $1_fc_rpm_4igg_ind.txt, get corr($1_rip1_rpm, $1_rip2_rpm) for every pair
    # Output order: r1,2; r1,3; ... r1,n; r2,3; r2,4;... r2,n; ... r(n-1),n

if [ $rip_file_num -gt 1 ]
then
    # there is more than 1 replicate
    
    for i in $(seq 1 $(($rip_file_num - 1))) # for each ripi with i < number(rip)
    do
        ((ripi_col = 10 + rip_file_num + i)) # column number of ripi_rpm

        for j in $(seq $(($i + 1)) $rip_file_num ) # for each ripj with i < j <= number(rip)
        do

            ((ripj_col = 10 + rip_file_num + j)) # column number of ripj_rpm

            ## correlation between ripi_rpm and ripj_rpm
            # get 2 lists
            cat $1_fc_rpm_4igg_ind.txt | cut -f${ripi_col} > $1_ripi_rpm.txt
            cat $1_fc_rpm_4igg_ind.txt | cut -f${ripj_col} > $1_ripj_rpm.txt
            # corr in python:
            pearson_r=$(python /proj/calabrlb/users/Zhiyue/23_sp/code/230313_corr.py $1_ripi_rpm.txt $1_ripj_rpm.txt)

            ## print i,j correlation to output
            if [ ${i} -eq 1 -a ${j} -eq 2 ] # for the format of output
            then
                if [ ${rip_file_num} -eq 2 ]
                then
                    # first also last; print with \n
                    echo $pearson_r >> $1_vet_peaks.txt
                else
                    # first not last; print without \n
                    echo -n $pearson_r >> $1_vet_peaks.txt
                fi
            else
                rip_file_num_minus1=$(($rip_file_num - 1))
                if [ ${i} -eq ${rip_file_num_minus1} -a ${j} -eq ${rip_file_num} ]
                then
                    # last; print with \n
                    echo ';'$pearson_r >> $1_vet_peaks.txt
                else
                    # middle; print without \n
                    echo -n ';'$pearson_r >> $1_vet_peaks.txt
                fi
            fi

        done
    done

else

    # there is only one replicate, print NA
    echo NA >> $1_vet_peaks.txt

fi