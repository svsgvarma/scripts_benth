#!/bin/bash


#Separate reads
##############################################################################################################

#./Run-RNA-workflow-4.0_Features-count_eXpress.sh > Run-RNA-workflow-4.0_Features-count_eXpress.sh.log 2>&1

echo "run script for rna-seq-analysis"

######################
work_dir=/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/
ref=/home/gala0002/proj/proj_Ramesh/Ref_Nicotiana_benthamiana/Supplemental_dataset_3_NbC_CDS_clean.fa
express="/bioinfo/express-1.5.1-linux_x86_64/express"

mkdir -p ${work_dir}SL-2400_4.0_express_trained/
out_dir=${work_dir}SL-2400_4.0_express_trained/

in_dir=${work_dir}SL-2400_3.0_Align-STAR_CDS/


cd ${in_dir}

#for nbr in NG-14833_Pas_E_1_lib233711_5747_1/ NG-14833_Tit_25dpa_3_lib233710_Merge/; do
#for nbr in NG-14833_Pas_E_1_lib233711_5747_1/ NG-14833_Pas_E_2_lib233712_5747_1/ NG-14833_Pas_E_3_lib233713_5767_1/ NG-14833_Pas_M_1_lib233714_5767_1/ NG-14833_Pas_M_2_lib233715_5767_1/ NG-14833_Pas_M_3_lib233716_5767_1/ NG-14833_Reg_15dpa_1_lib233717_5767_1/ NG-14833_Reg_25dpa_2_lib233721_5747_6/ NG-14833_Reg_25dpa_3_lib233722_5747_6/ NG-14833_Tit_15dpa_1_lib233705_5747_1/ NG-14833_Tit_15dpa_3_lib233707_5747_1/ NG-14833_Tit_25dpa_1_lib233708_5747_1/ NG-14833_Tit_25dpa_2_lib233709_5747_1/ NG-14833_Reg_15dpa_2_lib233718_Merge/ NG-14833_Reg_15dpa_3_lib233719_Merge/ NG-14833_Reg_25dpa_1_lib233720_Merge/ NG-14833_Tit_15dpa_2_lib233706_Merge/ NG-14833_Tit_25dpa_3_lib233710_Merge/; do

for nbr in `ls $in_dir`
do

nbr1=$(echo $nbr | sed 's=/[^/]*$==;s/\.$//')
mkdir -p ${out_dir}${nbr}/

temp_dir=${out_dir}${nbr}/
cd ${temp_dir}

#All-reads
echo "Processing sample: ${in_dir}${nbr}/${nbr1}" 

#Samtools and eXpress feature count (FPKM/TPM)

#express [options]* <target_seqs.fasta> <aligned_reads.(sam/bam)>
#nice -n 5 samtools sort -n Sample_484-10-1-sortmerna-trimmomatic-STARAligned.out.bam -o Sample_484-10-1-sortmerna-trimmomatic-STARAligned_sort.bam
#nice -n 5 /bioinfo/express-1.5.1-linux_x86_64/express ${ref} ${workdir}QE-1467-Merge_Align_STAR/Sample_484-10-1/Sample_484-10-1-sortmerna-trimmomatic-STARAligned_sort.bam -o ${express_out}Sample_484-10-1_noSAlign.express.count

#nice -n 5 samtools sort -n ${in_dir}${nbr}/${nbr1}-sortmerna-trimmomatic-STARAligned.sortedByCoord.out.bam -o ${in_dir}${nbr}/${nbr1}-sortmerna-trimmomatic-STARAligned_sort.bam
#--rf-stranded

nice -n 5 ${express} --rf-stranded --no-update-check ${ref} ${in_dir}${nbr}/${nbr1}-sortmerna-trimmomatic-STARAligned_sort.bam -o ${out_dir}${nbr1}

done

echo "Script done all...."
