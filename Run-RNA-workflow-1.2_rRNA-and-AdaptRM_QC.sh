#!/bin/bash -l


echo "run script for rna-seq-analysis Cleaning and QC"

#./Run-RNA-workflow-1.2_rRNA-and-AdaptRM_QC.sh > Run-RNA-workflow-1.2_rRNA-and-AdaptRM_QC.sh.log 2>&1

#1. Raw Data QC Assessment

work_dir=/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/
#ref=${work_dir}ta_IWGSC_MIPSv2.2_HighConf_CDS_2014Jul18.fa
#gff=${work_dir}wheatcdsta_IWGSC_MIPSv2.2_HighConf_2014Jul18.gtf

mkdir -p ${work_dir}SL-2400_1.1_sort-trim_Merge_QC/
out_dir=${work_dir}SL-2400_1.1_sort-trim_Merge_QC/

cd ${work_dir}SL-2400_1.1_sort-trim/

#for nbr in Sample_484-10-1/; do
<<COMM
NG-14833_Pas_E_1_lib233711_5747_1/
NG-14833_Pas_E_2_lib233712_5747_1/
NG-14833_Pas_E_3_lib233713_5767_1/
NG-14833_Pas_M_1_lib233714_5767_1/
NG-14833_Pas_M_2_lib233715_5767_1/
NG-14833_Pas_M_3_lib233716_5767_1/
NG-14833_Reg_15dpa_1_lib233717_5767_1/
NG-14833_Reg_25dpa_2_lib233721_5747_6/
NG-14833_Reg_25dpa_3_lib233722_5747_6/
NG-14833_Tit_15dpa_1_lib233705_5747_1/
NG-14833_Tit_15dpa_3_lib233707_5747_1/
NG-14833_Tit_25dpa_1_lib233708_5747_1/
NG-14833_Tit_25dpa_2_lib233709_5747_1/

COMM


for nbr in `ls /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/SL-2400_1.1_sort-trim/`

do
echo "Sample_DIR $nbr"/

#for nbr in NG-14833_Pas_E_2_lib233712_5747_1/ NG-14833_Pas_E_3_lib233713_5767_1/ NG-14833_Pas_M_1_lib233714_5767_1/ NG-14833_Pas_M_2_lib233715_5767_1/ NG-14833_Pas_M_3_lib233716_5767_1/ NG-14833_Reg_15dpa_1_lib233717_5767_1/ NG-14833_Reg_25dpa_2_lib233721_5747_6/ NG-14833_Reg_25dpa_3_lib233722_5747_6/ NG-14833_Tit_15dpa_1_lib233705_5747_1/ NG-14833_Tit_15dpa_3_lib233707_5747_1/ NG-14833_Tit_25dpa_1_lib233708_5747_1/ NG-14833_Tit_25dpa_2_lib233709_5747_1/; do
nbr1=$(echo $nbr | sed 's=/[^/]*$==;s/\.$//')
mkdir -p ${out_dir}${nbr}/

temp_dir=${out_dir}${nbr}/
cd ${temp_dir}
echo "Processing sample: ${nbr1}"

#Command to Merge two lanes
gunzip -c ${work_dir}SL-2400_1.1_sort-trim/${nbr}"/"${nbr1}-sortmerna-trimmomatic_1.fq.gz > ${temp_dir}${nbr1}-L1_1.fq
gunzip -c ${work_dir}SL-2400_1.1_sort-trim/${nbr}"/"${nbr1}-sortmerna-trimmomatic_2.fq.gz > ${temp_dir}${nbr1}-L1_2.fq

#Command to run fastqc
nice -n 5 fastqc -o ${temp_dir} -t 15 --noextract ${temp_dir}${nbr1}-L1_1.fq ${temp_dir}${nbr1}-L1_2.fq

#remove temp files

rm ${temp_dir}${nbr1}-L*_*.fq
#rm ${temp_dir}${nbr1}-Merge_*.fq
done

echo "Done script..."


