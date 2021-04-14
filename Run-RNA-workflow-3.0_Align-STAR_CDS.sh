#!/bin/bash


#Rna-seq using STAR
#./Run-RNA-workflow-3.0_Align-STAR_CDS.sh > Run-RNA-workflow-3.0_Align-STAR_CDS.sh.log 2>&1


echo "run script for rna-seq-analysis"
##########################################
#STAR-Aligner
##########################################

#make REF index
<<COMMENT
ref=/home/gala0002/proj/proj_Ramesh/Ref_Nicotiana_benthamiana/
nice -n 5 STAR --runThreadN 35 \
--runMode genomeGenerate \
--genomeDir ${ref} \
--genomeFastaFiles ${ref}/Supplemental_dataset_3_NbC_CDS.fasta \
--genomeSAindexNbases 8 \
--limitGenomeGenerateRAM 300000000000

COMMENT

ulimit -n 65535
ulimit -c unlimited

#Make STAR alignment
######################
work_dir=/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/
ref=/home/gala0002/proj/proj_Ramesh/Ref_Nicotiana_benthamiana/

mkdir -p ${work_dir}SL-2400_3.0_Align-STAR_CDS/
out_dir=${work_dir}SL-2400_3.0_Align-STAR_CDS/

######################################################
#####---STAR with 2 Lanes:  QE-1467 and QE-1467_180206 
######################################################

cd ${work_dir}SL-2400_1.1_sort-trim/

#for nbr in Sample_484-10-1/; do
#for nbr in Sample*/; do

for nbr in `ls /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/SL-2400_1.1_sort-trim/`

do
echo "Sample_DIR $nbr"/

nbr1=$(echo $nbr | sed 's=/[^/]*$==;s/\.$//')
nbr2=$(echo $nbr | cut -d"_" -f1)

echo $nbr1
echo $nbr2

mkdir -p ${out_dir}${nbr2}/

temp_dir=${out_dir}${nbr2}/
cd ${temp_dir}

#All-reads
echo "Processing sample: ${nbr1}"

#5.2.2) Make START alignment
#/SortedByCoordinate /Unsorted

nice -n 5 STAR --genomeDir ${ref} \
--runMode alignReads \
--twopassMode Basic \
--alignSJDBoverhangMin 10 \
--readFilesIn ${work_dir}SL-2400_1.1_sort-trim/${nbr}"/"${nbr1}-sortmerna-trimmomatic_1.fq.gz ${work_dir}SL-2400_1.1_sort-trim/${nbr}"/"${nbr1}-sortmerna-trimmomatic_2.fq.gz \
--runThreadN 40 \
--outFileNamePrefix ${temp_dir}${nbr2}-sortmerna-trimmomatic-STAR \
--outSAMtype BAM SortedByCoordinate \
--readFilesCommand zcat
done

echo "Done for independent samples..."


echo "Script done...."

