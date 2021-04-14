#!/bin/bash -l

echo "run script for rna-seq-analysis Cleaning and QC"

#./Run-RNA-workflow-1.1_rRNA-and-AdaptRM.sh > Run-RNA-workflow-1.1_rRNA-and-AdaptRM.sh.log 2>&1
#1. Raw Data rRNA and AdaptRM Assessment

work_dir=/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/

mkdir -p ${work_dir}SL-2400_1.1_sort-trim/
out_dir1=${work_dir}SL-2400_1.1_sort-trim/

################################
#####---Lane: NG-14833
################################

cd ${work_dir}200130_A00605_0102_AHJ75GDRXX/
#for nbr in Sample_484-10-1/; do
#for nbr in NG-14833_*/; do
#nbr=$(echo $nbr | sed 's=/[^/]*$==;s/\.$//')
nbr=$1
mkdir ${out_dir1}${nbr}

#----
SORTMERNADIR=/bioinfo/sortmerna-2.1b
TRIMMOMATIC=/bioinfo/Trimmomatic-0.36

data_dir=${work_dir}200130_A00605_0102_AHJ75GDRXX/
scripts=/bioinfo/sortmerna-2.1b/scripts/

temp_dir=${out_dir1}${nbr}"/"
cd ${temp_dir}

gunzip -c ${data_dir}${nbr}"_R1_001.fastq.gz" > ${temp_dir}${nbr}_R1.fq
gunzip -c ${data_dir}${nbr}"_R2_001.fastq.gz" > ${temp_dir}${nbr}_R2.fq

${scripts}merge-paired-reads.sh ${temp_dir}${nbr}_R1.fq ${temp_dir}${nbr}_R2.fq ${temp_dir}${nbr}-interleaved.fq

nice -n 5 ${SORTMERNADIR}/sortmerna --ref ${SORTMERNADIR}/rRNA_databases/silva-bac-16s-id90.fasta,${SORTMERNADIR}/index/silva-bac-16s-db:\
${SORTMERNADIR}/rRNA_databases/silva-bac-23s-id98.fasta,${SORTMERNADIR}/index/silva-bac-23s-db:\
${SORTMERNADIR}/rRNA_databases/silva-arc-16s-id95.fasta,${SORTMERNADIR}/index/silva-arc-16s-db:\
${SORTMERNADIR}/rRNA_databases/silva-arc-23s-id98.fasta,${SORTMERNADIR}/index/silva-arc-23s-db:\
${SORTMERNADIR}/rRNA_databases/silva-euk-18s-id95.fasta,${SORTMERNADIR}/index/silva-euk-18s-db:\
${SORTMERNADIR}/rRNA_databases/silva-euk-28s-id98.fasta,${SORTMERNADIR}/index/silva-euk-28s:\
${SORTMERNADIR}/rRNA_databases/rfam-5s-database-id98.fasta,${SORTMERNADIR}/index/rfam-5s-db:\
${SORTMERNADIR}/rRNA_databases/rfam-5.8s-database-id98.fasta,${SORTMERNADIR}/index/rfam-5.8s-db \
--reads ${temp_dir}${nbr}-interleaved.fq --aligned ${temp_dir}${nbr}-sortmerna_aligned-rRNA \
--other ${temp_dir}${nbr}-sortmerna --log -v --paired_in -a 15 --fastx

${scripts}unmerge-paired-reads.sh ${temp_dir}${nbr}-sortmerna.fq ${temp_dir}${nbr}-sortmerna_1.fq ${temp_dir}${nbr}-sortmerna_2.fq

find . -name "${nbr}-sortmerna_[1,2].fq" | xargs -P 2 -I {} gzip {} ; cd -

#4. Quality trimming and Adapter removal
nice -n 5 java -jar $TRIMMOMATIC/trimmomatic-0.36.jar PE -threads 15 -phred33 -trimlog ${temp_dir}trimmomatic_${nbr}_logFile \
${temp_dir}${nbr}-sortmerna_1.fq.gz ${temp_dir}${nbr}-sortmerna_2.fq.gz \
${temp_dir}${nbr}-sortmerna-trimmomatic_1.fq.gz ${temp_dir}${nbr}-sortmerna-trimmomatic-unpaired_1.fq.gz \
${temp_dir}${nbr}-sortmerna-trimmomatic_2.fq.gz ${temp_dir}${nbr}-sortmerna-trimmomatic-unpaired_2.fq.gz \
ILLUMINACLIP:$TRIMMOMATIC/adapters/"TruSeq3-PE-2.fa":2:30:10 SLIDINGWINDOW:5:20 MINLEN:20

echo "Done sortmerna, trimommatic...."

rm ${temp_dir}${nbr}_R*.fq
rm ${temp_dir}${nbr}-interleaved.fq
rm ${temp_dir}trimmomatic_${nbr}_logFile
rm ${temp_dir}${nbr}-sortmerna.fq
rm ${temp_dir}${nbr}-sortmerna-trimmomatic-unpaired_*
rm ${temp_dir}${nbr}-sortmerna_*

echo "Done script sample > "${nbr}
#done
echo "script all done...."
