#!/bin/bash


# ./00_04_Run-all_ncbiblast.sh > 00_04_Run-all_ncbiblast.sh.log 2>&1


# Run all scripts 
cd /data/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v3/

cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_delta8K.DESeq2.DE_results_anno.tsv | cut -f1 > DESeq2.DE_results.all_Mk-v-Delta.txt
cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_PMTVWT.DESeq2.DE_results_anno.tsv | cut -f1 > DESeq2.DE_results.all_Mk-v-PMTVWT.txt

grep -w -A 1 -f DESeq2.DE_results.all_Mk-v-Delta.txt /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/Ref_Nicotiana_benthamiana/Supplemental_dataset_4_NbC_AA.fasta | sed 's/--//g' > /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/all_DESeq2.DE_results.all_Mock_vs_delta8K_AA.fa
grep -w -A 1 -f DESeq2.DE_results.all_Mk-v-PMTVWT.txt /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/Ref_Nicotiana_benthamiana/Supplemental_dataset_4_NbC_AA.fasta | sed 's/--//g' > /home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/all_DESeq2.DE_results.all_Mock_vs_PMTVWT_AA.fa


####################################################################
# Ref_N.b ----vs---  # Ref_N.t
####################################################################
NCBIBLAST="/bioinfo/BLAST/ncbi-blast-2.9.0+/bin"

workdir=/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/Ref_Ntab-TN90/
#REF_Q_Delta8K="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_AA.fa"
#REF_Q_PMTVWT="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_AA.fa"

REF_Q_Delta8K="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/all_DESeq2.DE_results.all_Mock_vs_delta8K_AA.fa"
REF_Q_PMTVWT="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/all_DESeq2.DE_results.all_Mock_vs_PMTVWT_AA.fa"

REF_A="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/KEGG_KofamKOALA/Ref_Ntab-TN90/GCF_000715135.1_Ntab-TN90_protein.faa"
cd $workdir

#$NCBIBLAST/makeblastdb -in $REF_A -out $workdir/GCF_000715135.1_Ntab-TN90_protein -input_type fasta -dbtype prot -parse_seqids

nice -n 5 $NCBIBLAST/blastp -db $workdir/GCF_000715135.1_Ntab-TN90_protein \
-query $REF_Q_Delta8K -num_threads 40 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-05 -out $workdir/N.benQ-Delta8K_vs_N.tab_blastp-max1_proteins.tsv


nice -n 5 $NCBIBLAST/blastp -db $workdir/GCF_000715135.1_Ntab-TN90_protein \
-query $REF_Q_PMTVWT -num_threads 40 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-05 -out $workdir/N.benQ-PMTVWT_vs_N.tab_blastp-max1_proteins.tsv


#######

cat N.benQ-Delta8K_vs_N.tab_blastp-max1_proteins.tsv | cut -f2 > N.benQ-Delta8K_vs_N.tab_blastp-max1_proteins_XP.tsv
cat N.benQ-PMTVWT_vs_N.tab_blastp-max1_proteins.tsv | cut -f2 > N.benQ-PMTVWT_vs_N.tab_blastp-max1_proteins_XP.tsv


rm N.benQ-Delta8K_vs_N.tab_blastp-max1_proteins_XP-GID.tsv

while read i
do
#echo "$i"

LANG=C grep -m 1 -w "$i" GCF_000715135.1_Ntab-TN90_genomic.gff | awk '{split($9,a,";"); print a[3]}' >> N.benQ-Delta8K_vs_N.tab_blastp-max1_proteins_XP-GID.tsv

done < N.benQ-Delta8K_vs_N.tab_blastp-max1_proteins_XP.tsv


#######
rm N.benQ-Delta8K_vs_N.tab_blastp-max1_proteins_XP-GID.tsv

while read i
do
#echo "$i"

LANG=C grep -m 1 -w "$i" GCF_000715135.1_Ntab-TN90_genomic.gff | awk '{split($9,a,";"); print a[3]}' >> N.benQ-PMTVWT_vs_N.tab_blastp-max1_proteins_XP-GID.tsv

done < N.benQ-PMTVWT_vs_N.tab_blastp-max1_proteins_XP.tsv


#######


echo "done script..."

