#!/bin/bash


#Pars DEG

echo "run script "

###############################################
#all samples Uniq
#########

workdir="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/DESeq2_genes/"

python Run-RNA-workflow-7.0_Pars-annotation-DEG_CDS-DEseq2.py Express_counts_all_Mock-PMTVWT-delta8k.matrix.Delta8K_vs_QMock.DESeq2.DE_results.P5e-2_C1.0.DE.subset all_DESeq2.DE_results.P5e-2_C1.0_Delta8K-vs-Mock_anno.tsv $workdir
python Run-RNA-workflow-7.0_Pars-annotation-DEG_CDS-DEseq2.py Express_counts_all_Mock-PMTVWT-delta8k.matrix.PMTVWT_vs_QMock.DESeq2.DE_results.P5e-2_C1.0.DE.subset all_DESeq2.DE_results.P5e-2_C1.0_PMTVWT-vs-Mock_anno.tsv $workdir

python Run-RNA-workflow-7.0_Pars-annotation-DEG_CDS-DEseq2.py Express_counts_all_Mock-PMTVWT-delta8k.matrix.Delta8K_vs_QMock.DESeq2.DE_results all_DESeq2.DE_results.Delta8K-vs-Mock_anno.tsv $workdir
python Run-RNA-workflow-7.0_Pars-annotation-DEG_CDS-DEseq2.py Express_counts_all_Mock-PMTVWT-delta8k.matrix.PMTVWT_vs_QMock.DESeq2.DE_results all_DESeq2.DE_results.PMTVWT-vs-Mock_anno.tsv $workdir


#########
# filter DE list 
###
cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_delta8K.DESeq2.DE_results | wc -l
cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_PMTVWT.DESeq2.DE_results | wc -l 

###

cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_delta8K.DESeq2.DE_results | cut -f7,11 | awk '($1 >= 1 || $1 <= -1) && ($2 <=0.01)' | wc -l
cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_PMTVWT.DESeq2.DE_results | cut -f7,11 | awk '($1 >= 1 || $1 <= -1) && ($2 <=0.01)' | wc -l

###
cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_delta8K.DESeq2.DE_results | cut -f7,11 | awk '($1 >= 2 || $1 <= -2) && ($2 <=0.01)' | wc -l
cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_PMTVWT.DESeq2.DE_results | cut -f7,11 | awk '($1 >= 2 || $1 <= -2) && ($2 <=0.01)' | wc -l

#####

###
cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_delta8K.DESeq2.DE_results | cut -f7,11 | awk '($1 >= 1 || $1 <= -1) && ($2 <=0.05)' | wc -l
cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_PMTVWT.DESeq2.DE_results | cut -f7,11 | awk '($1 >= 1 || $1 <= -1) && ($2 <=0.05)' | wc -l

###
cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_delta8K.DESeq2.DE_results | cut -f7,11 | awk '($1 >= 2 || $1 <= -2) && ($2 <=0.05)' | wc -l
cat Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_PMTVWT.DESeq2.DE_results | cut -f7,11 | awk '($1 >= 2 || $1 <= -2) && ($2 <=0.05)' | wc -l


echo "Script done...."



