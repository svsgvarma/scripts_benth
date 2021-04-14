#!/bin/bash


#Rna-seq using STAR
#./Run-RNA-workflow-5.0_DE_CDS.sh > Run-RNA-workflow-5.0_DE_CDS.sh.log 2>&1


echo "run script for rna-seq-analysis"

############################################
# 2.2. Identifying differentially expressed (DE) transcripts
############################################
#Extracting differentially expressed transcripts and generating heatmaps

Trinity="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5/Trinity"
Trinity_path="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5"
work_dir=/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/
ref=/home/gala0002/proj/proj_Ramesh/Ref_Nicotiana_benthamiana/

mkdir -p ${work_dir}SL-2400_3.0_Align-STAR_CDS/
out_dir=${work_dir}SL-2400_3.0_Align-STAR_CDS/

cd $work_dir

mkdir -p ${work_dir}"DESeq2_genes"


############################################
#Build Transcript and Gene Expression Matrices
#https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#express-output
# perl $Trinity_path/util/align_and_estimate_abundance.pl \
# --transcripts Trinity.fasta \
# --seqType fq \
# --samples_file Metadata_CDS_Mock-PMTVWT-delta8k_merge.txt \
# --est_method salmon \
# --trinity_mode --prep_reference \
# --output_dir outdir_estimate-ab > estimate-ab.log 2>&1


############################################
#Run the DE analysis at the gene level

#DESeq2
$Trinity_path/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix Express_CDS_Mock-PMTVWT-delta8k.gene.counts.matrix \
--samples_file Metadata_CDS_Mock-PMTVWT-delta8k_merge.txt \
--method DESeq2 \
--output DESeq2_genes > Trinity_DGE-run.log 2>&1

#Extracting differentially expressed transcripts and generating heatmaps
#Extract those differentially expressed (DE) transcripts that are at least 4-fold (C is set to 2^(2) ) differentially expressed at a significance of <= 0.001 (-P 1e-3) in any of the pairwise sample comparisons
cd DESeq2_genes/
nice -n 5 $Trinity_path/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix ../Express_CDS_Mock-PMTVWT-delta8k.gene.counts.matrix \
--samples ../Metadata_CDS_Mock-PMTVWT-delta8k_merge.txt -P 5e-2 -C 1.0 > Trinity_DGE-analyze.log 2>&1

############################################

echo "Script done...."