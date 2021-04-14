#!/bin/bash


# ./Run-RNA-workflow-6.0_ncbiblast.sh > Run-RNA-workflow-6.0_ncbiblast.sh.log 2>&1


####################################################################
# Ref_Nicotiana_benthamiana
####################################################################

##cat trinity_genes_magicblast-000001735.4_cds.tsv | awk '$3 != 0' | wc -l

############################################################################################################################################################################################################
# Run all scripts 

####################################################################
# Ref_Nicotiana_benthamiana ----vs---  # iTAK , Nicotiana_benthamiana, transcription factors
####################################################################

#cat ITAK_Tobacco-transcription_factor.txt | awk -F'\t' '{ print ">"$1"|"$2"\n"$4 }' > ITAKv18.12_Ext-Nbe_cds.fa

NCBIBLAST="/bioinfo/BLAST/ncbi-blast-2.9.0+/bin"
REF_Q="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/SL-2400_5.0_Annotation/TF_DB/"
query="/home/gala0002/proj/proj_Ramesh/Ref_Nicotiana_benthamiana/Supplemental_dataset_3_NbC_CDS_clean.fa"
outfl="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/SL-2400_5.0_Annotation/NbC-cds_blastn-max1_ITAKv18.12_Ext-NbC-cds.tsv"
workdir=/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/SL-2400_5.0_Annotation/TF_DB/
cd $workdir

#$NCBIBLAST/makeblastdb -in $REF_Q/ITAKv18.12_Ext-Nbe_cds.fa -out $REF_Q/ITAKv18.12_Ext-Nbe_cds -input_type fasta -dbtype nucl -parse_seqids

nice -n 5 $NCBIBLAST/blastn -db $REF_Q/ITAKv18.12_Ext-Nbe_cds \
-query $query -num_threads 30 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-05 -out $outfl


# Run all scripts 

####################################################################
# Ref_Chenopodium_Quinoa ----vs---  # PlantTFDBv5.0_Ext-Cqu_cds.fa, Nicotiana_benthamiana, transcription factors
####################################################################

# PlantTFDBv5.0_Nbe_TF_list.txt, PlantTFDBv5.0_Nbe_pep.fa
#sed 's/ /_/g' PlantTFDBv5.0_Nbe_pep.fa > PlantTFDBv5.0_Nbe_pep_pars.fa 
#cat PlantTFDBv5.0_Nbe_pep.fa | sed 's/ /_/g' | sed 's/_Nicotiana_benthamiana|/_/g' | sed 's/|/ /g' > PlantTFDBv5.0_Nbe_pep_pars.fa 

#$BLAST/makeblastdb -in $REF/GCF_000001735.4_TAIR10.1_cds_from_genomic_clean.fa -out $REF/GCF_000001735.4_TAIR10.1_cds_from_genomic_clean -input_type fasta -dbtype nucl -parse_seqids

NCBIBLAST="/bioinfo/BLAST/ncbi-blast-2.9.0+/bin"
REF_Q="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/SL-2400_5.0_Annotation/TF_DB/"
query="/home/gala0002/proj/proj_Ramesh/Ref_Nicotiana_benthamiana/Supplemental_dataset_4_NbC_AA.fasta"
outfl="/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/SL-2400_5.0_Annotation/NbC-prot_blastp-max1_PlantTFDBv5.0_Ext-NbC-pep.tsv"
workdir=/home/gala0002/proj/proj_Ramesh/RNA-seq_2020/SL-2400_5.0_Annotation/TF_DB/
cd $workdir

#$NCBIBLAST/makeblastdb -in $REF_Q/PlantTFDBv5.0_Nbe_pep_pars.fa -out $REF_Q/PlantTFDBv5.0_Nbe_pep_pars -input_type fasta -dbtype prot -parse_seqids

nice -n 5 $NCBIBLAST/blastp -db $REF_Q/PlantTFDBv5.0_Nbe_pep_pars \
-query $query -num_threads 30 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-05 -out $outfl


echo "done script..."

