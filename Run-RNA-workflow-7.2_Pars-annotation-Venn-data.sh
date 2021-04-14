#!/bin/bash


#Pars DEG

echo "run script "


##########################################
##########################################
work_dir=/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/SL-2400_6.0_Venn-diagram_P5e-2_FC1/
DESeq2_genes=/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/DESeq2_genes/

cd $work_dir

grep -f venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911.txt -w ${DESeq2_genes}all_DESeq2.DE_results.P5e-2_C1.0_Delta8K-vs-Mock_anno.tsv > venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911_anno.tsv
grep -f venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500.txt -w ${DESeq2_genes}all_DESeq2.DE_results.P5e-2_C1.0_PMTVWT-vs-Mock_anno.tsv> venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500_anno.tsv
#grep -f venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791.txt -w ${DESeq2_genes}Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_PMTVWT.DESeq2.DE_results.P5e-2_C1.0.DE.subset_anno.tsv > venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_anno.tsv

grep -f venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791.txt -w ${DESeq2_genes}all_DESeq2.DE_results.P5e-2_C1.0_*-vs-Mock_anno.tsv > Delta8K_Inter_PMTVWT_791.tsv
cat Delta8K_Inter_PMTVWT_791.tsv | awk '{split($1,a,":"); print a[2],substr($0,length($1)+1)}' | awk -F"\t" '!_[$1]++' > venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_anno.tsv

rm Delta8K_Inter_PMTVWT_791.tsv


############
############
# Inser header for all files 

for f in *_anno.tsv; do
  cat venn_anno.head.tsv $f > $f.new
  mv $f.new $f
done

############

cat *_Inter*_anno.tsv > venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_Merge-REG-PAS-TIT_anno.tsv
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_Merge-REG-PAS-TIT_anno.tsv | awk '{split($1,a,":"); print a[2],substr($0,length($1)+1)}' | awk -F"\t" '!_[$32]++' | awk '$32 != "."' > venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_Merge-REG-PAS-TIT_anno_uniq.tsv

############
# Up and Down regulated 

# Delta
cat venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911_REG_anno_uniq.tsv | awk -F"\t" '{if($7 >= 1 ) print $2,$7 }' | wc -l ; echo "Pas"
cat venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911_REG_anno_uniq.tsv | awk -F"\t" '{if($7 <= -1 ) print $3,$7 }' | wc -l ; echo "Reg"

# PMTVWT
cat venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500_PAS_anno_uniq.tsv | awk -F"\t" '{if($7 >= 1 ) print $2,$7 }' | wc -l ; echo "all"
cat venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500_PAS_anno_uniq.tsv | awk -F"\t" '{if($7 <= -1 ) print $3,$7 }' | wc -l ; echo "PAS"

# TIT
cat venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_TIT_anno_uniq.tsv| awk -F"\t" '{if($7 >= 1 ) print $2,$7 }' | wc -l ; echo "all"
cat venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_TIT_anno_uniq.tsv | awk -F"\t" '{if($7 <= -1 ) print $3,$7 }' | wc -l ; echo "TIT"

# TIT, REG, PAS
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_Inter_RegvTit_PAS-TIT-REG_anno_uniq.tsv| awk -F"\t" '{if($7 >= 1 ) print $2,$7 }' | wc -l ; echo "all"
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_Inter_RegvTit_PAS-TIT-REG_anno_uniq.tsv | awk -F"\t" '{if($7 <= -1 ) print $3,$7 }' | wc -l ; echo "all"


#########################
# count TF database genes 
#PasvReg_1371 PasvTit_645 RegvTit_96

for i in PasvReg_1371 PasvTit_645 RegvTit_96
do
echo $i
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_${i}_anno_uniq_pars.tsv | cut -f47-49 | grep "|" | wc -l
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_${i}_anno_uniq_pars.tsv | cut -f50-53 | grep "|" | wc -l
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_${i}_anno_uniq_pars.tsv | cut -f42 | awk '$1 !="."' | wc -l 

done

####

cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_TIT_anno_uniq.tsv | cut -f47-49 | grep "|" | wc -l
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_TIT_anno_uniq.tsv | cut -f50-53 | grep "|" | wc -l

####

cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_anno_uniq.tsv | cut -f47-49 | grep "|" | wc -l
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_anno_uniq.tsv | cut -f50-53 | grep "|" | wc -l


### Count AT genes numbers in each sample 
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_PAS_anno_uniq.tsv | cut -f42 | awk '$1 !="."' | wc -l 
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_REG_anno_uniq.tsv | cut -f42 | awk '$1 !="."' | wc -l 
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_TIT_anno_uniq.tsv | cut -f42 | awk '$1 !="."' | wc -l 

cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*PAS-TIT-REG_anno_uniq.tsv | cut -f42 | awk '$1 !="."' | wc -l 

##########################
#SNP calling filter on CDS


########
for i in PasvReg_1371 PasvTit_645 RegvTit_96
do
echo $i
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_${i}_anno_uniq_pars.tsv | cut -f1 | grep -v "Quinoa_GeneID" | awk '{split($1,a,"|"); print a[2]}' > CDS_${i}_anno_uniq.tsv
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_${i}_anno_uniq_pars.tsv | cut -f31 | grep -v "Gene_Name" > Genes_${i}_anno_uniq.tsv

done

########
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_PAS_anno_uniq.tsv | cut -f1 | grep -v "Quinoa_GeneID" | awk '{split($1,a,"|"); print a[2]}' > CDS_Pas_anno_uniq.tsv
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_REG_anno_uniq.tsv | cut -f1 | grep -v "Quinoa_GeneID" | awk '{split($1,a,"|"); print a[2]}' > CDS_Reg_anno_uniq.tsv
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_TIT_anno_uniq.tsv | cut -f1 | grep -v "Quinoa_GeneID" | awk '{split($1,a,"|"); print a[2]}' > CDS_Tit_anno_uniq.tsv

cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*PAS-TIT-REG_anno_uniq.tsv | cut -f1 | grep -v "Quinoa_GeneID" | awk '{split($1,a,"|"); print a[2]}' > CDS_Pas-Reg-Tit_anno_uniq.tsv
# total no of SNPs
#for i in Pas Reg Tit; do  echo $i; cat ${SNP_DIR}variants_${i}/filtered_output_sID.vcf | grep -v "^#" | wc -l; done

for i in Pas Reg Tit
do 
echo $i

##### 
#SNP_DIR=/Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_2.1_STAR_Genome-CDS-ass/
SNP_DIR=/home/gala0002/proj/proj_Asa/NG-14833_2.1_STAR_Genome-CDS-ass/

grep -w -f CDS_${i}_anno_uniq.tsv ${SNP_DIR}variants_${i}/filtered_output_sID.vcf | grep -v "##contig=" > CDS_${i}_anno_uniq_SNPsel.tsv

###### make VCF 
cat ${SNP_DIR}variants_${i}/filtered_output_sID.vcf | grep "^#" > head_${i}

## merge vcf header
cat head_${i} CDS_${i}_anno_uniq_SNPsel.tsv > CDS_${i}_anno_uniq_SNPsel.vcf

##
rm head_*

done

echo "done..."

######


##########################
#SNP calling filter on Genes (whole genome SNP calling)

cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_PAS_anno_uniq.tsv | cut -f31 | grep -v "Gene_Name" > Genes_Pas_anno_uniq.tsv
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_REG_anno_uniq.tsv | cut -f31 | grep -v "Gene_Name" > Genes_Reg_anno_uniq.tsv
cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*Inter*_TIT_anno_uniq.tsv | cut -f31 | grep -v "Gene_Name" > Genes_Tit_anno_uniq.tsv

cat venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_*PAS-TIT-REG_anno_uniq.tsv| cut -f31 | grep -v "Gene_Name" > Genes_Pas-Reg-Tit_anno_uniq.tsv


# total SNP count 
#for i in Pas Reg Tit; do  echo $i; cat ${SNP_DIR}variants_${i}/filtered_output_sID_flr2_snpeff_Prot.vcf | grep -v "^#" | wc -l; done

for i in Pas Reg Tit
do 
echo $i

##### 
#SNP_DIR=/Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_2.1_STAR_Genome-CDS-ass/
SNP_DIR=/home/gala0002/proj/proj_Asa/NG-14833_2.1_STAR_Genome-ass/

grep -w -f Genes_${i}_anno_uniq.tsv ${SNP_DIR}variants_${i}/filtered_output_sID_flr2_snpeff_Prot.vcf | grep -v "#" > Genes_${i}_anno_uniq_SNPsel.tsv

###### make VCF 
cat ${SNP_DIR}variants_${i}/filtered_output_sID_flr2_snpeff_Prot.vcf | grep "^#" > head_${i}

## merge vcf header
cat head_${i} Genes_${i}_anno_uniq_SNPsel.tsv > Genes_${i}_anno_uniq_SNPsel.vcf


##
rm head_*

done

echo "done..."

######

for i in Pas Reg Tit Pas-Reg-Tit
do 
echo $i
echo "$i splice_acceptor_variant: "
grep -w "splice_acceptor_variant" Genes_${i}_anno_uniq_SNPsel.vcf | wc -l; 
echo "splice_donor_variant:"  
grep -w "splice_donor_variant" Genes_${i}_anno_uniq_SNPsel.vcf | wc -l;
echo "stop_gained:" 
grep -w "stop_gained" Genes_${i}_anno_uniq_SNPsel.vcf | wc -l;
echo "frameshift_variant:"
grep -w "frameshift_variant" Genes_${i}_anno_uniq_SNPsel.vcf | wc -l;
echo "stop_lost:" 
grep -w "stop_lost" Genes_${i}_anno_uniq_SNPsel.vcf | wc -l; 
echo "start_lost:"; 
grep -w "start_lost" Genes_${i}_anno_uniq_SNPsel.vcf | wc -l; 


done

######



### Merge using GATK 

REF_Quinoa=/home/gala0002/proj/proj_Asa/REF_Genome/Ref_Chenopodium_Quinoa/
work_dir_STAR=/home/gala0002/proj/proj_Asa/NG-14833_6.0_Venn-diagram/
#GATK4=/bioinfo/GATK/gatk-4.1.4.0/gatk/
#GATK37=/home/gala0002/softwares/GATK/GenomeAnalysisTK-3.7-0/
GATK38=/bioinfo/GATK/GenomeAnalysisTK-3.8-1/

#java -Xmx8g -jar /home/gala0002/softwares/GATK/GenomeAnalysisTK-3.7-0/GenomeAnalysisTK.jar
#/usr/lib/jvm/java-8-openjdk-amd64/bin/java -Xmx200g -jar 
# --minimumN 3 -genotypeMergeOptions UNIQUIFY \
/usr/lib/jvm/java-8-openjdk-amd64/bin/java -Xmx10g -jar ${GATK38}GenomeAnalysisTK.jar \
-T CombineVariants -nt 8 \
-genotypeMergeOptions UNIQUIFY \
-R ${REF_Quinoa}GCF_001683475.1_ASM168347v1_genomic.fna \
--variant ${work_dir_STAR}Genes_Pas_anno_uniq_SNPsel.vcf \
--variant ${work_dir_STAR}Genes_Reg_anno_uniq_SNPsel.vcf \
--variant ${work_dir_STAR}Genes_Tit_anno_uniq_SNPsel.vcf \
-o ${work_dir_STAR}Genes_all_anno_uniq_SNPsel_combine_unfy.vcf 

## common SNPs

grep -w -f Genes_all_anno_uniq.tsv /home/gala0002/proj/proj_Asa/NG-14833_2.1_STAR_Genome-ass/Variants_snpeff_Prot_allsamp-combine_min3-uniq.vcf | grep -v "#" > Genes_all_anno_uniq_SNPsel_common.tsv

################################################
#Pars with Python script
python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvReg_1371_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvReg_1371_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_6.0_Venn-diagram/
python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_645_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_645_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_6.0_Venn-diagram/
python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_RegvTit_96_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_RegvTit_96_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_6.0_Venn-diagram/

python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_TIT_anno_uniq.tsv venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_TIT_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_6.0_Venn-diagram/
python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500_PAS_anno_uniq.tsv venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500_PAS_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_6.0_Venn-diagram/
python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911_REG_anno_uniq.tsv venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911_REG_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_6.0_Venn-diagram/

python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_Inter_RegvTit_PAS-TIT-REG_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_Inter_RegvTit_PAS-TIT-REG_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_6.0_Venn-diagram/


###############################################
#all samples Uniq
#########

workdir="/Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_2.1_STAR_Genome-CDS-ass/DESeq2_genes/"

cat ${workdir}Express_CDS_Pas-Reg-Tit.gene.counts.matrix.Pas_vs_Reg.DESeq2.DE_results.P1e-2_C1.5.DE.subset_anno.tsv | awk -F"\t" '!_[$31]++' | awk '$31 != "."' > ${workdir}Express_CDS_Pas-Reg-Tit.gene.counts.matrix.Pas_vs_Reg.DESeq2.DE_results.P1e-2_C1.5.DE.subset_anno_uniq.tsv
cat ${workdir}Express_CDS_Pas-Reg-Tit.gene.counts.matrix.Pas_vs_Tit.DESeq2.DE_results.P1e-2_C1.5.DE.subset_anno.tsv | awk -F"\t" '!_[$31]++' | awk '$31 != "."' > ${workdir}Express_CDS_Pas-Reg-Tit.gene.counts.matrix.Pas_vs_Tit.DESeq2.DE_results.P1e-2_C1.5.DE.subset_anno_uniq.tsv
cat ${workdir}Express_CDS_Pas-Reg-Tit.gene.counts.matrix.Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.5.DE.subset_anno.tsv | awk -F"\t" '!_[$31]++' | awk '$31 != "."' > ${workdir}Express_CDS_Pas-Reg-Tit.gene.counts.matrix.Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.5.DE.subset_anno_uniq.tsv

python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911_anno.tsv venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911_anno_pars.tsv /Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/SL-2400_6.0_Venn-diagram_P5e-2_FC1/
python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_anno.tsv venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_anno_pars.tsv /Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/SL-2400_6.0_Venn-diagram_P5e-2_FC1/
python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500_anno.tsv venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500_anno_pars.tsv /Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/SL-2400_6.0_Venn-diagram_P5e-2_FC1/

#####

python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py all_DESeq2.DE_results.P5e-2_C1.0_PMTVWT-vs-Mock_anno.tsv all_DESeq2.DE_results.P5e-2_C1.0_PMTVWT-vs-Mock_anno_pars.tsv /Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/DESeq2_genes/
python Run-RNA-workflow-7.3_Pars-annotation-Venn-data-pars.py all_DESeq2.DE_results.P5e-2_C1.0_Delta8K-vs-Mock_anno.tsv all_DESeq2.DE_results.P5e-2_C1.0_Delta8K-vs-Mock_anno_pars.tsv /Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/DESeq2_genes/


#########
#Annotation with GO 
grep -w -A 1 -f venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911.txt /Volumes/Seagate/Backup-MAC_HD2/proj_Ramesh/RNA-seq_2020/Ref_Nicotiana_benthamiana/Supplemental_dataset_4_NbC_AA.fasta | sed 's/--//g' > venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911_AA.fa
grep -w -A 1 -f venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500.txt /Volumes/Seagate/Backup-MAC_HD2/proj_Ramesh/RNA-seq_2020/Ref_Nicotiana_benthamiana/Supplemental_dataset_4_NbC_AA.fasta | sed 's/--//g' > venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500_AA.fa
grep -w -A 1 -f venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791.txt /Volumes/Seagate/Backup-MAC_HD2/proj_Ramesh/RNA-seq_2020/Ref_Nicotiana_benthamiana/Supplemental_dataset_4_NbC_AA.fasta | sed 's/--//g' > venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_AA.fa
####

grep -w -A 1 -f all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Ramesh/RNA-seq_2020/Ref_Nicotiana_benthamiana/Supplemental_dataset_4_NbC_AA.fasta | sed 's/--//g' > all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_AA.fa
grep -w -A 1 -f all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_Ramesh/RNA-seq_2020/Ref_Nicotiana_benthamiana/Supplemental_dataset_4_NbC_AA.fasta | sed 's/--//g' > all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_AA.fa


#Annotation with KEGG
#

cat venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911_anno_pars.tsv | cut -f1,40 | sed 's/"//'g | sed 's/	/;/'g | awk -F";" '{ for (i = 2; i <= NF; ++i) { print $1, $i } }' OFS="\t" | awk '$2!="."' | awk '$1!="NB_GeneID"'> venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911_GOSLIM.tsv
cat venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500_anno_pars.tsv | cut -f1,40 | sed 's/"//'g | sed 's/	/;/'g | awk -F";" '{ for (i = 2; i <= NF; ++i) { print $1, $i } }' OFS="\t" | awk '$2!="."' | awk '$1!="NB_GeneID"'> venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500_GOSLIM.tsv
cat venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_anno_pars.tsv | cut -f1,40 | sed 's/"//'g | sed 's/	/;/'g | awk -F";" '{ for (i = 2; i <= NF; ++i) { print $1, $i } }' OFS="\t" | awk '$2!="."' | awk '$1!="NB_GeneID"' > venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_GOSLIM.tsv

####

# GOSLIM
cat all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars.tsv | cut -f1,40 | sed 's/"//'g | sed 's/	/;/'g | awk -F";" '{ for (i = 2; i <= NF; ++i) { print $1, $i } }' OFS="\t" | awk '$2!="."' | awk '$1!="NB_GeneID"'> all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_GOSLIM.tsv
cat all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars.tsv | cut -f1,40 | sed 's/"//'g | sed 's/	/;/'g | awk -F";" '{ for (i = 2; i <= NF; ++i) { print $1, $i } }' OFS="\t" | awk '$2!="."' | awk '$1!="NB_GeneID"'> all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_GOSLIM.tsv

# GO_P
cat all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars.tsv | cut -f1,32 | sed 's/"//'g | sed 's/	/;/'g | awk -F";" '{ for (i = 2; i <= NF; ++i) { print $1, $i } }' OFS="\t" | awk '$2!="."' | awk '$1!="NB_GeneID"'> all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_GO-P.tsv
cat all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars.tsv | cut -f1,32 | sed 's/"//'g | sed 's/	/;/'g | awk -F";" '{ for (i = 2; i <= NF; ++i) { print $1, $i } }' OFS="\t" | awk '$2!="."' | awk '$1!="NB_GeneID"'> all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_GO-P.tsv

# GO_F
cat all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars.tsv | cut -f1,34 | sed 's/"//'g | sed 's/	/;/'g | awk -F";" '{ for (i = 2; i <= NF; ++i) { print $1, $i } }' OFS="\t" | awk '$2!="."' | awk '$1!="NB_GeneID"'> all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_GO-F.tsv
cat all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars.tsv | cut -f1,34 | sed 's/"//'g | sed 's/	/;/'g | awk -F";" '{ for (i = 2; i <= NF; ++i) { print $1, $i } }' OFS="\t" | awk '$2!="."' | awk '$1!="NB_GeneID"'> all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_GO-F.tsv

# GO_C
cat all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars.tsv | cut -f1,36 | sed 's/"//'g | sed 's/	/;/'g | awk -F";" '{ for (i = 2; i <= NF; ++i) { print $1, $i } }' OFS="\t" | awk '$2!="."' | awk '$1!="NB_GeneID"'> all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_GO-C.tsv
cat all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars.tsv | cut -f1,36 | sed 's/"//'g | sed 's/	/;/'g | awk -F";" '{ for (i = 2; i <= NF; ++i) { print $1, $i } }' OFS="\t" | awk '$2!="."' | awk '$1!="NB_GeneID"'> all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_GO-C.tsv



#####

python Run-RNA-workflow-7.4_Pars-annotation-GO_KEGG.py all_DESeq2.DE_results.P5e-2_C1.0_Delta8K-vs-Mock_anno_pars.tsv all_DESeq2.DE_results.P5e-2_C1.0_Delta8K-vs-Mock_anno_pars-GO.tsv /Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/DESeq2_genes/
python Run-RNA-workflow-7.4_Pars-annotation-GO_KEGG.py all_DESeq2.DE_results.P5e-2_C1.0_PMTVWT-vs-Mock_anno_pars.tsv all_DESeq2.DE_results.P5e-2_C1.0_PMTVWT-vs-Mock_anno_pars-GO.tsv /Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/DESeq2_genes/


##############################
#| sed 's/NB_GeneID	GOSLIM//'g 
cat CDS_Pas_anno.tsv | awk '{split($1,a,"_"); print "XP_"a[5]}' | sort -u > CDS_Pas_anno_protID.tsv

grep -f CDS_Pas_anno_protID.tsv -A 1 /home/gala0002/proj/proj_Asa/REF_Genome/Ref_Chenopodium_Quinoa/GCF_001683475.1_ASM168347v1_protein_clean.fa | sed 's/--//g' > CDS_Pas_anno_protID.fa


echo "Script done...."



