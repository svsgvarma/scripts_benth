library(edgeR)
library(DESeq2)

library("pheatmap")

#https://molbiocloud.com/help/tiki-index.php?page=Differential-Expression-Analysis-Using-DESeq2#Differential_Expression_Analysis_Using_DESeq2

#res = read.table("/Volumes/Seagate/Backup-MAC_HD2/proj_Santosh/scripts/Trans.isoform.counts.matrix.AtgamFem_vs_AtgamMale.DESeq2.DE_results",sep="\t", header = TRUE, row.names = 1)
#res[ , !(names(res) %in% c("sampleA","sampleB","baseMeanA","baseMeanB") )]

#Importing the files
#outdir_plot="/Volumes/Seagate/Backup-MAC_HD2/proj_Asa/NG-14833_4.3_DEG_results/Results_plots/"
outdir_plot="/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/SL-2400_6.0_Venn-diagram_P5e-2_FC1/Results_plots_CDS/"

dir.create(file.path(outdir_plot), showWarnings = F, recursive = FALSE)

##################################################################################################
#allsamp_UP = c("AtgamFem-UP", "AtgamLarv-UP", "AtgamMale-UP")
#for (sampID_UP in allsamp_UP) {

#allsamp = c("Pas_vs_Reg", "Pas_vs_Tit", "Reg_vs_Tit")
#allsamp = c("Mk-v-Delta")
#allsamp = c("Mk-v-PMTVWT")
allsamp = c("Mk-v-Delta_Inter_Mk-v-PMTVWT")

#all count table:
for (sampID in allsamp) {
  #dgeFull_ann <- "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/SL-2400_6.0_Venn-diagram_P5e-2_FC1/venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_911_anno_pars.tsv"
  #dgeFull_ann <- "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/SL-2400_6.0_Venn-diagram_P5e-2_FC1/venn_DESeq2.DE_results.P5e-2_C1.0_PMTVWT_1500_anno_pars.tsv"
  dgeFull_ann <- "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/SL-2400_6.0_Venn-diagram_P5e-2_FC1/venn_DESeq2.DE_results.P5e-2_C1.0_Delta8K_Inter_PMTVWT_791_anno_pars.tsv"

  dgeFull_mat <- "/Volumes/Seagate/Backup-MAC_HD2/proj_Ramesh/RNA-seq_2020/Express_counts_all_Mock-PMTVWT-delta8k.matrix"
  #dgeFull_mat <- "/Volumes/Seagate/Backup-MAC_HD2/proj_Ramesh/RNA-seq_2020/DESeq2_genes/Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_delta8K.DESeq2.count_matrix"
  #dgeFull_mat <- "/Volumes/Seagate/Backup-MAC_HD2/proj_Ramesh/RNA-seq_2020/DESeq2_genes/Express_counts_all_Mock-PMTVWT-delta8k.matrix.Mock_vs_PMTVWT.DESeq2.count_matrix"
  #dgeFull_mat <- "/Volumes/Seagate/Backup-MAC_HD2/proj_Ramesh/RNA-seq_2020/DESeq2_genes/Express_counts_all_Mock-PMTVWT-delta8k.matrix.PMTVWT_vs_delta8K.DESeq2.count_matrix"
  
  #res_sel = read.table(dgeFull_ann,sep="\t", header = TRUE, row.names = 1, fill=FALSE)
  res_sel = read.table(dgeFull_ann,sep="\t", header = TRUE, row.names = 1, comment.char="#", na.strings=".", stringsAsFactors=FALSE, quote="", fill=FALSE)
  
  read.mat = read.table(dgeFull_mat,sep="\t", header = TRUE, row.names = 1)
  
  #coldata = data.frame(row.names = colnames(read.mat),group = rep(c("gt1","gt2"),each = 4),treatment = rep(c("Mock","delta8K"),each = 4))
  #coldata = data.frame(row.names = colnames(read.mat),group = rep(c("gt1","gt2"),each = 4),treatment = rep(c("Mock","PMTVWT"),each = 4))
  coldata = data.frame(row.names = colnames(read.mat),group = rep(c("gt1","gt2","gt3"),each = 4),treatment = rep(c("delta8K","Mock","PMTVWT"),each = 4))
  
  ####
  # Define the annotation
  annotation_row = data.frame(res_sel$Gene.Anno)
  #4. Run DESeqDataSetFromMatrix
  dds <- DESeqDataSetFromMatrix(countData = read.mat, colData = coldata, design = ~ treatment) 
  de.genes <- rownames(res_sel)
  vsd <- vst(dds, blind=FALSE)
  vsd_sel <- vsd[de.genes,]
  mat = assay(vsd_sel)[ order(res_sel$padj), ] # select the all genes with the lowest padj
  mat = mat - rowMeans(mat) # Subtract the row means from each value
  
  ###################
  # replace row names 
  rownames(mat) <- res_sel$Gene.Anno
  # Optional, but to make the plot nicer:
  df = as.data.frame(colData(vsd_sel)[,c("treatment")]) # Create a dataframe with a column of the conditions
  colnames(df) = "condition" # Rename the column header
  rownames(df) = colnames(mat) # add rownames
  
  ####################
  ## and plot the actual heatmap
  pheatmap(mat, cellwidth = 30, cellheight = 12, fontsize = 10, filename = paste(outdir_plot,sampID,"_Heatmap-Diff-exp-genes_selgenes_pairs.pdf", sep=""), annotation_col=df)
  #########################################################################
}



