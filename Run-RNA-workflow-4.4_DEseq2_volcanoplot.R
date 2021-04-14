##########################
#  devtools::install_github('kevinblighe/EnhancedVolcano')
#https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

#install.packages("pasilla")
#cat  Express_CDS_Pas-Reg-Tit.gene.counts.matrix | awk '{split($1,a,"_"); print a[4]"_"a[5]"_"a[6]"\t"$0}' | awk '{$2=""; print}' | sed s'/__/Pas_rep1/'g > Express_CDS_Pas-Reg-Tit.gene.counts.matrix_pars 
####
library(DESeq2)
library(edgeR)
library(EnhancedVolcano)
####
data = read.table("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/Express_counts_all_Mock-PMTVWT-delta8k.matrix", header=T, row.names=1, com='')

#### sample 1
col_ordering = c(1,2,3,4,5,6,7,8)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("Delta8K", 4), rep("QMock", 4))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","Delta8K","QMock")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Delta8K"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "QMock"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="Delta8K", sampleB="QMock", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])

#### sample 2
col_ordering = c(9,10,11,12,5,6,7,8)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("PMTVWT", 4), rep("QMock", 4))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","PMTVWT","QMock")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "PMTVWT"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "QMock"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="PMTVWT", sampleB="QMock", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])

#write.table(res, file='Express_CDS_Pas-Reg-Tit.gene.counts.matrix.Pas_vs_Reg.DESeq2.DE_results', sep='	', quote=FALSE)
#write.table(rnaseqMatrix, file='Express_CDS_Pas-Reg-Tit.gene.counts.matrix.Pas_vs_Reg.DESeq2.count_matrix', sep='	', quote=FALSE)
#source("/data/bioinfo/trinityrnaseq-Trinity-v2.6.5/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")

#### plot
#pdf("Trans.isoform.matrix.Larv_Male_Fem.DESeq2.DE_results.Volcano.pdf")
#plot_MA_and_Volcano(rownames(res), log2(res$baseMean+1), res$log2FoldChange, res$padj)
#dev.off()

#######################
#sampl1
featureNames <- rownames(res)
logFoldChange <- res$log2FoldChange
FDR <- res$padj

#######################
#library(EnhancedVolcano)
#######################
#selectLab = rownames(res)[1:40],
######
#### plot
#pdf("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/volcano_plots/Volcano_DESeq2_DE.CDS_Delta8K-vs-Mock.pdf")
pdf("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/volcano_plots/Volcano_DESeq2_DE.CDS_PMTVWT-vs-Mock.pdf")

EnhancedVolcano(res,
                title = "Volcano plot",
                #subtitle = "Delta8K-vs-Mock",
                subtitle = "PMTVWT-vs-Mock",
                selectLab = rownames(res)[0],
                lab = featureNames,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                labSize = 2.5,
                colAlpha = 1,
                legend=c('NS','Log2 FC','Adj p-value',
                         'Adj p-value & Log2FC'),
                legendPosition = 'top',
                legendLabSize = 13,
                legendIconSize = 4.0)
dev.off()

######

# plot_Volcano = function(featureNames, logFoldChange, FDR, xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch=20, top_gene_labels_show=20) {
#   plot(logFoldChange, -1*log10(FDR), col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);
#   text(logFoldChange[1:top_gene_labels_show], (-1*log10(FDR))[1:top_gene_labels_show], labels=featureNames[1:top_gene_labels_show], cex= 0.7, pos=3)
# }
# plot_Volcano(featureNames, logFoldChange, FDR);

####


