#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")
#install.packages("clusterProfiler")

######################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("KEGGprofile")
BiocManager::install("enrichKEGG")

library(clusterProfiler)
library(dplyr)

library(enrichKEGG)
library(KEGGprofile)
library(KEGGREST)

###############

# gseMKEGG or enrichMKEGG

#https://yulab-smu.github.io/clusterProfiler-book/chapter6.html#kegg-gene-set-enrichment-analysis

#Chapter 6 KEGG analysis
#https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

#######
#3.1 Input data
#######
#d <- read.csv(your_csv_file)
#DE.table <- read.delim("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v3/all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars_KEGG-InFC.tsv")
#DE.table <- read.delim("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v3/all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars_KEGG-InFC.tsv")

DE.table <- read.delim("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.PMTVWT-vs-Mock_anno_KEGG-In.tsv")
#DE.table <- read.delim("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.Delta8K-vs-Mock_anno_KEGG-In.tsv")


#########################################################
#Prepare Input
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = DE.table[!duplicated(DE.table[c("Gene")]),]

# Create a vector of the gene unuiverse
kegg_gene_list <- dedup_ids$FC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- dedup_ids$Gene

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)


#Create gseKEGG object


kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'nta',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05, 
               pAdjustMethod = "none")

head(kk2, n = 20L)



#Merge DF

outdat_PM <- kk2
outdat_delta <- kk2

joined_df <- merge(outdat_PM, outdat_delta, by.x = "ID", 
                   by.y = "ID", all.x = T, all.y = T)

###plot

require(DOSE)
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

dotplot(outdat_PM , showCategory = 24, title = "Enriched Pathways in PMTVWT" , split=".sign") + facet_grid(.~.sign)
dotplot(outdat_delta, showCategory = 24, title = "Enriched Pathways in Delta8K" , split=".sign") + facet_grid(.~.sign)



#dataframe3 <- merge(outdat_PM, outdat_delta, by=c("pathway.code"),all=TRUE)

write.table(joined_df, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.delta8K-PMTVWT-vs-Mock_gseKEGG_join_all-0.05.tsv", row.names=T, sep="\t")



#########################################################

DE.table <- read.delim("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.PMTVWT-vs-Mock_anno_KEGG-In.tsv")
#DE.table <- read.delim("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.Delta8K-vs-Mock_anno_KEGG-In.tsv")

## assume that 1st column is ID
## 2nd column is fold change
## feature 1: numeric vector
geneList <- DE.table[,2]

## feature 2: named vector
names(geneList) <- as.character(DE.table[,1])

## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)

#######

#6.1 KEGG over-representation test
#data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 1]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'nta',
                 pvalueCutoff = 0.05)
head(kk)

### Merge

outdat_PM <- kk
outdat_delta <- kk

joined_df <- merge(outdat_PM, outdat_delta, by.x = "ID", 
                   by.y = "ID", all.x = T, all.y = T)

#dataframe3 <- merge(outdat_PM, outdat_delta, by=c("pathway.code"),all=TRUE)

write.table(joined_df, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.delta8K-PMTVWT-vs-Mock_enrichKEGG_join_all-0.05.tsv", row.names=T, sep="\t")


#6.3 KEGG Module over-representation test
mkk <- enrichMKEGG(gene = gene,
                   organism = 'nta')

head(mkk)
#6.4 KEGG Module Gene Set Enrichment Analysis

mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'nta')
head(mkk2)

#write.table(kk, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars_KEGG-In_enrichKEGG.tsv", row.names=T, sep="\t")
#write.table(kk, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars_KEGG-In_enrichKEGG.tsv", row.names=T, sep="\t")

write.table(kk, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.all_Mock_vs_PMTVWT_KEGG-In_enrichKEGG.tsv", row.names=T, sep="\t")
#write.table(kk, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.all_Mock_vs_delta8K_KEGG-In_enrichKEGG.tsv", row.names=T, sep="\t")


#########
#https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#browsekegg
#Chapter 12 Visualization of Functional
#12.11 browseKEGG

browseKEGG(kk, 'nta04712')


#12.12 pathview from pathview package
library("pathview")

nta04712 <- pathview(gene.data  = geneList,
                     pathway.id = 'nta04712',
                     species    = "nta",
                     limit      = list(gene=max(abs(geneList)), cpd=1))


####################################################################################
#https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#clusterprofiler

#KEGG enrichment analysis

#####
geneList <- DE.table %>% filter(P.Value < "0.05", !is.na(P.Value), abs(FC) < "1" ) 


## feature 1: numeric vector
genes <- geneList[,2]

## feature 2: named vector
names(genes) <- as.character(geneList[,1])

## feature 3: decreasing order
genes <- sort(genes, decreasing = TRUE)

##
kk <- enrichKEGG(gene = names(genes), organism = 'nta')
head(kk, n=10)



########
#https://www.infoworld.com/article/3454356/how-to-merge-data-in-r-using-r-merge-dplyr-or-datatable.html
#Merge DF

outdat_PM <- outdat
outdat_delta <- outdat

joined_df <- merge(outdat_PM, outdat_delta, by.x = "pathway.code", 
                   by.y = "pathway.code", all.x = T, all.y = T)

#dataframe3 <- merge(outdat_PM, outdat_delta, by=c("pathway.code"),all=TRUE)

write.table(joined_df, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v3/all_DESeq2.DE_results.all_Mock_vs_delta8K-PMTVWT_KEGG-In_KEGGREST_join_all.tsv", row.names=T, sep="\t")



#############################################

