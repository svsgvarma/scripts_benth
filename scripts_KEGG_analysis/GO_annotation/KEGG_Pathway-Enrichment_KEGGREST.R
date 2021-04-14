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
#4. KEGG Pathway Enrichment Testing With KEGGREST
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html
## Pull all pathways for AT  
#pathways.list <- keggList("pathway", "ath")
pathways.list <- keggList("pathway", "nta")

head(pathways.list)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
head(genes.by.pathway)

##
#DE.table <- read.delim("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/master/friday/I5_v_C_time6.txt")

DE.table <- read.delim("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.P5e-2_C1.0_PMTVWT-vs-Mock_anno_KEGG-In.tsv")
#DE.table <- read.delim("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.P5e-2_C1.0_Delta8K-vs-Mock_anno_KEGG-In.tsv")

#DE.table <- read.delim("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.PMTVWT-vs-Mock_anno_KEGG-In.tsv")
#DE.table <- read.delim("/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.Delta8K-vs-Mock_anno_KEGG-In.tsv")

DE.table_rmdup = DE.table[!duplicated(DE.table[c("Gene")]),]

geneList <- DE.table_rmdup$P.Value
#geneList <- DE.table$FDR
names(geneList) <- DE.table_rmdup$Gene
head(geneList)
##
# GeneList filtered 
#geneList <- (DE.table %>% filter(P.Value < "0.05") %>% select(P.Value))[,1]
#names(geneList) <-(DE.table %>% filter(P.Value < "0.05") %>% select(Gene))[,1]
#head(geneList)

#Apply Wilcoxon rank-sum test to each pathway, testing if “in” p-values are smaller than “out” p-values:

# Wilcoxon test for each pathway
pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                             function(pathway) {
                               pathway.genes <- genes.by.pathway[[pathway]]
                               list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                               list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                               scores.in.pathway <- geneList[list.genes.in.pathway]
                               scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                               myvector <- paste(list.genes.in.pathway, collapse = '/')
                               
                               if (length(scores.in.pathway) > 0){
                                 p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
                                 #print(p.value)
                                 #print(scores.in.pathway)
                                 #print(p.value)
                                 #print(scores.not.in.pathway)
                               } else{
                                 p.value <- NA
                               }
                               return(c(p.value = p.value, Annotated = length(list.genes.in.pathway), Annotated_tot = length(pathway.genes), Annotated_genelist = myvector ))
                             }
))

# Assemble output table
outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- pathways.list[outdat$pathway.code]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat$Annotated_tot <- pVals.by.pathway[,"Annotated_tot"]
outdat$Annotated_genelist <- pVals.by.pathway[,"Annotated_genelist"]
outdat <- outdat[order(outdat$p.value),]
head(outdat)

#write.table(outdat, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars_KEGG-In_KEGGREST.tsv", row.names=T, sep="\t")
#write.table(outdat, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars_KEGG-In_KEGGREST.tsv", row.names=T, sep="\t")

#write.table(outdat, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.PMTVWT-vs-Mock_KEGGREST.tsv", row.names=T, sep="\t")
#write.table(outdat, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.delta8K-vs-Mock_KEGGREST.tsv", row.names=T, sep="\t")


########
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



########
#https://www.infoworld.com/article/3454356/how-to-merge-data-in-r-using-r-merge-dplyr-or-datatable.html
#Merge DF

outdat_PM <- outdat
outdat_delta <- outdat

joined_df <- merge(outdat_PM, outdat_delta, by.x = "pathway.code", 
                   by.y = "pathway.code", all.x = T, all.y = T)

#dataframe3 <- merge(outdat_PM, outdat_delta, by=c("pathway.code"),all=TRUE)

write.table(joined_df, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_Delta8K-and-PMTVWT_vs_Mock_KEGGREST_join_all.tsv", row.names=T, sep="\t")
write.table(joined_df, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_Delta8K-and-PMTVWT_vs_Mock_KEGGREST_join_DE-pval-0.05.tsv", row.names=T, sep="\t")



#############################################

