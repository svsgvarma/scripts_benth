source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
install.packages("clusterProfiler")
library(clusterProfiler)

gda <- read.delim("GOannotation.txt", header=TRUE)
dim(gda)
head(gda)

term2gene=gda[, c("description", "geneID")]

#gene=read.table('common_24hr-48hr-72hr.txt')
gene=read.table('common_M-all.txt')
gene <- as.vector(gene[,1])
gene

x = enricher(gene, TERM2GENE=term2gene)
summary(x)

dotplot(x)
barplot(x, showCategory = 12)
enrichMap(x)

x
#######################


######################

bitr_kegg("K01489", "kegg", "Path", "ko") -> x
ko2name(x$Path) -> y
merge(x, y, by.x='Path', by.y='nta')

nic <- search_kegg_organism('Nicotiana tabacum', by='scientific_name')



#########

##4.3 KEGG pathway enrichment analysis

kk <- enrichKEGG(gene = gene, organism = "human", pvalueCutoff = 0.01,readable = TRUE)
head(summary(kk))

######################


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGGprofile")

library(enrichKEGG)
library(KEGGprofile)
library(KEGGREST)

###############
#4. KEGG Pathway Enrichment Testing With KEGGREST
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html
## Pull all pathways for AT  
pathways.list_nta <- keggList("pathway", "nta")
head(pathways.list_nta)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list_nta)) 
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
DE.table <- read.delim("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/master/friday/I5_v_C_time6.txt")

geneList <- DE.table$P.Value
names(geneList) <- DE.table$Gene
head(geneList)

