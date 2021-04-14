#!/usr/bin/RScript

# create_metadata-table.R
## create a metadata table that will be used later

# this table will contain information required to split the samples
# into groups and define contrasts for the differential expression analysis.
# You may build such table using your favorite spreadsheet editor or
# let R make it for you as shown below
# this code is saved under 'build-metadata.R'

# R CMD BATCH Run-RNA-workflow-4.2.1_create-metadata-table.R

#library("stringr")

#########################################
##--PO_PO-PP----

# Read CSV into R
MyData <- read.csv(file="/home/gala0002/proj/proj_Asa/NG-14833_4.2_express_combine-all/Express_counts_all_Pas-Reg-Tit_summary.csv", header=TRUE,check.names=FALSE, sep=",")
metadata <- as.data.frame(t(MyData[1,-1]))
# column name Rename 

library(data.table)
setDT(metadata, keep.rownames = TRUE)[]
names(metadata)[1:2]<- c("accession","read_count")

###
#PPD <- "SRP033351"
#ena.url <- paste("http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=",
#                 PPD,
#                 "&result=read_run",
#                 "&fields=run_accession,read_count",
#                 sep="")

# assemble into a data.frame
#metadata <- as.data.frame(read.table(url(metadata), header=TRUE, sep="\t"))

# add missing information about cells & treatment

# samples_all=c("484-10-S-1-O",
#               "484-10-S-2-O",
#               "484-10-S-3-O",
#               "484-12-S-1-O",
#               "484-12-S-2-O",
#               "484-12-S-3-O",
#               "485-10-S-1-O",
#               "485-10-S-2-O",
#               "485-10-S-3-O",
#               "485-12-S-1-O",
#               "485-12-S-2-O",
#               "485-12-S-3-O",
#               "506-10-R-1-O",
#               "506-10-R-2-O",
#               "506-10-R-3-O",
#               "506-12-R-1-O",
#               "506-12-R-2-O",
#               "506-12-R-3-O",
#               "565-10-R-1-O",
#               "565-10-R-2-O",
#               "565-10-R-3-O",
#               "565-12-R-1-O",
#               "565-12-R-2-O",
#               "565-12-R-3-O",
#               "Nimbus-10-S-1-P",
#               "Nimbus-10-S-2-P",
#               "Nimbus-10-S-3-P",
#               "Nimbus-12-S-1-P",
#               "Nimbus-12-S-2-P",
#               "Nimbus-12-S-3-P",
#               "Stigg-10-R-1-P",
#               "Stigg-10-R-2-P",
#               "Stigg-10-R-3-P",
#               "Stigg-12-R-1-P",
#               "Stigg-12-R-2-P",
#               "Stigg-12-R-3-P")

samples_all=c("Pas-15dpa-1",
"Pas-15dpa-2",
"Pas-15dpa-3",
"Pas-25dpa-1",
"Pas-25dpa-2",
"Pas-25dpa-3",
"Reg-15dpa-1",
"Reg-15dpa-2",
"Reg-15dpa-3",
"Reg-25dpa-1",
"Reg-25dpa-2",
"Reg-25dpa-3",
"Tit-15dpa-1",
"Tit-15dpa-2",
"Tit-15dpa-3",
"Tit-25dpa-1",
"Tit-25dpa-2",
"Tit-25dpa-3")

metadata$samples_all <- samples_all

# sample the first line
t(metadata[1,])

# add new columns
#metadata$cells <- str_extract(metadata$sample, "\\b[A-Za-z0-9]+")
#metadata$treatment <- mapply(sub,paste(metadata$cells,"_",sep=""), "", metadata$sample)

###
out_all <- strsplit(as.character(metadata$samples_all),'-')
#metadata_ad <- with(metadata, data.frame(attr = attr))
metadata_ad <- cbind(metadata, data.frame(t(sapply(out_all, `[`))))
names(metadata_ad)[4:6] <- c("sampID", "days", "dup")

# save to file for reuse
write.table(metadata_ad, file="/home/gala0002/proj/proj_Asa/NG-14833_4.2_express_combine-all/Express_counts_all_Pas-Reg-Tit_metadata_all.txt", row.names=F, sep="\t", quote=FALSE)

# view table
metadata_ad
#########################################
