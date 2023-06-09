#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
library(data.table)
library(reshape)

df=fread(args[1])

for(i in unique(df$Metadata_pert_iname)){
    print(length(which(df$Metadata_pert_iname==i)))
}