#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)
library(data.table)
library(reshape)

df<-fread(args[1])
col_start<-grep("Cells",colnames(df))[[1]]
    #get mean and mad and rename the grouping variables
df_dat<-df[,col_start:ncol(df)]
df<-cbind(Metadata_Well=df$Metadata_Well,df_dat)
df<-melt(df)
max(df$value,na.rm=T)
min(df$value,na.rm=T)