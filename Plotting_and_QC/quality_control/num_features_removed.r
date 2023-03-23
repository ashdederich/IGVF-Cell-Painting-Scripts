#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

library(data.table)
norm<-fread(args[1])
featsel<-fread(args[2])

norm_new<-norm[,grep("Cells",colnames(norm))[[1]]:ncol(norm)]
featsel_new<-featsel[,grep("Cells",colnames(featsel))[[1]]:ncol(featsel)]

norm<-cbind(Metadata_Well=norm$Metadata_Well,norm_new)
featsel<-cbind(Metadata_Well=featsel$Metadata_Well,featsel_new)

norm<-melt(norm,id.vars="Metadata_Well")
featsel<-melt(featsel,id.vars="Metadata_Well")

norm_feat<-unique(norm$variable)
featsel_feat<-unique(featsel$variable)

length(featsel_feat)
length(norm_feat)-length(featsel_feat)