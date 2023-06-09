#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)

rownum<-read.table(args[1])

let<-data.frame(Row=sapply(rownum,function(i) LETTERS[i]))
write.table(let,"metadata_rowforcol.txt",row.names=FALSE,col.names=FALSE)
