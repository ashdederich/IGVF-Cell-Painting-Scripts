#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)

row<-read.table(args[1])

myLetters<-LETTERS[1:26]
let<-data.frame(Row=match(row$V1,myLetters))
write.table(let,"metadata_row.txt",row.names=FALSE,col.names=FALSE)
