#!/usr/bin/env Rscript
#Merge each file outputted from CellProfiler (e.g. Cells.csv, Nucleus.csv, and Cytoplasm.csv) into one jointed plate csv file for input into PyCytominer.

library(data.table)
library(reshape)
library(R.utils)

args=commandArgs(trailingOnly=TRUE)
file_1=args[1]
file_2=args[2]
file_3=args[3]


prep_files<-function(input_file){
    cell_location<-gsub(".*IGVF|.csv.gz*","",input_file)
    input_file=fread(input_file)
    # only keep Metadata_Plate and Metadata_Well and the measurements
    file<-input_file[,which(colnames(input_file)=='AreaShape_Area'):ncol(input_file)]
    #add the cellular location to the beginning of the measurement for each input file
    colnames(file)<-paste(cell_location,colnames(file),sep="_")
    file<-cbind(Metadata_Plate=input_file$Metadata_Plate,Metadata_Well=input_file$Metadata_Well,ObjectNumber=input_file$ObjectNumber,ImageNumber=input_file$ImageNumber,Metadata_Site=input_file$Metadata_Site,file)
}

file_1=prep_files(file_1)
file_2=prep_files(file_2)
file_3=prep_files(file_3)

merged<-merge(file_1,file_2,by=c("Metadata_Plate","Metadata_Well","ObjectNumber","ImageNumber","Metadata_Site"))
merged<-merge(merged,file_3,by=c("Metadata_Plate","Metadata_Well","ObjectNumber","ImageNumber","Metadata_Site"))
merged=subset(merged,select=-c(ObjectNumber,ImageNumber,Metadata_Site))
write.csv(merged,gzfile(paste0(unique(merged$Metadata_Plate),".csv.gz")))