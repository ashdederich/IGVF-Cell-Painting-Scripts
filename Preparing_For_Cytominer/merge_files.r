#!/usr/bin/env Rscript
#Merge each file outputted from CellProfiler (e.g. Cells.csv, Nucleus.csv, and Cytoplasm.csv) into one jointed plate csv file for input into PyCytominer.

library(data.table)
library(reshape)
library(R.utils)

args=commandArgs(trailingOnly=TRUE)
file_1=args[1]
file_2=args[2]
file_3=args[3]

#function to make files back into short and wide format
reshape_files<-function(a_file){
    data<-fread(a_file)
    data_cast<-cast(data,Metadata_Plate+Metadata_Well~variable,value="Median",fun.aggregate=sum)
    return(data_cast)
}

merge_files<-function(file1,file2,file3){
    #create a list of file names and reshape them with the cast function
    list_o_files<-list(file1,file2,file3)
    reshaped<-lapply(list_o_files,reshape_files)

    #merge the files together into one and write is as a csv file
    merged<-merge(reshaped[[1]],reshaped[[2]],by=c("Metadata_Plate","Metadata_Well"),all.x=TRUE,all.y=TRUE)
    merged1<-merge(merged,reshaped[[3]],by=c("Metadata_Plate","Metadata_Well"),all.x=TRUE,all.y=TRUE)
    platename<-unique(merged1$Metadata_Plate)
    platename<-gsub('/','',platename)
    write.csv(merged1,gzfile(paste0(platename,".csv.gz")),row.names = FALSE)
}

merge_files(file_1,file_2,file_3)
