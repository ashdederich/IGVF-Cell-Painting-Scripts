#!/usr/bin/env Rscript
#Merge replicates separated across plates into one data sheet.

library(data.table)
library(reshape)
library(R.utils)

args=commandArgs(trailingOnly=TRUE)
file_1=args[1]
file_2=args[2]
file_3=args[3]
file_4=args[4]

#first need to reshape the dataframe to be tall and skinny for ease of merging the files

#change data frames from short and wide to tall and skinny - My Data
reshape_files<-function(a_file){
    data<-fread(a_file)
    data_new=data[,grep("Cells",colnames(data))[[1]]:ncol(data)]
    data_new<-cbind(Metadata_Plate=data$Metadata_Plate,Metadata_Well=data$Metadata_Well,data_new)
    data_melt<-melt(data_new)
    names(data_melt)[names(data_melt)=="variable"]<-"Measurement"
    names(data_melt)[names(data_melt)=="value"]<-"Median"
    return(data_melt)
}

#now binding the files
bind_and_cast<-function(file1,file2,file3,file4){
    #create a list of file names and reshape them with the reshape_files function
    list_o_files<-list(file1,file2,file3,file4)
    reshaped<-lapply(list_o_files,reshape_files)

    #rbind the files together into one, cast it, and write is as a csv file
    binded<-rbind(reshaped[[1]],reshaped[[2]],reshaped[[3]],reshaped[[4]])
    binded_cast<-cast(binded,Metadata_Plate+Metadata_Well~Measurement,value="Median",fun.aggregate=sum)
    write.csv(binded_cast,gzfile("CellProfiler-Output_merged.csv.gz"),row.names = FALSE)
}

#now need to cast the dataframe
bind_and_cast(file_1,file_2,file_3,file_4)