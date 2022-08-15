#!/usr/bin/env Rscript --vanilla
#Aggregate data output from Cell Painting.
#This file requires one input file input file to aggregate: this is the output file from CellProfiler.

args=commandArgs(trailingOnly=TRUE)
library(data.table)
library(reshape)
library(tools)

datafile=fread(args[1])
filename=args[1]
plateid=unique(datafile$Metadata_Plate)

aggregate_data<-function(input_file,file_name){
    # only keep Metadata_Plate and Metadata_Well and the measurements
    file<-input_file[,which(colnames(input_file)=='AreaShape_Area'):ncol(input_file)]
    #add the cellular location to the beginning of the measurement for each input file
    cell_location<-gsub("IGVF|.csv*","",filename)
    colnames(file)<-paste(cell_location,colnames(file),sep="_")
    file<-cbind(input_file$Metadata_Plate,input_file$Metadata_Well,file)
    colnames(file)[1]<-"Metadata_Plate"
    colnames(file)[2]<-"Metadata_Well"
    
    #find start location of measurements
    col_start<-grep(paste0(cell_location,"_","AreaShape_Area"),colnames(file))
    #get median and mad and rename the grouping variables
    median_ofdata=aggregate(file[,col_start:ncol(file)],by=list(file$Metadata_Well),median)
    mad_ofdata=aggregate(file[,col_start:ncol(file)],by=list(file$Metadata_Well),mad)
    names(median_ofdata)[names(median_ofdata)=="Group.1"]<-"Metadata_Well"
    names(mad_ofdata)[names(mad_ofdata)=="Group.1"]<-"Metadata_Well"

    #add plate info
    plate<-rep(plateid,nrow(median_ofdata))
    median_ofdata=cbind(plate,median_ofdata)
    mad_ofdata=cbind(plate,mad_ofdata)

    median_ofdata=melt(median_ofdata)
    mad_ofdata=melt(mad_ofdata)
    names(median_ofdata)[names(median_ofdata)=="value"]<-"Median"
    names(mad_ofdata)[names(mad_ofdata)=="value"]<-"MAD"

    data_med_mad<-merge(median_ofdata,mad_ofdata,by=c("Metadata_Well","variable","plate"))
    file_ext=file_path_sans_ext(file_name,compression=FALSE)
    write.csv(data_med_mad,paste0(file_ext,"_aggregated",".csv"))
}

aggregate_data(datafile,filename)