#!/usr/bin/env Rscript
#Aggregate data output from Cell Painting.
#This file requires one input file input file to aggregate: this is the output file from CellProfiler.

args=commandArgs(trailingOnly=TRUE)
library(data.table)
library(reshape)
library(tools)
library(R.utils)

#input file - i.e. IGVFCells.csv.gz
datafile=args[1]

aggregate_data<-function(input_file){
    #get the filename
    filename<-sub("\\.csv.gz.*","",input_file)
    #get the cell location to append to the feature name
    cell_location<-gsub(".*IGVF|.csv.gz*","",input_file)
    input_file=fread(input_file)
    plateid=unique(input_file$Metadata_Plate)
    # only keep Metadata_Plate and Metadata_Well and the measurements
    file<-input_file[,which(colnames(input_file)=='AreaShape_Area'):ncol(input_file)]
    #add the cellular location to the beginning of the measurement for each input file
    colnames(file)<-paste(cell_location,colnames(file),sep="_")
    file<-cbind(Metadata_Plate=input_file$Metadata_Plate,Metadata_Well=input_file$Metadata_Well,file)
    
    #find start location of measurements
    col_start<-grep(paste0(cell_location,"_","AreaShape_Area"),colnames(file))
    #get mean and mad and rename the grouping variables
    mean_ofdata=aggregate(file[,col_start:ncol(file)],by=list(file$Metadata_Well),mean)
    mad_ofdata=aggregate(file[,col_start:ncol(file)],by=list(file$Metadata_Well),mad)
    names(mean_ofdata)[names(mean_ofdata)=="Group.1"]<-"Metadata_Well"
    names(mad_ofdata)[names(mad_ofdata)=="Group.1"]<-"Metadata_Well"

    #add plate info
    plate<-rep(plateid,nrow(mean_ofdata))
    mean_ofdata=cbind(plate,mean_ofdata)
    mad_ofdata=cbind(plate,mad_ofdata)

    #reshaping the data frames into long and skinny format
    mean_ofdata=melt(mean_ofdata)
    mad_ofdata=melt(mad_ofdata)
    names(mean_ofdata)[names(mean_ofdata)=="value"]<-"Median"
    names(mad_ofdata)[names(mad_ofdata)=="value"]<-"MAD"

    #merging the median and mad data frames by their shared columns
    data_med_mad<-merge(mean_ofdata,mad_ofdata,by=c("Metadata_Well","variable","plate"))
    #outputting the file
    write.csv(data_med_mad,gzfile(paste0(filename,"_aggregated",".csv.gz")))
}

aggregate_data(datafile)