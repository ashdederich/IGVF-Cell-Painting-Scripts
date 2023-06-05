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
#should I add in an argument to read the platemap file and count the number of wells there should be?

aggregate_data<-function(input_file,fileheader,plateid){
    #get the filename
    filename<-sub("\\.csv.gz.*","",input_file)
    #get the cell location to append to the feature name
    cell_location<-gsub(paste0(".*",fileheader,"|.csv.gz*"),"",input_file)
    input_file=fread(input_file,blank.lines.skip=FALSE)
    wellcount=length(unique(input_file$Metadata_Well)) #set up a qc for wellcount to make sure we retain the correct number of wells
    #to count number of replicates, read in the platemap file, then exclude the DMSO wells. Then count then number of times each unique value occurs.
    # only keep Metadata_Plate and Metadata_Well and the measurements
    file<-input_file[,which(colnames(input_file)=='AreaShape_Area'):ncol(input_file)]
    #add the cellular location to the beginning of the measurement for each input file
    colnames(file)<-paste(cell_location,colnames(file),sep="_")
    file<-cbind(Metadata_Plate=plateid,Metadata_Well=input_file$Metadata_Well,file)
    
    #find start location of measurements
    col_start<-grep(paste0(cell_location,"_","AreaShape_Area"),colnames(file))
    #get mean and mad and rename the grouping variables
    mean_ofdata=aggregate(file[,col_start:ncol(file)],by=list(file$Metadata_Well),mean)
    names(mean_ofdata)[names(mean_ofdata)=="Group.1"]<-"Metadata_Well"

    #add plate info
    plate<-rep(plateid,nrow(mean_ofdata))
    mean_ofdata=cbind(Metadata_Plate=plate,mean_ofdata)

    #reshaping the data frames into long and skinny format
    mean_ofdata=melt(mean_ofdata)
    names(mean_ofdata)[names(mean_ofdata)=="value"]<-"Mean"
    agg_wellcount<-length(unique(mean_ofdata$Metadata_Well))
    if(wellcount=agg_wellcount){ # AND num reps are correct AND anything else???{ #should these be nested if statements?
        #outputting the file
        write.csv(mean_ofdata,gzfile(paste0(filename,"_aggregated",".csv.gz")))
    }else{
        print("After aggregation, there is not the correct number of wells and/or replicates.")
    }

}

aggregate_data(datafile,fileheader=args[2],plateid=args[3])