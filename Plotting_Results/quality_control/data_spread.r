#!/usr/bin/env Rscript
#Calculate robust MAD for all features and then find the correlation between median/mad for two files
#This file requires one input file input file to aggregate: this is the output file from CellProfiler.

args=commandArgs(trailingOnly=TRUE)
library(data.table)
library(reshape)
library(tools)
library(ggplot2)
library(zplyr)

plate1=args[1]
plate2=args[2]
filetype=args[3]

plate1=sub("\\/.*","",plate1)

if(grepl("cp",filetype,fixed=TRUE)==TRUE){
    filename="CellProfiler-Output"
    filename_sp=gsub("-"," ", filename,fixed=TRUE)
    mydf=paste0(plate1,"/",plate1,".csv.gz")
    compdf=paste0(plate2,"/",basename(plate2),".csv.gz")
} else if(grepl("py",filetype,fixed=TRUE)==TRUE && grepl("feat",filetype,fixed=TRUE)==TRUE) {
    filename="PyCyto-Feature-Normalized"
    filename_sp=gsub("-"," ", filename,fixed=TRUE)
    mydf=paste0(plate1,"/",plate1,"_normalized_feature_select_batch.csv.gz")
    compdf=paste0(plate2,"/",basename(plate2),"_normalized_feature_select_batch.csv.gz")
} else if(grepl("py",filetype,fixed=TRUE)==TRUE && grepl("neg",filetype,fixed=TRUE)==TRUE){
    filename="PyCyto-NegCon-Normalized"
    filename_sp=gsub("-"," ", filename,fixed=TRUE)
    mydf=paste0(plate1,"/",plate1,"_normalized_feature_select_negcon_batch.csv.gz")
    compdf=paste0(plate2,"/",basename(plate2),"_normalized_feature_select_negcon_batch.csv.gz")
} else {
    print("There is no matching file")
}

mydf=fread(mydf)
compdf=fread(compdf)

aggregate_data<-function(file,file_name){
    plate=unique(file$Metadata_Plate)
    #melt the df
    file<-melt(file)
    names(file)[names(file)=="variable"]<-"Measurement"
    names(file)[names(file)=="value"]<-"Median"
    #get median and mad and rename the grouping variables
    median_ofdata=aggregate(file$Median,by=list(file$Measurement),median)
    mad_ofdata=aggregate(file$Median,by=list(file$Measurement),mad)
    names(median_ofdata)[names(median_ofdata)=="Group.1"]<-"Measurement"
    names(median_ofdata)[names(median_ofdata)=="x"]<-"Measurement_Median"
    names(mad_ofdata)[names(mad_ofdata)=="Group.1"]<-"Measurement"
    names(mad_ofdata)[names(mad_ofdata)=="x"]<-"Measurement_MAD"

    data_med_mad<-merge(median_ofdata,mad_ofdata,by="Measurement")
    data_med_mad<-cbind(Metadata_Plate=rep(plate,nrow(data_med_mad)),data_med_mad)
    return(data_med_mad)
}

mydf<-aggregate_data(mydf)
names(mydf)[names(mydf)=="Measurement_Median"]<-"FirstDF_Measurement_Median"
names(mydf)[names(mydf)=="Measurement_MAD"]<-"FirstDF_Measurement_MAD"

compdf<-aggregate_data(compdf)
names(compdf)[names(compdf)=="Measurement_Median"]<-"CompDF_Measurement_Median"
names(compdf)[names(compdf)=="Measurement_MAD"]<-"CompDF_Measurement_MAD"
#find difference in spread of data between both datasets and the correlation between the two

calc_spread<-function(firstdf,compdf){
    #merge the two dataframes
    merged<-merge(firstdf,compdf,by=c("Metadata_Plate","Measurement"))
    #find the difference in median
    merged$Median_Difference=(merged$FirstDF_Measurement_Median - merged$CompDF_Measurement_Median)
    median_stats<-data.frame(Median=round(median(merged$Median_Difference,na.rm=TRUE),2),SD=round(sd(merged$Median_Difference,na.rm=TRUE),2))
    median_stats<-cbind(Median_Difference=1,median_stats)
    #find the difference in MAD
    merged$MAD_Difference=(merged$FirstDF_Measurement_MAD - merged$CompDF_Measurement_MAD)
    mad_stats<-data.frame(Median=round(median(merged$MAD_Difference,na.rm=TRUE),2),SD=round(sd(merged$MAD_Difference,na.rm=TRUE),2))
    mad_stats<-cbind(MAD_Difference=1,mad_stats)
    #plot both in same plot and ggsave it
    p_median<-ggplot(merged,aes(x=Median_Difference)) + geom_histogram(binwidth=1,color="black") + labs(title=paste0("Difference between Median Values of each Feature Comparing\nBroad to UTSW Data Set with the ",filename_sp,"\nPlate ",unique(merged$Metadata_Plate))) + theme(strip.text = element_text(size = 6.2)) + geom_abs_text(data=median_stats,mapping = aes(label = paste0("Median=",Median)),color="black",size=3,xpos=0.15,ypos=0.9) + geom_abs_text(data=median_stats,mapping = aes(label = paste0("Standard Deviation=",SD)),color="black",size=3,xpos=0.15,ypos=0.85)
    ggsave(paste0("MedianDifference_",filename,"_",unique(merged$Metadata_Plate),".png"), type = "cairo")
    p_mad<-ggplot(merged,aes(x=MAD_Difference)) + geom_histogram(binwidth=1,color="black") + labs(title=paste0("Difference between MAD of each Feature Comparing\nBroad to UTSW Data Set with the ", filename_sp, "\nPlate ",unique(merged$Metadata_Plate))) + theme(strip.text = element_text(size = 6.2)) + geom_abs_text(data=mad_stats,mapping = aes(label = paste0("Median=",Median)),color="black",size=3,xpos=0.15,ypos=0.9) + geom_abs_text(data=mad_stats,mapping = aes(label = paste0("Standard Deviation=",SD)),color="black",size=3,xpos=0.15,ypos=0.85)
    ggsave(paste0("MADDifference_",filename,"_",unique(merged$Metadata_Plate),".png"), type = "cairo")
}

calc_spread(mydf,compdf)
