#!/usr/bin/env Rscript
#Aggregate data output from Cell Painting.
#This file requires one input file input file to aggregate: this is the output file from CellProfiler.

#GOAL: TAKE MEAN AND STD OF EACH FEATURE AND CREATE CORRELATION PLOTS AND NORMAL CURVE
args=commandArgs(trailingOnly=TRUE)
library(data.table)
library(reshape)
library(tools)
library(ggplot2)
library(ggpmisc)
library(zplyr)
library(dplyr)

datafile=fread(args[1])
comparisondf=fread(args[2])
filegroup=args[3]
group1="UTSW"
group2="Broad"
group1_name=args[4]
group2_name=args[5]

calc_mean<-function(input_file){
    #take mean and std of each column
    col_start<-grep("Cells",colnames(input_file))[[1]]
    input_file<-data.table(input_file)
    setkey(input_file,Metadata_broad_sample)
    mean_ofdata=aggregate(input_file[,col_start:ncol(input_file)], list(input_file$Metadata_broad_sample), mean)
    names(mean_ofdata)[names(mean_ofdata)=="Group.1"]<-"Metadata_broad_sample"

    #reshape the files
    mean_ofdata<-melt(mean_ofdata)
    names(mean_ofdata)[names(mean_ofdata)=="variable"]<-"Measurement"
    names(mean_ofdata)[names(mean_ofdata)=="value"]<-"Mean"
    mean_ofdata$Metadata_broad_sample[which(mean_ofdata$Metadata_broad_sample=="")]<-"DMSO"
    return(mean_ofdata)
}

calc_sd<-function(input_file){
    #take mean and std of each column
    col_start<-grep("Cells",colnames(input_file))[[1]]
    sd_ofdata=apply(input_file[,col_start:ncol(input_file)],2,sd)

    #reshape the files
    sd_ofdata<-melt(sd_ofdata)
    names(sd_ofdata)[names(sd_ofdata)=="value"]<-"Standard_Deviation"
    sd_ofdata<-data.frame(Measurement=rownames(sd_ofdata),Standard_Deviation=sd_ofdata$Standard_Deviation)
    return(sd_ofdata)
}

calc_and_merge<-function(input_file,groupname){
    input_mean<-calc_mean(input_file)
    input_sd<-calc_sd(input_file)
    names(input_mean)[names(input_mean)=="Mean"]<-paste0(groupname,"_Mean")
    names(input_sd)[names(input_sd)=="Standard_Deviation"]<-paste0(groupname,"_Standard_Deviation")
    mean_sd<-merge(input_mean,input_sd,by="Measurement")
    return(mean_sd)
}

my_mean_sd<-calc_and_merge(datafile,group1)
comp_mean_sd<-calc_and_merge(comparisondf,group2)
merged<-merge(my_mean_sd,comp_mean_sd,by=c("Measurement","Metadata_broad_sample"))
merged<-na.omit(merged)

correlate<-function(input_file){
    mean_lm<-lm(Broad_Mean~UTSW_Mean,input_file)
    mean_lmcoef<-coef(mean_lm)
    mean_lmcoef<-data.frame(Mean_Slope=round(mean_lmcoef[[2]],3))
    mean_lmcoef$Mean_Slope<-paste0("Mean Slope=",mean_lmcoef$Mean_Slope)

    sd_lm<-lm(Broad_Standard_Deviation~UTSW_Standard_Deviation,input_file)
    sd_lmcoef<-coef(sd_lm)
    sd_lmcoef<-data.frame(SD_Slope=round(sd_lmcoef[[2]],3))
    sd_lmcoef$SD_Slope<-paste0("SD Slope=",sd_lmcoef$SD_Slope)

    lmcoef=cbind(Metadata_broad_sample=1,UTSW_Mean=1,Broad_Mean=1,UTSW_Standard_Deviation=1,Broad_Standard_Deviation=1,mean_lmcoef,sd_lmcoef)
    return(lmcoef)
}

correl<-correlate(merged)

#plotting the correlation results
ggplot(merged, aes(x=UTSW_Mean,Broad_Mean,group=Metadata_broad_sample,color=Metadata_broad_sample)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x, colour="black") + facet_wrap(~Metadata_broad_sample,scales="free") + labs(title=paste0("Correlation between Broad and UTSW Mean Feature Values\nUTSW Outliers Removed"),x=paste0(group1_name," Mean"), y=paste0(group2_name," Mean")) + geom_abs_text(data=correl,mapping = aes(label = Mean_Slope),color="black",size=3.8,xpos=0.155,ypos=0.89) + stat_poly_eq(label.x="left",label.y="top",size=2.5,color="black")
ggsave(paste0("MeanCorrelationAcross_",filegroup,"_Features.png"), type = "cairo")

ggplot(merged,aes(UTSW_Standard_Deviation,Broad_Standard_Deviation,group=Metadata_broad_sample,color=Metadata_broad_sample)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x,colour="black") + facet_wrap(~Metadata_broad_sample,scales="free") + stat_poly_eq(label.x="left",label.y="top",size=2.5,color="black") + geom_abs_text(data=correl,mapping = aes(label = SD_Slope),color="black",size=3.8,xpos=0.143,ypos=0.89) + labs(title=paste0("Correlation of Standard Deviations between\nBroad and UTSW Feature Values"))
ggsave(paste0("SDCorrelationAcross_",filegroup,"_Features.png"), type = "cairo")