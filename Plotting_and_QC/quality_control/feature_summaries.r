#!/usr/bin/env Rscript

library(data.table)
library(reshape)
library(dplyr)

mydf_10<-fread("Plate1_IGVF_Test_10_normalized_feature_select_batch.csv.gz")
mydf_3=fread("../Plate1_IGVF_Test_7_5/Plate1_IGVF_Test_7_5_normalized_feature_select_batch.csv")
compdf=fread("../../2020_11_04_CPJUMP1/Feature_Normalized_merged.csv.gz")

calc_mean<-function(input_file){
    #take mean and std of each column
    col_start<-grep("Cells",colnames(input_file))[[1]]
    mean_ofdata=colMeans(input_file[,col_start:ncol(input_file)])

    #reshape the files
    mean_ofdata<-melt(mean_ofdata)
    names(mean_ofdata)[names(mean_ofdata)=="value"]<-"Mean"
    mean_ofdata<-data.frame(Measurement=rownames(mean_ofdata),Mean=mean_ofdata$Mean)
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

calc_min<-function(input_file){
    col_start<-grep("Cells",colnames(input_file))[[1]]
    min_ofdata<-apply(input_file[,col_start:ncol(input_file)],2,min)

    min_ofdata<-melt(min_ofdata)
    names(min_ofdata)[names(min_ofdata)=="value"]<-"Minimum"
    min_ofdata<-data.frame(Measurement=rownames(min_ofdata),Minimum=min_ofdata$Minimum)
    return(min_ofdata)
}

calc_max<-function(input_file){
    col_start<-grep("Cells",colnames(input_file))[[1]]
    max_ofdata<-apply(input_file[,col_start:ncol(input_file)],2,max)

    max_ofdata<-melt(max_ofdata)
    names(max_ofdata)[names(max_ofdata)=="value"]<-"Maximum"
    max_ofdata<-data.frame(Measurement=rownames(max_ofdata),Maximum=max_ofdata$Maximum)
    return(max_ofdata)
}

mydf_10_mean<-calc_mean(mydf_10)
mydf_10_sd<-calc_sd(mydf_10)
mydf_10_min<-calc_min(mydf_10)
mydf_10_max<-calc_max(mydf_10)
names(mydf_10_mean)[names(mydf_10_mean)=="Mean"]<-"Mean_10uM"
names(mydf_10_sd)[names(mydf_10_sd)=="Standard_Deviation"]<-"SD_10uM"
names(mydf_10_min)[names(mydf_10_min)=="Minimum"]<-"Min_10uM"
names(mydf_10_max)[names(mydf_10_max)=="Maximum"]<-"Max_10uM"

mydf_3_mean<-calc_mean(mydf_3)
mydf_3_sd<-calc_sd(mydf_3)
mydf_3_min<-calc_min(mydf_3)
mydf_3_max<-calc_max(mydf_3)
names(mydf_3_mean)[names(mydf_3_mean)=="Mean"]<-"Mean_3.33uM"
names(mydf_3_sd)[names(mydf_3_sd)=="Standard_Deviation"]<-"SD_3.33uM"
names(mydf_3_min)[names(mydf_3_min)=="Minimum"]<-"Min_3.33uM"
names(mydf_3_max)[names(mydf_3_max)=="Maximum"]<-"Max_3.33uM"

compdf_mean<-calc_mean(compdf)
compdf_sd<-calc_sd(compdf)
compdf_min<-calc_min(compdf)
compdf_max<-calc_max(compdf)
names(compdf_mean)[names(compdf_mean)=="Mean"]<-"Mean_Broad"
names(compdf_sd)[names(compdf_sd)=="Standard_Deviation"]<-"SD_Broad"
names(compdf_min)[names(compdf_min)=="Minimum"]<-"Min_Broad"
names(compdf_max)[names(compdf_max)=="Maximum"]<-"Max_Broad"

merged<-inner_join(mydf_10_mean,mydf_10_sd,by="Measurement")
merged<-inner_join(merged,mydf_10_min,by="Measurement")
merged<-inner_join(merged,mydf_10_max,by="Measurement")
merged<-inner_join(merged,mydf_3_mean,by="Measurement")
merged<-inner_join(merged,mydf_3_sd,by="Measurement")
merged<-inner_join(merged,mydf_3_min,by="Measurement")
merged<-inner_join(merged,mydf_3_max,by="Measurement")
merged<-inner_join(merged,compdf_mean,by="Measurement")
merged<-inner_join(merged,compdf_sd,by="Measurement")
merged<-inner_join(merged,compdf_max,by="Measurement")
merged<-inner_join(merged,compdf_min,by="Measurement")
write.csv(merged,"feature_summary.csv")