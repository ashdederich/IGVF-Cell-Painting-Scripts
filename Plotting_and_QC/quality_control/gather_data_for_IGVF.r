#!/usr/bin/env Rscript
#This file takes the feature-selected, normalized file from CellProfiler and, using the subset_name argument, searches for all features that include the subset_name. It then finds the maximum feature value for each location, subsets the input file for those features, and reshapes the dataframe for input into GeneData.
#This file requires two arguments: the feature-selected, normalized file and the subset_name argument, which can currently only take on four values: AreaShape, Texture, Intensity, and Radial Distribution.

args=commandArgs(trailingOnly=TRUE)
library(data.table)
library(reshape)
library(tools)
library(R.utils)
library(dplyr)
library(berryFunctions)
library(stringr)

df<-fread(args[1])
subset_name<-args[2]
if(subset_name=="Radial_Distribution"){
    subset_name_search="RadialDistribution"
}else{
    subset_name_search=subset_name
}

df_sub<-df %>% select(contains(subset_name_search))
df_sub<-cbind(Metadata_Well=df$Metadata_Well,df_sub)
df_sub<-melt(df_sub)
names(df_sub)[names(df_sub)=="variable"]<-"Measurement"
meas<-as.character(df_sub$Measurement)
if(subset_name=="AreaShape"){
    locations<-sapply(strsplit(meas,"_"),"[",1)
    df_sub$Locations<-locations
    DT<-data.table(df_sub)
    max_meas<-DT[,.SD[which.max(value)],by=Locations] %>% select(Measurement)
    max_df<-merge(DT,max_meas,by="Measurement")
}else if(subset_name=="Texture" | subset_name=="Intensity"){
    #taking a guess hear that the locations occur after the third '_' - COULD BE WRONG
    locations<-t(sapply(strsplit(meas,"_"),"["))[,4]
    df_sub$Locations<-locations
    DT<-data.table(df_sub)
    max_meas<-DT[,.SD[which.max(value)],by=Locations] %>% select(Measurement)
    max_df<-merge(DT,max_meas,by="Measurement")
}else if(subset_name=="Radial_Distribution"){
    #CCN = cells,cytoplasm,nuclei
    locationsCCN<-sapply(strsplit(meas,"_"),"[",1)
    #DEAM=DNA,ER,AGP,Mito
    locationsDEAM<-t(sapply(strsplit(meas,"_"),"["))[,4]
    df_sub$LocationsCCN<-locationsCCN
    df_sub$LocationsDEAM<-locationsDEAM
    DT<-data.table(df_sub)
    max_measCCN<-DT[,.SD[which.max(value)],by=LocationsCCN] %>% select(Measurement)
    max_measDEAM<-DT[,.SD[which.max(value)],by=LocationsDEAM] %>% select(Measurement)
    max_dfCCN<-merge(DT,max_measCCN,by="Measurement")
    max_dfCCN=max_dfCCN[,-5]
    names(max_dfCCN)[names(max_dfCCN)=="LocationsCCN"]<-"Locations"
    max_dfDEAM<-merge(DT,max_measDEAM,by="Measurement")
    max_dfDEAM<-max_dfDEAM[,-4]
    names(max_dfDEAM)[names(max_dfDEAM)=="LocationsDEAM"]<-"Locations"
    max_df<-rbind(max_dfCCN,max_dfDEAM)   
}else{
    "There were no measurements within the subset_name context that contained Cells, Cytoplasm, Nuclei, DNA, Mito, AGP, or ER."
}

letter_to_number<-function(df){
    myLetters<-LETTERS[1:26]
    let<-data.frame(Row=match(df$Row,myLetters))
    return(let)
}

reshape_df<-function(df,location){
    df<-df%>%filter(Locations==location)
    df$Row<-substr(df$Metadata_Well,1,1)
    df$Row<-letter_to_number(df)
    df$Col<-substr(df$Metadata_Well,2,3)
    df<-cast(df,Row~Col,value="value")
    num_cols_missing<-24-ncol(df[,2:ncol(df)])
    if(num_cols_missing>0){
        col_start=ncol(df)+1
        df[,col_start:24]<-0
    }
    df<-setnafill(df,fill=0)
    colnames(df)<-c(" ",colnames(df)[2]:24)
    #add blank rows and colnames
    df<-insertRows(df,1,new=" ")
    df<-insertRows(df,2,new=" ")
    df[2,1]<-location
    df<-insertRows(df,3,new=colnames(df))
    return(df)
}

grab_measurements<-function(df){
    meas<-unique(df$Measurement)
    write.csv(meas,paste0(subset_name,"_Measurements.csv"),row.names=FALSE)
}

if(subset_name=="AreaShape"){
    cells=reshape_df(max_df,"Cells")
    cells=cells[-1,]
    cytoplasm=reshape_df(max_df,"Cytoplasm")
    nuclei=reshape_df(max_df,"Nuclei")
    final_df<-rbind(cells,cytoplasm,nuclei)
    write.table(final_df,sep=',',gzfile(paste0(subset_name,"_NormalizedValues.csv.gz")),row.names=FALSE,col.names=FALSE)
}else if(subset_name=="Texture" | subset_name=="Intensity"){
    DNA=reshape_df(max_df,"DNA")
    DNA=DNA[-1,]
    Mito=reshape_df(max_df,"Mito")
    ER=reshape_df(max_df,"ER")
    AGP=reshape_df(max_df,"AGP")
    final_df=rbind(DNA,Mito,ER,AGP)
    write.table(final_df,sep=',',gzfile(paste0(subset_name,"_NormalizedValues.csv.gz")),row.names=FALSE,col.names=FALSE)
}else if(subset_name=="Radial_Distribution"){
    cells=reshape_df(max_df,"Cells")
    cells=cells[-1,]
    cytoplasm=reshape_df(max_df,"Cytoplasm")
    nuclei=reshape_df(max_df,"Nuclei")
    Mito=reshape_df(max_df,"Mito")
    ER=reshape_df(max_df,"ER")
    AGP=reshape_df(max_df,"AGP")
    if(any(grepl("DNA",meas,fixed=TRUE))==TRUE){
        DNA=reshape_df(max_df,"DNA")
        DNA=DNA[-1,]
        final_df=rbind(cells,cytoplasm,nuclei,DNA,Mito,ER,AGP)
    }else{
        final_df=rbind(cells,cytoplasm,nuclei,Mito,ER,AGP)
    }
    write.table(final_df,sep=',',gzfile(paste0(subset_name,"_NormalizedValues.csv.gz")),row.names=FALSE,col.names=FALSE)
}else{
    print("I have not yet set up this function to search for any CellProfiler measurements that do not include AreaShape, Texture, Intensity, or Radial Distribution.")
}

grab_measurements(max_df)