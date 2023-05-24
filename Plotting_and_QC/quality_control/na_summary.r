#!/usr/bin/env Rscript
#Aggregate data output from Cell Painting.
#This file requires one input file input file to aggregate: this is the output file from CellProfiler.

args=commandArgs(trailingOnly=TRUE)
library(data.table)
library(reshape)
library(tools)
library(R.utils)
library(ggplot2)

#cells=args[1]
#cyto=args[2]
#nuclei=args[3]
#image=fread(args[4])

input_file<-args[1]

na_summary<-function(data){
    cell_location<-gsub(".*IGVF|.csv.gz*","",data)
    data<-fread(data)
    #plate<-unique(data$Metadata_Plate)
    data_dd<-data[,grep("AreaShape_Area",colnames(data)):ncol(data)]
    data=cbind(Metadata_Plate=data$Metadata_Plate,Metadata_Well=data$Metadata_Well,data_dd)
    data_melt<-melt(data,id.vars=c("Metadata_Plate","Metadata_Well"))
    na_count<-aggregate(value~variable,data_melt,function(x) { sum(is.na(x)) }, na.action=NULL)
    num_cells=length(which(data_melt$variable=="AreaShape_Area"))
    na_count$value=na_count$value/num_cells
    names(na_count)[names(na_count)=="value"]<-paste0(cell_location,"_NA_Count")
    na_count$CellCompartment=cell_location
    return(na_count)
}

na_count_by_row<-function(df){
    df<-transform(df,Row=substr(Metadata_Well,1,1),Column=substr(Metadata_Well,2,3))
    myLetters<-LETTERS[1:26]
    df$Row<-match(df$Row,myLetters)
    df<-cbind(plate=rep(plate,nrow(df)),df)
    df_melt<-melt(df,id.vars=c("Metadata_Well","Row","Column","plate"))
    df_melt$Row<-as.character(df_melt$Row)
    df_melt$Column<-sub("^0+","",df_melt$Column)
    df_melt$Column<-as.character(df_melt$Column)
    names(df_melt)[names(df_melt)=="value"]<-"NA_Count"
    return(df)
}

file_na<-na_summary(input_file)

ggplot(file_na,aes(x=variable,y=NA_Count))+geom_bar(stat='identity')+labs(title=paste0("Percentage of NAs Output by CellProfiler Per Feature, ",unique(file_na$CellCompartment)," File"),y="Percent NA",x="CellProfiler Feature in Alphabetical Order")+theme_bw()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())