#!/usr/bin/env Rscript

library(data.table)
library(reshape)
library(dplyr)
library(sqldf)
library(ggplot2)

args=commandArgs(trailingOnly=TRUE)

plate1=args[1]
plate2=args[2]
filetype=args[3]

#getting the correct file
if(grepl("feat",filetype,fixed=TRUE)==TRUE){
    filename="Feature-Normalized"
    filename_sp=gsub("-"," ", filename,fixed=TRUE)
    mydf=paste0(plate1,"/",plate1,"_normalized_feature_select_batch.csv.gz")
    compdf=paste0(plate2,"/",basename(plate2),"_normalized_feature_select_batch.csv.gz")
} else if(grepl("neg",filetype,fixed=TRUE)==TRUE){
    filename="NegCon-Normalized"
    filename_sp=gsub("-"," ", filename,fixed=TRUE)
    mydf=paste0(plate1,"/",plate1,"_normalized_feature_select_negcon_batch.csv.gz")
    compdf=paste0(plate2,"/",basename(plate2),"_normalized_feature_select_negcon_batch.csv.gz")
} else {
    print("There is no matching filetype")
}

#reading file and getting plate ID
mydf=fread(mydf)
compdf=fread(compdf)

#melt the data and get a vector of measurements
reshape_data<-function(df){
    df_new<-df[,grep("Cells",colnames(df))[[1]]:ncol(df)]
    df<-cbind(Metadata_pert_iname=df$Metadata_pert_iname,df_new)
    df_melt<-melt(df)
    df_feat<-data.frame(Measurement=unique(df_melt$variable))
    return(df_feat)
}
mydf_meas<-reshape_data(mydf)
mydf_length<-length(mydf_meas$Measurement)
compdf_meas<-reshape_data(compdf)
compdf_length<-length(compdf_meas$Measurement)

#finding each length of features features
shared<-intersect(mydf_meas$Measurement,compdf_meas$Measurement)
shared_length<-length(shared)
mydf_unique<-sqldf('SELECT * FROM mydf_meas EXCEPT SELECT * FROM compdf_meas')
mydf_unique_length<-length(mydf_unique$Measurement)
compdf_unique<-sqldf('SELECT * FROM compdf_meas EXCEPT SELECT * FROM mydf_meas')
compdf_unique_length<-length(compdf_unique$Measurement)

meas_df<-data.frame(X=c(as.vector(mydf_meas$Measurement),as.vector(mydf_unique$Measurement),shared,as.vector(compdf_unique$Measurement),as.vector(compdf_meas$Measurement)),Category=rep(c("UTSW All","UTSW Unique","Intersection of Features","Broad Unique","Broad All"),times=c(mydf_length,mydf_unique_length,shared_length,compdf_unique_length,compdf_length)))
meas_df$Category<-factor(meas_df$Category,levels=c("Broad All","Broad Unique","Intersection of Features","UTSW Unique","UTSW All"))

ggplot(meas_df,aes(x=Category,group=Category)) + geom_bar(color="black",stat="count") + labs(title=paste0("Comparison of Number of Features Retained\n",filename_sp," Data\nPlate ",plateid)) + geom_text(aes(label = ..count..), stat = "count", vjust = -0.2,color="black")
ggsave(paste0("VennPlot_",filename,"_",plateid,".png"), type = "cairo")

#find df with max # of rows 
#df<-data.frame(Vec=c(deparse(substitute(shared)),deparse(substitute(mydf_unique)),deparse(substitute(compdf_unique))),NumRow=c(nrow(shared),nrow(mydf_unique),nrow(compdf_unique)))
#maxdf<-get(df[which.max(df$NumRow),1])

#features<-data.frame(plate1=c(mydf_unique$Measurement,rep(NA,nrow(maxdf)-nrow(mydf_unique))),Shared=c(shared$Measurement,rep(NA,nrow(maxdf)-length(shared$Measurement))),plate2=c(compdf_unique$Measurement,rep(NA,nrow(maxdf)-length(compdf_unique$Measurement))))
#names(features)[names(features)=="plate1"]<-paste0("Unique To ",plate1)
#names(features)[names(features)=="plate2"]<-paste0("Unique To ",plate2)

#write.csv(features,file=gzfile(paste0("FeatureVenn_",mydf_plateid,"vs",compdf_plateid,"_",filename,".csv.gz")),row.names=FALSE)