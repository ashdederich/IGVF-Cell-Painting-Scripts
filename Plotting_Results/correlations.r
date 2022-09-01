#!/usr/bin/env Rscript
#Create a correlation plot between two dataframes

#load in data
library(data.table)
library(reshape)
library(ggplot2)
library(ggpmisc)

args=commandArgs(trailingOnly=TRUE)

#read in the data frames
mydf=fread(args[1])
comparisondf=fread(args[2])
filetype=args[3]

#change data frames from short and wide to tall and skinny - My Data
mydf_new<-mydf[,grep("Cells",colnames(mydf))[[1]]:ncol(mydf)]
mydf_new<-cbind(mydf$Metadata_broad_sample,mydf_new)
colnames(mydf_new)[1]<-"Metadata_broad_sample"
mydf_new$Metadata_broad_sample[which(mydf_new$Metadata_broad_sample=="")]<-"DMSO"
mydf_new<-melt(mydf_new)
names(mydf_new)[names(mydf_new)=="variable"]<-"Measurement"
names(mydf_new)[names(mydf_new)=="value"]<-"Ashley_Median"

#change data frames from short and wide to tall and skinny - Broad Data
comparisondf_new<-comparisondf[,grep("Cells",colnames(comparisondf))[[1]]:ncol(comparisondf)]
comparisondf_new<-cbind(comparisondf$Metadata_broad_sample,comparisondf_new)
colnames(comparisondf_new)[1]<-"Metadata_broad_sample"
comparisondf_new$Metadata_broad_sample[which(comparisondf_new$Metadata_broad_sample=="")]<-"DMSO"
comparisondf_new<-melt(comparisondf_new)
names(comparisondf_new)[names(comparisondf_new)=="variable"]<-"Measurement"
names(comparisondf_new)[names(comparisondf_new)=="value"]<-"Broad_Median"

#merge data sheets by type
df_all<-merge(mydf_new,comparisondf_new,by=c("Metadata_broad_sample","Measurement"))

#creating correlation plot by Metadata_Compound
ggplot(df_all, aes(x=Ashley_Median,Broad_Median)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x, colour="black") + stat_poly_eq() + labs(x="Ashley's Measurements",y="Broad's Measurements",title=paste0("Correlation between Broad Data Set and Ours\nAfter PyCytominer Using the ",filetype," File"))
ggsave(paste0("CorrelationPlot_AllData_",filetype,".pdf"))

#summarizing by compound
ggplot(df_all, aes(x=Ashley_Median,y=Broad_Median,group=Metadata_broad_sample,color=Metadata_broad_sample)) + geom_point() + geom_smooth(method='lm',formula=y~x) + facet_wrap(~Metadata_broad_sample,scales="free") + labs(title=paste0("Correlation between Broad Data Set and Ours\nAfter PyCytominer Using the ",filetype," File")) + stat_poly_eq(label.x="right",label.y="bottom",size=2.5,color="black") + theme(strip.text = element_text(size = 6.2))
ggsave(paste0("CorrelationPlot_ByCmpd_",filetype,".pdf"))
