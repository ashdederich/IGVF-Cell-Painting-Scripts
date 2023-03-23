#!/usr/bin/env Rscript
#Create a correlation plot between two dataframes

#load in data
library(data.table)
library(reshape)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(plyr)
library(zplyr)
library(dplyr)

args=commandArgs(trailingOnly=TRUE)

#read in the data frames
mydf=args[1]
comparisondf=args[2]
metadata=fread(args[3])

mydf=fread(paste0(mydf,"/",mydf,".csv.gz"))
comparisondf=paste0(comparisondf,"/",basename(comparisondf),".csv.gz")
plateid<-unique(mydf$Metadata_Plate)

#for broad
broad_batchid=basename(dirname(dirname(comparisondf)))
comparisondf=fread(comparisondf)
broad_plateid=unique(comparisondf$Metadata_Plate)[1]
broad_barcode=fread(paste0("../../metadata/platemaps/",broad_batchid,"/barcode_platemap.csv"))
broad_barcode=broad_barcode %>% filter(Assay_Plate_Barcode==broad_plateid) %>% pull(var=Plate_Map_Name)
broad_platemap=fread(paste0("../../metadata/platemaps/",broad_batchid,"/platemap/",broad_barcode,".txt"))
broad_platemap=data.frame(Metadata_Well=broad_platemap$well_position,broad_sample=broad_platemap$broad_sample) #only get well_position and broad_sample information
broad_metadata=inner_join(broad_platemap,metadata,by="broad_sample") #join metadata and platemap by broad_sample ID to get metadata pert iname
broad_metadata=data.frame(Metadata_Well=broad_metadata$Metadata_Well,Metadata_pert_iname=broad_metadata$pert_iname)
comparisondf<-merge(broad_metadata,comparisondf,by="Metadata_Well")

#change data frames from short and wide to tall and skinny - My Data
mydf_new=merge(broad_metadata,mydf,by="Metadata_Well")
mydf_subset=mydf_new[mydf_new$Metadata_pert_iname=="DMSO",] # only get DMSO wells
mydf_melt=melt(mydf_subset)
names(mydf_melt)[names(mydf_melt)=="variable"]<-"Measurement"
names(mydf_melt)[names(mydf_melt)=="value"]<-"UTSW_Median"

#change data frames from short and wide to tall and skinny - Broad Data
comparisondf_subset=comparisondf[comparisondf$Metadata_pert_iname=="DMSO",] # only get JUMP cmpds
comparisondf_new<-melt(comparisondf_subset)
names(comparisondf_new)[names(comparisondf_new)=="variable"]<-"Measurement"
names(comparisondf_new)[names(comparisondf_new)=="value"]<-"Broad_Median"

#merge data sheets by type
df_all<-inner_join(x=mydf_melt,y=comparisondf_new,by=c("Metadata_pert_iname","Measurement"))
#calculate the slope and add to the plot
df_all_lm<-lm(UTSW_Median~Broad_Median,df_all)
df_all_lmcoef<-coef(df_all_lm)
df_all_lmcoef<-data.frame(Slope=round(df_all_lmcoef[[2]],3))
df_all_lmcoef$Slope<-paste0("Slope=",df_all_lmcoef$Slope)
df_all_lmcoef=cbind(UTSW_Median=1,Broad_Median=1,df_all_lmcoef)

#creating correlation plot
ggplot(df_all, aes(x=UTSW_Median,y=Broad_Median)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x,color="black") + stat_regline_equation(aes(label = ..rr.label..)) + labs(title=paste0("Correlation between Broad and UTSW Data Sets - DMSO Only\nPlate ",plateid),x="UTSW Median", y="Broad Median") + geom_abs_text(data=df_all_lmcoef,mapping = aes(label = Slope),color="black",size=3.8,xpos=0.08,ypos=0.88)
ggsave(paste0("CorrelationPlot_",plateid,"-DMSOonly.png"), type = "cairo",width=10,height=7)