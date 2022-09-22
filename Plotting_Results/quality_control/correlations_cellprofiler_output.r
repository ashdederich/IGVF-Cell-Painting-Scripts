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
mydf=fread(args[1])
comparisondf=args[2]
filetype=args[3]
cmpd_df<-fread(args[4])

#getting metadata information for each file and merging with each respective dataframe
#for utsw
wd=getwd()
utsw_batchid=basename(wd)
utsw_plateid=unique(mydf$Metadata_Plate)[[1]]
utsw_barcode=fread(paste0("../../metadata/platemaps/",utsw_batchid,"/barcode_platemap.csv"))
utsw_barcode=utsw_barcode %>% filter(Assay_Plate_Barcode==utsw_plateid) %>% pull(var=Plate_Map_Name)
utsw_platemap=fread(paste0("../../metadata/platemaps/",utsw_batchid,"/platemap/",utsw_barcode,".txt"))
names(utsw_platemap)[names(utsw_platemap)=="well_position"]<-"Metadata_Well"
utsw_platemap=utsw_platemap[,1:3] #only get well_position and broad_sample information
mydf<-merge(utsw_platemap,mydf,by="Metadata_Well")
names(mydf)[names(mydf)=="broad_sample"]<-"Metadata_broad_sample"

#for broad
broad_batchid=basename(dirname(comparisondf))
comparisondf=fread(comparisondf)
broad_plateid=unique(comparisondf$Metadata_Plate)[1]
broad_barcode=fread(paste0("../../metadata/platemaps/",broad_batchid,"/barcode_platemap.csv"))
broad_barcode=broad_barcode %>% filter(Assay_Plate_Barcode==broad_plateid) %>% pull(var=Plate_Map_Name)
broad_platemap=fread(paste0("../../metadata/platemaps/",broad_batchid,"/platemap/",broad_barcode,".txt"))
names(broad_platemap)[names(broad_platemap)=="well_position"]<-"Metadata_Well"
broad_platemap=broad_platemap[,1:2] #only get well_position and broad_sample information
comparisondf<-merge(broad_platemap,comparisondf,by="Metadata_Well")
names(comparisondf)[names(comparisondf)=="broad_sample"]<-"Metadata_broad_sample"

#change data frames from short and wide to tall and skinny - My Data
mydf_new<-mydf[,grep("Cells",colnames(mydf))[[1]]:ncol(mydf)]
mydf_new<-cbind(mydf$Metadata_broad_sample,mydf_new)
colnames(mydf_new)[1]<-"Metadata_broad_sample"
mydf_new$Metadata_broad_sample[which(mydf_new$Metadata_broad_sample=="")]<-"DMSO"
compounds=unique(cmpd_df$Metadata_broad_sample)
mydf_subset=mydf_new[mydf_new$Metadata_broad_sample %in% compounds,]
mydf_new<-melt(mydf_subset)
names(mydf_new)[names(mydf_new)=="variable"]<-"Measurement"
names(mydf_new)[names(mydf_new)=="value"]<-"UTSW_Median"

#change data frames from short and wide to tall and skinny - Broad Data
comparisondf_new<-comparisondf[,grep("Cells",colnames(comparisondf))[[1]]:ncol(comparisondf)]
comparisondf_new<-cbind(comparisondf$Metadata_broad_sample,comparisondf_new)
colnames(comparisondf_new)[1]<-"Metadata_broad_sample"
comparisondf_new$Metadata_broad_sample[which(comparisondf_new$Metadata_broad_sample=="")]<-"DMSO"
comparisondf_subset=comparisondf_new[comparisondf_new$Metadata_broad_sample %in% compounds,]
comparisondf_new<-melt(comparisondf_subset)
names(comparisondf_new)[names(comparisondf_new)=="variable"]<-"Measurement"
names(comparisondf_new)[names(comparisondf_new)=="value"]<-"Broad_Median"

#merge data sheets by type
df_all<-inner_join(x=mydf_new,y=comparisondf_new,by=c("Metadata_broad_sample","Measurement"))
#calculate the slope and add to the plot
df_all_lm<-lm(UTSW_Median~Broad_Median,df_all)
df_all_lmcoef<-coef(df_all_lm)
df_all_lmcoef<-data.frame(Slope=round(df_all_lmcoef[[2]],3))
df_all_lmcoef$Slope<-paste0("Slope=",df_all_lmcoef$Slope)
df_all_lmcoef=cbind(UTSW_Median=1,Broad_Median=1,df_all_lmcoef)

#creating correlation plot
ggplot(df_all, aes(x=UTSW_Median,Broad_Median)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x,color="black") + stat_regline_equation(aes(label = ..rr.label..)) + labs(title=paste0("Correlation between Broad Data Set and UTSW\nUsing the ",filetype," File")) + geom_abs_text(data=df_all_lmcoef,mapping = aes(label = Slope),color="black",size=3.8,xpos=0.08,ypos=0.88)
ggsave(paste0("CorrelationPlot_AllData_",filetype,".png"), type = "cairo",width=10,height=7)

#summarizing by compound
#calculate lm for each broad sample
cmpd_lm<-dlply(df_all,"Metadata_broad_sample",function(df) lm(UTSW_Median~Broad_Median,data=df))
cmpd_lm_coef<-ldply(cmpd_lm,coef)
names(cmpd_lm_coef)[names(cmpd_lm_coef)=="Broad_Median"]<-"Slope"
cmpd_lm_coef[,2]=c(1:nrow(cmpd_lm_coef))
names(cmpd_lm_coef)[names(cmpd_lm_coef)=="(Intercept)"]<-"UTSW_Median"
cmpd_lm_coef=cbind(cmpd_lm_coef,Broad_Median=c(1:nrow(cmpd_lm_coef)))
cmpd_lm_coef[,3]=round(cmpd_lm_coef$Slope,3)
cmpd_lm_coef$Slope<-ldply(paste0("Slope=",cmpd_lm_coef$Slope))
colnames(cmpd_lm_coef[,3])<-"Slope"

ggplot(df_all, aes(x=UTSW_Median,y=Broad_Median,group=Metadata_broad_sample,color=Metadata_broad_sample)) + geom_point() + geom_smooth(method='lm',formula=y~x) + facet_wrap(~Metadata_broad_sample,scales="free") + labs(title=paste0("Correlation between Broad Data Set and Ours\nUsing the ",filetype," File")) + stat_poly_eq(label.x="left",label.y="top",size=2.5,color="black") + theme(strip.text = element_text(size = 6.2),axis.text=element_text(size=5)) + geom_abs_text(data=cmpd_lm_coef,mapping = aes(label = Slope),color="black",size=2.5,xpos=0.7,ypos=0.18)
ggsave(paste0("CorrelationPlot_ByCmpd_",filetype,".png"), type = "cairo")