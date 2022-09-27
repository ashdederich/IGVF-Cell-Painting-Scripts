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
filetype_sp<-gsub("-", " ", filetype, fixed=TRUE)
cmpd_df<-fread(args[4])
comp1<-args[5]
comp1<-gsub("_"," ", comp1,fixed=TRUE)
comp2<-args[6]
comp2<-gsub("_"," ", comp2,fixed=TRUE)

#take cmpd_df, melt it, and get a list of the measurements kept
cmpd_df_new<-cmpd_df[,grep("Cells",colnames(cmpd_df))[[1]]:ncol(cmpd_df)]
cmpd_df_new<-cbind(Metadata_pert_iname=cmpd_df$Metadata_pert_iname,Metadata_Well=cmpd_df$Metadata_Well,cmpd_df_new)
cmpd_df_new<-melt(cmpd_df_new,id.vars=c("Metadata_pert_iname","Metadata_Well"))
cmpd_df_new<-data.frame(Metadata_pert_name=cmpd_df_new$Metadata_pert_iname,Metadata_Well=cmpd_df_new$Metadata_Well)
cmpd_df_new<-cmpd_df_new[!duplicated(cmpd_df_new),]
cmpd_meas<-as.character(unique(cmpd_df_new$variable))
compounds=unique(cmpd_df$Metadata_pert_iname) 

#getting metadata information for each file and merging with each respective dataframe
#for utsw
#wd=getwd()
#utsw_batchid=basename(dirname(wd))
#utsw_plateid=unique(mydf$Metadata_Plate)[[1]]
#utsw_barcode=fread(paste0("../../../metadata/platemaps/",utsw_batchid,"/barcode_platemap.csv"))
#utsw_barcode=utsw_barcode %>% filter(Assay_Plate_Barcode==utsw_plateid) %>% pull(var=Plate_Map_Name)
#utsw_platemap=fread(paste0("../../../metadata/platemaps/",utsw_batchid,"/platemap/",utsw_barcode,".txt"))
#names(utsw_platemap)[names(utsw_platemap)=="well_position"]<-"Metadata_Well"
#utsw_platemap=utsw_platemap[,1:3] #only get well_position and broad_sample information
#mydf<-merge(utsw_platemap,mydf,by="Metadata_Well")
#names(mydf)[names(mydf)=="broad_sample"]<-"Metadata_pert_iname"
mydf=inner_join(cmpd_df_new,mydf,by="Metadata_Well")

#for broad
broad_batchid=basename(dirname(comparisondf))
comparisondf=fread(comparisondf)
broad_plateid=unique(comparisondf$Metadata_Plate)[1]
broad_barcode=fread(paste0("../../../metadata/platemaps/",broad_batchid,"/barcode_platemap.csv"))
broad_barcode=broad_barcode %>% filter(Assay_Plate_Barcode==broad_plateid) %>% pull(var=Plate_Map_Name)
broad_platemap=fread(paste0("../../../metadata/platemaps/",broad_batchid,"/platemap/",broad_barcode,".txt"))
names(broad_platemap)[names(broad_platemap)=="well_position"]<-"Metadata_Well"
broad_platemap=broad_platemap[,1:2] #only get well_position and broad_sample information
comparisondf<-merge(broad_platemap,comparisondf,by="Metadata_Well")
names(comparisondf)[names(comparisondf)=="broad_sample"]<-"Metadata_pert_iname"

#change data frames from short and wide to tall and skinny - My Data
mydf_new<-mydf[,grep("Cells",colnames(mydf))[[1]]:ncol(mydf)]
mydf_new<-cbind(mydf$Metadata_pert_iname,mydf_new)
colnames(mydf_new)[1]<-"Metadata_pert_iname"
mydf_new$Metadata_pert_iname[which(mydf_new$Metadata_pert_iname=="")]<-"DMSO"
mydf_subset=mydf_new[mydf_new$Metadata_pert_iname %in% compounds,] # only get JUMP cmpds
mydf_new<-melt(mydf_subset)
names(mydf_new)[names(mydf_new)=="variable"]<-"Measurement"
names(mydf_new)[names(mydf_new)=="value"]<-"UTSW_Median"
mydf_new=mydf_new[mydf_new$Measurement %in% cmpd_meas,] #only get measurements that pycytominer selected as variable

#change data frames from short and wide to tall and skinny - Broad Data
comparisondf_new<-comparisondf[,grep("Cells",colnames(comparisondf))[[1]]:ncol(comparisondf)]
comparisondf_new<-cbind(comparisondf$Metadata_pert_iname,comparisondf_new)
colnames(comparisondf_new)[1]<-"Metadata_pert_iname"
comparisondf_new$Metadata_pert_iname[which(comparisondf_new$Metadata_pert_iname=="")]<-"DMSO"
comparisondf_subset=comparisondf_new[comparisondf_new$Metadata_pert_iname %in% compounds,] # only get JUMP cmpds
comparisondf_new<-melt(comparisondf_subset)
names(comparisondf_new)[names(comparisondf_new)=="variable"]<-"Measurement"
names(comparisondf_new)[names(comparisondf_new)=="value"]<-"Broad_Median"
comparisondf_new=comparisondf_new[comparisondf_new$Measurement %in% cmpd_meas,]

#merge data sheets by type
df_all<-inner_join(x=mydf_new,y=comparisondf_new,by=c("Metadata_pert_iname","Measurement"))
#calculate the slope and add to the plot
df_all_lm<-lm(UTSW_Median~Broad_Median,df_all)
df_all_lmcoef<-coef(df_all_lm)
df_all_lmcoef<-data.frame(Slope=round(df_all_lmcoef[[2]],3))
df_all_lmcoef$Slope<-paste0("Slope=",df_all_lmcoef$Slope)
df_all_lmcoef=cbind(UTSW_Median=1,Broad_Median=1,df_all_lmcoef)

#creating correlation plot
ggplot(df_all, aes(x=UTSW_Median,Broad_Median)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x,color="black") + stat_regline_equation(aes(label = ..rr.label..)) + labs(title=paste0("Correlation between Broad Data Set and UTSW Using the\n",filetype_sp," File"),x=paste0(comp1," Median"), y=paste0(comp2," Median")) + geom_abs_text(data=df_all_lmcoef,mapping = aes(label = Slope),color="black",size=3.8,xpos=0.08,ypos=0.88)
ggsave(paste0("CorrelationPlot_AllData_",filetype,".png"), type = "cairo",width=10,height=7)

#summarizing by compound
#calculate lm for each broad sample
cmpd_lm<-dlply(df_all,"Metadata_pert_iname",function(df) lm(UTSW_Median~Broad_Median,data=df))
cmpd_lm_coef<-ldply(cmpd_lm,coef)
names(cmpd_lm_coef)[names(cmpd_lm_coef)=="Broad_Median"]<-"Slope"
cmpd_lm_coef[,2]=c(1:nrow(cmpd_lm_coef))
names(cmpd_lm_coef)[names(cmpd_lm_coef)=="(Intercept)"]<-"UTSW_Median"
cmpd_lm_coef=cbind(cmpd_lm_coef,Broad_Median=c(1:nrow(cmpd_lm_coef)))
cmpd_lm_coef[,3]=round(cmpd_lm_coef$Slope,3)
cmpd_lm_coef$Slope<-ldply(paste0("Slope=",cmpd_lm_coef$Slope))
colnames(cmpd_lm_coef[,3])<-"Slope"

ggplot(df_all, aes(x=UTSW_Median,y=Broad_Median,group=Metadata_pert_iname,color=Metadata_pert_iname)) + geom_point() + geom_smooth(method='lm',formula=y~x) + facet_wrap(~Metadata_pert_iname,scales="free") + labs(title=paste0("Correlation between Broad Data Set and Ours Using the\n",filetype_sp," File"),x=paste0(comp1," Median"), y=paste0(comp2," Median")) + stat_poly_eq(label.x="left",label.y="top",size=2.5,color="black") + theme(strip.text = element_text(size = 6.2),axis.text=element_text(size=5)) + geom_abs_text(data=cmpd_lm_coef,mapping = aes(label = Slope),color="black",size=2.5,xpos=0.7,ypos=0.18)
ggsave(paste0("CorrelationPlot_ByCmpd_",filetype,".png"), type = "cairo")