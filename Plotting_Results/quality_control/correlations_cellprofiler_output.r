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
filetype=args[4]
comp1<-args[5]
comp2<-args[6]

#set up title and measurements files to detect feature-summarized or negcon summarized
mydf_meas<-sub("\\.csv.gz.*","",mydf)

if(grepl("feature",filetype,fixed=TRUE)==TRUE){
    title="CP-Output-with-Feature-Normalized-Features"
    title_sp=gsub("-"," ",title,fixed=TRUE)
    plottitle="CP-Ouput-Feature-Normalized"
    mydf_meas_file=paste0(mydf_meas,"/",mydf_meas,"_normalized_feature_select_batch.csv.gz")
    compdf_meas_file=paste0(comparisondf,"/",basename(comparisondf),"_normalized_feature_select_batch.csv.gz")
    mydf=paste0(mydf,"/",mydf,".csv.gz")
    compdf=paste0(comparisondf,"/",basename(comparisondf),".csv.gz")
} else if(grepl("neg",filetype,fixed=TRUE)==TRUE){
    title="CP-Output-with-NegCon-Normalized-Features"
    title_sp=gsub("-"," ",title,fixed=TRUE)
    plottitle="CP-Output-NegCon-Normalized"
    mydf_meas_file=paste0(mydf_meas,"_normalized_feature_select_negcon_batch.csv.gz")
    compdf_meas_file=paste0(comparisondf,"/",basename(comparisondf),"_normalized_feature_select_negcon_batch.csv.gz")
    mydf=paste0(mydf,"/",mydf,".csv.gz")
    compdf=paste0(comparisondf,"/",basename(comparisondf),".csv.gz")
} else {
    print("There is no matching filetype")
}

comp1=gsub("-"," ", comp1,fixed=TRUE)
comp2=gsub("-"," ", comp2,fixed=TRUE)

#getting intersection of features
reshape_data<-function(file){
    df<-fread(file)
    df_new<-df[,grep("Cells",colnames(df))[[1]]:ncol(df)]
    df<-cbind(Metadata_pert_iname=df$Metadata_pert_iname,df_new)
    df_melt<-melt(df)
    df_feat<-df_melt$variable
    return(df_feat)
}

mydf_feat<-reshape_data(mydf_meas_file)
comp_feat<-reshape_data(compdf_meas_file)
features<-intersect(mydf_feat,comp_feat)
compounds<-c("NVS-PAK1-1","aloxistatin","FK-866","AMG900","LY2109761","dexamethasone","quinidine","TC-S-7004","DMSO")

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
#cmpd_pert_well<-data.frame(Metadata_pert_iname=cmpd_df_new$Metadata_pert_iname,Metadata_Well=cmpd_df_new$Metadata_Well)
#mydf=inner_join(cmpd_pert_well,mydf,by="Metadata_Well")

#getting metadata for broad
broad_batchid=basename(dirname(dirname(comparisondf)))
comparisondf=fread(comparisondf)
broad_plateid=unique(comparisondf$Metadata_Plate)[1]
broad_barcode=fread(paste0("../../metadata/platemaps/",broad_batchid,"/barcode_platemap.csv"))
broad_barcode=broad_barcode %>% filter(Assay_Plate_Barcode==broad_plateid) %>% pull(var=Plate_Map_Name)
broad_platemap=fread(paste0("../../metadata/platemaps/",broad_batchid,"/platemap/",broad_barcode,".txt"))
names(broad_platemap)[names(broad_platemap)=="well_position"]<-"Metadata_Well"
broad_platemap=broad_platemap[,1:2] #only get well_position and broad_sample information
broad_metadata=inner_join(broad_platemap,metadata,by="broad_sample")
broad_metadata=data.frame(Metadata_Well=broad_metadata$Metadata_Well,Metadata_pert_iname=broad_metadata$pert_iname,Metadata_broad_sample=broad_metadata$broad_sample)
comparisondf<-merge(broad_metadata,comparisondf,by="Metadata_Well")

#change data frames from short and wide to tall and skinny - My Data
mydf=fread(mydf)
plateid=unique(mydf$Metadata_Plate)
mydf_new=merge(broad_metadata,mydf,by="Metadata_Well")

subset_and_melt<-function(df,compounds=compounds,features=features){
    df$Metadata_broad_sample[which(df$Metadata_pert_iname=="")]<-"DMSO"
    df_subset=df[df$Metadata_pert_iname %in% compounds,] # only get JUMP cmpds
    df_melt<-melt(df_subset)
    names(df_melt)[names(df_melt)=="variable"]<-"Measurement"
    names(df_melt)[names(df_melt)=="value"]<-"Median"
    df_melt=df_melt[df_melt$Measurement %in% features,] #only get measurements that pycytominer selected as variable
    return(df_melt)
}

mydf_formerge<-subset_and_melt(mydf_new,compounds,features)
names(mydf_formerge)[names(mydf_formerge)=="Median"]<-"UTSW_Median"
compdf_formerge<-subset_and_melt(comparisondf,compounds,features)
names(compdf_formerge)[names(compdf_formerge)=="Median"]<-"Broad_Median"

#merge data sheets by type
df_all<-inner_join(x=mydf_formerge,y=compdf_formerge,by=c("Metadata_pert_iname","Measurement"))
#calculate the slope and add to the plot
df_all_lm<-lm(UTSW_Median~Broad_Median,df_all)
df_all_lmcoef<-coef(df_all_lm)
df_all_lmcoef<-data.frame(Slope=round(df_all_lmcoef[[2]],3))
df_all_lmcoef$Slope<-paste0("Slope=",df_all_lmcoef$Slope)
df_all_lmcoef=cbind(UTSW_Median=1,Broad_Median=1,df_all_lmcoef)

#creating correlation plot
ggplot(df_all, aes(x=UTSW_Median,Broad_Median)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x,color="black") + stat_regline_equation(aes(label = ..rr.label..)) + labs(title=paste0("Correlation between Broad and UTSW Data Sets Using\n",title_sp,"\nPlate ", plateid),x=paste0(comp1," Normalized Feature Value"), y=paste0(comp2," Normalized Feature Value")) + geom_abs_text(data=df_all_lmcoef,mapping = aes(label = Slope),color="black",size=3.8,xpos=0.08,ypos=0.84)
ggsave(paste0("CorrelationAllData_",plateid,"_",plottitle,".png"), type = "cairo",width=10,height=7)

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

ggplot(df_all, aes(x=UTSW_Median,y=Broad_Median,group=Metadata_pert_iname,color=Metadata_pert_iname)) + geom_point() + geom_smooth(method='lm',formula=y~x) + facet_wrap(~Metadata_pert_iname,scales="free") + labs(title=paste0("Correlation between Broad and UTSW Data Sets Using\n",title_sp,"\nPlate ", plateid),x=paste0(comp1," Normalized Feature Value"), y=paste0(comp2," Normalized Feature Value")) + stat_poly_eq(label.x="left",label.y="top",size=2.5,color="black") + theme(strip.text = element_text(size = 6.2),axis.text=element_text(size=5)) + geom_abs_text(data=cmpd_lm_coef,mapping = aes(label = Slope),color="black",size=2.5,xpos=0.7,ypos=0.18)
ggsave(paste0("CorrelationByCmpd_",plateid,"_",plottitle,".png"), type = "cairo")