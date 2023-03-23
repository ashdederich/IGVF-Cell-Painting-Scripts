#!/usr/bin/env Rscript
#Create a correlation plot between two dataframes

#load in data
library(plyr)
library(dplyr)
library(data.table)
library(reshape)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(zplyr)

args=commandArgs(trailingOnly=TRUE)

#read in the data frames
mydf=args[1]
comparisondf=args[2]
metadata=fread(args[3])
filetype=args[4]
batch=args[5]

#set up title and measurements files to detect feature-summarized or negcon summarized
if(grepl("feat",filetype,fixed=TRUE)==TRUE){
    title="CP-Output-with-Feature-Normalized-Features"
    title_sp=gsub("-"," ",title,fixed=TRUE)
    plottitle="CP-Output-with-Feature-Normalized-Features"
    mydf_meas_file=paste0(mydf,"/",mydf,"_normalized_feature_select_batch.csv.gz")
    #compdf_meas_file=paste0(comparisondf,"/",basename(comparisondf),"_normalized_feature_select_batch.csv.gz")
    mydf=paste0(mydf,"/",mydf,".csv.gz")
    #compdf=paste0(comparisondf,"/",basename(comparisondf),".csv.gz")
} else if(grepl("neg",filetype,fixed=TRUE)==TRUE){
    title="CP-Output-with-NegCon-Normalized-Features"
    title_sp=gsub("-"," ",title,fixed=TRUE)
    plottitle="CP-Output-with-NegCon-Normalized-Features"
    mydf_meas_file=paste0(mydf,"/",mydf,"_normalized_feature_select_negcon_batch.csv.gz")
    #compdf_meas_file=paste0(comparisondf,"/",basename(comparisondf),"_normalized_feature_select_negcon_batch.csv.gz")
    mydf=paste0(mydf,"/",mydf,".csv.gz")
    #compdf=paste0(comparisondf,"/",basename(comparisondf),".csv.gz")
} else if(grepl("cp",filetype,fixed=TRUE)==TRUE){
    title="CP-Output"
    title_sp=gsub("-"," ",title,fixed=TRUE)
    plottitle=title
    mydf=paste0(mydf,"/",basename(mydf),".csv.gz")
    #compdf=paste0(comparisondf,"/",basename(comparisondf),".csv.gz")
}else{
    print("There is no matching filetype")
}

compdf=comparisondf
compdf_meas_file=args[6]

#getting intersection of features
reshape_data<-function(file){
    df<-fread(file)
    df_new<-df[,grep("Cells",colnames(df))[[1]]:ncol(df)]
    df<-cbind(Metadata_pert_iname=df$Metadata_pert_iname,df_new)
    df_melt<-melt(df)
    df_feat<-df_melt$variable
    return(df_feat)
}

if(grepl("feat",filetype,fixed=TRUE)==TRUE | grepl("neg",filetype,fixed=TRUE)==TRUE){
    mydf_feat<-reshape_data(mydf_meas_file)
    comp_feat<-reshape_data(compdf_meas_file)
    features<-intersect(mydf_feat,comp_feat)
}

compounds<-c("NVS-PAK1-1","aloxistatin","FK-866","AMG900","LY2109761","dexamethasone","quinidine","TC-S-7004","DMSO")

get_metadata<-function(df,metadata_file,batchid){
    df<-fread(df)
    plateid=as.list(unique(df$Metadata_Plate))
    barcode=fread(paste0("../../../../metadata/platemaps/",batchid,"/barcode_platemap.csv"))
    barcode=barcode %>% filter(Assay_Plate_Barcode==plateid) %>% pull(var=Plate_Map_Name)
    platemap=fread(paste0("../../../../metadata/platemaps/",batchid,"/platemap/",barcode,".txt"))
    names(platemap)[names(platemap)=="well_position"]<-"Metadata_Well"
    platemap=data.frame(Metadata_Well=platemap$Metadata_Well,pert_iname=platemap$pert_iname) #only get well_position and pert_iname information
    metadata=inner_join(platemap,metadata_file,by="pert_iname")
    metadata=data.frame(Metadata_Well=metadata$Metadata_Well,Metadata_pert_iname=metadata$pert_iname)
    df<-merge(metadata,df,by="Metadata_Well")
    return(df)
}

platemap_df_reshape<-function(df){
    df<-data.frame(Metadata_Well=df$well_position,pert_iname=df$pert_iname)
    return(df)
}

#a version of get_metadata where each plate has its own unique platemap
#here I will need to read each platemap and join them row-wise
#I will need to makesure that they have the same col names in the same order
get_metadata_multipleplates<-function(df,metadata_file,batchid){
    df<-fread(df)
    plateid=as.list(unique(df$Metadata_Plate))
    barcode=fread(paste0("../../../../metadata/platemaps/",batchid,"/barcode_platemap.csv"))
    barcode=lapply(plateid,function(x) barcode %>% filter(Assay_Plate_Barcode==x) %>% pull(var=Plate_Map_Name))
    platemap=lapply(barcode, function(x) fread(paste0("../../../../metadata/platemaps/",my_batch,"/platemap/",x,".txt")))
    platemap=lapply(platemap,platemap_df_reshape)
    platemap<-rbind(platemap[[1]],platemap[[2]],platemap[[3]],platemap[[4]])
    metadata=inner_join(platemap,metadata_file,by="pert_iname")
    metadata=data.frame(Metadata_Well=metadata$Metadata_Well,Metadata_pert_iname=metadata$pert_iname)
    df<-merge(metadata,df,by="Metadata_Well")
    return(df)
}

#add metadata information - if each plate has a different metadata file, change the code to supply different metadata files
mydf_new=get_metadata(mydf,metadata,batch)
plateid=unique(mydf_new$Metadata_Plate)
compdf_new=get_metadata(compdf,metadata,batch)

#change data frames from short and wide to tall and skinny - My Data
subset_and_melt<-function(df,cmpds,feats=features){
    df_subset=df[df$Metadata_pert_iname %in% cmpds,] # only get JUMP cmpds
    df_melt<-melt(df_subset)
    names(df_melt)[names(df_melt)=="variable"]<-"Measurement"
    names(df_melt)[names(df_melt)=="value"]<-"Median"
    if(grepl("feat",filetype,fixed=TRUE)==TRUE | grepl("neg",filetype,fixed=TRUE)==TRUE){
        df_melt_sub=df_melt[df_melt$Measurement %in% feats,] #only get measurements that pycytominer selected as variable
        return(df_melt_sub)
    }else{
        return(df_melt)
    }
}

mydf_formerge<-subset_and_melt(mydf_new,compounds)
names(mydf_formerge)[names(mydf_formerge)=="Median"]<-"UTSW_Median"
compdf_formerge<-subset_and_melt(compdf_new,compounds)
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
ggplot(df_all, aes(x=UTSW_Median,Broad_Median)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x,color="black") + stat_regline_equation(aes(label = ..rr.label..)) + labs(title=paste0("Correlation between Broad and UTSW Data Sets Using\n",title_sp,"\nPlate ", plateid),x="UTSW Feature Value", y="Broad Feature Value") + geom_abs_text(data=df_all_lmcoef,mapping = aes(label = Slope),color="black",size=3.8,xpos=0.08,ypos=0.84)
ggsave(paste0("CorrelationAllData_20220914_v_jump",plottitle,".png"), type = "cairo",width=10,height=7)

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

ggplot(df_all, aes(x=UTSW_Median,y=Broad_Median,group=Metadata_pert_iname,color=Metadata_pert_iname)) + geom_point() + geom_smooth(method='lm',formula=y~x) + facet_wrap(~Metadata_pert_iname,scales="free") + labs(title=paste0("Correlation between Broad and UTSW Data Sets Using\n",title_sp,"\nPlate ", plateid),x="UTSW Feature Value", y="Broad Feature Value") + stat_poly_eq(label.x="left",label.y="top",size=2.5,color="black") + theme(strip.text = element_text(size = 6.2),axis.text=element_text(size=5)) + geom_abs_text(data=cmpd_lm_coef,mapping = aes(label = Slope),color="black",size=2.5,xpos=0.7,ypos=0.18)
ggsave(paste0("CorrelationByCmpd_20220914_v_jump","_",plottitle,".png"), type = "cairo")