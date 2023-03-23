#!/usr/bin/env Rscript
#Create a correlation plot between two dataframes

#load in data
library(data.table)
library(reshape)
library(ggplot2)
library(ggpmisc)
library(plyr)
library(zplyr)

args=commandArgs(trailingOnly=TRUE)

#read in the data frames

plate1=args[1]
plate2=args[2]
filetype=args[3]
jump_compounds=args[4]

plate1=sub("\\/.*","",plate1)

if(grepl("norm",filetype,fixed=TRUE)==TRUE){
    if(grepl("feat",filetype,fixed=TRUE)==TRUE){
        filename="Feature-Normalized"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
        mydf=paste0(plate1,"/",plate1,"_normalized.csv.gz")
        compdf=paste0(plate2,"/",basename(plate2),"_normalized.csv.gz")
    } else if(grepl("neg",filetype,fixed=TRUE)==TRUE){
        filename="Negative-Control-Normalized"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
        mydf=paste0(plate1,"/",plate1,"_normalized_negcon.csv.gz")
        compdf=paste0(plate2,"/",basename(plate2),"_normalized_negcon.csv.gz")
    }
} else if(grepl("sel",filetype,fixed=TRUE)==TRUE){
    if(grepl("feat",filetype,fixed=TRUE)==TRUE){
        filename="Feature-Normalized-Selected-Features"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
        #mydf=paste0(plate1,"/",plate1,"_normalized_feature_select_batch.csv.gz")
        #compdf=paste0(plate2,"/",basename(plate2),"_normalized_feature_select_batch.csv.gz")
    } else if(grepl("neg",filetype,fixed=TRUE)==TRUE){
        filename="Negative-Control-Normalized-Selected-Features"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
        #mydf=paste0(plate1,"/",plate1,"_normalized_feature_select_negcon_batch.csv.gz")
        #compdf=paste0(plate2,"/",basename(plate2),"_normalized_feature_select_negcon_batch.csv.gz")
    }
} else {
    print("There is no matching file")
}

mydf=plate1
compdf=plate2
mydf=fread(mydf)
compdf=fread(compdf)
mydf_plateid=unique(mydf$Metadata_Plate)
mydf_plateid=args[5]
mydf_plateid_nosp=gsub(" ","-",mydf_plateid,fixed=TRUE)
compdf_plateid=unique(compdf$Metadata_Plate)
compounds<-c("NVS-PAK1-1","aloxistatin","FK-866","AMG900","LY2109761","dexamethasone","quinidine","TC-S-7004","DMSO")

reshape_data<-function(df){
    df_new<-df[,grep("Cells",colnames(df))[[1]]:ncol(df)]
    meta=df[,1:6]
    df<-cbind(meta,df_new)
    df_melt<-melt(df)
    names(df_melt)[names(df_melt)=="variable"]<-"Measurement"
    names(df_melt)[names(df_melt)=="value"]<-"Median"
    return(df_melt)
}

mydf<-reshape_data(mydf)
mydf_feat<-unique(mydf$Measurement)
names(mydf)[names(mydf)=="Median"]<-"UTSW_Median"
compdf<-reshape_data(compdf)
compdf_feat<-unique(compdf$Measurement)
names(compdf)[names(compdf)=="Median"]<-"Broad_Median"

feats<-intersect(mydf_feat,compdf_feat)

subset_data<-function(df,compounds,jump_cmpds,features){
    #subset data for features
    df_featsub=df[df$Measurement %in% features,]
    #if jump_compounds is set to true, subset for compounds, else don't do it
    if(jump_cmpds==TRUE){
        df_sub=df_featsub[df_featsub$Metadata_pert_iname %in% compounds,]
        return(df_sub)
    } else {
        return(df_featsub)
    }
}

mydf_t<-cast(df_sub,Metadata_pert_iname+Metadata_solvent+Metadata_Plate+Metadata_Well+Metadata_pert_type+Metadata_control_type~Measurement,value="UTSW_Median")

mydf_sub<-subset_data(mydf,compounds,jump_compounds,feats)
compdf_sub<-subset_data(compdf,compounds,jump_compounds,feats)

#merge data sheets by type
df_all<-merge(mydf_sub,compdf_sub,by=c("Metadata_pert_iname","Measurement"))
#calculate the slope and add to the plot
#df_all_lm<-lm(UTSW_Median~Broad_Median,df_all)
#df_all_lmcoef<-coef(df_all_lm)
#df_all_lmcoef<-data.frame(Slope=round(df_all_lmcoef[[2]],3))
#df_all_lmcoef$Slope<-paste0("Slope=",df_all_lmcoef$Slope)
#df_all_lmcoef=cbind(UTSW_Median=1,Broad_Median=1,df_all_lmcoef)

#creating correlation plot by Metadata_Compound
ggplot(df_all, aes(x=UTSW_Median,Broad_Median)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x, colour="black") + stat_poly_eq() + labs(title=paste0("Correlation between Broad and UTSW Data Sets\nAfter Data Normalization with ",filename_sp," File\nPlate ",mydf_plateid),x="UTSW Normalized Feature Value",y="Broad Normalized Feature Value")

# + geom_abs_text(data=df_all_lmcoef,mapping = aes(label = Slope),color="black",size=3.8,xpos=0.1,ypos=0.9)
ggsave(paste0("PyCyto_CorrelationAll_",mydf_plateid_nosp,"_",filename,".png"), type = "cairo")

#summarizing by compound
#calculate lm for each broad sample
if(jump_compounds==TRUE){
    cmpd_lm<-dlply(df_all,"Metadata_pert_iname",function(df) lm(UTSW_Median~Broad_Median,data=df))
    cmpd_lm_coef<-ldply(cmpd_lm,coef)
    names(cmpd_lm_coef)[names(cmpd_lm_coef)=="Broad_Median"]<-"Slope"
    cmpd_lm_coef[,2]=c(1:nrow(cmpd_lm_coef))
    names(cmpd_lm_coef)[names(cmpd_lm_coef)=="(Intercept)"]<-"UTSW_Median"
    cmpd_lm_coef=cbind(cmpd_lm_coef,Broad_Median=c(1:nrow(cmpd_lm_coef)))
    cmpd_lm_coef[,3]=round(cmpd_lm_coef$Slope,3)
    cmpd_lm_coef$Slope<-ldply(paste0("Slope=",cmpd_lm_coef$Slope))
    colnames(cmpd_lm_coef[,3])<-"Slope"

    ggplot(df_all, aes(x=UTSW_Median,y=Broad_Median,group=Metadata_pert_iname,color=Metadata_pert_iname)) + geom_point() + geom_smooth(method='lm',formula=y~x) + facet_wrap(~Metadata_pert_iname,scales="free") + labs(title=paste0("Correlation between Broad and UTSW Data Sets\nAfter Data Normalization with ",filename_sp," File\nPlate ",mydf_plateid),x="UTSW Normalized Feature Value",y="Broad Normalized Feature Value") + stat_poly_eq(label.x="left",label.y="top",size=2.5,color="black") + theme(strip.text = element_text(size = 6.2)) + geom_abs_text(data=cmpd_lm_coef,mapping = aes(label = Slope),color="black",size=2.5,xpos=0.7,ypos=0.18)
    ggsave(paste0("PyCyto_CorrelationByCmpd_",mydf_plateid_nosp,"_",filename,".png"), type = "cairo")
}