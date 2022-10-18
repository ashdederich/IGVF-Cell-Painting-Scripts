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
mydf=fread(args[1])
comparisondf=fread(args[2])
filetype=args[3]
filetype_sp=gsub("-"," ", filetype,fixed=TRUE)
jump_compounds=args[4]
compounds<-c("NVS-PAK1-1","aloxistatin","FK-866","AMG900","LY2109761","dexamethasone","quinidine","TC-S-7004","DMSO")

melt_df<-function(df,jump_compounds=TRUE) {
    df_new<-df[,grep("Cells",colnames(df))[[1]]:ncol(df)]
    df_new<-cbind(df$Metadata_pert_iname,df_new)
    colnames(df_new)[1]<-"Metadata_pert_iname"
    df_new_melt<-melt(df_new)
    names(df_new_melt)[names(df_new_melt)=="variable"]<-"Measurement"
    names(df_new_melt)[names(df_new_melt)=="value"]<-"Median"
    if(jump_compounds==TRUE) {
        df_new_melt=df_new_melt[df_new_melt$Metadata_pert_iname %in% compounds,] # only get JUMP cmpds
    }
    return(df_new_melt)
}

mydf_new_melt<-melt_df(mydf)
names(mydf_new_melt)[names(mydf_new_melt)=="Median"]<-"UTSW_Median"

comparisondf_new_melt<-melt_df(comparisondf)
names(comparisondf_new_melt)[names(comparisondf_new_melt)=="Median"]<-"Broad_Median"

#merge data sheets by type
df_all<-merge(mydf_new_melt,comparisondf_new_melt,by=c("Metadata_pert_iname","Measurement"))
#calculate the slope and add to the plot
df_all_lm<-lm(UTSW_Median~Broad_Median,df_all)
df_all_lmcoef<-coef(df_all_lm)
df_all_lmcoef<-data.frame(Slope=round(df_all_lmcoef[[2]],3))
df_all_lmcoef$Slope<-paste0("Slope=",df_all_lmcoef$Slope)
df_all_lmcoef=cbind(UTSW_Median=1,Broad_Median=1,df_all_lmcoef)

#creating correlation plot by Metadata_Compound
ggplot(df_all, aes(x=UTSW_Median,Broad_Median)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x, colour="black") + stat_poly_eq() + labs(title=paste0("Correlation between Broad Data Set and UTSW\nAfter PyCytominer Using the ",filetype_sp," File")) + geom_abs_text(data=df_all_lmcoef,mapping = aes(label = Slope),color="black",size=3.8,xpos=0.1,ypos=0.9)
ggsave(paste0("CorrelationPlot_AllData_",filetype,".png"), type = "cairo")

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

ggplot(df_all, aes(x=UTSW_Median,y=Broad_Median,group=Metadata_pert_iname,color=Metadata_pert_iname)) + geom_point() + geom_smooth(method='lm',formula=y~x) + facet_wrap(~Metadata_pert_iname,scales="free") + labs(title=paste0("Correlation between Broad Data Set and Ours\nAfter PyCytominer Using the ",filetype_sp," File")) + stat_poly_eq(label.x="left",label.y="top",size=2.5,color="black") + theme(strip.text = element_text(size = 6.2)) + geom_abs_text(data=cmpd_lm_coef,mapping = aes(label = Slope),color="black",size=2.5,xpos=0.7,ypos=0.18)
ggsave(paste0("CorrelationPlot_ByCmpd_",filetype,".png"), type = "cairo")