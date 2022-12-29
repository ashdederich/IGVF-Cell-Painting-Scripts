#!/usr/bin/env Rscript
#Calculate robust MAD for all features and then find the correlation between median/mad for two files
#This file requires one input file input file to aggregate: this is the output file from CellProfiler.

args=commandArgs(trailingOnly=TRUE)
library(data.table)
library(reshape)
library(tools)
library(ggplot2)
library(zplyr)

plate1=args[1]
plate2=args[2]
filetype=args[3]
dmso=args[4]
outliers=args[5]

#adding default option for dmso boolean if argument not given
if(dmso==""){
    dmso=FALSE
}

if(outliers==""){
    outliers=FALSE
}

plate1=sub("\\/.*","",plate1)

#create an if outlier or if cp ==true, then join with metadata
if(grepl("cp",filetype,fixed=TRUE)==TRUE | outliers==TRUE){
    if(dmso==TRUE){
        filename="CellProfiler-DMSO-Only"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
    }else{
        filename="CellProfiler-Output"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
    }
    #read in dataframes
    my_cpdf=fread(paste0(plate1,"/",plate1,".csv.gz"))
    comp_cpdf=fread(paste0(plate2,"/",basename(plate2),".csv.gz"))
    #read in metadata
    metadata=fread("../../metadata/external_metadata/JUMP-Target-1_compound_metadata.tsv")
    broad_batchid=basename(dirname(plate2))
    broad_plateid=unique(comp_cpdf$Metadata_Plate)[1]
    broad_barcode=fread(paste0("../../metadata/platemaps/",broad_batchid,"/barcode_platemap.csv"))
    broad_barcode=broad_barcode %>% filter(Assay_Plate_Barcode==broad_plateid) %>% pull(var=Plate_Map_Name)
    broad_platemap=fread(paste0("../../metadata/platemaps/",broad_batchid,"/platemap/",broad_barcode,".txt"))
    names(broad_platemap)[names(broad_platemap)=="well_position"]<-"Metadata_Well"
    broad_platemap=broad_platemap[,1:2] #only get well_position and broad_sample information
    broad_metadata=inner_join(broad_platemap,metadata,by="broad_sample")
    broad_metadata=data.frame(Metadata_Well=broad_metadata$Metadata_Well,Metadata_pert_iname=broad_metadata$pert_iname,Metadata_broad_sample=broad_metadata$broad_sample)
    broad_metadata=subset(broad_metadata,select=c("Metadata_Well","Metadata_pert_iname"))
    #merge both dataframes with metadata
    my_cpdf<-merge(broad_metadata,my_cpdf,by="Metadata_Well")
    comp_cpdf<-merge(broad_metadata,comp_cpdf,by="Metadata_Well")
}

#separate from cp/outlier if statement, becuase we need both files if outlier is TRUE
#if pycytominer is needed, then select these files
if(grepl("py",filetype,fixed=TRUE)==TRUE && grepl("feat",filetype,fixed=TRUE)==TRUE){
    if(dmso==TRUE){
        filename="PyCytominer-Feature-Normalized-DMSO-Only"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
    }else{
        filename="PyCytominer-Feature-Normalized"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
    }
    mydf=paste0(plate1,"/",plate1,"_normalized_feature_select_batch.csv.gz")
    compdf=paste0(plate2,"/",basename(plate2),"_normalized_feature_select_batch.csv.gz")
} else if(grepl("py",filetype,fixed=TRUE)==TRUE && grepl("neg",filetype,fixed=TRUE)==TRUE){
    if(dmso==TRUE){
        filename="PyCytominer-NegCon-Normalized-DMSO-Only"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
    }else{
        filename="PyCytominer-NegCon-Normalized"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
    }
    mydf=paste0(plate1,"/",plate1,"_normalized_feature_select_negcon_batch.csv.gz")
    compdf=paste0(plate2,"/",basename(plate2),"_normalized_feature_select_negcon_batch.csv.gz")
} else {
    print("There is no matching file")
}

if(grepl("py",filetype,fixed=TRUE)==TRUE) {
    mydf=fread(mydf)
    compdf=fread(compdf)
}

#function to aggregate data, with a built-in if dmso==TRUE statement (no need to have if dmso==TRUE statement in later functions, then.
aggregate_data<-function(file){
    if(dmso==TRUE){
        if(grepl("py",filetype,fixed=TRUE)==TRUE){
            #get metadata wanted and join with measurements, then melt
            file_meta<-data.frame(Metadata_Plate=file$Metadata_Plate,Metadata_Well=file$Metadata_Well,Metadata_pert_iname=file$Metadata_pert_iname)
            col_start<-grep("Cells",colnames(file))[[1]]
            file_data<-file[,col_start:ncol(file)]
            file<-melt(data.frame(c(file_meta,file_data)),id.vars=c("Metadata_Plate","Metadata_Well","Metadata_pert_iname"))
        }else{
            file<-melt(file,id.vars=c("Metadata_Plate","Metadata_Well","Metadata_pert_iname"))
        }
        #subset the data to only get the DMSO data
        file<-file[which(file$Metadata_pert_iname=="DMSO"),]
    }else{
        file_meta<-data.frame(Metadata_Plate=file$Metadata_Plate,Metadata_Well=file$Metadata_Well)
        col_start<-grep("Cells",colnames(file))[[1]]
        file_data<-file[,col_start:ncol(file)]
        file<-melt(data.frame(c(file_meta,file_data)),id.vars=c("Metadata_Plate","Metadata_Well"))
    }
    #melt the df
    names(file)[names(file)=="variable"]<-"Measurement"
    names(file)[names(file)=="value"]<-"Median"
    plate<-unique(file$Metadata_Plate)
    #get median and mad and rename the grouping variables
    median_ofdata=aggregate(file$Median,by=list(file$Measurement),median)
    mad_ofdata=aggregate(file$Median,by=list(file$Measurement),mad)
    names(median_ofdata)[names(median_ofdata)=="Group.1"]<-"Measurement"
    names(median_ofdata)[names(median_ofdata)=="x"]<-"Measurement_Median"
    names(mad_ofdata)[names(mad_ofdata)=="Group.1"]<-"Measurement"
    names(mad_ofdata)[names(mad_ofdata)=="x"]<-"Measurement_MAD"

    data_med_mad<-merge(median_ofdata,mad_ofdata,by="Measurement")
    data_med_mad<-cbind(Metadata_Plate=rep(plate,nrow(data_med_mad)),data_med_mad)
    return(data_med_mad)
}

#for both cp or pycyto data
if(dmso==FALSE){
    mydf_agg<-aggregate_data(mydf)
    names(mydf_agg)[names(mydf_agg)=="Measurement_Median"]<-"UTSW_Median"
    names(mydf_agg)[names(mydf_agg)=="Measurement_MAD"]<-"UTSW_MAD"
}

if(dmso==FALSE){
    compdf_agg<-aggregate_data(compdf)
    names(compdf_agg)[names(compdf_agg)=="Measurement_Median"]<-"Broad_Median"
    names(compdf_agg)[names(compdf_agg)=="Measurement_MAD"]<-"Broad_MAD"
}

#take aggregated file
#find outliers
#subset aggregated file to remove outliers
#return the subsetted file
outlier_removal<-function(cp_df,py_df=NULL){
    cp_agg=aggregate_data(cp_df)
    py_agg=aggregate_data(py_df)
    #only get dmso values first if dmso is set to TRUE
    #remove outliers from CP file for each col - Median and MAD
    median_outliers<-cp_agg[,"Measurement_Median"] %in% boxplot.stats(cp_agg[,"Measurement_Median"])$out
    mad_outliers<-cp_agg[,"Measurement_MAD"] %in% boxplot.stats(cp_agg[,"Measurement_MAD"])$out
    #replace outliers with NAs
    cp_agg$Measurement_Median[median_outliers==TRUE]<-NA
    cp_agg$Measurement_MAD[mad_outliers==TRUE]<-NA
    #then, if filetype = py, we match those up to the wells and measurements retained and replace outliers with NA
    if(grepl("py",filetype,fixed=TRUE)==TRUE){
        #which measurements have an NA in cp_agg for median?
        cp_median_na<-as.vector(cp_agg$Measurement[which(is.na(cp_agg$Measurement_Median))])
        #which measurements have an NA in cp_agg for MAD?
        cp_mad_na<-as.vector(cp_agg$Measurement[which(is.na(cp_agg$Measurement_MAD))])
        #for those measurements in cp_agg with an NA, find those that match the measurements in py_agg and replace that data with an NA
        py_agg$Measurement_Median[which(py_agg$Measurement%in%cp_median_na)]<-NA
        py_agg$Measurement_MAD[which(py_agg$Measurement%in%cp_mad_na)]<-NA      
        return(py_agg)
    }else{
        return(cp_agg)
    }
}

#now set up if, else statement for cp or pycytominer files
#if dmso is false and outlier is true, then perform outlier_removal function
if(dmso==FALSE && outliers==TRUE){
    if(grepl("cp",filetype,fixed=TRUE)==TRUE){
        my_cpdf=outlier_removal(cp_df=my_cpdf)
        comp_cpdf=outlier_removal(cp_df=comp_cpdf)
    }else if(grepl("py",filetype,fixed=TRUE)==TRUE){
        mydf=outlier_removal(cp_df=my_cpdf,py_df=mydf)
        compdf=outlier_removal(cp_df=comp_cpdf,py_df=compdf)
    }
}

calc_spread<-function(firstdf,compdf){
    #merge the two dataframes
    merged<-merge(firstdf,compdf,by=c("Metadata_Plate","Measurement"))
    #median plots
    if(grepl("cp",filetype,fixed=TRUE)==TRUE){
        ggplot(merged) + geom_point(aes(x=1:nrow(merged),y=UTSW_Median,color="UTSW_Median"),show.legend=TRUE) + geom_point(aes(x=1:nrow(merged),y=Broad_Median,color="Broad_Median"),show.legend=TRUE) + labs(title=paste0("Median DMSO Values per Well, CellProfiler Output\nPlate ",unique(merged$Metadata_Plate)),x="All DMSO Measurements",y="Median Feature Value",fill="Group") + guides(color=guide_legend("Group")) + geom_line(aes(x=1:nrow(merged),y=quantile(merged[,"UTSW_Median"],prob=0.95,na.rm=TRUE),color="UTSW_Median"),size=0.8) + geom_line(aes(x=1:nrow(merged),y=quantile(merged[,"Broad_Median"],prob=0.95,na.rm=TRUE),color="Broad_Median"),size=0.8) + scale_y_continuous(breaks=seq(round(min(merged$UTSW_Median,na.rm=TRUE),0),round(max(merged$UTSW_Median,na.rm=TRUE)+700,0),by=500))
    }else{
        max_utsw<-max(merged$UTSW_Median,na.rm=TRUE)
        max_broad<-max(merged$Broad_Median)
        if(max_utsw>max_broad){
            plot_max=max_utsw
        }else{
            plot_max=max_broad
        }
        min_utsw=min(merged$UTSW_Median,na.rm=TRUE)
        min_broad=min(merged$Broad_Median,na.rm=TRUE)
        if(min_utsw<min_broad){
            plot_min=min_utsw
        }else{
            plot_min=min_broad
        }
        ggplot(merged) + geom_point(aes(x=1:nrow(merged),y=UTSW_Median,color="UTSW_Median"),show.legend=TRUE) + geom_point(aes(x=1:nrow(merged),y=Broad_Median,color="Broad_Median"),show.legend=TRUE) + labs(title=paste0("Median Values per Well\n",filename_sp,"\nPlate ",unique(merged$Metadata_Plate)),x="All DMSO Measurements",y="Median Feature Value",fill="Group") + guides(color=guide_legend("Group")) + geom_line(aes(x=1:nrow(merged),y=quantile(merged[,"UTSW_Median"],prob=0.95,na.rm=TRUE),color="UTSW_Median"),size=0.8) + geom_line(aes(x=1:nrow(merged),y=quantile(merged[,"Broad_Median"],prob=0.95,na.rm=TRUE),color="Broad_Median"),size=0.8) + scale_y_continuous(breaks=seq(round(plot_min-0.1,1),round(plot_max+0.25,1),by=0.1))
    }
    ggsave(paste0("MedianValues_",filename,"_",unique(merged$Metadata_Plate),".png"), type = "cairo")
    #mad plot
    if(grepl("cp",filetype,fixed=TRUE)==TRUE){
        ggplot(merged) + geom_point(aes(x=1:nrow(merged),y=UTSW_MAD,color="UTSW_MAD"),show.legend=TRUE) + geom_point(aes(x=1:nrow(merged),y=Broad_MAD,color="Broad_MAD"),show.legend=TRUE) + labs(title=paste0("MAD DMSO Values per Well, CellProfiler Output\nPlate ",unique(merged$Metadata_Plate)),x="All DMSO Measurements",y="MAD Feature Value",fill="Group") + guides(color=guide_legend("Group")) + geom_line(aes(x=1:nrow(merged),y=quantile(merged[,"UTSW_MAD"],prob=0.95,na.rm=TRUE),color="UTSW_MAD"),size=0.8) + geom_line(aes(x=1:nrow(merged),y=quantile(merged[,"Broad_MAD"],prob=0.95,na.rm=TRUE),color="Broad_MAD"),size=0.8) + scale_y_continuous(breaks=seq(round(min(merged$UTSW_MAD,na.rm=TRUE),0),round(max(merged$UTSW_MAD,na.rm=TRUE)+300,0),by=50))
    }else{
        max_utsw<-max(merged$UTSW_MAD,na.rm=TRUE)
        max_broad<-max(merged$Broad_MAD)
        if(max_utsw>max_broad){
            plot_max=max_utsw
        }else{
            plot_max=max_broad
        }
        min_utsw=min(merged$UTSW_MAD,na.rm=TRUE)
        min_broad=min(merged$Broad_MAD,na.rm=TRUE)
        if(min_utsw<min_broad){
            plot_min=min_utsw
        }else{
            plot_min=min_broad
        }        
        ggplot(merged) + geom_point(aes(x=1:nrow(merged),y=UTSW_MAD,color="UTSW_MAD"),show.legend=TRUE) + geom_point(aes(x=1:nrow(merged),y=Broad_MAD,color="Broad_MAD"),show.legend=TRUE) + labs(title=paste0("MAD Values per Well\n",filename_sp,"\nPlate ",unique(merged$Metadata_Plate)),x="All DMSO Measurements",y="MAD Feature Value",fill="Group") + guides(color=guide_legend("Group")) + geom_line(aes(x=1:nrow(merged),y=quantile(merged[,"UTSW_MAD"],prob=0.95,na.rm=TRUE),color="UTSW_MAD"),size=0.8) + geom_line(aes(x=1:nrow(merged),y=quantile(merged[,"Broad_MAD"],prob=0.95,na.rm=TRUE),color="Broad_MAD"),size=0.8) + scale_y_continuous(breaks=seq(round(plot_min-0.1,1),round(plot_max+0.75,1),by=0.1))
    }
    ggsave(paste0("MADValues_",filename,"_",unique(merged$Metadata_Plate),".png"), type = "cairo")
}

if(dmso==FALSE){
    calc_spread(mydf_agg,compdf_agg)
}

#will need to remove this
if(dmso==TRUE){
    #may get rid of this OR add a conditional to the if statement - if outliers==FALSE - or add an if,else statement here so if outliers==FALSE, then it subsets the data for dmso vals only.
    compdf<-compdf[which(compdf$Metadata_pert_iname=="DMSO"),]
    mydf<-mydf[which(mydf$Metadata_pert_iname=="DMSO"),]
    mydf<-aggregate_data(mydf)
    names(mydf)[names(mydf)=="Measurement_Median"]<-"UTSW_Median"
    names(mydf)[names(mydf)=="Measurement_MAD"]<-"UTSW_MAD"
    compdf<-aggregate_data(compdf)
    names(compdf)[names(compdf)=="Measurement_Median"]<-"Broad_Median"
    names(compdf)[names(compdf)=="Measurement_MAD"]<-"Broad_MAD"
    #if outlier remove is true, perform outlier_removal function
    calc_spread(mydf,compdf)
}