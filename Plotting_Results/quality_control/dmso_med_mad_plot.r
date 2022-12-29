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

if(grepl("cp",filetype,fixed=TRUE)==TRUE){
    if(dmso==TRUE){
        filename="CellProfiler-DMSO-Only"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
    }else{
        filename="CellProfiler-Output"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
    }
    mydf=paste0(plate1,"/",plate1,".csv.gz")
    compdf=paste0(plate2,"/",basename(plate2),".csv.gz")
} else if(grepl("py",filetype,fixed=TRUE)==TRUE && grepl("feat",filetype,fixed=TRUE)==TRUE) {
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


#create an if outlier or if cp ==true, then join with metadata
if(grepl("cp",filetype,fixed=TRUE)==TRUE | outliers==TRUE){
    if(dmso==TRUE){
        filename="CellProfiler-DMSO-Only"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
    }else{
        filename="CellProfiler-Output"
        filename_sp=gsub("-"," ", filename,fixed=TRUE)
    }
    mydf=fread(paste0(plate1,"/",plate1,".csv.gz"))
    compdf=fread(paste0(plate2,"/",basename(plate2),".csv.gz"))
    metadata=fread("../../metadata/external_metadata/JUMP-Target-1_compound_metadata.tsv")
    broad_batchid=basename(dirname(plate2))
    broad_plateid=unique(compdf$Metadata_Plate)[1]
    broad_barcode=fread(paste0("../../metadata/platemaps/",broad_batchid,"/barcode_platemap.csv"))
    broad_barcode=broad_barcode %>% filter(Assay_Plate_Barcode==broad_plateid) %>% pull(var=Plate_Map_Name)
    broad_platemap=fread(paste0("../../metadata/platemaps/",broad_batchid,"/platemap/",broad_barcode,".txt"))
    names(broad_platemap)[names(broad_platemap)=="well_position"]<-"Metadata_Well"
    broad_platemap=broad_platemap[,1:2] #only get well_position and broad_sample information
    broad_metadata=inner_join(broad_platemap,metadata,by="broad_sample")
    broad_metadata=data.frame(Metadata_Well=broad_metadata$Metadata_Well,Metadata_pert_iname=broad_metadata$pert_iname,Metadata_broad_sample=broad_metadata$broad_sample)
    broad_metadata=subset(broad_metadata,select=c("Metadata_Well","Metadata_pert_iname"))
    compdf<-merge(broad_metadata,compdf,by="Metadata_Well")
}else{
    mydf=fread(mydf)
    compdf=fread(compdf)
}

aggregate_data<-function(file){
    if(dmso==TRUE){
        if(grepl("py",filetype,fixed=TRUE)==TRUE){
            file_meta<-data.frame(Metadata_Plate=file$Metadata_Plate,Metadata_Well=file$Metadata_Well,Metadata_pert_iname=file$Metadata_pert_iname)
            col_start<-grep("Cells",colnames(file))[[1]]
            file_data<-file[,col_start:ncol(file)]
            file<-melt(data.frame(c(file_meta,file_data)),id.vars=c("Metadata_Plate","Metadata_Well","Metadata_pert_iname"))
        }else{
            file<-melt(file,id.vars=c("Metadata_Plate","Metadata_Well","Metadata_pert_iname"))
        }
    }else{
        file<-melt(file,id.vars=c("Metadata_Plate","Metadata_Well"))
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
outlier_removal<-function(plate){
    file=paste0(plate,"/",plate,".csv.gz")
    file=fread(file)
    file_agg=aggregate_data(file)
    #need to remove outliers from CP file
    #then, if filetype = py, we match those up to the wells and measurements retained and replace outliers with NA
    #find outliers and replace with NA for each column

}

#if dmso is false and outlier is true, then perform outlier_removal function

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

if(dmso==TRUE){
    if(grepl("cp",filetype,fixed=TRUE)==TRUE){
        metadata=fread("../../metadata/external_metadata/JUMP-Target-1_compound_metadata.tsv")
        broad_batchid=basename(dirname(plate2))
        broad_plateid=unique(compdf$Metadata_Plate)[1]
        broad_barcode=fread(paste0("../../metadata/platemaps/",broad_batchid,"/barcode_platemap.csv"))
        broad_barcode=broad_barcode %>% filter(Assay_Plate_Barcode==broad_plateid) %>% pull(var=Plate_Map_Name)
        broad_platemap=fread(paste0("../../metadata/platemaps/",broad_batchid,"/platemap/",broad_barcode,".txt"))
        names(broad_platemap)[names(broad_platemap)=="well_position"]<-"Metadata_Well"
        broad_platemap=broad_platemap[,1:2] #only get well_position and broad_sample information
        broad_metadata=inner_join(broad_platemap,metadata,by="broad_sample")
        broad_metadata=data.frame(Metadata_Well=broad_metadata$Metadata_Well,Metadata_pert_iname=broad_metadata$pert_iname,Metadata_broad_sample=broad_metadata$broad_sample)
        broad_metadata=subset(broad_metadata,select=c("Metadata_Well","Metadata_pert_iname"))
        compdf<-merge(broad_metadata,compdf,by="Metadata_Well")
        mydf<-merge(broad_metadata,mydf,by="Metadata_Well")
    }
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