#!/usr/bin/env Rscript
#Aggregate data output from Cell Painting.
#This file requires one input file input file to aggregate: this is the output file from CellProfiler.

#GOAL: TAKE MEAN AND STD OF EACH FEATURE AND CREATE CORRELATION PLOTS AND NORMAL CURVE
args=commandArgs(trailingOnly=TRUE)
library(data.table)
library(reshape)
library(tools)
library(ggplot2)
library(ggpmisc)
library(zplyr)
library(dplyr)

datafile=fread(args[1])
comparisondf=fread(args[2])
filegroup=args[3]

calc_mean<-function(input_file){
    #take mean and std of each column
    col_start<-grep("Cells",colnames(input_file))[[1]]
    mean_ofdata=colMeans(input_file[,col_start:ncol(input_file)])

    #reshape the files
    mean_ofdata<-melt(mean_ofdata)
    names(mean_ofdata)[names(mean_ofdata)=="value"]<-"Mean"
    mean_ofdata<-data.frame(Measurement=rownames(mean_ofdata),Mean=mean_ofdata$Mean)
    return(mean_ofdata)
}

calc_sd<-function(input_file){
    #take mean and std of each column
    col_start<-grep("Cells",colnames(input_file))[[1]]
    sd_ofdata=apply(input_file[,col_start:ncol(input_file)],2,sd)

    #reshape the files
    sd_ofdata<-melt(sd_ofdata)
    names(sd_ofdata)[names(sd_ofdata)=="value"]<-"Standard_Deviation"
    sd_ofdata<-data.frame(Measurement=rownames(sd_ofdata),Standard_Deviation=sd_ofdata$Standard_Deviation)
    return(sd_ofdata)
}


mymean=calc_mean(datafile)
my_out_mean<-boxplot(mymean$Mean, plot=FALSE)$out
my_NOout_mean<-mymean[-which(mymean$Mean %in% my_out_mean),]
mysd=calc_sd(datafile)
my_out_sd<-boxplot(mysd$Standard_Deviation, plot=FALSE)$out
my_NOout_sd<-mysd[-which(mysd$Standard_Deviation %in% my_out_sd),]
names(my_NOout_mean)[names(my_NOout_mean)=="Mean"]<-"UTSW_Mean"
names(my_NOout_sd)[names(my_NOout_sd)=="Standard_Deviation"]<-"UTSW_Standard_Deviation"
my_mean_sd<-merge(my_NOout_mean,my_NOout_sd,by="Measurement")

compmean=calc_mean(comparisondf)
#comp_out_mean<-boxplot(compmean$Mean, plot=FALSE)$out
#comp_NOout_mean<-compmean[-which(compmean$Mean %in% comp_out_mean),]
compsd=calc_sd(comparisondf)
#comp_out_sd<-boxplot(compsd$Standard_Deviation, plot=FALSE)$out
#comp_NOout_sd<-compsd[-which(compsd$Standard_Deviation %in% comp_out_sd),]
names(compmean)[names(compmean)=="Mean"]<-"Broad_Mean"
names(compsd)[names(compsd)=="Standard_Deviation"]<-"Broad_Standard_Deviation"
comp_mean_sd<-merge(compmean,compsd,by="Measurement")

merged<-merge(my_mean_sd,comp_mean_sd,by="Measurement")
merged<-na.omit(merged)

correlate<-function(input_file){
    mean_lm<-lm(Broad_Mean~UTSW_Mean,input_file)
    mean_lmcoef<-coef(mean_lm)
    mean_lmcoef<-data.frame(Mean_Slope=round(mean_lmcoef[[2]],3))
    mean_lmcoef$Mean_Slope<-paste0("Mean Slope=",mean_lmcoef$Mean_Slope)

    sd_lm<-lm(Broad_Standard_Deviation~UTSW_Standard_Deviation,input_file)
    sd_lmcoef<-coef(sd_lm)
    sd_lmcoef<-data.frame(SD_Slope=round(sd_lmcoef[[2]],3))
    sd_lmcoef$SD_Slope<-paste0("SD Slope=",sd_lmcoef$SD_Slope)

    lmcoef=cbind(UTSW_Mean=1,Broad_Mean=1,UTSW_Standard_Deviation=1,Broad_Standard_Deviation=1,mean_lmcoef,sd_lmcoef)
    return(lmcoef)
}

correl<-correlate(merged)

#plotting the correlation results
ggplot(merged, aes(x=UTSW_Mean,Broad_Mean)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x, colour="black") + stat_poly_eq() + labs(title=paste0("Correlation between Broad and UTSW Mean Feature Values\nUTSW Outliers Removed")) + geom_abs_text(data=correl,mapping = aes(label = Mean_Slope),color="black",size=3.8,xpos=0.155,ypos=0.89)
ggsave(paste0("MeanCorrelationAcross_",filegroup,"_Features.png"), type = "cairo")

ggplot(merged,aes(UTSW_Standard_Deviation,Broad_Standard_Deviation)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x,colour="black") + stat_poly_eq() + geom_abs_text(data=correl,mapping = aes(label = SD_Slope),color="black",size=3.8,xpos=0.143,ypos=0.89) + labs(title=paste0("Correlation of Standard Deviations between\nBroad and UTSW Feature Values"))
ggsave(paste0("SDCorrelationAcross_",filegroup,"_Features.png"), type = "cairo")

#plot bell curve of both data sets, too
origdf<-datafile[,grep("Cells",colnames(datafile))[[1]]:ncol(datafile)]
origdf<-cbind(Metadata_Plate=datafile$Metadata_Plate,origdf)
origdf=melt(origdf)
names(origdf)[names(origdf)=="variable"]<-"Measurement"
names(origdf)[names(origdf)=="value"]<-"UTSW_Value"

compdf<-comparisondf[,grep("Cells",colnames(comparisondf))[[1]]:ncol(comparisondf)]
compdf<-cbind(Metadata_Plate=comparisondf$Metadata_Plate,compdf)
compdf<-melt(compdf)
names(compdf)[names(compdf)=="variable"]<-"Measurement"
names(compdf)[names(compdf)=="value"]<-"Broad_Value"

merg_dfs<-inner_join(x=origdf,y=compdf,by="Measurement")

ggplot(merg_dfs) + stat_function(aes(x=Broad_Value,color="Broad"),fun=dnorm,args=list(mean=mean(merg_dfs$Broad_Value,na.rm=TRUE),sd=sd(merg_dfs$Broad_Value,na.rm=TRUE))) + stat_function(aes(x=UTSW_Value,color="UTSW"),fun=dnorm,args=list(mean=mean(merg_dfs$UTSW_Value,na.rm=TRUE),sd=sd(merg_dfs$UTSW_Value,na.rm=TRUE))) + ylab("") + xlab("") + theme(legend.position = 'bottom',legend.text = element_text(color='black',face='bold'),legend.title = element_text(color='black',face='bold')) + labs(color='Group',y='') + scale_color_manual(values=c('blue','black')) + labs(title=paste0("Distribution of All ",filegroup," Feature Values"))
ggsave(paste0("BellCurveAcross_",filegroup,"_Features.png"), type = "cairo")



producers <- data.frame(   
    surname =  c("Spielberg","Scorsese","Hitchcock","Tarantino","Polanski"),    
    nationality = c("US","US","UK","US","Poland"),    
    stringsAsFactors=FALSE)

# Create destination dataframe
movies <- data.frame(    
    surname = c("Spielberg",
		"Scorsese",
                "Hitchcock",
              	"Hitchcock",
                "Spielberg",
                "Tarantino",
                "Polanski",
                "Fake",
                "Fake"),    
    title = c("Super 8",
    		"Taxi Driver",
    		"Psycho",
    		"North by Northwest",
    		"Catch Me If You Can",
    		"Reservoir Dogs","Chinatown",
            "This is a Fake Movie",
            "This is a Second Fake Movie"),                
     		stringsAsFactors=FALSE)

# Merge two datasets
m1 <- merge(producers, movies, by.x = "surname")
m1