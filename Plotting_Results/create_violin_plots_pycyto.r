#!/usr/bin/env Rscript

library(data.table)
library(reshape)
library(stringr)
library(ggplot2)
library(dplyr)


args=commandArgs(trailingOnly=TRUE)

#read in the file and change blanks to DMSO
datafile=fread(args[1])
filetype=args[2]

datafile$Metadata_broad_sample[which(datafile$Metadata_broad_sample=="")]<-"DMSO"
plateid=unique(datafile$Metadata_Plate)

#create replicate information
broad_well<-paste0(datafile$Metadata_broad_sample,"_",datafile$Metadata_Well)
datafile<-cbind(broad_well,datafile)
datafile<-datafile %>% group_by(Metadata_broad_sample) %>% mutate(Replicate = cumsum(!duplicated(broad_well)))
datafile$Replicate[which(datafile$Metadata_broad_sample=="DMSO")]<-c(rep(1:4,length(datafile$Metadata_broad_sample[which(datafile$Metadata_broad_sample=="DMSO")])),1)
broad_replicate<-paste0(datafile$Metadata_broad_sample,"_",datafile$Replicate)
datafile<-cbind(broad_replicate=broad_replicate,datafile)
datafile = subset(datafile, select = -c(Replicate,broad_well) )

#melt the data frame with the replicate information
df_new<-datafile[,grep("Cells",colnames(datafile))[[1]]:ncol(datafile)]
df_new<-cbind(datafile$broad_replicate,datafile$Metadata_Well,df_new)
colnames(df_new)[1]<-"broad_replicate"
colnames(df_new)[2]<-"Metadata_Well"
df_new<-melt(df_new)
names(df_new)[names(df_new)=="variable"]<-"Measurement"
names(df_new)[names(df_new)=="value"]<-"Median"

#create a randomized dataframe
df_new_rand<-data.frame(Metadata_Well=df_new$Metadata_Well,broad_replicate=sample(df_new$broad_replicate),Measurement=df_new$Measurement,Median=df_new$Median)
df_new_rand[c('broad_sample','replicate')]<-str_split_fixed(df_new_rand$broad_replicate,'_',2)
df_new_rand = subset(df_new_rand, select = -c(broad_replicate) )

#separate broad_sample into two columns and remove that column
df_new[c('broad_sample','replicate')]<-str_split_fixed(df_new$broad_replicate,'_',2)
df_new = subset(df_new, select = -c(broad_replicate) )

#cast the data frames so that each replicate has its own median and find the correlations. This is done for the normal dataframe and the randomized dataframe.
df_new_cast<-cast(df_new,Measurement+broad_sample~replicate,value="Median",fun.aggregate=sum)
names(df_new_cast)[names(df_new_cast)==1]<-"A"
names(df_new_cast)[names(df_new_cast)==2]<-"B"
names(df_new_cast)[names(df_new_cast)==3]<-"C"
names(df_new_cast)[names(df_new_cast)==4]<-"D"

df_new_rand_cast<-cast(df_new_rand,Measurement+broad_sample~replicate,value="Median",fun.aggregate=sum)
names(df_new_rand_cast)[names(df_new_rand_cast)==1]<-"A"
names(df_new_rand_cast)[names(df_new_rand_cast)==2]<-"B"
names(df_new_rand_cast)[names(df_new_rand_cast)==3]<-"C"
names(df_new_rand_cast)[names(df_new_rand_cast)==4]<-"D"

#create a key based on broad sample and find correlations between replicates
df_new_cast<-data.table(df_new_cast)
setkey(df_new_cast,broad_sample)
df_cor<-df_new_cast[,list(repA.B=cor(A,B, use="complete.obs"),repA.C=cor(A,C, use="complete.obs"),repA.D=cor(A,D, use="complete.obs"),repB.C=cor(B,C, use="complete.obs"),repB.C=cor(B,D, use="complete.obs"),repB.D=cor(B,D, use="complete.obs")),by=key(df_new_cast)]
df_cor_melt<-melt(df_cor)
names(df_cor_melt)[names(df_cor_melt)=="variable"]<-"Measurement"
names(df_cor_melt)[names(df_cor_melt)=="value"]<-"Correlation"
df_cor_melt<-cbind(fake=rep("Non-Randomized",nrow(df_cor_melt)),df_cor_melt)


df_new_rand_cast<-data.table(df_new_rand_cast)
setkey(df_new_rand_cast,broad_sample)
df_cor_rand<-df_new_rand_cast[,list(repA.B=cor(A,B, use="complete.obs"),repA.C=cor(A,C, use="complete.obs"),repA.D=cor(A,D, use="complete.obs"),repB.C=cor(B,C, use="complete.obs"),repB.C=cor(B,D, use="complete.obs"),repB.D=cor(B,D, use="complete.obs")),by=key(df_new_rand_cast)]
df_cor_rand_melt<-melt(df_cor_rand)
names(df_cor_rand_melt)[names(df_cor_rand_melt)=="variable"]<-"Measurement"
names(df_cor_rand_melt)[names(df_cor_rand_melt)=="value"]<-"Correlation"
df_cor_rand_melt<-cbind(fake=rep("Randomized",nrow(df_cor_rand_melt)),df_cor_rand_melt)

#combine the two dataframes and create graph
df_cor_all<-rbind(df_cor_melt,df_cor_rand_melt)
quantile_randomized<-quantile(df_cor_all$Correlation[which(df_cor_all$fake=="Randomized")],0.95,na.rm=TRUE)
quantile_randomized=quantile_randomized[[1]]
ggplot(df_cor_all,aes(factor(fake),Correlation)) + geom_violin(scale='width') + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2) + theme(axis.text.x=element_text(angle=45,hjust=1)) + ggtitle(paste0(filetype,"\nReplicate Correlations")) +xlab("Comparison Group") + ylab("Correlation Value") + geom_hline(yintercept=quantile_randomized, color='red',linetype='dashed')
ggsave(paste0(plateid,"_",filetype,"_","violin_plot.pdf"))