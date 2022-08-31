#load in data
library(data.table)
library(reshape)
library(ggplot2)
library(ggpmisc)

args=commandArgs(trailingOnly=TRUE)

#read in the data frames
mydf=fread(args[1])
comparisondf=fread(args[2])

#change data frames from short and wide to tall and skinny - My Data
mydf_new<-mydf[,grep("Cells",colnames(mydf))[[1]]:ncol(mydf)]
mydf_new<-cbind(mydf$Metadata_broad_sample,mydf_new)
colnames(mydf_new)[1]<-"Metadata_broad_sample"
mydf_new$Metadata_broad_sample[which(mydf_new$Metadata_broad_sample=="")]<-"DMSO"
mydf_new<-melt(mydf_new)
names(mydf_new)[names(mydf_new)=="variable"]<-"Measurement"
names(mydf_new)[names(mydf_new)=="value"]<-"Ashley_Median"

#change data frames from short and wide to tall and skinny - Broad Data
comparisondf_new<-comparisondf[,grep("Cells",colnames(comparisondf))[[1]]:ncol(comparisondf)]
comparisondf_new<-cbind(comparisondf$Metadata_broad_sample,comparisondf_new)
colnames(comparisondf_new)[1]<-"Metadata_broad_sample"
comparisondf_new$Metadata_broad_sample[which(comparisondf_new$Metadata_broad_sample=="")]<-"DMSO"
comparisondf_new<-melt(comparisondf_new)
names(comparisondf_new)[names(comparisondf_new)=="variable"]<-"Measurement"
names(comparisondf_new)[names(comparisondf_new)=="value"]<-"Broad_Median"

#merge data sheets by type
df_all<-merge(mydf_new,comparisondf_new,by=c("Metadata_broad_sample","Measurement"))

#creating correlation plot by Metadata_Compound
ggplot(df_all, aes(Ashley_Median,Broad_Median)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x, colour="black") + stat_poly_eq() + labs(x="Ashley's Measurements",y="Broad's Measurements",title="Correlation between Broad Data Set and Ours\nAfter PyCytominer")
___________________________________________________________________________________

#summarizing by compound

#will need to change this to a facet plot, one plot per compound
ggplot(cyt_all, aes(x=Well_Median,y=Median_Ashley,group=Metadata_Compound,color=Metadata_Compound)) + geom_point() + geom_smooth(method='lm',formula=y~x) + facet_wrap(~Metadata_Compound) + labs(title="Correlation of Ashley and Melissa\nRaw Median Cytoplasm Data by Control") + stat_poly_eq()

ggplot(cell_all, aes(x=Well_Median,y=Median_Ashley,group=Metadata_Cmpd,color=Metadata_Cmpd)) + geom_point() + geom_smooth(method='lm',formula=y~x) + facet_wrap(~Metadata_Cmpd) + labs(title="Correlation of Ashley and Melissa\nRaw Median Cellular Data by Control") + stat_poly_eq()

ggplot(nucl_all, aes(x=Well_Median,y=Median_Ashley,group=Metadata_Compound,color=Metadata_Compound)) + geom_point() + geom_smooth(method='lm',formula=y~x) + facet_wrap(~Metadata_Compound) + labs(title="Correlation of Ashley and Melissa\nRaw Median Nuclei Data by Control") + stat_poly_eq()
___________________________________________________________________________________
#BROAD DATA (ANALYZED BY OUR CP PIPE) AND MINE COMPARISON

cyt_broad=read.csv("cyto_broad_combined.csv")

cyt_ashley=read.csv("../../Mellissa_JUMP_cmpd_test/analysis_results/Plate1_IGVF_Test_10/cyto_rawdata_summarized.csv")
cyt_ashley$broad_sample[which(cyt_ashley$broad_sample=="")]<-"DMSO"

#combining broad_sample and replicate column, randomizing by broad-replicate assignments, and then separating that column
cyt_ashley_rand<-data.frame(Well=cyt_ashley$Well,broad_replicate=paste0(cyt_ashley$broad_sample,'_',cyt_ashley$replicate),variable=cyt_ashley$variable,Median_Ashley=cyt_ashley$Median)
cyt_ashley_rand<-data.frame(Well=cyt_ashley_rand$Well,broad_replicate=sample(cyt_ashley_rand$broad_replicate),variable=cyt_ashley$variable,Median_Ashley=cyt_ashley_rand$Median_Ashley)
library('stringr')
cyt_ashley_rand[c('broad_sample','replicate')]<-str_split_fixed(cyt_ashley_rand$broad_replicate,'_',2)

#read in metadata for melissa experiment
library(data.table)
metadata<-fread('../../Mellissa_JUMP_cmpd_test/external_metadata/IGVF_test_metadata.txt')
platemap<-fread('../../Mellissa_JUMP_cmpd_test/external_metadata/IGVF_Test_platemap.txt')
ext_meta<-merge(metadata,platemap,by="broad_sample")
names(ext_meta)[names(ext_meta) == 'well_position'] <- 'Metadata_Well'
ext_meta=data.frame(broad_sample=ext_meta$broad_sample,Well=ext_meta$Metadata_Well)

cyt_ashley=merge(ext_meta,cyt_ashley,by="Well")

#aggregate ashley data by compound, as there is 4 data points per cmpd per measurement in ashley data and 1 data point per cmpd per measurement in broad data.
cyt_ashley=aggregate(cyt_ashley, by=list(cyt_ashley$broad_sample,cyt_ashley$variable), mean)

names(cyt_ashley)[names(cyt_ashley) == 'Group.1'] <- 'broad_sample'

names(cyt_ashley)[names(cyt_ashley) == 'Group.2'] <- 'variable'

cyt_broad=aggregate(cyt_broad, by=list(cyt_broad$broad_sample,cyt_broad$variable), mean)

names(cyt_broad)[names(cyt_broad) == 'Group.1'] <- 'broad_sample'

names(cyt_broad)[names(cyt_broad) == 'Group.2'] <- 'variable'

cyt_ashley=data.frame(broad_sample=cyt_ashley[,1],variable=cyt_ashley[,2],Median_Ashley=cyt_ashley$Median_Ashley,MAD_Ashley=cyt_ashley$MAD_Ashley)

cyt_broad=data.frame(broad_sample=cyt_broad[,1],variable=cyt_broad[,2],Median=cyt_broad$Median,MAD=cyt_broad$MAD)

#merge data sheets by compound and variable
cyt_all<-merge(cyt_ashley,cyt_broad,by=c("broad_sample","variable"))

#running correlation by Metadata_Compound
library(ggplot2)
library(ggpmisc)

cell_lm1<-lm(Median~Median_Ashley,data=cell_all)

cell_lm2<-lm(MAD~MAD_Ashley,data=cell_all)

cell_rsq<-data.frame(model=c("Median","MAD"),r2=c(summary(cell_lm1)$adj.r.squared,summary(cell_lm2)$adj.r.squared),slope=c(coef(cell_lm1)[2],coef(cell_lm2)[2]))

ggplot(cyt_all, aes(Median,Median_Ashley)) + geom_point(colour="black") + geom_smooth(method='lm',formula=y~x, colour="black") + geom_point(aes(MAD,MAD_Ashley),colour="red") + geom_smooth(aes(MAD, MAD_Ashley),method='lm', formula=y~x, colour='red') + labs(x="Broad Measurement",y="Ashley Measurement",title="Correlation of Cytoplasm Raw Data Points\nBroad Data Analyzed by HTSCore CP Pipeline") + annotate("text",x=1190,y=15000,label=paste("Median r2 = ", round(cyt_rsq$r2[1],2)),colour="black") + annotate("text",x=1100,y=14000,label=paste("MAD r2 = ", round(cyt_rsq$r2[2],2)),colour="red") + annotate("text",x=7900,y=15200,label=paste("slope = ", round(cyt_rsq$slope[1],2)),colour="black") + annotate("text",x=3900,y=13000,label=paste("slope = ", round(cyt_rsq$slope[2],2)),colour="red")

___________________________________________________________________________________
#summarizing by compound

cell_lm<-dlply(cell_all,"broad_sample",function(df) lm(Median~Median_Ashley,data=df))
cell_lm_coef<-ldply(cell_lm,coef)

cell_dat<-data.frame(label=c())

#will need to change this to a facet plot, one plot per compound
ggplot(cyt_all, aes(x=Median,y=Median_Ashley,group=broad_sample,color=broad_sample)) + geom_point() + geom_smooth(method='lm',formula=y~x, colour="black") + facet_wrap(~broad_sample) + labs(title="Correlation of Ashley and Broad\nRaw Median Cytoplasm Data\nBroad Data Analyzed by HTS Core CP Pipe",x="Broad Median",y="HTS Core Median") + stat_poly_eq(size=3) + theme(strip.text=element_text(size=7))