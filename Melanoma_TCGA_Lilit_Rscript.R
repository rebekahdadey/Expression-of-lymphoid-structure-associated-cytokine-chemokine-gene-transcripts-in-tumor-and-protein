#Help Lilit with melanoma TCGA samples to determine if gene expression correlates with survival
#Melanoma data download

#note, to run TCGA biolinks, must be on R4.1.0
#followed instructions on https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html
#https://www.biostars.org/p/153013/
setwd("/ix/jluke/Projects/becky/For/Lilit")

#issues with Biocmanager downloading TCGAbiolinks so had to install it from source.
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("dplyr")
library("tidyverse")
library("ensembldb")
library(tidyverse)
library(hrbrthemes)
library(viridis)

suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Hs.eg.db"))
BiocManager::install("ensembldb")
#first steps are to specify to R what data set you want to identify
query<- GDCquery(project="TCGA-SKCM",
                 data.category="Transcriptome Profiling",
                 data.type = "Gene Expression Quantification",
                 sample.type=c("Primary Tumor", "Metastatic"), #note, filtered out 1 solid tissue normal and 1 "additional metastatic". not sure what that means 
                 experimental.strategy="RNA-Seq")

#make results into table
lihc_res = getResults(query)

summary(factor(lihc_res$sample_type)) #look at numbers of samples per tissue type

#Additional Metastatic            Metastatic         Primary Tumor   Solid Tissue Normal 
###########1                   ######368                   103                     1 

#Download data for the query
GDCdownload(query = query)
#load the downloaded data into R
tcga_data = GDCprepare(query) 
dim(tcga_data) #60,660  #473



#take a look at column data, note, in this TCGA package use colData(), rowData(), assays()
colnames(colData(tcga_data))

#save tcga data as R obj
saveRDS(object = tcga_data,
        file = "tcga_data.RDS",
        compress = FALSE)

#look at what data is avilable from this data set
colnames(colData(tcga_data))
TCGAbiolinks:::getProjectSummary("TCGA-SKCM")


head(assay(tcga_data)[,1:10]) # expression of first 6 genes and first 10 samples
head(rowData(tcga_data))     # ensembl id and gene id of the first 6 genes.

#edgeR::DGEList converts the count matrix into an edgeR object.
dge= DGEList(counts=assay(tcga_data),
             samples=colData(tcga_data),
             genes=as.data.frame(rowData(tcga_data)))

#computes counts per million  (prior counts adds 3 so there is no 0)
#prior count=average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE.
dge.cpm = cpm(dge, log=F, prior.count=3)

write.csv(dge.cpm,
          file = gsub('.csv$',paste0(sample.ex.string,'.CPM.sm',sm.total,'.csv'), 
                      data.file))

dge.logcpm = cpm(dge, log=T, prior.count=3)

write.csv(dge.logcpm,
          file = gsub('.csv$',paste0(sample.ex.string,'.logCPM.sm',sm.total,'.csv'), 
                      data.file))


dge.cpm = cpm(dge, log=F, prior.count=3)


cpm.thres<-3

#create new column called total. Apply function to the rows of the dge cpm matrix
#function is to sum up the cpms per gene
dge.cpm.stats = data.frame(total = 
                             apply(dge.cpm, 1, 
                                   function(x) sum(x>cpm.thres, na.rm = T)))

#find genes which have more total CPM than the number of genes/2
keep = which(dge.cpm.stats$total >= as.integer(ncol(dge.cpm) / 2))
print(length(keep)) #12014


## recalculate lib size after filtering out low expr genes (Mathew Stephens)
dim(dge) ## 60660
dge = dge[keep,,keep.lib.sizes=FALSE]
dim(dge) ## 12014

## norm AFTER filtering genes and before voom
dge = calcNormFactors(dge, method="TMM") 

## produce new CPM and logCPM matrix after filtering
write.csv(cpm(dge, log=F, prior.count=0.5),
          file = "Filtered_TCGA_CPM.csv")
#log=T returns log2 values
write.csv(cpm(dge, log=T, prior.count=3),
          file = "Filtered_TCGA_logT.csv")

saveRDS(dge, file="Filtered_DGE_object.rds") 
Filtered_TCGA_logT<-read.csv("Filtered_TCGA_logT.csv")
#get Filtered logT TCGA into right format
colnames(Filtered_TCGA_logT) = gsub(".", "-", colnames(Filtered_TCGA_logT), fixed=T)
Filtered_TCGA_logT<- Filtered_TCGA_logT %>% column_to_rownames(var="X")
write.csv(Filtered_TCGA_logT,"Filtered_TCGA_logT.csv")
Filtered_TCGA_logT<-read.csv("Filtered_TCGA_logT.csv")
#access clnical data
clinical = tcga_data@colData

dim(clinical)

#this dataset does not contain days until last folllow up 
clin_df = clinical[
                   c("barcode",
                     "patient",
                     "vital_status",
                     "days_to_death",
                     "gender",
                     "days_to_last_follow_up",
                     "sample_type")]

# create a new boolean variable that has TRUE for dead patients
# and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# show first 10 samples
head(clin_df)

Surv(clin_df$overall_survival, clin_df$deceased) #adds censoring data

#try to find out if gender has any impact on survival (just an example)
Surv(clin_df$overall_survival, clin_df$deceased) ~ clin_df$gender

fit = survfit(Surv(overall_survival, deceased) ~ gender, data=clin_df)

print(fit)

a<-ggsurvplot(fit, data=clin_df, pval=T) #no statistically significant difference between genders

#export figure as high res tiff
tiff("Lilit_Melanoma_all_samples_gender_survival_plot.tiff", res=300, width=8, height=6, units='in')
a
dev.off()

write.csv(clin_df, "clin_df.csv")

#Gene expression, need to change ENSEMBL id to gene id for downstream purposes
#check df doesn't have any duplicates
any(duplicated(rownames(Filtered_TCGA_logT))) #false
any(duplicated(colnames(Filtered_TCGA_logT))) #false
#just to note that mapIDs doesn't like the the .version on the Filtered_TCGA ensembl Id's
rownames(Filtered_TCGA_logT)=gsub("\\..*","", rownames(Filtered_TCGA_logT))
write.csv(Filtered_TCGA_logT,"Filtered_TCGA_logT.csv")
Filtered_TCGA_logT<-read.csv("Filtered_TCGA_logT.csv")
Filtered_TCGA_logT$SYMBOL <- mapIds(org.Hs.eg.db, keys=rownames(Filtered_TCGA_logT),column="SYMBOL",keytype="ENSEMBL",multiVals="first")


#notice there are 2 genes that have two rows. delete these so we can put Gene name as row names


Filtered_TCGA_logT<-Filtered_TCGA_logT %>% distinct(SYMBOL, .keep_all = TRUE)
Filtered_TCGA_logT <- subset(Filtered_TCGA_logT, !is.na(Filtered_TCGA_logT$SYMBOL))
rownames(Filtered_TCGA_logT)=Filtered_TCGA_logT$SYMBOL
#can now remove the symbol column since its now in the rownames
Filtered_TCGA_logT = subset(Filtered_TCGA_logT, select = -c(SYMBOL) )
write.csv(Filtered_TCGA_logT,"Filtered_TCGA_logT.csv")

#make a smaller data frame with only data from APRIL. ***Lilit, here you would write in whatever gene you want in this space***
APRIL.expression<-Filtered_TCGA_logT[c("TNFSF13"),] 
str(APRIL.expression)
APRIL.expression<-as.data.frame(APRIL.expression)
APRIL.expression$TNFSF13=as.numeric(APRIL.expression$TNFSF13)
APRIL.expression<-rename(APRIL.expression, APRIL.expression=TNFSF13)
median_value = median(APRIL.expression$TNFSF13)
print(median_value) # 2.399596 #this is the medium expression so split high and low based off this

rownames(APRIL.expression) = gsub(".", "-", rownames(APRIL.expression), fixed=T)
Filtered_TCGA_logT.w.APRIL<-cbind(APRIL.expression, clin_df)
saveRDS(Filtered_TCGA_logT.w.APRIL, "Filtered_TCGA_logT.w.APRIL.rds")

Filtered_TCGA_logT.w.APRIL<-as.data.frame(Filtered_TCGA_logT.w.APRIL)

Filtered_TCGA_logT.w.APRIL$gene = NA
Filtered_TCGA_logT.w.APRIL$gene=ifelse(Filtered_TCGA_logT.w.APRIL$APRIL.expression >= median_value, "UP", "DOWN")


# we can fit a survival model, like we did in the previous section
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=Filtered_TCGA_logT.w.APRIL)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=Filtered_TCGA_logT.w.APRIL)$pval
print(pval)

b<-ggsurvplot(fit, data=Filtered_TCGA_logT.w.APRIL, pval=T, risk.table=T, title=paste("Overall Survival with TCGA APRIL expression"))
tiff("Lilit_Melanoma_all_samples_APRIL_Expression_OSS.tiff",res=300, width=8, height=6, units='in')
b
dev.off()

saveRDS(Filtered_TCGA_logT.w.APRIL, "Filtered_TCGA_logT.w.APRIL.rds")


##Primary vs metastastic split########----------------------------------
######################################################


clin_df_split = clinical[
  c("barcode",
    "patient",
    "vital_status",
    "days_to_death",
    "gender",
    "days_to_last_follow_up",
    "sample_type")]



# create a new boolean variable that has TRUE for dead patients
# and FALSE for live patients
clin_df_split$deceased = clin_df_split$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin_df_split$overall_survival = ifelse(clin_df_split$deceased,
                                  clin_df_split$days_to_death,
                                  clin_df_split$days_to_last_follow_up)

# show first 10 samples
head(clin_df_split)

Surv(clin_df_split$overall_survival, clin_df_split$deceased) #adds censoring data

write.csv(clin_df_split, "clin_df_split.csv")
clin_df_split<-read.csv("clin_df_split.csv")

Filtered_TCGA_logT.w.APRIL.met<-cbind(APRIL.expression, clin_df_split)
Filtered_TCGA_logT.w.APRIL.met<-Filtered_TCGA_logT.w.APRIL.met %>% dplyr::filter(sample_type=="Metastatic")

saveRDS(Filtered_TCGA_logT.w.APRIL.met, "Filtered_TCGA_logT.w.APRIL.met.rds")

Filtered_TCGA_logT.w.APRIL.met<-as.data.frame(Filtered_TCGA_logT.w.APRIL.met)

Filtered_TCGA_logT.w.APRIL.met$gene = NA
Filtered_TCGA_logT.w.APRIL.met$gene=ifelse(Filtered_TCGA_logT.w.APRIL.met$APRIL.expression >= median_value, "UP", "DOWN")


# we can fit a survival model, like we did in the previous section
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=Filtered_TCGA_logT.w.APRIL.met)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=Filtered_TCGA_logT.w.APRIL.met)$pval
print(pval)

c<-ggsurvplot(fit, data=Filtered_TCGA_logT.w.APRIL.met, pval=T, risk.table=T, title=paste("Overall Survival with Metastatic Melanoma TCGA APRIL expression"))
tiff("Lilit_Metastatic_Melanoma_APRIL_Expression_OSS.tiff", res=300, width=8, height=6, units='in')
c
dev.off()


saveRDS(Filtered_TCGA_logT.w.APRIL.met, "Filtered_TCGA_logT.w.APRIL.met.rds")


##Primary samples!!##########################
Filtered_TCGA_logT.w.APRIL.met<-cbind(APRIL.expression, clin_df_split)
Filtered_TCGA_logT.w.APRIL.primary<-Filtered_TCGA_logT.w.APRIL.met %>% dplyr::filter(sample_type=="Primary Tumor")

Filtered_TCGA_logT.w.APRIL.primary$gene = NA
Filtered_TCGA_logT.w.APRIL.primary$gene=ifelse(Filtered_TCGA_logT.w.APRIL.primary$APRIL.expression >= median_value, "UP", "DOWN")


# we can fit a survival model, like we did in the previous section
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=Filtered_TCGA_logT.w.APRIL.primary)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=Filtered_TCGA_logT.w.APRIL.primary)$pval
print(pval)

d<-ggsurvplot(fit, data=Filtered_TCGA_logT.w.APRIL.primary, pval=T, risk.table=T, title=paste("Overall Survival with Primary Melanoma TCGA APRIL expression"))
tiff("Lilit_Primary_Melanoma_APRIL_Expression_OSS.tiff", res=300, width=8, height=6, units='in')
d
dev.off()

saveRDS(Filtered_TCGA_logT.w.APRIL.primary, "Filtered_TCGA_logT.w.APRIL.primary.rds")

##lilit wants to analyze more cytokines--------------------------------------
#need to merge met vs primary 
clin_df<-read.csv("clin_df.csv")
Filtered_TCGA_logT.w.cytokines<-cbind(Filtered_TCGA_logT.w.APRIL, clin_df)
Filtered_TCGA_logT<-read.csv("Filtered_TCGA_logT.csv")
Filtered_TCGA_logT<- Filtered_TCGA_logT %>% column_to_rownames(var="X")

cytokines<-Filtered_TCGA_logT[c("CXCL13","CCL19","CCL21","CXCL10","TNFSF13"),] 
str(cytokines)
cytokines<-as.data.frame(cytokines)
cytokines<-t(cytokines)
for (i in colnames(cytokines)){
  print(median(cytokines[,i]))
}

#[1] 2.714015 #CXCL13
#[1] 2.386499 #CCl19
#[1] 2.366086 #CCL21
#[1] 3.992584 #CXCL10
#[1] 2.399596 #TNFSF13
#check my loop is good
str(cytokines)
cytokines<-as.data.frame(cytokines)
str(cytokines)
median(cytokines$CXCL13) #yeah the loop is good
median(cytokines$CCL19)
2.714015->CXCL13
2.386499->CCl19
2.366086->CCL21
3.992584->CXCL10
2.399596->TNFSF13


Filtered_TCGA_logT.w.APRIL<-cbind(cytokines, Filtered_TCGA_logT.w.APRIL)

Filtered_TCGA_logT.w.APRIL$CXCL13_level=ifelse(Filtered_TCGA_logT.w.APRIL$CXCL13 >= CXCL13, "UP", "DOWN")
saveRDS(Filtered_TCGA_logT.w.APRIL, "Filtered_TCGA_logT.w.APRIL.rds")

Filtered_TCGA_logT.w.APRIL$CCl19_level=ifelse(Filtered_TCGA_logT.w.APRIL$CCL19 >= CCl19, "UP", "DOWN")
Filtered_TCGA_logT.w.APRIL$CCL21_level=ifelse(Filtered_TCGA_logT.w.APRIL$CCL21 >= CCL21, "UP", "DOWN")
Filtered_TCGA_logT.w.APRIL$CXCL10_level=ifelse(Filtered_TCGA_logT.w.APRIL$CXCL10 >= CXCL10, "UP", "DOWN")
Filtered_TCGA_logT.w.APRIL$TNFSF13_level=ifelse(Filtered_TCGA_logT.w.APRIL$TNFSF13 >= TNFSF13, "UP", "DOWN")
saveRDS(Filtered_TCGA_logT.w.APRIL, "Filtered_TCGA_logT.w.APRIL.rds")


#loop for overall all patients--------------

plot_list=list()
for (i in colnames(Filtered_TCGA_logT.w.APRIL[,17:21])){
  fit = survfit(Surv(overall_survival, deceased) ~ Filtered_TCGA_logT.w.APRIL[,i], data=Filtered_TCGA_logT.w.APRIL)
  p=ggsurvplot(fit, data=Filtered_TCGA_logT.w.APRIL, pval=T, risk.table=T, title=paste("Overall Survival with TCGA", i, "expression"))
  plot_list[[i]]=p
}

#make tiff
for (i in colnames(Filtered_TCGA_logT.w.APRIL[,17:21])){
  file_name=paste("Overall Survival with TCGA", i, "expression")
  tiff(file_name, units="in", width=10, height=5, res=300)
  print(plot_list[[i]])
  dev.off()
}

#------------------------------------------------------------
#met from APRIL dataframe
cols.keep<-c("X","sample_type")
clin_df_split$X<- gsub("-", ".", clin_df_split$X)
clin_df_split<-clin_df_split[,cols.keep]
clin_df_split<-clin_df_split %>% column_to_rownames(var="X")
Filtered_TCGA_logT.w.APRIL<-cbind(clin_df_split, Filtered_TCGA_logT.w.APRIL)
saveRDS(Filtered_TCGA_logT.w.APRIL, "Filtered_TCGA_logT.w.APRIL.rds" )

met.1<-Filtered_TCGA_logT.w.APRIL %>% dplyr::filter(sample_type=="Metastatic")
primary.1<-Filtered_TCGA_logT.w.APRIL %>% dplyr::filter(sample_type=="Primary Tumor")
saveRDS(met.1,"met.1.rds")
saveRDS(primary.1, "primary.1.rds")
#metastatic levels of cytokines

plot_list=list()
for (i in colnames(met.1[,17:21])){
  fit = survfit(Surv(overall_survival, deceased) ~ met.1[,i], data=met.1)
  p=ggsurvplot(fit, data=met.1, pval=T, risk.table=T, title=paste("Overall Metastatic Survival with TCGA", i, "expression"))
  plot_list[[i]]=p
}

#make pdf
for (i in colnames(met.1[,17:21])){
  file_name=paste("Overall Metastatic Survival with TCGA", i, "expression")
  tiff(file_name, units="in", width=7, height=5, res=300)
  print(plot_list[[i]])
  dev.off()
}

#primary------------------------------------

plot_list=list()
for (i in colnames(primary.1[,17:21])){
  fit = survfit(Surv(overall_survival, deceased) ~primary.1[,i], data=primary.1)
  p=ggsurvplot(fit, data=primary.1, pval=T, risk.table=T, title=paste("Overall Primary Survival with TCGA", i, "expression"))
  plot_list[[i]]=p
}

#make pdf
for (i in colnames(primary.1[,17:21])){
  file_name=paste("Overall Primary Survival with TCGA", i, "expression")
  tiff(file_name, units="in", width=7, height=5, res=300)
  print(plot_list[[i]])
  dev.off()
}

#---bar plot of expression betwen primary and met sites for each cytokine------

Filtered_TCGA_logT.w.APRIL<-Filtered_TCGA_logT.w.APRIL %>% group_by(sample_type)
str(Filtered_TCGA_logT.w.APRIL)

plot_list=list()
for (i in colnames(Filtered_TCGA_logT.w.APRIL[,2:6])){
    p=ggplot(Filtered_TCGA_logT.w.APRIL, aes(x=sample_type, y=Filtered_TCGA_logT.w.APRIL[,i], fill=Filtered_TCGA_logT.w.APRIL[,i], col=sample_type)) +
          geom_boxplot() +
          geom_jitter(color="black", size=0.4, alpha=0.9) +
          theme(legend.position="none", plot.title = element_text(size=11)) +
          ggtitle(paste(i, "expression in Metastatic vs Primary Melanoma TCGA"))+
          stat_compare_means(test='wilcox.test', comparisons=list(c("Metastatic", "Primary Tumor")))
    plot_list[[i]]=p}

str(Filtered_TCGA_logT.w.APRIL)
Filtered_TCGA_logT.w.APRIL<-ungroup(Filtered_TCGA_logT.w.APRIL)
Filtered_TCGA_logT.w.APRIL<-as.data.frame(Filtered_TCGA_logT.w.APRIL)

for (i in colnames(Filtered_TCGA_logT.w.APRIL[,2:6])){
  file_name=paste("Overall Expression in TCGA of", i, "in Primary vs Metastatic Melanoma")
  tiff(file_name,units="in", width=7, height=5, res=300)
  print(plot_list[[i]])
  dev.off()
}

