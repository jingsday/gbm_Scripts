library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library(maftools)
library(ggsurvplot)

setwd('/Users/lidiayung/PhD_project/project_gbm/gbm_DATA/gbm_DATA_tcga_pan_can_atlas_2018')
gbm.maf <-read.delim("data_mutations.txt", sep = '\t')

head(gbm.maf$Tumor_Sample_Barcode)
head(gbm.maf$Matched_Norm_Sample_Barcode)


clin_raw<- read.delim("data_clinical_patient.txt", sep = '\t',skip = 4,row.names = 'PATIENT_ID')
table(clin_raw$SUBTYPE)

clin_sample <- read.delim("data_clinical_sample.txt", sep = '\t',skip = 4)
table(clin_sample$SAMPLE_TYPE)
recurrent <-clin_sample[clin_sample$SAMPLE_TYPE=="Recurrence",]
recurrent$PATIENT_ID

clin_raw$SAMPLE_TYPE ="Primary"
clin_raw[recurrent$PATIENT_ID,]$SAMPLE_TYPE='Recurrence'
table(clin_raw$SAMPLE_TYPE)

clin_raw <-clin_raw[clin_raw$SAMPLE_TYPE=='Primary',]
clin_raw <-clin_raw[clin_raw$SUBTYPE=='GBM_IDHwt',]
clin_raw <- clin_raw[clin_raw$SEX != "", ]

RNA_raw <- read.delim("data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",check.names = FALSE)
RNA_raw[is.na(RNA_raw)] <- 0
RNA_raw <- RNA_raw[RNA_raw$Hugo_Symbol!='',]
RNA_raw <- RNA_raw[!duplicated(RNA_raw$Hugo_Symbol),]
rownames(RNA_raw) <- RNA_raw$Hugo_Symbol
RNA <- as.data.frame(t(RNA_raw[-1:-2]))

# Align clinical data:
clin <- clin_raw[row.names(clin_raw)%in% str_sub(row.names(RNA), end = -4),]
RNA <- RNA[str_sub(row.names(RNA), end = -4) %in% row.names(clin_raw), ]

clin$Tumor_Sample_Barcode= RNA[str_sub(row.names(RNA), end = -4) %in% row.names(clin_raw), ]
clin$Tumor_Sample_Barcode =paste0(rownames(clin),'-01')

#maf_data_view 

gbm <-read.maf(maf=gbm.maf,clinicalData= clin)

# Filter MAF data to include only matching samples
matching_barcodes <- intersect(gbm.maf$Tumor_Sample_Barcode, clin$Tumor_Sample_Barcode)

gbm<- subsetMaf(maf = gbm, tsb = matching_barcodes)


maf_data_view <- gbm@data
maf_gene.summary <-gbm@gene.summary

getSampleSummary(gbm)

plotmafSummary(maf=gbm,rmOutlier = TRUE,addStat = 'median',dashboard = TRUE,
               titvRaw=FALSE)

mafbarplot(maf=gbm)

oncoplot(maf=gbm,top=30)

gbm.titv = titv(maf=gbm,plot=FALSE,useSyn=TRUE)
plotTiTv(res=gbm.titv)  

#individual genes TP53, ERCC1,FGFR3
lollipopPlot(maf=gbm,
             gene='FGFR3',
             AACol = 'HGVSp_Short',
             showMutationRate = TRUE)

lollipopPlot(maf = gbm,
             gene = 'FGFR3',
             AACol = 'HGVSp_Short',
             showMutationRate = TRUE,
             refSeqID = 'NM_022965') # Specify the transcript you want to use

plotProtein(gene="ERCC1")

rainfallPlot(maf = gbm, detectChangePoints = TRUE, pointSize = 0.4)

plotVaf(maf=gbm)

getClinicalData(gbm)

# prepare for survival analysis

gbm@clinical.data$deceased <- ifelse(gbm@clinical.data$OS_STATUS == "0:LIVING", 1, 0)

gbm@clinical.data[,OS_MONTHS :=  as.numeric(OS_MONTHS)]

par(mar = c(1, 1, 1, 1))
plot(1:30)
#Survival analysis
p <- mafSurvival(maf=gbm, genes='OLG',time ='OS_MONTHS',Status = 'deceased')
title(main = "Survival Analysis for TP53")


#
fit<- survfit(Surv(clin$OS_MONTHS, clin$OS_STATUS == '1:DECEASED') ~ SEX, data = clin)
ggsurvplot(fit, data = clin,
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE    ,pval = TRUE)

km_fit <-survfit(Surv(clin$OS_MONTHS, clin$OS_STATUS == '1:DECEASED') ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = clin)
ggsurvplot(km_fit,data=clin, risk.table = TRUE,pval = TRUE)
km_fit
fit <- survfit(surv_obj ~ 1, data = clin)
ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival")
fit

##Predict genesets with survivals
prog_geneset = survGroup(maf = gbm, top = 2, geneSetSize = 2, time = 'OS_MONTHS', Status = 'deceased', verbose = FALSE)

print(head(prog_geneset))


mafSurvGroup(maf = gbm, geneSet = c("", ""), time = "days_to_last_followup", Status = "Overall_Survival_Status")
