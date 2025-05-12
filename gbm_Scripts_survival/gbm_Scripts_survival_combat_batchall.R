library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)

library(sva)
wkdir<- "/home/jing/Phd_project/project_GBM/gbm_DATA/"

clin_raw <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_clinical_patient.txt"), 
             sep = '\t',skip = 4,row.names = 'PATIENT_ID')

#Filtered by stages of interest
table(clin_raw$SUBTYPE)

clin_raw <-clin_raw[clin_raw$SUBTYPE=='GBM_IDHwt',]


RNA_raw <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_mrna_seq_v2_rsem.txt"),check.names = FALSE)
RNA_raw[is.na(RNA_raw)] <- 0
RNA_raw <- RNA_raw[RNA_raw$Hugo_Symbol!='',]
RNA_raw <- RNA_raw[!duplicated(RNA_raw$Hugo_Symbol),]
rownames(RNA_raw) <- RNA_raw$Hugo_Symbol
RNA <- as.data.frame(t(RNA_raw[-1:-2]))

clin <- clin_raw[str_sub(row.names(RNA), end = -4),]
clin <- clin[!is.na(clin$OS_MONTHS) & !is.na(clin$OS_STATUS), ]

rows_to_remove <- c("TCGA-06-0125-02", "TCGA-06-0190-02", 
                    "TCGA-06-0210-02", "TCGA-14-1034-02", 
                    "TCGA-06-0211-02", 'TCGA-06-0125', 'TCGA-06-0171', 'TCGA-06-0190', 
                    'TCGA-06-0210', 'TCGA-06-0211', 'TCGA-14-1034' )

clin <- clin[!rownames(clin) %in% rows_to_remove,]
clin <- clin[!grepl("\\.1$", rownames(clin)), ]

RNA <- RNA[str_sub(row.names(RNA), end = -4) %in% row.names(clin), ]
table(clin$SUBTYPE)

### 


clin_693_raw <- read.delim(paste0(wkdir,"gbm_DATA_CGGA/CGGA.mRNAseq_693_clinical.20200506.txt"), 
                           sep = '\t')
table(clin_693_raw$IDH_mutation_status)
clin_693 <- clin_693_raw[clin_693_raw$IDH_mutation_status %in% c('Wildtype'),]
clin_693 <- clin_693[clin_693$Grade %in% c('WHO IV'),]

cgga_693RNA_raw <- read.delim(paste0(wkdir,"gbm_DATA_CGGA/CGGA.mRNAseq_693.RSEM-genes.20200506.txt"),check.names = FALSE)
cgga_693RNA_raw[is.na(cgga_693RNA_raw)] <- 0

cgga_693RNA_raw <- cgga_693RNA_raw[cgga_693RNA_raw$Gene_Name!='',]
cgga_693RNA_raw <- cgga_693RNA_raw[!duplicated(cgga_693RNA_raw$Gene_Name),]

rownames(cgga_693RNA_raw) <- cgga_693RNA_raw$Gene_Name

cgga_693RNA <- as.data.frame(t(cgga_693RNA_raw[-1]))
row.names(clin_693) <- clin_693$CGGA_ID
# Align clinical data:
cgga_693RNA <- cgga_693RNA[str_sub(row.names(cgga_693RNA)) %in% row.names(clin_693), ]


#clin
clin_325_raw <- read.delim(paste0(wkdir,"gbm_DATA_CGGA/CGGA.mRNAseq_325_clinical.20200506.txt"), sep = '\t')
clin_325 <- clin_325_raw[clin_325_raw$PRS_type !="Secondary",]
clin_325 <- clin_325[clin_325$IDH_mutation_status %in% c('Wildtype'),]
clin_325 <- clin_325[clin_325$Grade %in% c('WHO IV'),]
clin_325 <- clin_325[clin_325$Age != 11,] #Remove pediatric patients

table(clin_325$IDH_mutation_status)
#RNA
cgga_325RNA_raw <- read.delim(paste0(wkdir,"gbm_DATA_CGGA/CGGA.mRNAseq_325.RSEM-genes.20200506.txt"),check.names = FALSE)
cgga_325RNA_raw[is.na(cgga_325RNA_raw)] <- 0

cgga_325RNA_raw <- cgga_325RNA_raw[cgga_325RNA_raw$Gene_Name!='',]
cgga_325RNA_raw <- cgga_325RNA_raw[!duplicated(cgga_325RNA_raw$Gene_Name),]

rownames(cgga_325RNA_raw) <- cgga_325RNA_raw$Gene_Name

cgga_325RNA <- as.data.frame(t(cgga_325RNA_raw[-1]))
row.names(clin_325) <- clin_325$CGGA_ID

# Align clinical data:
cgga_325RNA <- cgga_325RNA[row.names(cgga_325RNA) %in% row.names(clin_325), ]
intersect(rownames(cgga_325RNA), rownames(cgga_693RNA)) # no overlapping row names


comb_clin_cgga <- rbind(clin_325, clin_693)
table(comb_clin_cgga$IDH_mutation_status) # wildtype
comb_clin_cgga$OS_MONTHS <- comb_clin_cgga$OS/30 # Days /30 
### 


####

# Ensure gene names match — intersect common genes
common_genes <- intersect(colnames(cgga_325RNA), colnames(cgga_693RNA))#23271 shared genes
common_genes <- intersect(common_genes,colnames(RNA))#15892 shared genes

# Transpose so that genes are rows and samples are columns
expr_693 <- as.matrix(t(cgga_693RNA[, common_genes]))
expr_325 <- as.matrix(t(cgga_325RNA[,common_genes]))
expr_tcga <- as.matrix(t(RNA[,common_genes]))

# Merge into a single matrix
combined_expr <- cbind(expr_693, expr_325, expr_tcga)
batch <- c(rep(1, ncol(expr_693)), rep(2, ncol(expr_325)), rep(3, ncol(expr_tcga)))
adjusted_counts <- ComBat_seq(counts = combined_expr, batch = batch, group = NULL)

###
cgga325_adjusted_counts <- t(adjusted_counts[,colnames(adjusted_counts)%in% colnames(expr_325)])
cgga693_adjusted_counts <- t(adjusted_counts[,colnames(adjusted_counts)%in% colnames(expr_693)])
tcga_adjusted_counts <- t(adjusted_counts[,colnames(adjusted_counts)%in% colnames(expr_tcga)])


cgga_325RNA_log2 <- apply(cgga325_adjusted_counts, c(1, 2), function(value) log2(value + 1))
cgga_693RNA_log2 <- apply(cgga693_adjusted_counts, c(1, 2), function(value) log2(value + 1))
tcga_RNA_log2 <- apply(tcga_adjusted_counts, c(1, 2), function(value) log2(value + 1))

RNA_combined_log2 <- rbind(cgga_693RNA_log2, cgga_325RNA_log2, tcga_RNA_log2)

RNA_combined <- scale(RNA_combined_log2, center = TRUE, scale = TRUE)
RNA_combined <- as.data.frame(RNA_combined)

commmon_cols <- c('PRS_type','Gender','Age','OS','Censor..alive.0..dead.1.','OS_MONTHS',"Radio_status..treated.1.un.treated.0." ,
                  "Chemo_status..TMZ.treated.1.un.treated.0.")

# comb_clin_cgga <- rbind(clin_325, clin_693)
# table(comb_clin_cgga$IDH_mutation_status) # wildtype
# comb_clin_cgga$OS_MONTHS <- comb_clin_cgga$OS/30 # Days /30 


commmon_cols <- c('PRS_type','Gender','Age','OS','Censor..alive.0..dead.1.','OS_MONTHS',"Radio_status..treated.1.un.treated.0.")
                  #"Chemo_status..TMZ.treated.1.un.treated.0.")
comb_clin <- comb_clin_cgga[,commmon_cols]

comb_clin$OS_STATUS <- ifelse(comb_clin$Censor..alive.0..dead.1. == 0, "0:LIVING",'1:DECEASED')
table(comb_clin$OS_STATUS)
comb_clin$RADIATION_THERAPY <- ifelse(comb_clin$Radio_status..treated.1.un.treated.0. == 1, "Yes",'NO')
table(comb_clin$Radio_status..treated.1.un.treated.0.)
###Change colnames 
paste0(colnames(comb_clin),collapse = ',')

comb_clin_fi <- comb_clin[,c('PRS_type','Gender','Age','OS_STATUS','OS_MONTHS','RADIATION_THERAPY')]
colnames(comb_clin_fi) <- c('PRS_type','SEX','AGE','OS_STATUS','OS_MONTHS','RADIATION_THERAPY')

#TCGA
sample <-  read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_clinical_sample.txt"), sep = '\t',skip = 4)
# table(sample[sample$PATIENT_ID %in% rownames(clin),'TUMOR_TYPE'])

treatment <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_timeline_treatment.txt"), sep = '\t')

clin_dup <- clin[,c('SEX','AGE','OS_STATUS', "OS_MONTHS" ,'RADIATION_THERAPY')]
clin_dup$PATIENT_ID <- rownames(clin_dup)
clin_dup <- merge(clin_dup,sample[,c('PATIENT_ID','TUMOR_TYPE','SAMPLE_TYPE')],by = 'PATIENT_ID',all.x=T)
rownames(clin_dup) <- clin_dup$PATIENT_ID

# clin_dup <- merge(clin_dup,treatment[,c('PATIENT_ID','TREATMENT_TYPE')],by = 'PATIENT_ID',all.x=T)
# rownames(clin_dup) <- clin_dup$PATIENT_ID

clin_dup <- clin_dup[,c('SAMPLE_TYPE','SEX','AGE','OS_STATUS', "OS_MONTHS")]

# clin_dup_fi <- clin_dup[,c('SAMPLE_TYPE','SEX','AGE','OS_STATUS','OS_MONTHS','RADIATION_THERAPY')]
# rownames(clin_dup) <- clin_dup$PATIENT_ID
clin_dup$RADIATION_THERAPY <-NA
colnames(clin_dup) <- c('PRS_type','SEX','AGE','OS_STATUS','OS_MONTHS','RADIATION_THERAPY')

colnames(comb_clin_fi)
colnames(clin_dup)
clin_comb_final <- rbind(comb_clin_fi,clin_dup)
colnames(clin_comb_final)


clin_comb_final$SEX <- trimws(clin_comb_final$SEX)

surv_obj <- Surv(time = clin_comb_final$OS_MONTHS, 
                 event = clin_comb_final$OS_STATUS=="1:DECEASED")

###
#surv_obj 
fit <- survfit(surv_obj ~ 1, data = clin_comb_final)
ggsurvplot(fit, data = comb_clin, xlab = "Days", ylab = "Overall survival", surv.median.line = 'hv',risk.table = TRUE)


summary(fit)$table["median"]
ggsurvplot(fit, data = clin_comb_final, xlab = "Months", ylab = "Overall survival", 
           surv.median.line = 'hv', risk.table = TRUE)


fit <- survfit(surv_obj ~ SEX , data = clin_comb_final)
ggsurvplot(fit, data = clin_comb_final, xlab = "Month", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)
table(clin_comb_final$SEX)


library(progeny)
zscores = as.matrix(t(RNA_combined))
pathways <- progeny(zscores, scale=TRUE, organism="Human")  #, top = 100, perm = 1)
path_df = as.data.frame(pathways)
# fit multivariate model (COX proportional hazard) 
colnames(path_df)


comb_clin_filt <- clin_comb_final[clin_comb_final$OS_MONTHS > 0,]
RNA_combined_filt <- RNA_combined[clin_comb_final$OS_MONTHS > 0,]
path_filt <- path_df[!is.na(clin_comb_final$OS_MONTHS) & clin_comb_final$OS_MONTHS > 0, ]


# create a survival object consisting of times & censoring
surv_filt <- Surv(time = clin_comb_final$OS_MONTHS, 
                  event = clin_comb_final$OS_STATUS=="1:DECEASED")
fit <- survfit(surv_filt ~ 1, data = clin_comb_final)
ggsurvplot(fit, data = clin_comb_final, xlab = "Days", ylab = "Overall survival")

fit <- survfit(surv_filt ~ 1, data = clin_comb_final)
ggsurvplot(fit, data = clin_comb_final, xlab = "Days", ylab = "Overall survival")

which(is.na(surv_filt))
paste(which(is.na(surv_filt)),collapse = ' , ')


surv_filt <-surv_filt[-c(108 , 145 , 174 , 176 , 185 , 240 , 249 , 271),]
RNA_combined_filt <- RNA_combined_filt[-c(108 , 145 , 174 , 176 , 185 , 240 , 249 , 271),]
comb_clin_filt<-comb_clin_filt[-c(108 , 145 , 174 , 176 , 185 , 240 , 249 , 271),]
path_filt <-path_df[-c(108 , 145 , 174 , 176 , 185 , 240 , 249 , 271),]


fit.coxph <- coxph(surv_filt ~  EGFR + Hypoxia   + TGFb  + `JAK-STAT` + TNFa, 
                   data = path_filt)
ggforest(fit.coxph, data = path_filt)


fit.coxph <- coxph(surv_filt ~  MAPK+NFkB+p53+PI3K+TGFb, 
                   data = path_filt)
ggforest(fit.coxph, data = path_filt)


library("glmpath")
library("glmnet")
library("penalized")


# glmnet over progeny
fit_glm <- glmnet(path_filt,surv_filt,family="cox") # , alpha = 1, standardize = TRUE, maxit = 1000
print(fit_glm)

# analysing results
cfs = coef(fit_glm,s=0.001175)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]


coef_data <-data.frame(variable=meaning_coefs,coefficient=meaning_vals)
# Create bar plot
ggplot(coef_data, aes(x = reorder(variable, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Pathways", y = "Coefficient") +
  theme_minimal()

clin_raw$
fit_glm <- glmnet(RNA_combined_filt, surv_filt, family="cox",nlambda = 150)# , alpha = 1, standardize = TRUE, maxit = 1000)
print(fit_glm)

cvfit <- cv.glmnet(data.matrix(RNA_combined_filt),surv_filt,family="cox",type.measure = 'C')
plot(cvfit)
table(clin_raw$ETHNICITY)

cfs = coef(fit_glm,s= 0.001647)

meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]

### Lincs gene
stv <- read_excel('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_LINCS/ALL_DATA_2020_Jing_gbm_del.xlsx',
                  sheet = 'STVs')
common_lincs <- intersect(colnames(RNA_combined_filt), stv$Gene)
RNA_combined_lincs <- RNA_combined_filt[,common_lincs]


fit_glm <- glmnet(RNA_combined_lincs, surv_filt, family="cox",nlambda = 200)# , alpha = 1, standardize = TRUE, maxit = 1000)
print(fit_glm)

cvfit <- cv.glmnet(data.matrix(RNA_combined_lincs),surv_filt,family="cox",type.measure = 'C',nlambda=200)
plot(cvfit)

cvfit$lambda.1se


cfs = coef(fit_glm,s=  0.001284)#184 306 58.22 0.002419 #200 407 42.16 0.001284
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]


length(meaning_vals)
coef_data <-data.frame(variable=meaning_coefs,coefficient=meaning_vals)

common_lincs_genes <- intersect(colnames(RNA_combined_log2),coef_data$variable)
RNA_combined_log2_lincs<- RNA_combined_log2[, common_lincs_genes]

coef_data_subset <- coef_data[rownames(coef_data) %in% common_lincs_genes,]

coef_data_subset <- coef_data[common_lincs_genes,-1]
result2 <- RNA_combined_log2_lincs %*% coef_data_subset

rownames(result2) <- str_remove(rownames(result2), "-\\d\\d$")

clin_comb_final$DPD_prognosis <- result2[rownames(result2) %in% rownames(clin_comb_final), ]

hist(clin_comb_final$DPD_prognosis)

clin_comb_final$outcome <- ifelse(clin_comb_final$DPD_prognosis > 0, "Positive",
                            ifelse(clin_comb_final$DPD_prognosis == 0, "0", "Negative"))

surv_obj_filt <- Surv(time = clin_comb_final$OS_MONTHS, 
                      event = clin_comb_final$OS_STATUS=="1:DECEASED")

fit <- survfit(surv_obj_filt ~ outcome , data = clin_comb_final)
ggsurvplot(fit, data = clin_comb_final, xlab = "Days", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)
