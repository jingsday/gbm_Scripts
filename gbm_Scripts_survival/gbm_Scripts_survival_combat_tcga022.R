library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)

library(sva)
wkdir<- "/home/jing/Phd_project/project_GBM/gbm_DATA/"


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


comb_clin <- rbind(clin_325, clin_693)
comb_clin$Gender <- trimws(comb_clin$Gender)

# Ensure gene names match — intersect common genes
common_genes <- intersect(colnames(cgga_325RNA), colnames(cgga_693RNA),colnames(tcga_rna))#23271 shared genes
# Transpose so that genes are rows and samples are columns
expr_693 <- as.matrix(t(cgga_693RNA[, common_genes]))
expr_325 <- as.matrix(t(cgga_325RNA[,common_genes]))

# Merge into a single matrix
combined_expr <- cbind(expr_693, expr_325)
batch <- c(rep(1, ncol(expr_693)), rep(2, ncol(expr_325)))
adjusted_counts <- ComBat_seq(counts = combined_expr, batch = batch, group = NULL)

###
cgga325_adjusted_counts <- t(adjusted_counts[,colnames(adjusted_counts)%in% colnames(expr_325)])
cgga693_adjusted_counts <- t(adjusted_counts[,colnames(adjusted_counts)%in% colnames(expr_693)])

cgga_325RNA_log2 <- apply(cgga325_adjusted_counts, c(1, 2), function(value) log2(value + 1))
cgga_693RNA_log2 <- apply(cgga693_adjusted_counts, c(1, 2), function(value) log2(value + 1))


#TCGA
tcga_clin_raw <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_clinical_patient.txt"), 
                            sep = '\t',skip = 4,row.names = 'PATIENT_ID')

#clin profiles
clin_sample <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_clinical_sample.txt"), sep = '\t',skip = 4)
table(clin_sample$SAMPLE_TYPE)
table(clin_sample$TUMOR_TYPE)
tcga_clin <-tcga_clin_raw[tcga_clin_raw$SUBTYPE!='GBM_IDHmut-non-codel',]


tcga_rna_raw <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_mrna_seq_v2_rsem.txt"),check.names = FALSE)
tcga_rna_raw[is.na(tcga_rna_raw)] <- 0
tcga_rna_raw <- tcga_rna_raw[tcga_rna_raw$Hugo_Symbol!='',]
tcga_rna_raw <- tcga_rna_raw[!duplicated(tcga_rna_raw$Hugo_Symbol),]
rownames(tcga_rna_raw) <- tcga_rna_raw$Hugo_Symbol

tcga_rna_raw <- as.data.frame(t(tcga_rna_raw[-1:-2]))

tcga_clin <- tcga_clin[rownames(tcga_clin) %in% str_sub(row.names(tcga_rna_raw),end=-4),]
tcga_rna_raw <- tcga_rna_raw[str_sub(rownames(tcga_rna_raw),end=-4) %in% row.names(tcga_clin),]

intersect(str_sub(rownames(tcga_rna_raw),end=-4),rownames(tcga_clin))




# Step 4: Merge datasets by rows
RNA_combined <- rbind(cgga_693RNA_log2, cgga_325RNA_log2)

# Step 5: Apply z-scale to the merged dataset
RNA_combined <- scale(RNA_combined, center = TRUE, scale = TRUE)

# Convert back to a data frame if needed
RNA_combined <- as.data.frame(RNA_combined)


surv_obj <- Surv(time = comb_clin$OS, 
                 event = comb_clin$Censor..alive.0..dead.1.=="1")

#surv_obj 
fit <- survfit(surv_obj ~ 1, data = comb_clin)
ggsurvplot(fit, data = comb_clin, xlab = "Days", ylab = "Overall survival", surv.median.line = 'hv',risk.table = TRUE)

summary(fit)$table["median"]
ggsurvplot(fit, data = comb_clin, xlab = "Days", ylab = "Overall survival", 
           surv.median.line = 'hv', risk.table = TRUE)
comb_clin$Gender <- trimws(comb_clin$Gender)

fit <- survfit(surv_obj ~ Gender , data = comb_clin)
ggsurvplot(fit, data = comb_clin, xlab = "Month", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)


library(progeny)
zscores = as.matrix(t(RNA_combined))
pathways <- progeny(zscores, scale=TRUE, organism="Human")  #, top = 100, perm = 1)
path_df = as.data.frame(pathways)
# fit multivariate model (COX proportional hazard) 
colnames(path_df)

comb_clin_filt <- comb_clin[comb_clin$OS > 0,]
RNA_combined_filt <- RNA_combined[comb_clin$OS > 0,]
path_filt <- path_df[!is.na(comb_clin$OS) & comb_clin$OS > 0, ]


# create a survival object consisting of times & censoring
surv_filt <- Surv(time = comb_clin$OS, 
                  event = comb_clin$Censor..alive.0..dead.1.=="1")
fit <- survfit(surv_filt ~ 1, data = comb_clin)
ggsurvplot(fit, data = comb_clin, xlab = "Days", ylab = "Overall survival")

which(is.na(surv_filt))

paste(which(is.na(surv_filt)),collapse = ' , ')


RNA_combined_filt[is.na(RNA_combined_filt)] <- 0

surv_filt <-surv_filt[-c(109 , 146 , 175 , 177 , 186 , 241 , 250 , 272),]
RNA_combined_filt <- RNA_combined_filt[-c(109 , 146 , 175 , 177 , 186 , 241 , 250 , 272),]
comb_clin_filt<-comb_clin_filt[-c(109 , 146 , 175 , 177 , 186 , 241 , 250 , 272),]
path_filt <-path_df[-c(109 , 146 , 175 , 177 , 186 , 241 , 250 , 272),]


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
cfs = coef(fit_glm,s=0.001523)
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


fit_glm <- glmnet(RNA_combined_filt, surv_filt, family="cox",nlambda = 150)# , alpha = 1, standardize = TRUE, maxit = 1000)
print(fit_glm)

cvfit <- cv.glmnet(data.matrix(RNA_combined_filt),surv_filt,family="cox",type.measure = 'C')
plot(cvfit)


cfs = coef(fit_glm,s= 0.002275)

meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
#write.csv(meaning_vals,file=paste0(outdir,'gbm_survival_coeffs_full.csv'))

length(meaning_vals)
coef_data <-data.frame(variable=meaning_coefs,coefficient=meaning_vals)


cgga_693RNA_subset <- cgga_693RNA_log2[row.names(cgga_693RNA_log2)%in% row.names(clin_693), ]
cgga_693RNA_subset<- cgga_693RNA_subset[, coef_data$variable]
coef_data_subset <- coef_data[,-1]
result2 <- cgga_693RNA_subset %*% coef_data_subset

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


cfs = coef(fit_glm,s=  0.002419)#184 306 58.22 0.002419

#By increasing gene number it doesn't help either
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
#write.csv(meaning_vals,file=paste0(outdir,'gbm_survival_final_coeffs_lincs.csv'))


length(meaning_vals)
coef_data <-data.frame(variable=meaning_coefs,coefficient=meaning_vals)


cgga_693RNA_subset <- cgga_693RNA_log2[row.names(cgga_693RNA_log2)%in% row.names(clin_693), ]
cgga_693RNA_subset<- cgga_693RNA_subset[, coef_data$variable]
coef_data_subset <- coef_data[,-1]
result2 <- cgga_693RNA_subset %*% coef_data_subset


cgga_325RNA_subset<- cgga_325RNA_log2[, coef_data$variable]
result3 <- cgga_325RNA_subset %*% coef_data_subset

result <- rbind(result2,result3)

#populate values back
comb_clin$DPD_prognosis <- result[str_sub(row.names(result)) %in% row.names(comb_clin), ]
hist(comb_clin$DPD_prognosis)
comb_clin$outcome <- ifelse(comb_clin$DPD_prognosis > 0, "Positive",
                            ifelse(comb_clin$DPD_prognosis == 0, "0", "Negative"))

surv_obj_filt <- Surv(time = comb_clin$OS, 
                      event = comb_clin$Censor..alive.0..dead.1.=="1")

fit <- survfit(surv_obj_filt ~ outcome , data = comb_clin)
ggsurvplot(fit, data = comb_clin, xlab = "Days", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)

#
outdir <- '/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_survival/'
surv_plot <- ggsurvplot(
  fit,
  data = comb_clin,
  xlab = "Days",
  ylab = "Overall survival",
  pval = TRUE,
  risk.table = TRUE
)

# Save the plot as a PNG
ggsave(
  filename = paste0(outdir,"survival_final_plot.png"),
  plot = surv_plot$plot,
  device = "png",
  width = 8,
  height = 4,
  dpi = 300
)

ggsave(
  filename = paste0(outdir,"risk_final_table.png"),
  plot = surv_plot$table,         # Use the risk table component
  device = "png",
  width = 8,
  height = 2,
  dpi = 300
)

#
#clin
tcga_clin_raw <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_clinical_patient.txt"), 
                            sep = '\t',skip = 4,row.names = 'PATIENT_ID')

#clin profiles
clin_sample <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_clinical_sample.txt"), sep = '\t',skip = 4)
table(clin_sample$SAMPLE_TYPE)
table(clin_sample$TUMOR_TYPE)
tcga_clin <-tcga_clin_raw[tcga_clin_raw$SUBTYPE!='GBM_IDHmut-non-codel',]


tcga_rna_raw <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_mrna_seq_v2_rsem.txt"),check.names = FALSE)
tcga_rna_raw[is.na(tcga_rna_raw)] <- 0
tcga_rna_raw <- tcga_rna_raw[tcga_rna_raw$Hugo_Symbol!='',]
tcga_rna_raw <- tcga_rna_raw[!duplicated(tcga_rna_raw$Hugo_Symbol),]
rownames(tcga_rna_raw) <- tcga_rna_raw$Hugo_Symbol

tcga_rna_raw <- as.data.frame(t(tcga_rna_raw[-1:-2]))

tcga_clin <- tcga_clin[rownames(tcga_clin) %in% str_sub(row.names(tcga_rna_raw),end=-4),]
tcga_rna_raw <- tcga_rna_raw[str_sub(rownames(tcga_rna_raw),end=-4) %in% row.names(tcga_clin),]

intersect(str_sub(rownames(tcga_rna_raw),end=-4),rownames(tcga_clin))




#surv_obj 
surv_obj <- Surv(time = tcga_clin$OS_MONTHS, 
                 event = tcga_clin$OS_STATUS=="1:DECEASED")
fit <- survfit(surv_obj ~ 1, data = tcga_clin)
ggsurvplot(fit, data = tcga_clin, xlab = "Month", ylab = "Overall survival", surv.median.line = 'hv',risk.table = TRUE)


tcga_rna <- tcga_rna_raw
rownames(tcga_rna) <- rownames(tcga_rna)

tcga_clin <- tcga_clin[tcga_clin$SUBTYPE == 'GBM_IDHwt',]
table(tcga_clin$SUBTYPE)
#retrieve RNAs of interest
tcga_rna <- tcga_rna[str_sub(row.names(tcga_rna), end = -4) %in% row.names(tcga_clin), ]

# Remove rows whose rownames end with "-2"
rownames(tcga_rna)[grepl("-02$", rownames(tcga_rna))]
#"TCGA-06-0125-02" "TCGA-06-0171-02" "TCGA-06-0190-02" 
#"TCGA-06-0210-02" "TCGA-06-0211-02" "TCGA-14-1034-02"

check <- c("TCGA-06-0125-02", "TCGA-06-0171-02", "TCGA-06-0190-02", "TCGA-06-0210-02" ,
           "TCGA-06-0211-02" ,"TCGA-14-1034-02")
tcga_clin[rownames(tcga_clin) %in% str_sub(check, end = -4),] #all exist in clin
#
check_01 <- str_replace(check, "-02$", "-01")
rows_to_remove <- c("TCGA-06-0125-02", "TCGA-06-0190-02", 
                    "TCGA-06-0210-02", "TCGA-14-1034-02", 
                    "TCGA-06-0211-02")

tcga_rna <- tcga_rna[!rownames(tcga_rna) %in% rows_to_remove, ]
max(tcga_rna)
min(tcga_rna)

tcga_rna_log2 <- apply(tcga_rna, c(1, 2), function(value) log2(value + 1))
min(tcga_rna_log2) #0
max(tcga_rna_log2) #19.87406

#

tcga_common_genes <- intersect(colnames(tcga_rna_log2),coef_data$variable)
tcga_rna_log2_subset <- tcga_rna_log2[, tcga_common_genes]

intersect(rownames(coef_data),tcga_common_genes)

result1 <- tcga_rna_log2_subset %*% as.numeric(coef_data[tcga_common_genes,]$coefficient)

rownames(result1) <- str_sub(row.names(result1), end = -4)

tcga_clin$DPD_prognosis <- result1[str_sub(row.names(result1)) %in% row.names(tcga_clin), ]
hist(tcga_clin$DPD_prognosis)

#assign to positive,negative (I only have positive and decided to select a positive threshold)
tcga_clin$outcome <- ifelse(tcga_clin$DPD_prognosis > 0, "Positive",
                            ifelse(tcga_clin$DPD_prognosis == 0, "0", "Negative"))


rownames(clin_sample) <- clin_sample$SAMPLE_ID
clin_sample[rownames(clin_sample) %in% rownames(tcga_clin)]


surv_obj <- Surv(time = tcga_clin$OS_MONTHS, 
                 
                 event = tcga_clin$OS_STATUS=="1:DECEASED")
fit <- survfit(surv_obj ~ outcome , data = tcga_clin)

ggsurvplot(fit, data = tcga_clin, xlab = "Days", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)


fit.coxph <- coxph(surv_obj ~ SEX + AGE +outcome +RADIATION_THERAPY, data = tcga_clin)

ggforest(fit.coxph, data = tcga_clin)           
