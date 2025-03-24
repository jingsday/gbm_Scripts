library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)


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
cgga_693RNA_log2 <- apply(cgga_693RNA, c(1, 2), function(value) log2(value + 1))
#cgga_693RNA <- as.data.frame(scale(cgga_693RNA_log2, center = TRUE, scale = TRUE))


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

cgga_325RNA_log2 <- apply(cgga_325RNA, c(1, 2), function(value) log2(value + 1))
#cgga_325RNA <- as.data.frame(scale(cgga_325RNA_log2, center = TRUE, scale = TRUE))


common_genes <- intersect(colnames(cgga_325RNA), colnames(cgga_693RNA))
cgga_693RNA[,common_genes]

comb_clin <- rbind(clin_325, clin_693)
comb_clin$Gender <- trimws(comb_clin$Gender)

#RNA
common_genes <- intersect(colnames(cgga_693RNA_log2), colnames(cgga_325RNA_log2))

# Step 3: Subset to common genes
cgga_693RNA_subset <- cgga_693RNA_log2[, common_genes, drop = FALSE]
cgga_325RNA_subset <- cgga_325RNA_log2[, common_genes, drop = FALSE]

# Step 4: Merge datasets by rows
RNA_combined <- rbind(cgga_693RNA_subset, cgga_325RNA_subset)

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
cfs = coef(fit_glm,s=0.000925)
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

#Full gene set

fit_glm <- glmnet(RNA_combined_filt, surv_filt, family="cox")# , alpha = 1, standardize = TRUE, maxit = 1000)
print(fit_glm)

cvfit <- cv.glmnet(data.matrix(RNA_combined_filt),surv_filt,family="cox",type.measure = 'C')
plot(cvfit)

cfs = coef(fit_glm,s= 0.002168)

meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
#write.csv(meaning_vals,file=paste0(outdir,'gbm_survival_coeffs_full.csv'))

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
#assign to positive,negative (I only have positive and decided to select a positive threshold)
comb_clin$outcome <- ifelse(comb_clin$DPD_prognosis > 0, "Positive",
                            ifelse(comb_clin$DPD_prognosis == 0, "0", "Negative"))


### KP and hazard ratio
surv_obj_filt <- Surv(time = comb_clin$OS, 
                      event = comb_clin$Censor..alive.0..dead.1.=="1")

fit <- survfit(surv_obj_filt ~ outcome , data = comb_clin)
ggsurvplot(fit, data = comb_clin, xlab = "Days", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)



fit.coxph <- coxph(surv_obj_filt ~ Gender + Age +Histology+outcome , data = comb_clin)
ggforest(fit.coxph, data = comb_clin)


forest_plot <- ggforest(fit.coxph, data = comb_clin)

#Conclusions: full gene set not necessarily better


#Using previous STV from TCGA, to check if it applies to CGGA
common_genes <- intersect(colnames(RNA_combined_filt), stv$Gene)

cgga_693RNA_subset <- cgga_693RNA_log2[row.names(cgga_693RNA_log2)%in% row.names(clin_693), ]
cgga_693RNA_subset<- cgga_693RNA_subset[, common_genes]

result2 <- cgga_693RNA_subset %*% vector


cgga_325RNA_subset<- cgga_325RNA_log2[, common_genes]
result3 <- cgga_325RNA_subset %*% vector

result <- rbind(result2,result3)
#populate values back

comb_clin$DPD_prognosis <- result[str_sub(row.names(result)) %in% row.names(comb_clin), ]


hist(comb_clin$DPD_prognosis)
#assign to positive,negative (I only have positive and decided to select a positive threshold)
comb_clin$outcome <- ifelse(comb_clin$DPD_prognosis > 0, "Positive",
                            ifelse(comb_clin$DPD_prognosis == 0, "0", "Negative"))


### KP and hazard ratio
surv_obj_filt <- Surv(time = comb_clin$OS, 
                      event = comb_clin$Censor..alive.0..dead.1.=="1")

fit <- survfit(surv_obj_filt ~ outcome , data = comb_clin)
ggsurvplot(fit, data = comb_clin, xlab = "Days", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)



fit.coxph <- coxph(surv_obj_filt ~ Gender + Age +Histology+outcome , data = comb_clin)
ggforest(fit.coxph, data = comb_clin)


forest_plot <- ggforest(fit.coxph, data = comb_clin)

# Save the forest plot as an image
ggsave(paste0(outdir,"forest_plot.png"), plot = forest_plot, width = 10, height = 6, dpi = 300)

















#Lincs genes 

stv <- read_excel('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_LINCS/ALL_DATA_2020_Jing_gbm_del.xlsx',
                  sheet = 'STVs')

common_genes <- intersect(colnames(RNA_combined_filt), stv$Gene)

RNA_combined_lincs <- RNA_combined_filt[,common_genes]


fit_glm <- glmnet(RNA_combined_lincs, surv_filt, family="cox",nlambda = 200)# , alpha = 1, standardize = TRUE, maxit = 1000)
print(fit_glm)

cvfit <- cv.glmnet(data.matrix(RNA_combined_lincs),surv_filt,family="cox",type.measure = 'C',nlambda=200)
plot(cvfit)

cvfit$

cvfit$lambda.1se
print(lambda_1se)

cfs = coef(fit_glm,s=  0.001654)
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
result2 <- cgga_693RNA_subset %*% as.numeric(stv$GBM_survival)


cgga_325RNA_subset<- cgga_325RNA_log2[, coef_data$variable]
result3 <- cgga_325RNA_subset %*% coef_data_subset

result <- rbind(result2,result3)
#populate values back

comb_clin$DPD_prognosis <- result[str_sub(row.names(result)) %in% row.names(comb_clin), ]


hist(comb_clin$DPD_prognosis)
#assign to positive,negative (I only have positive and decided to select a positive threshold)
comb_clin$outcome <- ifelse(comb_clin$DPD_prognosis > 0, "Positive",
                            ifelse(comb_clin$DPD_prognosis == 0, "0", "Negative"))


### KP and hazard ratio
surv_obj_filt <- Surv(time = comb_clin$OS, 
                      event = comb_clin$Censor..alive.0..dead.1.=="1")

fit <- survfit(surv_obj_filt ~ outcome , data = comb_clin)
ggsurvplot(fit, data = comb_clin, xlab = "Days", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)



fit.coxph <- coxph(surv_obj_filt ~ Gender + Age +Histology+outcome , data = comb_clin)
ggforest(fit.coxph, data = comb_clin)


forest_plot <- ggforest(fit.coxph, data = comb_clin)

# Save the forest plot as an image
ggsave(paste0(outdir,"forest_plot.png"), plot = forest_plot, width = 10, height = 6, dpi = 300)


library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)

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


#testing on tcga 
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

surv_obj <- Surv(time = tcga_clin$OS_MONTHS, 
                 event = tcga_clin$OS_STATUS=="1:DECEASED")

#surv_obj 
fit <- survfit(surv_obj ~ 1, data = tcga_clin)
ggsurvplot(fit, data = tcga_clin, xlab = "Month", ylab = "Overall survival", surv.median.line = 'hv',risk.table = TRUE)




tcga_rna <- tcga_rna_raw[tcga_rna_raw$Hugo_Symbol!='',]
tcga_rna <- tcga_rna[!duplicated(tcga_rna$Hugo_Symbol),]
rownames(tcga_rna) <- tcga_rna$Hugo_Symbol
tcga_rna <- as.data.frame(t(tcga_rna[-1:-2]))


tcga_clin <- tcga_clin[tcga_clin$SUBTYPE == 'GBM_IDHwt',]
table(tcga_clin$SUBTYPE)
#retrieve RNAs of interest
tcga_rna <- tcga_rna[str_sub(row.names(tcga_rna), end = -4) %in% row.names(tcga_clin), ]

tcga_clin <- tcga_clin[str_sub(row.names(tcga_rna), end = -4),]

max(tcga_rna)
min(tcga_rna)

tcga_rna_log2 <- apply(tcga_rna, c(1, 2), function(value) log2(value + 1))
min(tcga_rna_log2) #0
max(tcga_rna_log2) #19.96911


#dot product

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
table(tcga_clin)

surv_obj <- Surv(time = tcga_clin$OS_MONTHS, 
                 event = tcga_clin$OS_STATUS=="1:DECEASED")
fit <- survfit(surv_obj ~ outcome , data = tcga_clin)
ggsurvplot(fit, data = tcga_clin, xlab = "Days", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)


fit.coxph <- coxph(surv_obj ~ SEX + AGE +outcome +RADIATION_THERAPY, data = tcga_clin)

ggforest(fit.coxph, data = tcga_clin)                                        

colnames(tcga_clin)


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



# Combine both source CGGA and TCGA patients profiles 

#clin columns
tcga_clin$Gender <- tcga_clin$SEX
tcga_clin$Age <- tcga_clin$AGE
tcga_clin$Censor..alive.0..dead.1. <- tcga_clin$OS_STATUS
tcga_clin$OS <- tcga_clin$OS_MONTHS*30
colnames(tcga_clin)
common_columns <- c('Gender','Age','Censor..alive.0..dead.1.','OS')

nano_clin <- list(tcga_clin, clin_325, clin_693)
nano_clin <- do.call(rbind, lapply(nano_clin, function(x) x[,common_columns, drop = FALSE]))

library(ggplot2)

# Age distribution by Gender
ggplot(nano_clin, aes(x = Age, fill = Gender)) +
  geom_histogram(binwidth = 5, alpha = 0.7, position = "dodge") +
  facet_wrap(~ Censor..alive.0..dead.1.) +
  labs(
    title = "Patient Nanogram",
    x = "Age",
    y = "Count"
  ) +
  theme_minimal()

#Check MAPK and p53 activities 
#whether there are groups of patients with single mutations that affect survival significantly
#Results: 1 with high MAPK activity and 269 with low
#Results: 2 with 23 high tp53 activity, but no significant group
#Results:  no significant dif among hypoxia groups

# plotting Kaplan-Mayer curves
pathway = 'Hypoxia'
pathway_data = path_filt$Hypoxia
# sort age 
uni_path = sort(unique(pathway_data))
# store results
results_path = matrix(1,length(uni_path))

# do a for loop for every unique value of age mat
for (i in 2:(length(uni_path)-1)){ # Starting from 2 because the first element would yield no elements higher than.
  path_i = 1*(pathway_data>uni_path[i])
  # survdiff is the function from the survival package 
  logrank = survdiff(surv_filt ~ path_i)
  # store in results_age
  results_path[i] = logrank$pvalue
}
# Plot unique elements of age against p-value
plot(uni_path, results_path, log = 'y')
# Select minimum P-value
min_p_path = which.min(results_path)
# here are 1 good thresholds, -1 
opt_thr = uni_path[min_p_path]
#opt_JAK = opt_thr
#opt_thr = 0.0
# I recalculated the P-value as a sanity check
pval = survdiff(surv_filt ~ pathway_data>opt_thr)$pvalue
nplus = sum(pathway_data>opt_thr)   # how many patients we have in high group
nminus = sum(pathway_data<opt_thr)   # how many patients we have in low group
# fit Kaplan Meier model
path_filt <- path_filt %>% mutate(Hypoxia = ifelse(Hypoxia >= 1, "high", ifelse(Hypoxia < 1.001*opt_thr,"low","intermediate")))
nhigh = sum(pathway_data>1)
ninter = sum((pathway_data<1) & (pathway_data > 1.001*opt_thr))
nlow = sum(pathway_data<1.001*opt_thr)
#KM = survfit(surv_filt ~ pathway_data>opt_thr,data = path_filt)
KM = survfit(surv_filt ~ Hypoxia,data = path_filt)
# Plot Kaplan Meier 
#plot(KM, lwd = 3, col = c(1,2), cex.axis = 1.5, xlab = 'Months', ylab = 'Survival Probability' , cex.lab = 1.5)
#ggsurvplot(KM, data = path_filt,pval = TRUE,xlab = 'Overall survival time, months',
#           legend.labs=c(paste('Low ',pathway,' activity, ',nminus,' patient(s)',sep = ''),paste('High ',pathway,' activity, ',nplus,' patient(s)',sep = '')),
#           palette = c('blue','red'),legend.title="")
p <- ggsurvplot(KM, data = path_filt,pval = TRUE,xlab = 'Overall survival time, months',
                legend.labs=c(paste("High MAPK activity,\n",nhigh," patients",sep=""),paste("Intermediate MAPK activity,\n",ninter," patients",sep=""),paste("Low MAPK activity,\n",nlow," patients",sep="")),
                legend.title=""
)

ggpar(p, 
      font.main = c(13, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.caption = c(16, "bold"), 
      font.legend = c(13, "bold"), 
      font.tickslab = c(14, "bold"))

