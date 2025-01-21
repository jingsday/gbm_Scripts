#Not this one either
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)


setwd("/Users/lidiayung/project/project_gbm/gbm_DATA/gbm_DATA_tcga_pan_can_atlas_2018")
paste0(wkdir,"gbm_DATA_TCGA/data_mrna_seq_v2_rsem.txt")

clin_raw<- read.delim("/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_TCGA/data_clinical_patient.txt", sep = '\t',skip = 4)
rownames(clin_raw) <- clin_raw$PATIENT_ID

RNA_raw <- read.delim("/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_TCGA/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",check.names = FALSE)
RNA_raw[is.na(RNA_raw)] <- 0
RNA_raw <- RNA_raw[RNA_raw$Hugo_Symbol!='',]
RNA_raw <- RNA_raw[!duplicated(RNA_raw$Hugo_Symbol),]
rownames(RNA_raw) <- RNA_raw$Hugo_Symbol
RNA <- as.data.frame(t(RNA_raw[-1:-2]))
# Align clinical data:
clin <- clin_raw[str_sub(row.names(RNA), end = -4),]

clin$SAMPLE_TYPE ="Primary"
#recurrent primary
clin_sample <- read.delim("/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_TCGA/data_clinical_sample.txt", sep = '\t',skip = 4)

recurrent <-clin_sample[clin_sample$SAMPLE_TYPE=="Recurrence",]
recurrent$PATIENT_ID

clin[recurrent$PATIENT_ID,]$SAMPLE_TYPE='Recurrence'
table(clin$SAMPLE_TYPE)


# create a survival object consisting of times & censoring
surv_obj <- Surv(time = clin$OS_MONTHS, 
                 event = clin$OS_STATUS=="1:DECEASED")


#surv_obj 
fit <- survfit(surv_obj ~ 1, data = clin)
p<-ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival", surv.median.line = 'hv',risk.table = TRUE)


surv_summary <- summary(fit)
# Calculate the median survival time
median_surv <- summary(fit)$table["median"]
median_surv

p <- ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival", 
                surv.median.line = 'hv', risk.table = TRUE)

# Customize the x-axis to include the value 13.61 along with other specified breaks
p$plot <- p$plot + 
  scale_x_continuous(breaks = c(0, 13.61, 20, 40, 60, 80))

# Print the plot
print(p)

# Fit a survival curve
fit <- survfit(surv_obj ~ 1, data = clin)


# Print the median survival time
print(paste("Median Overall Survival:", median_surv, "months"))

plot_title <- paste("KM Survival Curve Median OS:", round(median_surv,1), "months")

# Plot the survival curve with the median survival line
ggsurvplot(fit, data = clin, 
           xlab = "Month", 
           ylab = "Overall Survival", 
           surv.median.line = 'hv', # Add horizontal and vertical lines at the median
           risk.table = TRUE,
           title=plot_title)       # Add a risk table below the plot


#primary and recurrent
fit <- survfit(surv_obj ~ SAMPLE_TYPE , data = clin)
ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)


# fit multivariate model (COX proportional hazard) 
fit.coxph <- coxph(surv_obj ~ EGFR + ERBB2 + ERBB3 + ERBB4 + NRG1 + NRG2 + NRG3 + NRG4, 
                   data = RNA)

ggforest(fit.coxph, data = RNA)
#new.coxph <- predict(fit.coxph, newdata = RNA[,c("KRAS", "ERBB2")])
#new.coxph
fit.coxph <- coxph(surv_obj ~ KSR1 + KSR2 + IQGAP1 + GAB1 + GAB2, 
                   data = RNA)
ggforest(fit.coxph, data = RNA)

fit.coxph <- coxph(surv_obj ~ SEX + AGE + SAMPLE_TYPE, data = clin)
ggforest(fit.coxph, data = clin)

### Remove n_recurrent=7 ###

library(progeny)
zscores = as.matrix(t(RNA))
pathways <- progeny(zscores, scale=TRUE, organism="Human")  #, top = 100, perm = 1)
path_df = as.data.frame(pathways)
# fit multivariate model (COX proportional hazard) 
colnames(path_df)

fit.coxph <- coxph(surv_obj ~  EGFR + Hypoxia   + TGFb  + `JAK-STAT` + TNFa, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)

fit.coxph <- coxph(surv_obj ~ NFkB + Androgen +  VEGF + `JAK-STAT` + TGFb, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)

fit.coxph <- coxph(surv_obj ~ MAPK + PI3K  +  p53 + TGFb  + `JAK-STAT` + TNFa + 
                     WNT + Androgen +  Estrogen + Hypoxia +  Trail + VEGF, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)

fit.coxph <- coxph(surv_obj ~ Androgen + EGFR + Estrogen + Hypoxia + `JAK-STAT` + 
                     MAPK + NFkB + p53 + PI3K + TGFb + TNFa + Trail + VEGF + WNT,
                   data = path_df)

ggforest(fit.coxph, data = path_df)

fit.coxph$coefficients[order(fit.coxph$coefficients,decreasing = TRUE)]

fit.coxph <- coxph(surv_obj ~ MAPK+Trail+EGFR+Androgen+Estrogen+TNFa+ `JAK-STAT`+
                     TGFb+Hypoxia+WNT+VEGF+NFkB+p53+PI3K,data = path_df)

ggforest(fit.coxph, data = path_df)

library(corrplot)
corrplot(cor(path_df), type = "upper", order = "hclust",tl.col = "black", tl.srt = 45)


# analysing on what MAPK activity depend
MAPK_df = cbind(MAPK = path_df$MAPK, PI3K = path_df$PI3K,
                KSR1 = RNA$KSR1, KSR2 = RNA$KSR2, IQGAP1 = RNA$IQGAP1, 
                IQGAP2 = RNA$IQGAP2, IQGAP3 = RNA$IQGAP3, GAB1 = RNA$GAB1, GAB2 = RNA$GAB2,
                KRAS = RNA$KRAS, NRAS = RNA$NRAS, HRAS = RNA$HRAS, BRAF=RNA$BRAF, 
                CRAF=RNA$RAF1, ARAF=RNA$ARAF)
corrplot(cor(MAPK_df), type = "upper", order = "hclust",tl.col = "black", tl.srt = 45)

# now creating object without zero times
clin_filt <- clin[clin$OS_MONTHS > 0,]
RNA_filt <- RNA[clin$OS_MONTHS > 0,]
path_filt <- path_df[clin$OS_MONTHS > 0,]
# create a survival object consisting of times & censoring
surv_filt <- Surv(time = clin_filt$OS_MONTHS, 
                 event = clin_filt$OS_STATUS=="1:DECEASED")
fit <- survfit(surv_filt ~ 1, data = clin_filt)
ggsurvplot(fit, data = clin_filt, xlab = "Month", ylab = "Overall survival")

which(is.na(surv_filt))
RNA_filt[is.na(RNA_filt)] <- 0
surv_filt
RNA_filt <- RNA_filt[-c(127),]

surv_filt<- surv_filt[-c(127),]
clin_filt<-clin_filt[-c(127),]



recurrence <- clin_sample[clin_sample$SAMPLE_TYPE=="Recurrence",]

clin
# glmnet
library("glmpath")
library("glmnet")
library("penalized")

fit_glm <- glmnet(RNA_filt,surv_filt,family="cox") # , alpha = 1, standardize = TRUE, maxit = 1000
print(fit_glm)

# analysing results
cfs = coef(fit_glm,s= 0.00769)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
meaning_vals

coef_data <-data.frame(variable=meaning_coefs,coefficient=meaning_vals)

sorted_coef_abs <- coef_data[order(-coef_data$coefficient), ]
sorted_coef_abs
top_100 <- head(sorted_coef_abs,15)
list(top_100$variable)
#write.csv(meaning_vals,file='/Users/lidiayung/Downloads/gbm_survival_coeffs_glmnet.csv')
print(paste(meaning_coefs,collapse=" + "))
# cutting to find most important
ncut = 15
vals_surv = sort(abs(meaning_vals),decreasing = TRUE)[1:ncut]
print(paste(top_100$variable,collapse=" + "))
# if we want to run coxph we have to copy-paste the string
fit.coxph <- coxph(surv_filt ~ LCN9+FAHD2A+RTP2+TGM7+EN2+PCDHGB2+OR2AG1+C3orf31+GALE+LENEP+C2orf57+CTAGE6+TMEM233+HOXB2+GPR110, data = RNA_filt)
ggforest(fit.coxph, data = RNA_filt)


top_100 <- head(sorted_coef_abs,15)
list(top_100$variable)

# Assuming top_100 is a data frame containing the top 100 coefficients
paste(top_100$variable, collapse = "+")


# Create bar plot
ggplot(top_100, aes(x = reorder(variable, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Genes", y = "Coefficient") +
  theme_minimal()


#
rna_raw <- read.delim("data_mrna_seq_v2_rsem.txt",check.names = FALSE)


rna_raw[is.na(rna_raw)] <- 0
rna_raw <- rna_raw[rna_raw$Hugo_Symbol!='',]
rna_raw <- rna_raw[!duplicated(rna_raw$Hugo_Symbol),]
rownames(rna_raw) <- rna_raw$Hugo_Symbol
rna <- as.data.frame(t(rna_raw[-1:-2]))

#retrieve RNAs of interest
rna <- rna[str_sub(row.names(rna), end = -4) %in% row.names(clin_raw), ]
rna <- rna[-127,]


max(rna)
min(rna)

#rna <- rna[-261,]
rna_log2 <- apply(rna, c(1, 2), function(value) log2(value + 1))
min(rna_log2) #0
max(rna_log2) #19.96911


#dot product
rna_log2_subset <- rna_log2[, coef_data$variable]
coef_data_subset <-coef_data[,-1]

result <- rna_log2_subset %*% coef_data_subset


clin_filt$DPD_prognosis <- result[str_sub(row.names(result), end = -4) %in% row.names(clin_filt), ]
#assign to positive,negative (I only have positive and decided to select a positive threshold)
clin_filt$outcome <- ifelse(clin_filt$DPD_prognosis > 110, "High",
                       ifelse(clin_filt$DPD_prognosis == 0, "0", "Low"))

hist(clin_filt$DPD_prognosis)


### hazard ratio
surv_obj_filt <- Surv(time = clin_filt$OS_MONTHS, 
                 event = clin_filt$OS_STATUS=="1:DECEASED")


fit <- survfit(surv_obj_filt ~ outcome , data = clin_filt)
ggsurvplot(fit, data = clin_filt, xlab = "Month", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)


fit.coxph <- coxph(surv_obj_filt ~ SEX + AGE + SAMPLE_TYPE+outcome, data = clin_filt)
ggforest(fit.coxph, data = clin_filt)

# Create bar plot
ggplot(coef_data, aes(x = reorder(variable, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Genes", y = "Coefficient") +
  theme_minimal()

path_filt <- path_filt[-which(is.na(path_filt)),]
surv_obj_filt <-
which(is.na(surv_obj_filt))
# glmnet over progeny
fit_glm <- glmnet(path_filt,surv_obj_filt,family="cox") # , alpha = 1, standardize = TRUE, maxit = 1000
print(fit_glm)
# analysing results
cfs = coef(fit_glm,s=0.007320)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
#write.csv(meaning_vals,file='survival_coeffs_progeny_glmnet.csv')
print(paste(meaning_coefs,collapse=" + "))
# cutting to find most important
#ncut = 10
#vals_surv = sort(abs(meaning_vals),decreasing = TRUE)[1:ncut]
#print(paste(names(vals_surv),collapse=" + "))
# if we want to run coxph we have to copy-paste the string
fit.coxph <- coxph(surv_filt ~ Androgen + EGFR + Estrogen + Hypoxia + `JAK-STAT` + MAPK + p53 + PI3K + TGFb + Trail + VEGF + WNT, data = path_filt)
ggforest(fit.coxph, data = RNA_filt)


# plotting Kaplan-Mayer curves
pathway = 'MAPK'
pathway_data = path_filt$MAPK
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
# here are 2 good thresholds, -1 and 1
opt_thr = uni_path[min_p_path]
#opt_JAK = opt_thr
#opt_thr = 0.0
# I recalculated the P-value as a sanity check
pval = survdiff(surv_filt ~ pathway_data>opt_thr)$pvalue
nplus = sum(pathway_data>opt_thr)   # how many patients we have in high group
nminus = sum(pathway_data<opt_thr)   # how many patients we have in low group
# fit Kaplan Meier model
path_filt <- path_filt %>% mutate(MAPK_group = ifelse(MAPK >= 1, "high", ifelse(MAPK < 1.001*opt_thr,"low","intermediate")))
nhigh = sum(pathway_data>1)
ninter = sum((pathway_data<1) & (pathway_data > 1.001*opt_thr))
nlow = sum(pathway_data<1.001*opt_thr)
#KM = survfit(surv_filt ~ pathway_data>opt_thr,data = path_filt)
KM = survfit(surv_filt ~ MAPK_group,data = path_filt)
# Plot Kaplan Meier 
#plot(KM, lwd = 3, col = c(1,2), cex.axis = 1.5, xlab = 'Months', ylab = 'Survival Probability' , cex.lab = 1.5)
#ggsurvplot(KM, data = path_filt,pval = TRUE,xlab = 'Overall survival time, months',
#           legend.labs=c(paste('Low ',pathway,' activity, ',nminus,' patient(s)',sep = ''),paste('High ',pathway,' activity, ',nplus,' patient(s)',sep = '')),
#           palette = c('blue','red'),legend.title="")
p <- ggsurvplot(KM, data = path_filt,pval = TRUE,xlab = 'Overall survival time, months')#,
           #legend.labs=c(paste("High MAPK activity,\n",nhigh," patients",sep=""),paste("Intermediate MAPK activity,\n",ninter," patients",sep=""),paste("Low MAPK activity,\n",nlow," patients",sep="")),
           #legend.title=""
           #)

ggpar(p, 
      font.main = c(13, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.caption = c(16, "bold"), 
      font.legend = c(13, "bold"), 
      font.tickslab = c(14, "bold"))


# RF
library("randomForestSRC")
# building RF model
B <- 1000
# Building a RSF
status_surv <- clin_filt$OS_STATUS=="1:DECEASED"
time_surv <- clin_filt$OS_MONTHS
dataSetRF <- cbind(time_surv,status_surv,RNA_filt)

names(dataSetRF)<-make.names(names(dataSetRF))

RF_obj <- rfsrc(Surv(time_surv,status_surv)~., dataSetRF,  ntree = B,  membership = TRUE, importance=TRUE)
# Printing the RF object  
print(RF_obj)

# Vadiable importance
jk.obj <- subsample(RF_obj)
#pdf("VIMPsur.pdf", width = 15, height = 20)
#par(oma = c(0.5, 10, 0.5, 0.5))
#par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,17,1,1), mgp = c(4, 1, 0))
#plot(jk.obj, xlab = "Variable Importance (x 100)", cex = 1.2)
plot(jk.obj)
#dev.off()

fit.coxph <- coxph(surv_filt ~ TRIM67 + RASA1 + RHOF + NCAPD3 + NEU2 + ARNTL2 + CDK6, data = RNA_filt)
ggforest(fit.coxph, data = RNA_filt)
