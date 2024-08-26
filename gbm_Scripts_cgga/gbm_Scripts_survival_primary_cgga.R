library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)


setwd("/Users/lidiayung/PhD_project/project_GBM/gbm_DATA/gbm_DATA_CGGA")

clin_raw<- read.delim("CGGA.mRNAseq_693_clinical.20200506.txt", sep = '\t')#,skip = 4,row.names = 'PATIENT_ID')
table(clin_raw$PRS_type)

clin <- clin_raw[clin_raw$PRS_type =="Primary",]

table(clin$IDH_mutation_status)

clin <- clin[clin$IDH_mutation_status %in% c('Mutant'),]


RNA_raw <- read.delim("CGGA.mRNAseq_693.RSEM-genes.20200506.txt",check.names = FALSE)
RNA_raw[is.na(RNA_raw)] <- 0

RNA_raw <- RNA_raw[RNA_raw$Gene_Name!='',]
RNA_raw <- RNA_raw[!duplicated(RNA_raw$Gene_Name),]

rownames(RNA_raw) <- RNA_raw$Gene_Name
RNA <- as.data.frame(t(RNA_raw[-1]))
row.names(clin) <- clin$CGGA_ID
# Align clinical data:
#clin <- clin_raw[row.names(clin_raw)%in% str_sub(row.names(RNA), end = -4),]
RNA <- RNA[str_sub(row.names(RNA)) %in% row.names(clin), ]
RNA <- apply(RNA, c(1, 2), function(value) log2(value + 1))
RNA <- as.data.frame(scale(RNA, center = TRUE, scale = TRUE))

surv_obj <- Surv(time = clin$OS, 
                 event = clin$Censor..alive.0..dead.1.=="1")

#surv_obj 
fit <- survfit(surv_obj ~ 1, data = clin)
ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival", surv.median.line = 'hv',risk.table = TRUE)


surv_summary <- summary(fit)
# Calculate the median survival time
median_surv <- summary(fit)$table["median"]
median_surv

p <- ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival", 
                surv.median.line = 'hv', risk.table = TRUE)

# Customize the x-axis to include the value 13.61 along with other specified breaks
p$plot <- p$plot + 
  scale_x_continuous(breaks = c(0, 13.78, 20, 40, 60, 80))

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
fit <- survfit(surv_obj ~ Gender , data = clin)
ggsurvplot(fit, data = clin, xlab = "Days", ylab = "Overall survival",pval = TRUE,
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

fit.coxph <- coxph(surv_obj ~ Gender + Age , data = clin)
ggforest(fit.coxph, data = clin)

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

paste(names(fit.coxph$coefficients[order(fit.coxph$coefficients,decreasing = TRUE)]),collapse = '+')

fit.coxph <- coxph(surv_obj ~ `JAK-STAT`+Trail+MAPK+NFkB+Hypoxia+TGFb+EGFR+p53+WNT+PI3K+VEGF+Androgen+Estrogen+TNFa,data = path_df)

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
clin_filt <- clin[clin$OS > 0,]
RNA_filt <- RNA[clin$OS > 0,]
path_filt <- path_df[clin$OS > 0,]
# create a survival object consisting of times & censoring
surv_filt <- Surv(time = clin_filt$OS, 
                 event = clin_filt$Censor..alive.0..dead.1.=="1")
fit <- survfit(surv_filt ~ 1, data = clin_filt)
ggsurvplot(fit, data = clin_filt, xlab = "Day", ylab = "Overall survival")
which(is.na(surv_filt))
RNA_filt[is.na(RNA_filt)] <- 0
paste(which(is.na(surv_filt)),collapse = ' , ')
surv_filt <-surv_filt[-c(6 , 9 , 87 , 144 , 166 , 167 , 169 , 179 , 182 , 185 , 204),]
surv_filt
RNA_filt <- RNA_filt[-c(6 , 9 , 87 , 144 , 166 , 167 , 169 , 179 , 182 , 185 , 204),]

genes_info <- read.delim('/Users/lidiayung/PhD_project/project_UCD_blca/blca_DATA/blca_DATA_LINCS/geneinfo_beta.txt')

genes_lm <- genes_info[genes_info$feature_space %in%'landmark', ]$gene_symbol

RNA_filt <- RNA_filt[colnames(RNA_filt) %in%genes_lm]

clin_filt<-clin_filt[-c(6 , 9 , 87 , 144 , 166 , 167 , 169 , 179 , 182 , 185 , 204),]


# glmnet
library("glmpath")
library("glmnet")
library("penalized")

fit_glm <- glmnet(RNA_filt,surv_filt,family="cox") # , alpha = 1, standardize = TRUE, maxit = 1000
print(fit_glm)

cvfit <- cv.glmnet(data.matrix(RNA_filt),surv_filt,family="cox",type.measure = "C")
plot(cvfit)

#121 72.43 0.003615
# analysing results
cfs = coef(fit_glm,s=  0.003615)

meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
meaning_vals

coef_data <-data.frame(variable=meaning_coefs,coefficient=meaning_vals)

sorted_coef_abs <- coef_data[order(-coef_data$coefficient,decreasing = FALSE), ]
sorted_coef_abs
top_100 <- tail(sorted_coef_abs,15)
list(top_100$variable)
#write.csv(meaning_vals,file='/Users/lidiayung/PhD_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_cgga/gbm_cgga_survival_coeffs_glmnet.csv')
print(paste(meaning_coefs,collapse=" + "))
# cutting to find most important
ncut = 15
vals_surv = sort(abs(meaning_vals))[1:ncut]
print(paste(top_100$variable,collapse=" + "))
# if we want to run coxph we have to copy-paste the string
fit.coxph <- coxph(surv_filt ~ SPDEF + LYN + NMT1 + PRSS23 + KIT + DFFA + SIRT3 + ARL4C + ALDH7A1 + MEST + LAMA3 + INPP4B + RRP1B + SERPINE1 + ABCF1, data = RNA_filt)
ggforest(fit.coxph, data = RNA_filt)


top_100 <- tail(sorted_coef_abs,15)
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
rna_raw <- read.delim("CGGA.mRNAseq_693.RSEM-genes.20200506.txt",check.names = FALSE)

rna_raw[is.na(rna_raw)] <- 0
rna_raw <- rna_raw[rna_raw$Gene_Name!='',]
rna_raw <- rna_raw[!duplicated(rna_raw$Gene_Name),]
rownames(rna_raw) <- rna_raw$Gene_Name
rna <- as.data.frame(t(rna_raw[-1]))

#retrieve RNAs of interest
rna <- rna[str_sub(row.names(rna)) %in% row.names(clin_filt), ]
#rna <- rna[-143,]

clin_filt$Censor..alive.0..dead.1.
max(rna)
min(rna)

#rna <- rna[-261,]
rna_log2 <- apply(rna, c(1, 2), function(value) log2(value + 1))
min(rna_log2) #0
max(rna_log2) #19.96911


#dot product

coef_data

result <- rna_log2 %*% coef_data_subset$DPD_survival

clin_filt$DPD_prognosis <- result[row.names(clin_filt) %in% str_sub(row.names(result)), ]

hist(clin_filt$DPD_prognosis)
#assign to positive,negative (I only have positive and decided to select a positive threshold)
clin_filt$outcome <- ifelse(clin_filt$DPD_prognosis > -25, "Positive",
                       ifelse(clin_filt$DPD_prognosis == 0, "0", "Negative"))

### hazard ratio
surv_obj_filt <- Surv(time = clin_filt$OS, 
                 event = clin_filt$Censor..alive.0..dead.1.=="1")


fit <- survfit(surv_obj_filt ~ outcome , data = clin_filt)
ggsurvplot(fit, data = clin_filt, xlab = "Day", ylab = "Overall survival",pval = TRUE,
           risk.table =TRUE)

clin_filt$Gender <- trimws(clin_filt$Gender)

table(clin_filt$MGMTp_methylation_status)
fit.coxph <- coxph(surv_obj_filt ~ Gender+Age+MGMTp_methylation_status + outcome, data = clin_filt)
ggforest(fit.coxph, data = clin_filt)

# Create bar plot
ggplot(coef_data, aes(x = reorder(variable, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Genes", y = "Coefficient") +
  theme_minimal()
which(is.na(path_filt))
dim(path_filt)
surv_obj
path_filt <- path_filt[-c(143),]

# glmnet over progeny
fit_glm <- glmnet(path_filt,surv_obj,family="cox") # , alpha = 1, standardize = TRUE, maxit = 1000
print(fit_glm)

cv_fit <- cv.glmnet(data.matrix(path_filt),surv_filt,family="cox",type.measure = "C")
plot(cv_fit)
# analysing results
cfs = coef(fit_glm,s=0.011020)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
write.csv(meaning_vals,file='gbm_OUTPUT_glm_survival_coeffs.csv')
gbmprint(paste(meaning_coefs,collapse=" + "))
# cutting to find most important
ncut = 10
vals_surv = sort(abs(meaning_vals),decreasing = TRUE)[1:ncut]
print(paste(names(vals_surv),collapse=" + "))
# if we want to run coxph we have to copy-paste the string
fit.coxph <- coxph(surv_obj ~ EGFR + Trail + Estrogen + Androgen + VEGF + PI3K + `JAK-STAT`, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)


# plotting Kaplan-Mayer curves
pathway = '`JAK-STAT`'
pathway_data = path_filt$`JAK-STAT`

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
min_p_path
# here are 2 good thresholds, -1 and 1
opt_thr = uni_path[min_p_path]
opt_thr
#opt_JAK = opt_thr
#opt_thr = 0.0
# I recalculated the P-value as a sanity check
pval = survdiff(surv_filt ~ pathway_data>opt_thr)$pvalue
nplus = sum(pathway_data>opt_thr)   # how many patients we have in high group
nminus = sum(pathway_data<opt_thr)   # how many patients we have in low group
# fit Kaplan Meier model
path_filt <- path_filt %>% mutate(JAK_STAT_group = ifelse(`JAK-STAT` >= 2, "high", ifelse(`JAK-STAT` < 2.001*opt_thr,"low","intermediate")))
path_filt$`JAK-STAT`
nhigh = sum(pathway_data>2)
ninter = sum((pathway_data<2) & (pathway_data > 2.001*opt_thr))
nlow = sum(pathway_data<1.001*opt_thr)
#KM = survfit(surv_filt ~ pathway_data>opt_thr,data = path_filt)
KM = survfit(surv_filt ~ JAK_STAT_group,data = path_filt)

# Plot Kaplan Meier 
#plot(KM, lwd = 3, col = c(1,2), cex.axis = 1.5, xlab = 'Months', ylab = 'Survival Probability' , cex.lab = 1.5)
#ggsurvplot(KM, data = path_filt,pval = TRUE,xlab = 'Overall survival time, months',
#           legend.labs=c(paste('Low ',pathway,' activity, ',nminus,' patient(s)',sep = ''),paste('High ',pathway,' activity, ',nplus,' patient(s)',sep = '')),
#           palette = c('blue','red'),legend.title="")
ggsurvplot(KM, data = path_filt,pval = TRUE,xlab = 'Overall survival time, months')#,
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
