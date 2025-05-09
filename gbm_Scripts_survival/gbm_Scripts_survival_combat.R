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
common_genes <- intersect(colnames(cgga_325RNA), colnames(cgga_693RNA))#23271 shared genes
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
