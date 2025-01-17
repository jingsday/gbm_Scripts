
wkdir<- ("/home/jing/Phd_project/project_GBM/gbm_DATA/")
outdir<- ('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_survival/')

#clin
tcga_clin_raw <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_clinical_patient.txt"), 
                            sep = '\t',skip = 4,row.names = 'PATIENT_ID')

#clin profiles
clin_sample <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_clinical_sample.txt"), sep = '\t',skip = 4)
table(clin_sample$SAMPLE_TYPE)
# recurrent <-clin_sample[clin_sample$SAMPLE_TYPE=="Recurrence",]
# recurrent$PATIENT_ID
# 
# tcga_clin_raw$SAMPLE_TYPE ="Primary"
# tcga_clin_raw[recurrent$PATIENT_ID,]$SAMPLE_TYPE='Recurrence'
# table(tcga_clin_raw$SAMPLE_TYPE)
# 
# tcga_clin_raw <-tcga_clin_raw[tcga_clin_raw$SAMPLE_TYPE=='Primary',]
tcga_clin_raw <-tcga_clin_raw[tcga_clin_raw$SUBTYPE=='GBM_IDHwt',]
tcga_clin_raw <- tcga_clin_raw[tcga_clin_raw$SEX != "", ]

#RNA
tcga_RNA_raw <- read.delim(paste0(wkdir,"gbm_DATA_TCGA/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"),check.names = FALSE)
tcga_RNA_raw[is.na(tcga_RNA_raw)] <- 0
tcga_RNA_raw <- tcga_RNA_raw[tcga_RNA_raw$Hugo_Symbol!='',]
tcga_RNA_raw <- tcga_RNA_raw[!duplicated(tcga_RNA_raw$Hugo_Symbol),]
rownames(tcga_RNA_raw) <- tcga_RNA_raw$Hugo_Symbol
tcga_RNA <- as.data.frame(t(tcga_RNA_raw[-1:-2]))

# Align clinical data:
tcga_clin <- tcga_clin_raw[row.names(tcga_clin_raw)%in% str_sub(row.names(tcga_RNA), end = -4),]
tcga_RNA <- tcga_RNA[str_sub(row.names(tcga_RNA), end = -4) %in% row.names(tcga_clin_raw), ]
```