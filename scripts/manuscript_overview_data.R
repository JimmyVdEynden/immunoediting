###############################
# manuscript_overview_data.R
###############################

# Load data
###########
{
  TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
}

# Libraries
############
{
  library(svglite)
  library(WriteXLS)
}

# 0) Preprocess & Filtering
#############################
{
  # Filter 1) Restrict analysis to synonymous and non-synonymous mutations (indels were removed previously, also remove nonsense)
  TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"Variant_Classification"])&TCGA_maf[,"Variant_Classification"]%in%c("synonymous SNV","nonsynonymous SNV"),]
}

# Cancer abbreviations TCGA
############################
cancers_abbr<- c(
  'LAML',
  'ACC',
  'BLCA',
  'LGG',
  'BRCA',
  'CESC',
  'CHOL',
  'CRC',
  'ESCA',
  'GBM',
  'HNSC',
  'KICH',
  'KIRC',
  'KIRP',
  'LIHC',
  'LUAD',
  'LUSC',
  'DLBC',
  'MESO',
  'OV',
  'PAAD',
  'PCPG',
  'PRAD',
  'SARC',
  'SKCM',
  'STAD',
  'TGCT',
  'THYM',
  'THCA',
  'UCS',
  'UCEC',
  'UVM'
)

cancers_full<- c(
  'Acute Myeloid Leukemia',
  'Adrenocortical carcinoma',
  'Bladder Urothelial Carcinoma',
  'Brain Lower Grade Glioma',
  'Breast invasive carcinoma',
  'Cervical squamous cell carcinoma and endocervical adenocarcinoma',
  'Cholangiocarcinoma',
  'Colorectal cancer',
  'Esophageal carcinoma',
  'Glioblastoma multiforme',
  'Head and Neck squamous cell carcinoma',
  'Kidney Chromophobe',
  'Kidney renal clear cell carcinoma',
  'Kidney renal papillary cell carcinoma',
  'Liver hepatocellular carcinoma',
  'Lung adenocarcinoma',
  'Lung squamous cell carcinoma',
  'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma',
  'Mesothelioma',
  'Ovarian serous cystadenocarcinoma',
  'Pancreatic adenocarcinoma',
  'Pheochromocytoma and Paraganglioma',
  'Prostate adenocarcinoma',
  'Sarcoma',
  'Skin Cutaneous Melanoma',
  'Stomach adenocarcinoma',
  'Testicular Germ Cell Tumors',
  'Thymoma',
  'Thyroid carcinoma',
  'Uterine Carcinosarcoma',
  'Uterine Corpus Endometrial Carcinoma',
  'Uveal Melanoma'
)


# Table with cancer types
#################################
{
  # Get n samples and n mut per cancer
  mut_cancer_t<- table(TCGA_maf[,"Cancer"])
  sample_cancer_t<- table(TCGA_maf[!duplicated(TCGA_maf[,"Tumor_Sample_Barcode"]),"Cancer"])
  names(mut_cancer_t)[names(mut_cancer_t)=="COADREAD"]<- "CRC"
  names(sample_cancer_t)[names(sample_cancer_t)=="COADREAD"]<- "CRC"
  
  # Create table
  cancer_t<- cbind(cancers_abbr,cancers_full,sample_cancer_t[cancers_abbr],mut_cancer_t[cancers_abbr])
  cancer_t<- cancer_t[order(cancer_t[,"cancers_abbr"]),]
  colnames(cancer_t)<- c("TCGA abbreviation","cancer type","n samples","n mutations")
  cancer_t<- as.data.frame(cancer_t)
  
  # Save
  WriteXLS("cancer_t",ExcelFileName = paste0("results/tables/table_",fig_nr,".xls"),row.names = F)
}