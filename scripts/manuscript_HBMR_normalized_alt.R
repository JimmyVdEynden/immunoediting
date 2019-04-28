#######################################
# manuscript_HBMR_normalized_alt.R
#######################################

# Aim: calculate HBMRnorm for different parameter choices

# Specify whether observed or simulated data
if(useSim) sim<- "_sim"
if(!useSim) sim<- ""

# Cancers from study
cancers<- c('BLCA','BRCA','CESC','COADREAD','ESCA','GBM','HNSC','KIRP','LGG','LIHC','LUAD','LUSC','OV','PAAD','PRAD','SARC','SKCM','STAD','UCEC')

# Expected mut per subst type
subst_type_matrix_ls<- readRDS(file="data/subst_type_matrix_ls.rds")

# Load functions
source("scripts/functions/calculate_HBMR.R")
library(svglite)

# methods
methods<- c("Main","Kd50","best_Kd","expr","driver","excl_hypermut","MC3")

# summarize results
HBMR_norm<- matrix(NA,length(cancers),length(methods),dimnames = list(cancers, methods))
HBMR_norm_CI<- matrix(NA,length(cancers),length(methods),dimnames = list(cancers, methods))

for(method in methods){
  
  # "HLA specific" & "rank" not done because annotation not foreseen to normalize HBMR for this
  # cat(method,"\n")
  
  # Reference parameters
  Kd_cu<- 500
  subst_type_matrix<- subst_type_matrix_ls$triNT
  
  # Get data for specified method
  #############################
  if(method=="Kd50"){
    Kd_cu<- 50 
    subst_type_matrix<- subst_type_matrix_ls$triNT_Kd50
  }
  
  if(method=="MC3"){
    if(useSim) next
    TCGA_maf<- readRDS("data/maf_MC3_GP.rds") # in GPos format, no difference for calculation
    TCGA_maf$Variant_Classification<- TCGA_maf$variant # Change to format of other data
    TCGA_maf$Variant_Classification[TCGA_maf$Variant_Classification=="Missense_Mutation"]<- "nonsynonymous SNV" 
    TCGA_maf$Variant_Classification[TCGA_maf$Variant_Classification=="Silent"]<- "synonymous SNV"
  }
  else TCGA_maf<- readRDS(paste0("data/TCGA_maf",sim,"_mainHLAType.rds"))
  
  if(method=="excl_hypermut"){
    mut_load<- sort(table(TCGA_maf$Tumor_Sample_Barcode))
    samples_excl<- NULL
    for(cancer in cancers){
      samples_can<- unique(TCGA_maf[TCGA_maf$Cancer==cancer,"Tumor_Sample_Barcode"])
      mut_load_can<- sort(mut_load[samples_can])
      # max_can<- 0.05*sum(mut_load_can)
      max_can<- quantile(mut_load_can,0.95)
      samples_excl_temp<- names(mut_load_can)[max_can < mut_load_can]
      samples_excl<- c(samples_excl,samples_excl_temp)
    }
    # length(samples_excl) # 373
    TCGA_maf<- TCGA_maf[!TCGA_maf$Tumor_Sample_Barcode%in%samples_excl,]
  }
  
  if(method=="expr") TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"mRNA"])&TCGA_maf[,"mRNA"]>0,] 
  
  if(method=="driver"){
    driver_genes<- readRDS("data/CGC_v83.rds")
    TCGA_maf<- TCGA_maf[!TCGA_maf[,"Hugo_Symbol"]%in%driver_genes,] 
  }

  if(method=="best_Kd"){
    TCGA_maf$isHA<- rowMins(data.matrix(TCGA_maf[,grep("wt_HLA_.._aff",colnames(TCGA_maf))]))<500
    subst_type_matrix<- subst_type_matrix_ls$triNT_bestBinder
  }
  else if(method!="MC3") TCGA_maf$isHA<- TCGA_maf[,"wt_HLA_mean_aff"]<Kd_cu # Determined already for MC3
  
  # Calculate HBMR 
  #################
  HBMR<- calculate_HBMR(sample = TCGA_maf$Tumor_Sample_Barcode,
                        variant = TCGA_maf$Variant_Classification,
                        isHA = TCGA_maf$isHA,
                        subst_type3 = as.character(TCGA_maf$subst_type3),
                        exp_mut = subst_type_matrix,
                        cancer =  as.character(TCGA_maf$Cancer))[cancers,]
  
  # Combine results from different methods
  ########################################
  HBMR_norm[cancers,method]<- HBMR[cancers,"HBMR_norm"]
  HBMR_norm_CI[cancers,method]<- HBMR[cancers,"HBMR_norm_ciup"]
}

# Plot
#######
{
  # Sort
  HBMR_norm<- HBMR_norm[order(HBMR_norm[,"Main"]),]
  HBMR_norm_CI<- HBMR_norm_CI[rownames(HBMR_norm),]
  
  # Barplot
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"_HBMR_normalized_alt",sim,".pdf"),width = 14)
    else svglite(paste0("results/figs/svg/fig",fig_nr,"_HBMR_normalized_alt",sim,".svg"),width = 14)
    par(oma=c(6,2,2,2))
    # barplot
    bp<- barplot(t(HBMR_norm),las=2,ylab="Normalized HBMR",beside=T,ylim=c(0,1.3))
    # error bars
    arrows(bp,t(HBMR_norm),bp,t(HBMR_norm_CI),lwd = 1.5, angle = 90,code = 2, length = 0.03,xpd=NA)    
    abline(h=1,lty=2)
    legend("topleft",fill=grey.colors(length(methods)),border=NA,legend=methods,bty="n")
    dev.off()
  }
}

# Save results
###############
{
  save(HBMR_norm,HBMR_norm_CI,file = paste0("results/data/fig",fig_nr,"_HBMR_normalized_alt",sim,".RData"))
}
