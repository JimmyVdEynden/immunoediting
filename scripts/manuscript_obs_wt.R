###################################
# manuscript_obs_wt.R
###################################

# Set parameters
################
{
  if(!exists("HLA_method")) HLA_method<- "harmonic_mean"
  if(!exists("Kd_cu")) Kd_cu<- 500
  if(!exists("HLAType")) HLAType<- "main"
}

# Load data
###########
{
  if(HLAType=="main"){
    if(HLA_method=="harmonic_mean") TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
    else if(HLA_method=="rank_harmonic_mean") TCGA_maf<- readRDS("data/TCGA_maf_ranks.rds")
    else if(HLA_method=="bestBinder"){
      TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
      bestAff<- rowMins(data.matrix(TCGA_maf[,grep("wt_HLA_.._aff",colnames(TCGA_maf))]),na.rm=T)
      bestAff[is.infinite(bestAff)]<- NA
      TCGA_maf$wt_HLA_mean_aff<- bestAff
      bestAff<- rowMins(data.matrix(TCGA_maf[,grep("mut_HLA_.._aff",colnames(TCGA_maf))]),na.rm=T)
      bestAff[is.infinite(bestAff)]<- NA
      TCGA_maf$mut_HLA_mean_aff<- bestAff
    }
    # else TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType_6.rds") 
  }
  else TCGA_maf<- readRDS("data/TCGA_maf.rds")

  # color map
  source("scripts/functions/cancer_cols.R")
  
  # Expression filter
  if(exists("expr_cu")) TCGA_maf<- TCGA_maf[TCGA_maf[,"mRNA"]>expr_cu,] # +/- 25% lost for cu of 10

  # Driver gene filter
  if(exists("driver_genes")) TCGA_maf<- TCGA_maf[!TCGA_maf[,"Hugo_Symbol"]%in%driver_genes,] # +/- 25% lost for cu of 10
  
}

# Libraries
############
{
  library(svglite)
}

# 0) Preprocess & Filtering
#############################
{
  # Filter 1) Restrict analysis to synonymous and non-synonymous mutations (indels were removed previously, also remove nonsense)
  TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"Variant_Classification"])&TCGA_maf[,"Variant_Classification"]%in%c("synonymous SNV","nonsynonymous SNV"),]

  # Turn barcodes in factors
  TCGA_maf[,"Tumor_Sample_Barcode"]<- factor(TCGA_maf[,"Tumor_Sample_Barcode"])

  # Define high affinity HLA binders (default cut-off Kd<500nM)
  isHA_wt<- TCGA_maf[,"wt_HLA_mean_aff"]<Kd_cu
  isHA_mut<- TCGA_maf[,"mut_HLA_mean_aff"]<Kd_cu
  TCGA_maf<- cbind(TCGA_maf,isHA_wt=isHA_wt,isHA_mut=isHA_mut)
}

# 1) n/s Pan cancer analysis
#############################
{
  # General numbers
  TCGA_samples<- unique(as.character(TCGA_maf[,"Tumor_Sample_Barcode"]))
  TCGA_cancers<- unique(TCGA_maf[,"Cancer"])
  cat("n samples:",length(TCGA_samples),"  \n")
  cat("n cancers:",length(TCGA_cancers),"  \n")
  cat("n mutations:",nrow(TCGA_maf),"  \n")
  cat("prop HLA binding:",mean(TCGA_maf[,"isHA_wt"],na.rm=T),"  \n","  \n")

  # Nonsynonymous to synonymous? (ratio n/s) 
  obs_isHA_variant_t<- table(variant=TCGA_maf[,"Variant_Classification"],isHA=TCGA_maf[,"isHA_wt"])
  obs_isHA_variant_t<- obs_isHA_variant_t[c("synonymous SNV","nonsynonymous SNV"),]
  colnames(obs_isHA_variant_t)<- c("non-HLA binding","HLA binding")
  print.data.frame(as.data.frame.matrix(obs_isHA_variant_t))  
  # non-HLA binding HLA binding
  # synonymous SNV             408170      117507
  # nonsynonymous SNV         1043513      267179

  obs_ns<- obs_isHA_variant_t["nonsynonymous SNV",]/obs_isHA_variant_t["synonymous SNV",]
  obs_fisher_test<- fisher.test(obs_isHA_variant_t)
  obs_p<- format(obs_fisher_test$p.value,digits = 3)
  obs_HBMR<- format(obs_fisher_test$estimate,digits = 3)
  obs_prop_HLA<- prop.table(obs_isHA_variant_t,1)[,"HLA binding"] # sil 22.4%; miss 20.4% 

  cat("n/s HLA binding:",obs_ns["HLA binding"],"  \n")
  cat("n/s non-HLA binding:",obs_ns["non-HLA binding"],"  \n")
  cat("p=",obs_p,"  \n")
  cat("HBMR=",obs_HBMR,"  \n")
  cat("prop HLA synonymous:",obs_prop_HLA["synonymous SNV"],"  \n")
  cat("prop HLA non-synonymous:",obs_prop_HLA["nonsynonymous SNV"],"  \n")
}

# 2) n/s Per cancer analysis
#############################
{
  # Calculate
  cancer_t<- table(TCGA_maf[,"Cancer"])
  # cancers<- sort(names(cancer_t)[cancer_t>10000]) # Only the ones with a miniml number of mutations
  cancers<- sort(names(cancer_t)[prop.table(cancer_t)>0.005]) # Same as putting absolute nr of 10000, but like this same cancers used when each HLA allele considered as seperate mutation
  obs_HBMR_isHA<- matrix(NA,length(cancers),3,dimnames = list(cancers,c("p","HBMR","CI")))
  for(cancer in cancers){
    TCGA_maf_cancer<- TCGA_maf[TCGA_maf[,"Cancer"]==cancer,]
    sample_t<- table(as.character(TCGA_maf_cancer[,"Tumor_Sample_Barcode"]),TCGA_maf_cancer[,"isHA_wt"],TCGA_maf_cancer[,"Variant_Classification"])
    fish_res<- fisher.test(rbind(colSums(sample_t[,,"synonymous SNV"]),colSums(sample_t[,,"nonsynonymous SNV"])))
    OR_temp<- fish_res$estimate
    CI_temp<- fish_res$conf.int
    p_temp<- fish_res$p.value
    obs_HBMR_isHA[cancer,"p"]<- p_temp
    obs_HBMR_isHA[cancer,"HBMR"]<- OR_temp
    obs_HBMR_isHA[cancer,"CI"]<- CI_temp[2]
  }
  # sort
  obs_HBMR_isHA<- obs_HBMR_isHA[order(obs_HBMR_isHA[,"HBMR"]),]

  # Print results
  print.data.frame(as.data.frame(format(obs_HBMR_isHA,digits=3)))
                   
  # Plot
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_HBMR_ns.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_HBMR_ns.pdf"),width = 10)
    par(oma=c(6,2,2,2))
    # barplot
    bp<- barplot(obs_HBMR_isHA[,"HBMR"],las=2,ylab="HBMR",beside=T,ylim=c(0,1.2))
    # error bars
    arrows(bp,obs_HBMR_isHA[,"HBMR"],bp,obs_HBMR_isHA[,"CI"],lwd = 1.5, angle = 90,code = 2, length = 0.05,xpd=NA)    
    abline(h=1,lty=2)
    dev.off()
  }
}

# Save all data
################
{
  FULL_FILENAME <- parent.frame(2)$ofile # Only works when sourced, otherwise is.null
  FULL_FILENAME<- gsub("scripts/","",FULL_FILENAME)
  FULL_FILENAME<- paste0("fig",fig_nr,"_",FULL_FILENAME)
  save(obs_isHA_variant_t,obs_ns,obs_prop_HLA,obs_HBMR_isHA,file = paste0("results/data/",FULL_FILENAME,"Data"))
}
