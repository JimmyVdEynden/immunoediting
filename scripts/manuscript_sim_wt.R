###################################
# manuscript_sim_wt.R
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
    if(HLA_method=="harmonic_mean") TCGA_maf_sim<- readRDS("data/TCGA_maf_sim_mainHLAType.rds")
    else if(HLA_method=="rank_harmonic_mean") TCGA_maf_sim<- readRDS("data/TCGA_maf_sim_ranks.rds")
    else if(HLA_method=="bestBinder"){
      TCGA_maf_sim<- readRDS("data/TCGA_maf_sim_mainHLAType.rds")
      bestAff<- rowMins(data.matrix(TCGA_maf_sim[,grep("wt_HLA_.._aff",colnames(TCGA_maf_sim))]),na.rm=T)
      bestAff[is.infinite(bestAff)]<- NA
      TCGA_maf_sim$wt_HLA_mean_aff<- bestAff
      bestAff<- rowMins(data.matrix(TCGA_maf_sim[,grep("mut_HLA_.._aff",colnames(TCGA_maf_sim))]),na.rm=T)
      bestAff[is.infinite(bestAff)]<- NA
      TCGA_maf_sim$mut_HLA_mean_aff<- bestAff
    }
    # else TCGA_maf_sim<- readRDS("data/TCGA_maf_sim_mainHLAType_6.rds")
  }
  else TCGA_maf_sim<- readRDS("data/TCGA_maf_sim.rds")

  # color map
  source("scripts/functions/cancer_cols.R")
  
  # Expression filter
  if(exists("expr_cu")) TCGA_maf_sim<- TCGA_maf_sim[TCGA_maf_sim[,"mRNA"]>expr_cu,] # +/- 25% lost for cu of 10

  # Driver gene filter
  if(exists("driver_genes")) TCGA_maf_sim<- TCGA_maf_sim[!TCGA_maf_sim[,"Hugo_Symbol"]%in%driver_genes,] # +/- 25% lost for cu of 10
  
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
  TCGA_maf_sim<- TCGA_maf_sim[!is.na(TCGA_maf_sim[,"Variant_Classification"])&TCGA_maf_sim[,"Variant_Classification"]%in%c("synonymous SNV","nonsynonymous SNV"),]

  # Turn barcodes in factors
  TCGA_maf_sim[,"Tumor_Sample_Barcode"]<- factor(TCGA_maf_sim[,"Tumor_Sample_Barcode"])

  # Define high affinity HLA binders (default cut-off Kd<500nM)
  isHA_wt<- TCGA_maf_sim[,"wt_HLA_mean_aff"]<Kd_cu
  isHA_mut<- TCGA_maf_sim[,"mut_HLA_mean_aff"]<Kd_cu
  TCGA_maf_sim<- cbind(TCGA_maf_sim,isHA_wt=isHA_wt,isHA_mut=isHA_mut)
}

# 1) n/s Pan cancer analysis
#############################
{
  # General numbers
  TCGA_samples<- unique(as.character(TCGA_maf_sim[,"Tumor_Sample_Barcode"]))
  TCGA_cancers<- unique(TCGA_maf_sim[,"Cancer"])
  cat("n samples:",length(TCGA_samples),"  \n")
  cat("n cancers:",length(TCGA_cancers),"  \n")
  cat("n mutations:",nrow(TCGA_maf_sim),"  \n")
  cat("prop HLA binding:",mean(TCGA_maf_sim[,"isHA_wt"],na.rm=T),"  \n","  \n")

  # Nonsynonymous to synonymous? (ratio n/s) 
  sim_isHA_variant_t<- table(variant=TCGA_maf_sim[,"Variant_Classification"],isHA=TCGA_maf_sim[,"isHA_wt"])
  sim_isHA_variant_t<- sim_isHA_variant_t[c("synonymous SNV","nonsynonymous SNV"),]
  colnames(sim_isHA_variant_t)<- c("non-HLA binding","HLA binding")
  print.data.frame(as.data.frame.matrix(sim_isHA_variant_t))  
  # non-HLA binding HLA binding
  # synonymous SNV             408170      117507
  # nonsynonymous SNV         1043513      267179

  sim_ns<- sim_isHA_variant_t["nonsynonymous SNV",]/sim_isHA_variant_t["synonymous SNV",]
  sim_fisher_test<- fisher.test(sim_isHA_variant_t)
  sim_p<- format(sim_fisher_test$p.value,digits = 3)
  sim_HBMR<- format(sim_fisher_test$estimate,digits = 3)
  sim_prop_HLA<- prop.table(sim_isHA_variant_t,1)[,"HLA binding"] # sil 22.4%; miss 20.4% 

  cat("N/S HLA binding:",sim_ns["HLA binding"],"  \n")
  cat("N/S non-HLA binding:",sim_ns["non-HLA binding"],"  \n")
  cat("p=",sim_p,"  \n")
  cat("HBMR=",sim_HBMR,"  \n")
  cat("prop HLA synonymous:",sim_prop_HLA["synonymous SNV"],"  \n")
  cat("prop HLA non-synonymous:",sim_prop_HLA["nonsynonymous SNV"],"  \n")
}

# 2) n/s Per cancer analysis
#############################
{
  # Calculate
  cancer_t<- table(TCGA_maf_sim[,"Cancer"])
  # cancers<- sort(names(cancer_t)[cancer_t>10000]) # Only the ones with a miniml number of mutations
  cancers<- sort(names(cancer_t)[prop.table(cancer_t)>0.005]) # Same as putting absolute nr of 10000, but like this same cancers used when each HLA allele considered as seperate mutation
  sim_HBMR_isHA<- matrix(NA,length(cancers),3,dimnames = list(cancers,c("p","HBMR","CI")))
  for(cancer in cancers){
    TCGA_maf_sim_cancer<- TCGA_maf_sim[TCGA_maf_sim[,"Cancer"]==cancer,]
    sample_t<- table(as.character(TCGA_maf_sim_cancer[,"Tumor_Sample_Barcode"]),TCGA_maf_sim_cancer[,"isHA_wt"],TCGA_maf_sim_cancer[,"Variant_Classification"])
    fish_res<- fisher.test(rbind(colSums(sample_t[,,"synonymous SNV"]),colSums(sample_t[,,"nonsynonymous SNV"])))
    OR_temp<- fish_res$estimate
    CI_temp<- fish_res$conf.int
    p_temp<- fish_res$p.value
    sim_HBMR_isHA[cancer,"p"]<- p_temp
    sim_HBMR_isHA[cancer,"HBMR"]<- OR_temp
    sim_HBMR_isHA[cancer,"CI"]<- CI_temp[2]
  }
  # sort
  sim_HBMR_isHA<- sim_HBMR_isHA[order(sim_HBMR_isHA[,"HBMR"]),]

  # Print results
  print.data.frame(as.data.frame(format(sim_HBMR_isHA,digits=3)))
                   
  # Plot
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_HBMR_NS_sim.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_HBMR_NS_sim.pdf"),width = 10)
    par(oma=c(6,2,2,2))
    # barplot
    bp<- barplot(sim_HBMR_isHA[,"HBMR"],las=2,ylab="HBMR",beside=T,ylim=c(0,1.2))
    # error bars
    arrows(bp,sim_HBMR_isHA[,"HBMR"],bp,sim_HBMR_isHA[,"CI"],lwd = 1.5, angle = 90,code = 2, length = 0.05,xpd=NA)    
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
  # save(gene_ns_t_sim,sim_isHA_variant_t,sim_ns,sim_prop_HLA,sim_HBMR_isHA,sim_HBMR_isHA_expr,sim_isHA_expr,th_mRNA,sim_isHA_expr_t,sim_expr_ratio,sim_expr_prop_HLA,file = paste0("results/data/",FULL_FILENAME,"Data"))
  save(sim_isHA_variant_t,sim_ns,sim_prop_HLA,sim_HBMR_isHA,file = paste0("results/data/",FULL_FILENAME,"Data"))
}
