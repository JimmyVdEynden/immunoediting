#################################################
# manuscript_Remp.R
#################################################

# Calculate ratio of obs/exp (empirical pan-cancer approach nr of epitopes (Remp)
# Code was double checked using python code Alejandro: cfr scratch/check_Rooney_metric.R

# Load 
#######
{
  # Observed Mutation data
  TCGA_maf<- readRDS("data/TCGA_maf.rds")
  TCGA_maf$isHA<- TCGA_maf[,"mut_HLA_mean_aff"]<500

  # Exclude Non-expressed & driver genes
  TCGA_maf<- TCGA_maf[TCGA_maf[,"mRNA"]>0&!is.na(TCGA_maf[,"mRNA"]),] 
  driver_genes<- readRDS("data/CGC_v83.rds")
  TCGA_maf<- TCGA_maf[!TCGA_maf[,1]%in%driver_genes,]
  
  # Focus on cancers from study
  cancers<- sort(c('BRCA','COADREAD','LUAD','OV','HNSC','UCEC','LUSC','PRAD','LGG','STAD','SKCM','GBM','BLCA','LIHC','CESC','KIRP','SARC','PAAD','ESCA'))
}

# Load functions
################
{
  source("scripts/functions/calculate_R.R")
  source("scripts/functions/cancer_cols.R")
  library(svglite)
}

# Get metrics per sample
########################
{
  # Rsim (calculated using simulated mutations - observed data)
  R_ls<- calculate_R(TCGA_maf$Tumor_Sample_Barcode,
                     TCGA_maf$Variant_Classification,
                     TCGA_maf$isHA,
                     TCGA_maf$subst_type3
  )
  R<- R_ls[["R"]]
}

# Plot R
############
{
  # Get cancers for each sample 
  cancer_id<- TCGA_maf$Cancer
  names(cancer_id)<- TCGA_maf$Tumor_Sample_Barcode
  cancer_id<- cancer_id[!duplicated(names(cancer_id))]
  
  # Calculate median & p values per cancer
  R_med<- sort(tapply(R,cancer_id[names(R)],"median",na.rm=T)[cancers],decreasing = T)
  R_p<- tapply(R,cancer_id[names(R)],function(x) wilcox.test(x,mu = 1)$p.value)[names(R_med)]
  R_q<- p.adjust(R_p,"fdr")

  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_Remp.svg"),width=10)
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_Remp.pdf"),width=10)
    par(mfrow=c(2,3))
    boxplot(R~factor(cancer_id[names(R)],levels = names(R_med)),beside = T,las=2,outline=F,horizontal=T,frame.plot=F,col=c25[names(R_med)],staplewex=0,notch=T,xlab="R")
    if(sum(R_q<0.1)!=0) text(rep(2,sum(R_q<0.1)),which(R_q<0.1),paste0("P=",signif(R_p[which(R_q<0.1)],0.3)),xpd=NA)
    abline(v=1,lty=2)
    dev.off()
  }
  
  # Mut
  results_t<- cbind(median=R_med[cancers],p=R_p[cancers],q=p.adjust(R_p[cancers],"fdr"))
  results_t<- results_t[order(results_t[,"median"]),]
  print.data.frame(as.data.frame(format(results_t,digits=3)))
}

# Save results
###############
{
  save(R,file = paste0("results/data/fig",fig_nr,"Remp_.RData"))
}
