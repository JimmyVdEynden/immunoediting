###############################################
# manuscript_cor_cancerTC_HBMR.R
###############################################

# Correlate proportion of TC mutations for each cancer to its HBMR

# Load 
#######
{
  load("results/data/fig2_manuscript_obs_wt.RData")
  TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
  source("scripts/functions/cancer_cols.R")
}

# Plot
######
{
  # prop TC per cancer
  TCGA_mutProb<- prop.table(table(as.character(TCGA_maf[,"subst_type3"]),TCGA_maf[,"Cancer"]),2)
  cancer_TC<- colSums(TCGA_mutProb[grepl("TC.>T..",rownames(TCGA_mutProb)),])
  
  # HBMR per cancer
  cancer_HBMR<- obs_HBMR_isHA[,"HBMR"]
  
  # figure
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"_cor_cancerTC_HBMR.pdf"))
    else svglite(paste0("results/figs/svg/fig",fig_nr,"_cor_cancerTC_HBMR.svg"))
    plot(cancer_TC[names(cancer_HBMR)],cancer_HBMR,pch=16,frame.plot=F,xlab="%TCN>TNN",ylab="HBMR",col=c25[names(cancer_HBMR)],xlim=c(0,0.7),ylim=c(0.5,1.2))
    text(cancer_TC[names(cancer_HBMR)],cancer_HBMR-0.02,labels = names(cancer_HBMR),col = c25[names(cancer_HBMR)])
    cor_temp<- cor.test(cancer_TC[names(cancer_HBMR)],cancer_HBMR)
    p_temp<- format(cor_temp$p.value,digits=2)
    r_temp<- format(cor_temp$estimate,digits=2)
    legend("topleft",legend = c(paste0("r=",r_temp),paste0("P=",p_temp)),bty="n")
    dev.off()
  }
}

