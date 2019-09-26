###################################
# manuscript_sim_wt_validation.R
###################################

# Load data
#############
{
  TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
  TCGA_maf_sim<- readRDS("data/TCGA_maf_sim_mainHLAType.rds")
  NS_t <- readRDS("data/gene_NS_sites.rds")
  
  # CGC genes
  CGC_genes<- readRDS("data/CGC_v83.rds")
}

# Libraries & functions
##########################
{
  library(svglite)
  source("scripts/functions/calculate_dNdS.R")
}

# Calculate
###########
dNdS_obs<- calculate_dNdS(TCGA_maf$Hugo_Symbol,TCGA_maf$Variant_Classification,TCGA_maf$subst_type3,exp_mut = NS_t$triNT)
dNdS_sim<- calculate_dNdS(TCGA_maf_sim$Hugo_Symbol,TCGA_maf_sim$Variant_Classification,TCGA_maf_sim$subst_type3,exp_mut = NS_t$triNT)

# Get top 50 ranked CGC genes
dNdS_obs_top<- sort(dNdS_obs[!is.infinite(dNdS_obs[,"dNdS"])&rownames(dNdS_obs)%in%CGC_genes,"dNdS"],decreasing=T)[1:50]

# Plot
for(i in 1:2){
  if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_validate_sims.svg"),width = 20)
  else pdf(paste0("results/figs/pdf/fig",fig_nr,"_validate_sims.pdf"),width = 14)
  bp<- barplot(rbind(dNdS_obs_top,dNdS_sim[names(dNdS_obs_top),"dNdS"]),las=2,beside = T,ylab="dN/dS",col=c("black","grey"))
  # arrows(bp,rbind(dNdS_obs_top,dNdS_sim[names(dNdS_obs_top),"dNdS"]),bp,rbind(dNdS_obs[names(dNdS_obs_top),"ciup"],dNdS_sim[names(dNdS_obs_top),"ciup"]),lwd = 1.5, angle = 90,code = 2, length = 0.03,xpd=NA)    
  abline(h=1,lty=2)
  legend("topright",fill=c("black","grey"),legend=c("Observed","Simulated"),bty="n")
  dev.off()
}

# Compare variants
##################
{
  TCGA_maf_obs_var_t<- table(TCGA_maf$Variant_Classification)
  TCGA_maf_sim_var_t<- table(TCGA_maf_sim$Variant_Classification)
  TCGA_maf_obs_var_prop_t<- prop.table(TCGA_maf_obs_var_t)
  TCGA_maf_sim_var_prop_t<- prop.table(TCGA_maf_sim_var_t)
  
  # Plot
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_validate_sims_variants.svg"),width = 20)
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_validate_sims_variants.pdf"),width = 14)
    par(mfrow=c(1,2))
    pie(TCGA_maf_obs_var_prop_t,main="Observed mutations",col=c("blue","red","grey"),labels = paste0(names(TCGA_maf_obs_var_prop_t)," (",round(100*TCGA_maf_obs_var_prop_t,1),"%)"),xpd=NA)
    pie(TCGA_maf_sim_var_prop_t,main="Simulated mutations",col=c("blue","red","grey"),labels = paste0(names(TCGA_maf_sim_var_prop_t)," (",round(100*TCGA_maf_sim_var_prop_t,1),"%)"),xpd=NA)
    dev.off()
  }
}
