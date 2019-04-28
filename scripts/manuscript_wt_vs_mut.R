#####################################
# manuscript_wt_vs_mut.R
#####################################

# Compare wt vs mut affinities for each mutation

# Load
TCGA_maf<- readRDS("data/TCGA_maf.rds")
TCGA_maf_sim<- readRDS("data/TCGA_maf_sim.rds")

# 0) Preprocess & Filtering
#############################
{
  # Cancer with at least 10,000 mutations
  cancers<- names(table(TCGA_maf[,"Cancer"])[table(TCGA_maf[,"Cancer"])>10000])
  
  # Filter 1) Restrict analysis to non-synonymous mutations 
  TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"Variant_Classification"])&TCGA_maf[,"Variant_Classification"]=="nonsynonymous SNV",]
  TCGA_maf_sim<- TCGA_maf_sim[!is.na(TCGA_maf_sim[,"Variant_Classification"])&TCGA_maf_sim[,"Variant_Classification"]=="nonsynonymous SNV",]
  
  # Filter 2) Remove unknown or infinite HLA affinities
  TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"wt_HLA_mean_aff"])&!is.na(TCGA_maf[,"mut_HLA_mean_aff"])&!is.infinite(TCGA_maf[,"wt_HLA_mean_aff"])&!is.infinite(TCGA_maf[,"mut_HLA_mean_aff"]),]
  TCGA_maf_sim<- TCGA_maf_sim[!is.na(TCGA_maf_sim[,"wt_HLA_mean_aff"])&!is.na(TCGA_maf_sim[,"mut_HLA_mean_aff"])&!is.infinite(TCGA_maf_sim[,"mut_HLA_mean_aff"])&!is.infinite(TCGA_maf_sim[,"wt_HLA_mean_aff"]),]
  
  # Turn barcodes in factors
  TCGA_maf[,"Tumor_Sample_Barcode"]<- factor(TCGA_maf[,"Tumor_Sample_Barcode"])
  TCGA_maf_sim[,"Tumor_Sample_Barcode"]<- factor(TCGA_maf_sim[,"Tumor_Sample_Barcode"])
  
  # Define high affinity HLA binders (cut-off: Kd<500nM)
  # obs
  isHA_wt<- TCGA_maf[,"wt_HLA_mean_aff"]<500
  isHA_mut<- TCGA_maf[,"mut_HLA_mean_aff"]<500
  TCGA_maf<- cbind(TCGA_maf,isHA_wt=isHA_wt,isHA_mut=isHA_mut)
  # sim
  isHA_wt<- TCGA_maf_sim[,"wt_HLA_mean_aff"]<500
  isHA_mut<- TCGA_maf_sim[,"mut_HLA_mean_aff"]<500
  TCGA_maf_sim<- cbind(TCGA_maf_sim,isHA_wt=isHA_wt,isHA_mut=isHA_mut)
}

# Observations
###############

# Observations
for(i in 1:2){
  if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_cor_wt_mut_obs.svg"))
  else pdf(paste0("results/figs/pdf/fig",fig_nr,"_cor_wt_mut_obs.pdf"))
  idx_sel<- sample(1:nrow(TCGA_maf),10000)
  plot(log10(TCGA_maf[idx_sel,"wt_HLA_mean_aff"]),log10(TCGA_maf[idx_sel,"mut_HLA_mean_aff"]),xlab="HLA aff WT (log nM)",ylab="HLA aff MUT (log nM)",col=rgb(0.5,0.5,0.5,0.5),pch=16,cex=0.8,frame.plot=F,main="Observations")
  abline(h=log10(500),v=log10(500),lty=2,lwd=2)
  temp_cor<- cor.test(TCGA_maf[,"wt_HLA_mean_aff"],TCGA_maf[,"mut_HLA_mean_aff"],method="spearman") 
  legend("topleft",legend=c(paste0("P=",signif(temp_cor$p.value,2)),paste0("r=",signif(temp_cor$estimate,2))),bty="n")
  # Indicate aff changers
  idx_pt<- which(10*TCGA_maf[,"mut_HLA_mean_aff"]<TCGA_maf[,"wt_HLA_mean_aff"]&TCGA_maf[,"wt_HLA_mean_aff"]>=500&TCGA_maf[,"mut_HLA_mean_aff"]<500)
  idx_pt<- intersect(idx_pt,idx_sel)
  points(log10(TCGA_maf[idx_pt,"wt_HLA_mean_aff"]),log10(TCGA_maf[idx_pt,"mut_HLA_mean_aff"]),pch=16,col=rgb(1,0,0,1),cex=0.8)
  abline(c(-1,1))
  dev.off()
}

cat("HLA binding generating muts in observed data:",100*sum(10*TCGA_maf[,"mut_HLA_mean_aff"]<TCGA_maf[,"wt_HLA_mean_aff"]&TCGA_maf[,"wt_HLA_mean_aff"]>=500&TCGA_maf[,"mut_HLA_mean_aff"]<500)/nrow(TCGA_maf),"%","  \n")

# Simulations
for(i in 1:2){
  if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_cor_wt_mut_sim.svg"))
  else pdf(paste0("results/figs/pdf/fig",fig_nr,"_cor_wt_mut_sim.pdf"))
  idx_sel<- sample(1:nrow(TCGA_maf_sim),10000)
  plot(log10(TCGA_maf_sim[idx_sel,"wt_HLA_mean_aff"]),log10(TCGA_maf_sim[idx_sel,"mut_HLA_mean_aff"]),xlab="HLA aff WT (log nM)",ylab="HLA aff MUT (log nM)",col=rgb(0.5,0.5,0.5,0.5),pch=16,cex=0.8,frame.plot=F,main="Simulations")
  abline(h=log10(500),v=log10(500),lty=2,lwd=2)
  temp_cor<- cor.test(TCGA_maf_sim[,"wt_HLA_mean_aff"],TCGA_maf_sim[,"mut_HLA_mean_aff"],method="spearman") 
  legend("topleft",legend=c(paste0("P=",signif(temp_cor$p.value,2)),paste0("r=",signif(temp_cor$estimate,2))),bty="n")
  # Indicate aff changers
  idx_pt<- which(10*TCGA_maf_sim[,"mut_HLA_mean_aff"]<TCGA_maf_sim[,"wt_HLA_mean_aff"]&TCGA_maf_sim[,"wt_HLA_mean_aff"]>=500&TCGA_maf_sim[,"mut_HLA_mean_aff"]<500)
  idx_pt<- intersect(idx_pt,idx_sel)
  points(log10(TCGA_maf_sim[idx_pt,"wt_HLA_mean_aff"]),log10(TCGA_maf_sim[idx_pt,"mut_HLA_mean_aff"]),pch=16,col=rgb(1,0,0,1),cex=0.8)
  abline(c(-1,1))
  dev.off()
}

cat("HLA binding generating muts in simulated data:",100*sum(10*TCGA_maf_sim[,"mut_HLA_mean_aff"]<TCGA_maf_sim[,"wt_HLA_mean_aff"]&TCGA_maf_sim[,"wt_HLA_mean_aff"]>=500&TCGA_maf_sim[,"mut_HLA_mean_aff"]<500)/nrow(TCGA_maf_sim),"%","  \n")

# Compare mut/sim for different cut-offs of Kd and changes
cu_Kd<- c(5000,1000,500,100,50,10)
cu_affIncrease<- c(1,2,4,6,8,10,20,40,60,80,100)
obs_affChange_matrix<- matrix(NA,length(cu_Kd),length(cu_affIncrease),dimnames = list(cu_Kd,cu_affIncrease))
sim_affChange_matrix<- obs_affChange_matrix
for(cu_Kd_temp in cu_Kd){
  for(cu_affIncrease_temp in cu_affIncrease){
    affChange_temp<- cu_affIncrease_temp*TCGA_maf[,"mut_HLA_mean_aff"]<TCGA_maf[,"wt_HLA_mean_aff"]&TCGA_maf[,"wt_HLA_mean_aff"]>=cu_Kd_temp&TCGA_maf[,"mut_HLA_mean_aff"]<cu_Kd_temp
    obs_affChange_matrix[as.character(cu_Kd_temp),as.character(cu_affIncrease_temp)]<- 100*mean(affChange_temp,na.rm=T)
    affChange_sim_temp<- cu_affIncrease_temp*TCGA_maf_sim[,"mut_HLA_mean_aff"]<TCGA_maf_sim[,"wt_HLA_mean_aff"]&TCGA_maf_sim[,"wt_HLA_mean_aff"]>=cu_Kd_temp&TCGA_maf_sim[,"mut_HLA_mean_aff"]<cu_Kd_temp
    sim_affChange_matrix[as.character(cu_Kd_temp),as.character(cu_affIncrease_temp)]<- 100*mean(affChange_sim_temp,na.rm=T)
  }
}

for(i in 1:2){
  if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_probAffIncrease.svg"))
  else pdf(paste0("results/figs/pdf/fig",fig_nr,"_probAffIncrease.pdf"))
  plot(log10(cu_Kd),obs_affChange_matrix[,1],type="l",frame.plot=F,ylab="% HLA")
  points(log10(cu_Kd),sim_affChange_matrix[,1],type="l",lty=2)
  for(i in 2:length(cu_affIncrease)){
    points(log10(cu_Kd),obs_affChange_matrix[,i],type="l",col=i)
    points(log10(cu_Kd),sim_affChange_matrix[,i],type="l",col=i,lty=2)
  }
  abline(v=log10(500),lty=2)
  legend("topleft",legend = cu_affIncrease,lty = 1,col=1:length(cu_affIncrease),bty="n")
  dev.off()
}

