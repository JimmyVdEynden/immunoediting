########################################
# manuscript_cor_HBMR_aaProbs.R
########################################

# Load data
###########
{
  # triNT HMBR x aaProb matrix
  GPPM_triNT_aa_var_t <- readRDS("data/GPPM_triNT_aa_var_t.rds")
  # aa classes
  load(file = "data/aa_classes.RData")
  # HBBMR values
  if(byCancer){
    TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
    load("results/data/fig2_manuscript_obs_wt.RData") # Per cancer
    HBMR<- obs_HBMR_isHA[,1]
  }  
  else{
    load("results/data/fig3_manuscript_cor_triNT_HBMR.RData") # per triNT
    HBMR<- sort(triNT_matrix[,"HBMR"])
  }
}

# Libraries and functions
##########################
{
  library(plotrix)
  source("scripts/functions/plot_aa_matrix.R")
  library(svglite)
  library(corrplot)
  library(RColorBrewer)
}

# Analysis of N/S effect
############################################
{
  # Get aa matrices for triNT
  if(byCancer){
    TCGA_mutProb<- prop.table(table(TCGA_maf[,"Cancer"],as.character(TCGA_maf[,"subst_type3"])),1)
    # Normalize by sites
    n_triNT<- rowSums(GPPM_triNT_aa_var_t)
    TCGA_mutProb<- prop.table(t(t(TCGA_mutProb[,names(n_triNT)])/n_triNT),1)
    cancers<- sort(c('BRCA','COADREAD','LUAD','OV','HNSC','UCEC','LUSC','PRAD','LGG','STAD','SKCM','GBM','BLCA','LIHC','CESC','KIRP','SARC','PAAD','ESCA'))
    aa_S_matrix<- NULL
    aa_N_matrix<- NULL
    for(cancer in cancers){
      aa_matrix_temp<- colSums(GPPM_triNT_aa_var_t*TCGA_mutProb[cancer,])
      aa_S_matrix<- rbind(aa_S_matrix,aa_matrix_temp[,"synonymous SNV"])
      aa_N_matrix<- rbind(aa_N_matrix,aa_matrix_temp[,"nonsynonymous SNV"])
    }
    rownames(aa_S_matrix)<- cancers
    rownames(aa_N_matrix)<- cancers
  }
  else{
    aa_S_matrix<- GPPM_triNT_aa_var_t[,,"synonymous SNV"]
    aa_N_matrix<- GPPM_triNT_aa_var_t[,,"nonsynonymous SNV"]
  }
  aa_S_matrix<- prop.table(aa_S_matrix,1) 
  aa_N_matrix<- prop.table(aa_N_matrix,1) 
  aa_S_matrix<- t(aa_S_matrix)
  aa_N_matrix<- t(aa_N_matrix)
  
  # Synonymous versus non-synonymous mutations, grouped by aa class
  aa_S_matrix_grouped<- rbind(hp=colSums(aa_S_matrix[aa_hp,]),polar=colSums(aa_S_matrix[aa_polar,]),charged=colSums(aa_S_matrix[aa_charged,]))
  aa_N_matrix_grouped<- rbind(hp=colSums(aa_N_matrix[aa_hp,]),polar=colSums(aa_N_matrix[aa_polar,]),charged=colSums(aa_N_matrix[aa_charged,]))

  # Plot   
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"_aaProb_grouped.pdf"))
    else svglite(paste0("results/figs/svg/fig",fig_nr,"_aaProb_grouped.svg"))
    par(mfrow=c(2,1))
    plotted_matrix<- plot_aa_matrix(aa_S_matrix_grouped[,!is.na(aa_S_matrix_grouped[1,])],c("hp","polar","charged"),names(HBMR))
    plotted_matrix<- plot_aa_matrix(aa_N_matrix_grouped,c("hp","polar","charged"),names(HBMR))
    dev.off()  
  }
  
  # Get spearman correlations
  cat("Spearman HBMR vs ref aa class prob, synonymous:"," \n")
  for(aa in rownames(aa_S_matrix_grouped)){
    cor_temp<- cor.test(HBMR,aa_S_matrix_grouped[aa,names(HBMR)],method="spearman")
    cat(aa,"r=",cor_temp$estimate,"(p=",cor_temp$p.value,")","\n")
  }
  
  cat("\n","Spearman HBMR vs ref aa class prob, non-synonymous:"," \n")
  for(aa in rownames(aa_N_matrix_grouped)){
    cor_temp<- cor.test(HBMR,aa_N_matrix_grouped[aa,names(HBMR)],method="spearman")
    cat(aa,"r=",cor_temp$estimate,"(p=",cor_temp$p.value,")","\n")
  }
  
  # Synonymous versus non-synonymous mutations for individual aa
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"aaProb.pdf"))
    else svglite(paste0("results/figs/svg/fig",fig_nr,"aaProb.svg"))
    par(mfrow=c(2,1))
    plotted_matrix<- plot_aa_matrix(aa_S_matrix[,!is.na(aa_S_matrix[1,])],c(aa_hp,aa_polar,aa_charged),names(HBMR),add_loess = F,plot_title = "Synonymous")
    points(seq(0,10),rep(nrow(plotted_matrix)+1,11),col=color.scale(seq(0:10),cs1=c(1,0),cs2=c(1,0),cs3=c(1,0)),pch=15,cex=2,xpd=NA)
    plot_aa_matrix(aa_N_matrix[,!is.na(aa_N_matrix[1,])],c(aa_hp,aa_polar,aa_charged),names(HBMR),add_loess = F,plot_title = "Non-synonymous")
    dev.off()
  }
  
  # Get spearman correlations for ind amino acids
  cat("\n","Spearman HBMR vs ref aa prob, synonymous:"," \n")
  for(aa in c(aa_hp,aa_polar,aa_charged)){
    cor_temp<- cor.test(HBMR,aa_S_matrix[aa,names(HBMR)],method="spearman")
    cat(aa,"r=",cor_temp$estimate,"(p=",cor_temp$p.value,")","\n")
  }
  
  cat("\n","Spearman HBMR vs ref aa prob, non-synonymous:"," \n")
  for(aa in c(aa_hp,aa_polar,aa_charged)){
    cor_temp<- cor.test(HBMR,aa_N_matrix[aa,names(HBMR)],method="spearman")
    cat(aa,"r=",cor_temp$estimate,"(p=",cor_temp$p.value,")","\n")
  }
  
  # Plot correlations ind aaa
  cor_matrix<- cbind(
    t(cor(HBMR,t(aa_S_matrix[c(aa_hp,aa_polar,aa_charged),names(HBMR)]),method = "spearman")),
    t(cor(HBMR,t(aa_N_matrix[c(aa_hp,aa_polar,aa_charged),names(HBMR)]),method = "spearman"))
  )
  colnames(cor_matrix)<- c("S","N")
  cor_matrix_p<- matrix(NA,length(aa_all),2,dimnames = list(c(aa_hp,aa_polar,aa_charged),c("S","N")))
  for(aa in aa_all){
    cor_matrix_p[aa,"S"]<- cor.test(HBMR,aa_S_matrix[aa,names(HBMR)],method="spearman")$p.value
    cor_matrix_p[aa,"N"]<- cor.test(HBMR,aa_N_matrix[aa,names(HBMR)],method="spearman")$p.value
  }
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"aaProb_cor.pdf"))
    else svglite(paste0("results/figs/svg/fig",fig_nr,"aaProb_cor.svg"))
    par(mfrow=c(2,1))
    corrplot(cbind(cor_matrix,0,scale=seq(-1,1,by = 0.25)),p.mat = cbind(cor_matrix_p,1,0),sig.level=0.05,insig = "blank", col=brewer.pal(n = 9, name = "PuOr"),cl.pos = "b")
    dev.off()
  }
  
  # Plot correlations grouped aaa
  cor_matrix<- cbind(
    t(cor(HBMR,t(aa_S_matrix_grouped[,names(HBMR)]),method = "spearman")),
    t(cor(HBMR,t(aa_N_matrix_grouped[,names(HBMR)]),method = "spearman"))
  )
  colnames(cor_matrix)<- c("S","N")
  cor_matrix_p<- matrix(NA,3,2,dimnames = list(rownames(aa_S_matrix_grouped),c("S","N")))
  for(aa in rownames(aa_S_matrix_grouped)){
    cor_matrix_p[aa,"S"]<- cor.test(HBMR,aa_S_matrix_grouped[aa,names(HBMR)],method="spearman")$p.value
    cor_matrix_p[aa,"N"]<- cor.test(HBMR,aa_N_matrix_grouped[aa,names(HBMR)],method="spearman")$p.value
  }
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"aaGroupProb_cor.pdf"))
    else svglite(paste0("results/figs/svg/fig",fig_nr,"aaGroupProb_cor.svg"))
    par(mfrow=c(2,1))
    corrplot(cbind(cor_matrix,0,scale=seq(-1,1)),p.mat = cbind(cor_matrix_p,1,0),sig.level=0.05,insig = "blank", col=brewer.pal(n = 9, name = "PuOr"),cl.pos = "b")
    dev.off()
  }
}


