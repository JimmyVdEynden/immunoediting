#######################
# manuscript_cor_CYT.R
#######################

# Load data
############
{
  load("downloads/TCGA/TCGA_mRNA.RData")
  # load("data/TCGA_manifest.RData")
  TCGA_maf<- readRDS("data/TCGA_maf.rds")
}

# Libraries & functions
#######################
{
  source("scripts/functions/cancer_cols.R")
  # Geometric mean on columns
  geom_colMean<- function(mRNA) apply(X = mRNA, MARGIN = 2, function(x) exp(mean(log(x))))
}

# Create expr TIL matrix
#########################
{
  # Exclude normals from mRNA data
  TCGA_mRNA_T<- TCGA_mRNA[,!grepl("*.-11",colnames(TCGA_mRNA))]

  # Define expr markers
  expr_markers<- c(
    CYT=list(c("GZMA","PRF1")),
    CD8=list(c("CD8A")),
    CD3=list(c("CD3G"))
    )
  
  # Create expr matrix
  for(i in 1:length(expr_markers)){
    genes<- expr_markers[[i]]
    if(i==1) TIL_matrix<-  geom_colMean(TCGA_mRNA_T[genes,,drop=F])
    else TIL_matrix<- cbind(TIL_matrix, geom_colMean(TCGA_mRNA_T[genes,,drop=F]))
  }
  colnames(TIL_matrix)<- names(expr_markers)
  rownames(TIL_matrix)<- substr(rownames(TIL_matrix),1,12)
}

# Load sample-specific dNHLA/dNnonHLA data
##########################################
{
  load("results/data/summary_results.RData")
  Rsim<- results_t_sample$dNHLA_dNnonHLA
  names(Rsim)<- rownames(results_t_sample)
  Rsim<- Rsim[!is.infinite(Rsim)]
  names(TCGA_cancer_id)<- substr(names(TCGA_cancer_id),1,12)
}

# Fuse with TIL matrix
######################
{
  common_samples<- intersect(intersect(names(Rsim),rownames(TIL_matrix)),names(TCGA_cancer_id))
  sample_matrix<- cbind(TIL_matrix[common_samples,],Rsim=Rsim[common_samples],cancer_id=TCGA_cancer_id[common_samples])
}

# Plot CYT vs Rsim
###################
{
  # Get Spearman corr for all cancers & plot
  CYT_Rsim_cor<- matrix(NA,length(cancers)+1,2,dimnames = list(c("Pan",cancers),c("r","p")))
  for(cancer in c("Pan",cancers)){
    if(cancer=="Pan") sample_matrix_can<- sample_matrix
    else sample_matrix_can<- sample_matrix[sample_matrix[,"cancer_id"]==cancer,]
    p_cor<- cor.test(as.numeric(sample_matrix_can[,"CYT"]),as.numeric(sample_matrix_can[,"Rsim"]),method="spearman")$p.value
    r_cor<- cor.test(as.numeric(sample_matrix_can[,"CYT"]),as.numeric(sample_matrix_can[,"Rsim"]),method="spearman")$estimate
    CYT_Rsim_cor[cancer,]<- c(r_cor,p_cor) 
  }
  # CYT_Rsim_cor<- CYT_Rsim_cor[order(CYT_Rsim_cor[,"r"],decreasing = F),]
  # print.data.frame(as.data.frame(format(CYT_Rsim_cor,digits=3)))

  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_cor_CYT_coeff.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_cor_CYT_coeff.pdf"))
    plot(rep(1,nrow(CYT_Rsim_cor)),CYT_Rsim_cor[,"r"],pch="-",ylim=c(-0.15,0.15),axes=F,xlab=NA,ylab="Spearman correlation coeff.",xlim=c(0.8,1.5))
    axis(2)
    text(rep(1.2,nrow(CYT_Rsim_cor)),CYT_Rsim_cor[,"r"],labels = paste0(rownames(CYT_Rsim_cor)," (P=",format(CYT_Rsim_cor[,"p"],digits=2,scientific = T),")"),adj = 0)
    dev.off()
  }
  
  # Plot pan cancer analysis
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_cor_CYT.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_cor_CYT.pdf"))
    plot(as.numeric(sample_matrix[,"CYT"]),as.numeric(sample_matrix[,"Rsim"]),pch=16,col="grey",ylim=c(0,3),xlim=c(1,1000),lwd=2,xlab="Cytolytic activity (CYT)",ylab="dNHLA/dNnonHLA")
    points(loess.smooth(as.numeric(sample_matrix[,"CYT"]),as.numeric(sample_matrix[,"Rsim"])),type="l",col="black",lwd=2)
    abline(h=1,lty=2)
    cor_temp<-  cor.test(as.numeric(sample_matrix[,"CYT"]),as.numeric(sample_matrix[,"Rsim"]),method="spearman")
    p_temp<- format(cor_temp$p.value,digits=3,scientific = T)
    r_temp<- format(cor_temp$estimate,digits=2)
    legend("topright",legend=c(paste0("P=",p_temp),paste0("r=",r_temp)),bty="n")
    dev.off()
  }
}
  

