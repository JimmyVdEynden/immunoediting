###############################################
# manuscript_triNT_cancer_rankedHBMR.R
###############################################

# Show proportions of triNT subst types as a function of their HBMR rank for each cancer

# Load 
#######
{
  load("results/data/fig3_manuscript_cor_triNT_HBMR.RData")
  TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
}

# Plot
######
{
  library(svglite)
  # Color for plots
  substTypes<-c("C>A","C>G","C>T","T>A","T>C","T>G")
  substTypes_cols<-c("blue","black","red","grey","green","pink")
  names(substTypes_cols)<-substTypes
  
  # Data for plot
  HBMR<- sort(triNT_matrix[,"HBMR"])
  HBMR_names<- names(HBMR)
  HBMR_names_substTypes<- paste(substr(HBMR_names,2,2),substr(HBMR_names,6,6),sep=">")
  cancers<- c("BLCA","BRCA","CESC","COADREAD","ESCA","GBM","HNSC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SARC","SKCM","STAD","UCEC")
  
  # prop triNT subst type per cancer
  TCGA_mutProb<- prop.table(table(as.character(TCGA_maf[,"subst_type3"]),TCGA_maf[,"Cancer"]),2)
  
  # figure
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"_manuscript_triNT_cancer_rankedHBMR.pdf"))
    else svglite(paste0("results/figs/svg/fig",fig_nr,"_manuscript_triNT_cancer_rankedHBMR.svg"))
    par(mfrow=c(5,4),mar=c(2,4,2,2))
    for(cancer in cancers){
      bp<- barplot(TCGA_mutProb[,cancer][HBMR_names],las=2,col=substTypes_cols[HBMR_names_substTypes],names.arg = NA,main=cancer,border = NA,axes=F)
    }
    dev.off()
  }
}
