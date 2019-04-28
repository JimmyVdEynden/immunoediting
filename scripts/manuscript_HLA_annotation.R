###################################
# manuscript_HLA_annotation.R
###################################

# Load GPPM
############
{
  GPPM<- readRDS(file = "data/GPPM_inclHLAAlleles.rds")
}

# Filter GPPM
##############
{
  # No need for mutation info for this analysis!
  idx_mut<- which(!duplicated(paste(GPPM$gene,pos(GPPM),sep="_")))
  GPPM<- GPPM[idx_mut] 
}

# Libraries
###########
{
  library(corrplot)
  library(svglite)
}

# Get proportions of genome that are annotated as HLA binding
#############################################################
{
  for(j in 1:2){
    if(j==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_HLA_allele_pie.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_HLA_allele_pie.pdf"),width = 10)
    par(mfrow=c(2,4))
    for(i in 1:7){
      # cat(i," ")
      idx_allele<- grep("HLA_aff",names(mcols(GPPM)))[i] # A1, A2, B1, B2, C1, C2, mean
      allele_name<- c("HLA_A02_01","HLA_A01_01","HLA_B07_02","HLA_B08_01","HLA_C07_01","HLA_C07_02","HLA_Mean")[i] # Cfr read.table("temp/HLA_types_freq.txt",stringsAsFactors = F)
      allele_HLA_aff<- mcols(GPPM)[,idx_allele]
      allele_HLA_aff_t<- table(allele_HLA_aff<500)
      pie(allele_HLA_aff_t,labels = c(NA,"HLA-binding"),main=allele_name)
      text(0,0,paste0(100*signif(prop.table(allele_HLA_aff_t)["TRUE"],3),"%"),adj = c(0,0))
    }
    dev.off()
  }
}

# Get correlations between Kds
###############################
{
  allele_HLA_aff_t<- mcols(GPPM)[,grep("HLA_aff",names(mcols(GPPM)))]
  colnames(allele_HLA_aff_t)<- c("HLA_A02_01","HLA_A01_01","HLA_B07_02","HLA_B08_01","HLA_C07_01","HLA_C07_02","HLA_Mean")
  idx_sample<- sample(1:nrow(allele_HLA_aff_t),size = 1000000,replace = F) # Cannot do on all because of memory issues, sample 1milion
  allele_HLA_aff_cor<- cor(as.matrix(allele_HLA_aff_t[idx_sample,]),method = "spearman")

  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_HLA_allele_cor_plot.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_HLA_allele_cor_plot.pdf"),width = 10)
    corrplot(allele_HLA_aff_cor,type="lower",tl.srt=45,tl.col="black",diag=FALSE)
    dev.off()
  }
}


