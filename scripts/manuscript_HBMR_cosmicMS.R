#################################################
# manuscript_HBMR_cosmicMS.R
#################################################

# Aim: Show effect of different cosmic signatures on HBMR

# Load & process data
######################
{
  # Expectd mut per subst type
  subst_type_matrix_ls<- readRDS(file="data/subst_type_matrix_ls.rds")
}

# Show effect cosmic signatures
################################
{
  cosmic_ms<- readRDS("data/cosmic_ms.rds")
  rownames(cosmic_ms)<- paste(cosmic_ms[,"Trinucleotide"],">",substr(cosmic_ms[,"Trinucleotide"],1,1),substr(cosmic_ms[,"Substitution.Type"],3,3),substr(cosmic_ms[,"Trinucleotide"],3,3),sep="")
  cosmic_ms<- t(cosmic_ms[,-c(1:3,34:40)])

  # Compare with HBMR
  HBMR_exp<- NULL
  for(i in 1:nrow(cosmic_ms)){
    exp_muts<- subst_type_matrix_ls$triNT*cosmic_ms[i,rownames(subst_type_matrix_ls$triNT)]
    isHA_t<- colSums(exp_muts[,,"TRUE"])
    isNonHA_t<- colSums(exp_muts[,,"FALSE"])
    HBMR_exp<- c(HBMR_exp,(isHA_t["nonsynonymous SNV"]/isHA_t["synonymous SNV"])/(isNonHA_t["nonsynonymous SNV"]/isNonHA_t["synonymous SNV"]))
  }
  names(HBMR_exp)<- rownames(cosmic_ms)
  
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"_ExpHBMR_cosmic.pdf"))
    else svglite(paste0("results/figs/svg/fig",fig_nr,"_ExpHBMR_cosmic.svg"))
    barplot(sort(HBMR_exp,decreasing = T),las=2,ylab="HBMR",horiz = T,xlab="HBMR")
    dev.off()
  }
}


