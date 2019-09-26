###################################
# manuscript_plot_Rsim_all.R
###################################

# Aim: plot Rsim for different cancer types using prototypical geotype (wt & mut affinities) and HLA-specific genotype (mut affinities)

# load
#######
{
  if(exists("usePentaData")&&usePentaData){
    if(useSimulated_data) load("results/data/summary_results_pentaNT_sim.RData")
    else load("results/data/summary_results_pentaNT.RData")  
  }
  else{
    if(useSimulated_data) load("results/data/summary_results_sim.RData")
    else load("results/data/summary_results.RData")  
  }
  source("scripts/functions/cancer_cols.R")
}

# Plot
#######
{
  # Cancers
  cancers<- sort(c('BRCA','COADREAD','LUAD','OV','HNSC','UCEC','LUSC','PRAD','LGG','STAD','SKCM','GBM','BLCA','LIHC','CESC','KIRP','SARC','PAAD','ESCA'))
  
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_Rsim_mut_all.svg"),width=10)
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_Rsim_mut_all.pdf"),width=10)
    par(mfrow=c(2,3))
    # Get cancers for each sample 
    cancer_id<- TCGA_cancer_id

    # Main - wt
    Rsim<- results_t_sample_mainWt[,"dNHLA_dNnonHLA"]
    names(Rsim)<- rownames(results_t_sample_mainWt)
    Rsim_med<- sort(tapply(Rsim,cancer_id[names(Rsim)],"median",na.rm=T)[cancers],decreasing = T)
    Rsim_p<- tapply(Rsim,cancer_id[names(Rsim)],function(x) wilcox.test(x,mu = 1)$p.value)[names(Rsim_med)]
    Rsim_q<- p.adjust(Rsim_p,"fdr")
    boxplot(Rsim~factor(cancer_id[names(Rsim)],levels = names(Rsim_med)),beside = T,las=2,outline=F,horizontal=T,frame.plot=F,col=c25[names(Rsim_med)],staplewex=0,notch=T,xlab="Observed/Expected",main="Prototypical HLA - wt")
    if(sum(Rsim_q<0.1)!=0) text(rep(2,sum(Rsim_q<0.1)),which(Rsim_q<0.1),paste0("P=",signif(Rsim_p[which(Rsim_q<0.1)],0.3)),xpd=NA)
    abline(v=1,lty=2)
    # Main - mut
    Rsim<- results_t_sample_mainMut[,"dNHLA_dNnonHLA"]
    names(Rsim)<- rownames(results_t_sample_mainMut)
    Rsim_med<- sort(tapply(Rsim,cancer_id[names(Rsim)],"median",na.rm=T)[cancers],decreasing = T)
    Rsim_p<- tapply(Rsim,cancer_id[names(Rsim)],function(x) wilcox.test(x,mu = 1)$p.value)[names(Rsim_med)]
    Rsim_q<- p.adjust(Rsim_p,"fdr")
    boxplot(Rsim~factor(cancer_id[names(Rsim)],levels = names(Rsim_med)),beside = T,las=2,outline=F,horizontal=T,frame.plot=F,col=c25[names(Rsim_med)],staplewex=0,notch=T,xlab="Observed/Expected",main="Prototypical HLA - mut")
    if(sum(Rsim_q<0.1)!=0) text(rep(2,sum(Rsim_q<0.1)),which(Rsim_q<0.1),paste0("P=",signif(Rsim_p[which(Rsim_q<0.1)],0.3)),xpd=NA)
    abline(v=1,lty=2)
    # Spec - mut
    Rsim<- results_t_sample[,"dNHLA_dNnonHLA"]
    names(Rsim)<- rownames(results_t_sample)
    names(cancer_id)<- substr(names(cancer_id),1,12) 
    Rsim_med<- sort(tapply(Rsim,cancer_id[names(Rsim)],"median",na.rm=T)[cancers],decreasing = T)
    Rsim_p<- tapply(Rsim,cancer_id[names(Rsim)],function(x) wilcox.test(x,mu = 1)$p.value)[names(Rsim_med)]
    Rsim_q<- p.adjust(Rsim_p,"fdr")
    boxplot(Rsim~factor(cancer_id[names(Rsim)],levels = names(Rsim_med)),beside = T,las=2,outline=F,horizontal=T,frame.plot=F,col=c25[names(Rsim_med)],staplewex=0,notch=T,xlab="Observed/Expected",main="Specific HLA - mut")
    if(sum(Rsim_q<0.1)!=0) text(rep(2,sum(Rsim_q<0.1)),which(Rsim_q<0.1),paste0("P=",signif(Rsim_p[which(Rsim_q<0.1)],0.3)),xpd=NA)
    abline(v=1,lty=2)
    dev.off()
  }
  
}