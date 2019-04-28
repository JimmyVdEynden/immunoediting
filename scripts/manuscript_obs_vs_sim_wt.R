###################################
# manuscript_obs_vs_sim_wt.R
###################################

# Load data
###########
{
  # Observation data
  load(paste0("results/data/fig",fig_nr,"_manuscript_obs_wt.RData"))
  
  # Simulation data
  load(paste0("results/data/fig",fig_nr,"_manuscript_sim_wt.RData"))
}

# Libraries
############
{
  library(svglite)
}

# correlation obs - sim
###################################
{
  # Compare n/s
  obs_p<- format(fisher.test(obs_isHA_variant_t)$p.value,digits = 3)
  sim_p<- format(fisher.test(sim_isHA_variant_t)$p.value,digits = 3)
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_obs_sim_ns.svg"),width = 5)
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_obs_sim_ns.pdf"))
    barplot(rbind(obs_ns,sim_ns),ylab="n/s (N/S)",names.arg = c("Non HLA binding","HLA binding"),xpd=F,beside=T,col=c("black","grey"),ylim = c(1.8,2.8))
    legend("topright",fill=c("black","grey"),legend=c(paste0("Observed (P=",obs_p,")"),paste0("Simulated (P=",sim_p,")")),bty="n")
    dev.off()
  }
  
  # Compare HBMR
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_obs_sim_HBMR_correlation.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_obs_sim_HBMR_correlation.pdf"))
    par(oma=c(6,2,2,2))
    # barplot
    bp<- barplot(rbind(obs_HBMR_isHA[,"HBMR"],sim_HBMR_isHA[rownames(obs_HBMR_isHA),"HBMR"]),las=2,ylab="HBMR",beside=T,ylim=c(0,1.2),col=c("black","grey"))
    # error bars
    arrows(bp,rbind(obs_HBMR_isHA[,"HBMR"],sim_HBMR_isHA[rownames(obs_HBMR_isHA),"HBMR"]),bp,rbind(obs_HBMR_isHA[,"CI"],sim_HBMR_isHA[rownames(obs_HBMR_isHA),"CI"]),lwd = 1.5, angle = 90,code = 2, length = 0.05,xpd=NA)    
    abline(h=1,lty=2)
    # Pearson correlation    
    cor_temp<- cor.test(obs_HBMR_isHA[,"HBMR"],sim_HBMR_isHA[rownames(obs_HBMR_isHA),"HBMR"])
    p_temp<- format(cor_temp$p.value,digits=2)
    r_temp<- format(cor_temp$estimate,digits=2)
    legend("topleft",fill=c("black","grey",NA,NA),border=NA,legend=c("Observed mutations","Simulated mutations",paste0("r=",r_temp),paste0("P=",p_temp)),bty="n")
    dev.off()
  }
}

