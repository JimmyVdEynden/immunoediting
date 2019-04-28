####################################################
# manuscript_HBMR_powerAnalysis.R
####################################################

# Aim: HBMR power analysis: determine statistical power for n mutations available for each cancer type

# Load library to perform power analysis on Fisher test
# install.packages("exact2x2")
library(exact2x2)
library(svglite)

# Power analysis  
# p0 = HLA binding prop. without neg sel (control group, =synonymous mut!)
# n0 = total number of syn mut (control group)
# p1 = HLA binding prop. (simulated) negative selection (nonsyn mut!)
# n1 = total number of nonsyn mut (treatment group)

# Load cancer colors for figures
source("scripts/functions/cancer_cols.R")

# Get number of non-synonymous and synonymous mut for each cancer
TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"Variant_Classification"])&TCGA_maf[,"Variant_Classification"]%in%c("synonymous SNV","nonsynonymous SNV"),]
TCGA_maf_can_t<- table(TCGA_maf$Cancer,TCGA_maf$Variant_Classification)

# 1. Calculate power for 1) different reductions in HBMR (i.e. selection) and 2) different numbers of mutations
################################################################################################################
{
  negSelSim<- c(0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,0.99,1) # Simulated amount of remaining mutations after neg sel
  nMutSim<- c(100,300,1000,3000,10000,30000,100000,300000,1000000) # Simulated n mutations
  propHLA<- 0.221 # From HLA annotation (cfr fig.1)
  
  # Simulate n syn and non-syn mutations for these numbers
  prop_n<- prop.table(colSums(TCGA_maf_can_t))["nonsynonymous SNV"] # 71%, i.e Overall pan-cancer proportion non-syn mutations
  nMutSim_t<- round(prop_n*nMutSim)
  nMutSim_t<- cbind(nMutSim_t,nMutSim-nMutSim_t)
  rownames(nMutSim_t)<- nMutSim
  colnames(nMutSim_t)<- c("nonsynonymous SNV","synonymous SNV")
  
  # Calculate power given # mut (rows) and # negsel
  pwr_n_matrix<- matrix(NA,nrow(nMutSim_t),length(negSelSim),dimnames = list(rownames(nMutSim_t),negSelSim))
  for(all_total in rownames(pwr_n_matrix)){
    
    # Get n and s per cancer
    n_total<- nMutSim_t[all_total,"nonsynonymous SNV"]
    s_total<- nMutSim_t[all_total,"synonymous SNV"]
    
    # Power
    for(trueNegSel in colnames(pwr_n_matrix)){
      # Get number o non-syn after selection on HLA-binding non-syn (i.e. fill in new contingency)
      n_nonHLA<- (1-propHLA)*n_total 
      n_HLA<- n_total-n_nonHLA # Get tot non-syn after selection
      n_HLA_sim<- n_HLA*as.numeric(trueNegSel)
      n_total_sim<- n_nonHLA + n_HLA_sim
      propHLA_sim<- n_HLA_sim/(n_HLA_sim+n_nonHLA)
      pwr_tmp<- power2x2(p0 = propHLA,p1 = propHLA_sim,n0 = s_total,n1 = n_total_sim,approx = T)
      pwr_n_matrix[all_total,trueNegSel]<- pwr_tmp$power
    }
  }
}

# 1.2 Add lower positive predictive values
##########################################
{
  negSelSim<- c(0,0.2,0.4,0.6,0.8,0.85,0.9,0.95,0.99,1) # Simulated amount of remaining mutations after neg sel
  nMutSim<- c(10000,30000,100000,300000) # Simulated n mutations
  ppv<- c(1,0.8,0.6,0.4,0.2,0.1)
  propHLA<- 0.221 # From HLA annotation (cfr fig.1)
  
  
  # Simulate n syn and non-syn mutations for these numbers
  prop_n<- prop.table(colSums(TCGA_maf_can_t))["nonsynonymous SNV"] # 71%, i.e Overall pan-cancer proportion non-syn mutations
  nMutSim_t<- round(prop_n*nMutSim)
  nMutSim_t<- cbind(nMutSim_t,nMutSim-nMutSim_t)
  rownames(nMutSim_t)<- nMutSim
  colnames(nMutSim_t)<- c("nonsynonymous SNV","synonymous SNV")
  
  # Calculate power given # mut (rows) and # negsel
  pwr_ppv_matrix_ls<- list()
  pwr_ppv_HBMR_matrix_ls<- list()
  for(all_total in rownames(nMutSim_t)){
    # Get n and s per cancer
    n_total<- nMutSim_t[all_total,"nonsynonymous SNV"]
    s_total<- nMutSim_t[all_total,"synonymous SNV"]
    
    pwr_ppv_matrix<- matrix(NA,length(ppv),length(negSelSim),dimnames = list(ppv,negSelSim))
    pwr_ppv_HBMR_matrix<- matrix(NA,length(ppv),length(negSelSim),dimnames = list(ppv,negSelSim))
    for(ppv_temp in rownames(pwr_ppv_matrix)){
      
      
      # Power
      for(trueNegSel in colnames(pwr_ppv_matrix)){
        # Get number o non-syn after selection on HLA-binding non-syn (i.e. fill in new contingency)
        n_nonHLA<- (1-propHLA)*n_total 
        n_HLA<- n_total-n_nonHLA # Get tot non-syn after selection
        n_HLA_sim<- n_HLA-(n_HLA*(1-as.numeric(trueNegSel))*as.numeric(ppv_temp))
        # n_HLA_sim<- n_HLA*as.numeric(trueNegSel)
        n_total_sim<- n_nonHLA + n_HLA_sim
        propHLA_sim<- n_HLA_sim/(n_HLA_sim+n_nonHLA)
        pwr_tmp<- power2x2(p0 = propHLA,p1 = propHLA_sim,n0 = s_total,n1 = n_total_sim,approx = T)
        pwr_ppv_matrix[ppv_temp,trueNegSel]<- pwr_tmp$power
        # HBMR?
        s_nonHLA<- s_total*(1-propHLA)
        s_HLA<- s_total*(propHLA)
        pwr_ppv_HBMR_matrix[ppv_temp,trueNegSel]<- fisher.test(round(rbind(c(s_nonHLA,s_HLA),c(n_nonHLA,n_HLA_sim))))$estimate
      }
    }
    pwr_ppv_matrix_ls<- c(pwr_ppv_matrix_ls,list(pwr_ppv_matrix))
    pwr_ppv_HBMR_matrix_ls<- c(pwr_ppv_HBMR_matrix_ls,list(pwr_ppv_HBMR_matrix))
  }
  names(pwr_ppv_matrix_ls)<- rownames(nMutSim_t)
  names(pwr_ppv_HBMR_matrix_ls)<- rownames(nMutSim_t)
}

# 2. Calculate power for 1) different reductions in HBMR (i.e. selection) and 2) given numbers of mutation for each cancer
#########################################################################################################################
{
  pwr_can_matrix<- matrix(NA,nrow(TCGA_maf_can_t),length(negSelSim),dimnames = list(rownames(TCGA_maf_can_t),negSelSim))
  for(cancer in rownames(pwr_can_matrix)){
    
    # Get n and s per cancer
    n_total<- TCGA_maf_can_t[cancer,"nonsynonymous SNV"]
    s_total<- TCGA_maf_can_t[cancer,"synonymous SNV"]
    
    # Power
    for(trueNegSel in colnames(pwr_can_matrix)){
      # Get number o non-syn after selection on HLA-binding non-syn (i.e. fill in new contingency)
      n_nonHLA<- (1-propHLA)*n_total 
      n_HLA<- n_total-n_nonHLA # Get tot non-syn after selection
      n_HLA_sim<- n_HLA*as.numeric(trueNegSel)
      n_total_sim<- n_nonHLA + n_HLA_sim
      propHLA_sim<- n_HLA_sim/(n_HLA_sim+n_nonHLA)
      pwr_tmp<- power2x2(p0 = propHLA,p1 = propHLA_sim,n0 = s_total,n1 = n_total_sim,approx = T)
      pwr_can_matrix[cancer,trueNegSel]<- pwr_tmp$power
    }
  }
}

# Plot
######
{
  # Numbers
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_HBMR_Power.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_HBMR_Power.pdf"))
    par(mfrow=c(2,1))
    for(i in 1:nrow(pwr_n_matrix)){
      if(i==1) plot(negSelSim,pwr_n_matrix[i,],type="l",ylab="Power",xlab="Expected HBMR value (i.e. extent of negative selection)",col=c25[i],frame.plot=F,ylim=c(0,1.1))
      else points(negSelSim,pwr_n_matrix[i,],type="l",col=c25[i])
    }
    legend("bottomleft",legend = nMutSim,col = c25[1:length(nMutSim)],lty=1,cex=0.5)
    abline(v=0.8,h=0.8,lty=2)
    # 10000 muts means 0.98 power to detect reduction of 20%
    
    # Cancer numbers
    c32<- c(c25,rep("black",7))
    names(c32)<- c(names(c25),setdiff(rownames(pwr_can_matrix),names(c25)))
    cancer_sel<- rownames(TCGA_maf_can_t)[rowSums(TCGA_maf_can_t)>10000]
    c32[!names(c32)%in%cancer_sel]<- "black"
    lty32<- rep(1,32)
    names(lty32)<- names(c32)
    lty32[!names(lty32)%in%cancer_sel]<- 2
    
    for(i in 1:nrow(pwr_can_matrix)){
      if(i==1) plot(negSelSim,pwr_can_matrix[i,],type="l",ylab="Power",xlab="Expected HBMR value (i.e. extent of negative selection)",col=c32[rownames(pwr_can_matrix)[i]],frame.plot=F,ylim=c(0,1.1),lty=lty32[rownames(pwr_can_matrix)[i]])
      else points(negSelSim,pwr_can_matrix[i,],type="l",col=c32[rownames(pwr_can_matrix)[i]],lty=lty32[rownames(pwr_can_matrix)[i]])
    }
    legend("bottomleft",legend = paste0(cancer_sel, " (n=",rowSums(TCGA_maf_can_t[cancer_sel,]),")"),col = c32[cancer_sel],lty=1,cex=0.5)
    dev.off()
  }
}

# Plot ppv
##########
{
  # Numbers
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_HBMR_Power_ppv.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_HBMR_Power_ppv.pdf"))
    par(mfrow=c(2,2),oma=c(2,2,2,2))
    for(p in 1:length(pwr_ppv_matrix_ls)){
      # Power
      pwr_ppv_matrix<- pwr_ppv_matrix_ls[[p]]
      for(i in 1:nrow(pwr_ppv_matrix)){
        if(i==1){
          plot(colnames(pwr_ppv_matrix),pwr_ppv_matrix[i,],type="l",ylab="Power",xlab="Expected HBMR",col=i,frame.plot=F,axes=F,ylim=c(0,1.1),main=names(pwr_ppv_matrix_ls)[p])
          axis(2)
          axis(3)
        }
        else points(colnames(pwr_ppv_matrix),pwr_ppv_matrix[i,],type="l",col=i)
      }
      legend("bottomleft",legend = rownames(pwr_ppv_matrix),col = 1:nrow(pwr_ppv_matrix),lty=1,cex=0.5)
      abline(v=0.8,h=0.8,lty=2)
      # HBMR: independent of n, take once
      pwr_ppv_HBMR_matrix<- round(pwr_ppv_HBMR_matrix_ls[[p]],2)
      # plot.new()
      for(i in 1:nrow(pwr_ppv_HBMR_matrix)){
        # lines(x = c(0,1),y= c(0.2*i,0.2*i),col=i)
        text(x = as.numeric(colnames(pwr_ppv_HBMR_matrix)),y=((-0.1*i)),xpd=NA,labels = pwr_ppv_HBMR_matrix[i,],cex=0.5,col=i)
        text(x = -0.2,y=((-0.1*i)),xpd=NA,labels = rownames(pwr_ppv_HBMR_matrix)[i])
      }
    }
    dev.off()
  }
}

# Save results
##############
{
  save(pwr_can_matrix,pwr_n_matrix,pwr_ppv_matrix_ls,pwr_ppv_HBMR_matrix_ls,file = paste0("results/data/fig", fig_nr,"HBMR_power.RData"))
}
