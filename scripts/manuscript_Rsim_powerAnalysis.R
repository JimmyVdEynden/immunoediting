####################################################
# manuscript_Rsim_powerAnalysis.R
####################################################

# Aim: Rsim power analysis: determine statistical power for n mutations available for each cancer type

# Get expected ratio's
load("results/data/Rsim_HLASpecific_simData.RData")

# Get cancers to analyze
cancers<- sort(c('BRCA','COADREAD','LUAD','OV','HNSC','UCEC','LUSC','PRAD','LGG','STAD','SKCM','GBM','BLCA','LIHC','CESC','KIRP','SARC','PAAD','ESCA'))

# Load cancer colors for figures
source("scripts/functions/cancer_cols.R")

# Cancer IDs
load("data/TCGA_manifest.RData")
TCGA_cancer_id<- cbind(substr(rownames(TCGA_cancer_id),1,12),TCGA_cancer_id)
TCGA_cancer_id<- TCGA_cancer_id[!duplicated(TCGA_cancer_id[,1]),]
rownames(TCGA_cancer_id)<- TCGA_cancer_id[,1] 

#libraries
library(svglite)

# 1. Calculate power for 1) different reductions in Rsim (i.e. selection) and 2) different numbers of samples
################################################################################################################
{
  # Focus on ratio's from 19 analysed cancers
  Rsim_sel<- Rsim[TCGA_cancer_id[names(Rsim),2]%in%cancers]
  Rsim_sel<- Rsim_sel[!is.na(Rsim_sel)]
  
  # Distribution after logtransformation +/- normal
  Rsim_sel<- Rsim_sel[Rsim_sel!=0]
  Rsim_sel<- Rsim_sel[!is.infinite(Rsim_sel)]
  # hist(log10(Rsim_sel),breaks = 100,freq=F)
  sd_log<- sd(log10(Rsim_sel))
  # curve(dnorm(x, mean=0, sd = sd_log), col="darkblue", lwd=2, add=TRUE, yaxt="n")
  # qqplot?
  # qqnorm(log10(Rsim_sel), pch = 1, frame = FALSE)
  # qqline(log10(Rsim_sel), col = "steelblue", lwd = 2)
  
  negSelSim<- c(0.01,0.2,0.4,0.6,0.8,0.85,0.9,0.95,0.99,1) # Simulated amount of remaining mutations after neg sel
  nMutSim<- c(10,30,100,300,1000,3000) # Simulated n samples
  
  # Calculate power given # samples and # negsel
  pwr_n_matrix<- matrix(NA,length(nMutSim),length(negSelSim),dimnames = list(nMutSim,negSelSim))
  for(all_total in rownames(pwr_n_matrix)){
    for(trueNegSel in colnames(pwr_n_matrix)){
      counter<- 0
      for(i in 1:1000){
        R_temp<- rnorm(all_total,mean = 0,sd = sd_log) 
        R_temp_sel<- R_temp+log10(as.numeric(trueNegSel))
        if(wilcox.test(R_temp_sel,mu=0)$p.value<0.05) counter<- counter + 1
      }
      pwr_n_matrix[all_total,trueNegSel]<- counter/1000
    }
  }
}
  
# 1.2 Add lower positive predictive values
##########################################
{
  negSelSim<- c(0.01,0.2,0.4,0.6,0.8,0.85,0.9,0.95,0.99,1) # Simulated amount of remaining mutations after neg sel
  nMutSim<- c(100,300,1000,3000) # Simulated n mutations
  ppv<- c(1,0.8,0.6,0.4,0.2,0.1)

  # Calculate power given # mut (rows) and # negsel
  pwr_ppv_matrix_ls<- list()
  pwr_ppv_Rsim_matrix_ls<- list()
  for(all_total in nMutSim){
    # cat(all_total,"\n")
    pwr_ppv_matrix<- matrix(NA,length(ppv),length(negSelSim),dimnames = list(ppv,negSelSim))
    pwr_ppv_Rsim_matrix<- matrix(NA,length(ppv),length(negSelSim),dimnames = list(ppv,negSelSim))
    for(ppv_temp in rownames(pwr_ppv_matrix)){
      for(trueNegSel in colnames(pwr_ppv_matrix)){
        counter<- 0
        R_temp_all<- NULL
        for(i in 1:1000){
          R_temp<- rnorm(all_total,mean = 0,sd = sd_log) 
          R_temp_sel<- R_temp+log10(1-((1-as.numeric(trueNegSel))*as.numeric(ppv_temp)))
          if(wilcox.test(R_temp_sel,mu=0)$p.value<0.05) counter<- counter + 1
          R_temp_all<- c(R_temp_all,R_temp_sel)
        }
        pwr_ppv_matrix[ppv_temp,trueNegSel]<- counter/1000
        pwr_ppv_Rsim_matrix[ppv_temp,trueNegSel]<- median(R_temp_all,na.rm=T) 
      }
    }
    pwr_ppv_matrix_ls<- c(pwr_ppv_matrix_ls,list(pwr_ppv_matrix))
    pwr_ppv_Rsim_matrix_ls<- c(pwr_ppv_Rsim_matrix_ls,list(pwr_ppv_Rsim_matrix))
  }
  names(pwr_ppv_matrix_ls)<- nMutSim
  names(pwr_ppv_Rsim_matrix_ls)<- nMutSim
}
  
# Plot
######
{
  # Numbers
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_Rsim_Power.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_Rsim_Power.pdf"))
    par(mfrow=c(2,1))
    for(i in 1:nrow(pwr_n_matrix)){
      if(i==1) plot(negSelSim,pwr_n_matrix[i,],type="l",ylab="Power",xlab="Expected dNHLA/dNnonHLA value",col=i,frame.plot=F,ylim=c(0,1.1))
      else points(negSelSim,pwr_n_matrix[i,],type="l",col=i)
    }
    legend("bottomleft",legend = rownames(pwr_n_matrix),col = 1:nrow(pwr_n_matrix),lty=1,cex=0.5)
    abline(v=0.8,h=0.8,lty=2)
    dev.off()
  }
}

# Plot ppv
##########
{
  # Numbers
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_Rsim_Power_ppv.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_Rsim_Power_ppv.pdf"))
    par(mfrow=c(2,2),oma=c(2,2,2,2))
    for(p in 1:length(pwr_ppv_matrix_ls)){
      # Power
      pwr_ppv_matrix<- pwr_ppv_matrix_ls[[p]]
      for(i in 1:nrow(pwr_ppv_matrix)){
        if(i==1){
          plot(colnames(pwr_ppv_matrix),pwr_ppv_matrix[i,],type="l",ylab="Power",xlab="Expected dNHLA/dNnonHLA",col=i,frame.plot=F,axes=F,ylim=c(0,1.1),main=names(pwr_ppv_matrix_ls)[p])
          axis(2)
          axis(3)
        }
        else points(colnames(pwr_ppv_matrix),pwr_ppv_matrix[i,],type="l",col=i)
      }
      legend("bottomleft",legend = rownames(pwr_ppv_matrix),col = 1:nrow(pwr_ppv_matrix),lty=1,cex=0.5)
      abline(v=0.8,h=0.8,lty=2)
      pwr_ppv_Rsim_matrix<- round(pwr_ppv_Rsim_matrix_ls[[p]],2)
      pwr_ppv_Rsim_matrix<- format(10^pwr_ppv_Rsim_matrix,digits=2) # Convert back to non-log scale
      for(i in 1:nrow(pwr_ppv_Rsim_matrix)){
        text(x = as.numeric(colnames(pwr_ppv_Rsim_matrix)),y=((-0.1*i)),xpd=NA,labels = pwr_ppv_Rsim_matrix[i,],cex=0.5,col=i)
        text(x = -0.2,y=((-0.1*i)),xpd=NA,labels = rownames(pwr_ppv_Rsim_matrix)[i])
      }
    }
    dev.off()
  }
}

# Calculate amount of negative selection for given sample size (and cancer type)
####################################################################################
{
  # Get sample sizes
  n_samples<- c(
    BLCA=390,
    BRCA=904,
    CESC=281,
    COADREAD=426,
    ESCA=124,
    GBM=388,
    HNSC=488,
    KIRP=222,
    LGG=502,
    LIHC=312,
    LUAD=407,
    LUSC=306,
    OV=333,
    PAAD=145,
    PRAD=434,
    SARC=216,
    SKCM=466,
    STAD=390,
    UCEC=506
    )
  
  # Focus on ratio's from 19 analysed cancers
  Rsim_sel<- Rsim[TCGA_cancer_id[names(Rsim),2]%in%cancers]
  Rsim_sel<- Rsim_sel[!is.na(Rsim_sel)]
  
  # Distribution after logtransformation +/- normal
  Rsim_sel<- Rsim_sel[Rsim_sel!=0]
  Rsim_sel<- Rsim_sel[!is.infinite(Rsim_sel)]
  sd_log<- sd(log10(Rsim_sel))

  # Set simulated values
  negSelSim<- as.character(seq(0.8,1,0.001))
  nMutSim<- unique(sort(c(n_samples,seq(100,1000,100)))) # Simulated n samples
  
  # Function to calculate moving average
  ma<- function(a,n=5){
    ma_a<- rep(NA,length(a))
    for(i in 1:length(ma_a)){
      start_tmp<- i-((n-1)/2)
      stop_tmp<- i+((n-1)/2)
      if(start_tmp<1) start_tmp=1
      if(stop_tmp>length(ma_a)) stop_tmp=length(ma_a)
      ma_a[i]<- sum(as.numeric(a)[start_tmp:stop_tmp])/length(start_tmp:stop_tmp)
    }
    return(ma_a)
  }
  
  # Calculate power given # mut and # negsel
  purSel_n_matrix<- matrix(NA,length(nMutSim),length(negSelSim),dimnames = list(nMutSim,negSelSim))
  purSel_n<- rep(NA,length(nMutSim))
  names(purSel_n)<- nMutSim
  for(all_total in rownames(purSel_n_matrix)){
    cat(all_total,"\n")
    for(trueNegSel in colnames(purSel_n_matrix)){
      counter<- 0
      for(i in 1:1000){
        R_temp<- rnorm(all_total,mean = 0,sd = sd_log) 
        R_temp_sel<- R_temp+log10(as.numeric(trueNegSel))
        if(wilcox.test(R_temp_sel,mu=0)$p.value<0.05) counter<- counter + 1
      }
      pwr_tmp<- counter/1000
      purSel_n_matrix[all_total,trueNegSel]<- pwr_tmp
    }
  # Het detectabke amount at power 0.8
  # purSel_n[all_total]<- names(which(abs(purSel_n_matrix[all_total,]-0.8)==min(abs(purSel_n_matrix[all_total,]-0.8))))
  purSel_n[all_total]<- colnames(purSel_n_matrix)[which(abs(ma(purSel_n_matrix[all_total,])-0.8)==min(abs(ma(purSel_n_matrix[all_total,])-0.8)))] # Take moving average to avoid getting outlier values
  }

  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_Rsim_Power_n_vs_purSel.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_Rsim_Power_n_vs_purSel.pdf"))
    # plot(log10(as.numeric(names(purSel_n))),100*(1-as.numeric(purSel_n)),xlim=c(2,3),ylim=c(0,20),xlab="# samples",ylab="Detectable amount of immunogenic purifying selection(%)",type="l",axes=F)
    plot(as.numeric(names(purSel_n)),100*(1-as.numeric(purSel_n)),xlim=c(100,1000),ylim=c(0,20),xlab="# samples",ylab="Detectable amount of immunogenic purifying selection(%)",type="l",frame.plot=F)
    purSel_n_can_matrix<- matrix(NA,length(cancers),2,dimnames = list(cancers,c("n","purSel")))
    for(cancer in cancers){
      # Get sample size
      n_samples_tmp<- n_samples[cancer]
      # Add to matrix
      purSel_tmp<- 100*(1-as.numeric(purSel_n[as.character(n_samples_tmp)]))
      purSel_n_can_matrix[cancer,]<- c(n_samples_tmp,purSel_tmp)
      # Add point
      points(n_samples_tmp,purSel_tmp,pch=16,col=c25[cancer])
    }
    purSel_n_can_matrix<- purSel_n_can_matrix[order(purSel_n_can_matrix[,"n"]),]
    # Highest power
    if(i==1){
      cat("Highest power: ", rownames(purSel_n_can_matrix)[nrow(purSel_n_can_matrix)], purSel_n_can_matrix[nrow(purSel_n_can_matrix),],"\n")
      cat("Lowest power: ", rownames(purSel_n_can_matrix)[1], purSel_n_can_matrix[1,],"\n")
      # Highest power:  BRCA 904 4.2
      # Lowest power:  ESCA 124 11
    }
    dev.off()
  }
}

# Save results
##############
{
  save(sd_log,pwr_n_matrix,pwr_ppv_matrix_ls,pwr_ppv_Rsim_matrix_ls,file = paste0("results/data/fig", fig_nr,"Rsim_power.RData"))
}
