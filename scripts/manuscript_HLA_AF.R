#####################
# manuscript_HLA_AF.R
#####################

library(svglite)

# HLA AF in TCGA population
###########################
{
  # Load HLA types
  load("data/TCGA_HLA_types.RData")
  
  # Most frequent HLA alleles from TCGA
  HLA_A_t<- sort(table(c(TCGA_HLA_types_netMHC[,"HLA-A1"],TCGA_HLA_types_netMHC[,"HLA-A2"])),decreasing = T)
  HLA_B_t<- sort(table(c(TCGA_HLA_types_netMHC[,"HLA-B1"],TCGA_HLA_types_netMHC[,"HLA-B2"])),decreasing = T)
  HLA_C_t<- sort(table(c(TCGA_HLA_types_netMHC[,"HLA-C1"],TCGA_HLA_types_netMHC[,"HLA-C2"])),decreasing = T)
  HLA_freq<- c(HLA_A_t[1:2],HLA_B_t[1:2],HLA_C_t[1:2])
  HLA_freq<- HLA_freq/nrow(TCGA_HLA_types_netMHC)
  cat("HLA allele frequencies in TCGA:","  \n")
  for(i in 1:6){cat(names(HLA_freq)[i],signif(100*HLA_freq[i],3)," \n")}
}

# HLA AF in general (Caucasian) population
###################################################
{
  HLA_freq<- readRDS("data/AFN_HLA_freq.rds")
  HLA_freq<- HLA_freq[HLA_freq[,"HLA"]!="DQ",] # Only major alleles (i.e. A, B, C)
  rownames(HLA_freq)<- paste("HLA-",rownames(HLA_freq),sep="")
  rownames(HLA_freq)<- gsub("\\*","",rownames(HLA_freq))
  HLA_freq_A<- HLA_freq[grepl("HLA-A",rownames(HLA_freq)),]
  HLA_freq_A<- HLA_freq_A[order(HLA_freq_A[,"US_Bethesda_Cau"],decreasing = T),]
  HLA_freq_B<- HLA_freq[grepl("HLA-B",rownames(HLA_freq)),]
  HLA_freq_B<- HLA_freq_B[order(HLA_freq_B[,"US_Bethesda_Cau"],decreasing = T),]
  HLA_freq_C<- HLA_freq[grepl("HLA-C",rownames(HLA_freq)),]
  HLA_freq_C<- HLA_freq_C[order(HLA_freq_C[,"US_Bethesda_Cau"],decreasing = T),]
  cat( "\n")
  cat("HLA allele frequencies in general population:","  \n")
  for(HLA_freq_allele in list(HLA_freq_A,HLA_freq_B,HLA_freq_C)){
    cat(rownames(HLA_freq_allele)[1],HLA_freq_allele[1,"US_Bethesda_Cau"]," \n")
    cat(rownames(HLA_freq_allele)[2],HLA_freq_allele[2,"US_Bethesda_Cau"]," \n")
  }
}

# Correlation between TCGA and General population
#################################################
{
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_HLA_AFs.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_HLA_AFs.pdf"))
    par(mfrow=c(2,2))
    for(HLA_allele in c("A","B","C")){
      TCGA_freq_temp<- 200*prop.table(sort(table(c(TCGA_HLA_types_netMHC[,paste0("HLA-",HLA_allele,"1")],TCGA_HLA_types_netMHC[,paste0("HLA-",HLA_allele,"2")])),decreasing=T))
      Pop_freq_temp<- HLA_freq[grepl(HLA_allele,rownames(HLA_freq)),][names(TCGA_freq_temp),"US_Bethesda_Cau"]
      plot(as.numeric(TCGA_freq_temp),Pop_freq_temp,pch=16,col="blue",xlab="TCGA (%)",ylab = "Population (%)",main=paste("HLA",HLA_allele),frame.plot=F,cex=1.5)
      abline(c(0,1))
      cor_temp<- cor.test(as.numeric(TCGA_freq_temp),Pop_freq_temp)
      p_temp<- signif(cor_temp$p.value,3)
      r_temp<- signif(cor_temp$estimate,3)
      legend("topleft",legend = c(paste("p=",p_temp,sep=""),paste("r=",r_temp,sep="")),bty="n")  
      text(TCGA_freq_temp[1],Pop_freq_temp[1],labels = names(TCGA_freq_temp)[1],xpd=NA,adj = c(0.5,1.5))
      text(TCGA_freq_temp[2],Pop_freq_temp[2],labels = names(TCGA_freq_temp)[2],xpd=NA,adj = c(0.5,1.5))
    }
    dev.off()
  }
}



