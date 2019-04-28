#############################################
# manuscript_logreg_aa_HLABinder.R
#############################################

# Determine expected HLA binding for nonapeptide with n number of amino acids using a logistic regression model

# Functions and libraries
library(svglite)
library(corrplot)
library(RColorBrewer)

# Load aa classes
load("data/aa_classes.RData")

if(runIndividualPeptides==FALSE){
  # Load GPPM
  GPPM<- readRDS(file="data/GPPM_inclHLAAlleles.rds") # Takes some time to load
  
  # Take mean number of amino acids from 9 nonapeptides for 1 mill random locationq
  idx_pep<- sample(1:length(GPPM),size = 1000000,replace = F)
  pep<- GPPM$nonaPep_wt[idx_pep]
  pep_matrix<- matrix(NA,length(pep),length(aa_all),dimnames = list(NULL,aa_all))
  for(i in 1:nrow(pep_matrix)){
    # cat(i," ")
    pep_matrix[i,]<- table(factor(unlist(strsplit(paste0(unlist(pep[i]),collapse = ""),"")),levels=aa_all))
  }
  # Total from 9 peptides, calculate mean per peptide
  pep_matrix<- 9*(pep_matrix/rowSums(pep_matrix))
  
  # Get affinities for all peptides
  pep_aff<- GPPM$HLA_aff_mean[idx_pep]
}

if(runIndividualPeptides==TRUE){
  load("data/nonapep_random.RData")
  pep_aff<- nonapep_random_HLA_mean_aff
}

# Logistic regression
pep_isHA<- pep_aff<500
# aa_logreg_matrix<- matrix(NA,length(aa_all),2,dimnames = list(aa_all,c("estimate","p")))
# for(aa in aa_all){
#   # cat(aa," ")
#   n_aa<- as.numeric(pep_matrix[,aa])
#   glm_model<- glm(pep_isHA~n_aa,family=binomial(link = "logit"))
#   aa_logreg_matrix[aa,"estimate"]<- glm_model$coefficients["n_aa"]
#   aa_logreg_matrix[aa,"p"]<- summary(glm_model)$coefficients["n_aa",4]
# }
# cat("Logistic regression results:", " \n")
# print.data.frame(format(as.data.frame.matrix(aa_logreg_matrix[order(aa_logreg_matrix[,"estimate"]),]),digits=3))  
# 
# # Visualize regression coefficients
# gen_code_matrix<- rbind(matrix(GENETIC_CODE[1:16],4,4),matrix(GENETIC_CODE[17:32],4,4),matrix(GENETIC_CODE[33:48],4,4),matrix(GENETIC_CODE[49:64],4,4))
# gen_code_matrix_regrCoeff<- cbind(
#   aa_logreg_matrix[gen_code_matrix[,1],"estimate"],
#   aa_logreg_matrix[gen_code_matrix[,2],"estimate"],
#   c(aa_logreg_matrix[gen_code_matrix[1:2,3],"estimate"],NA,NA,aa_logreg_matrix[gen_code_matrix[5:16,3],"estimate"]),
#   c(aa_logreg_matrix[gen_code_matrix[1:2,4],"estimate"],NA,aa_logreg_matrix[gen_code_matrix[4:16,4],"estimate"])
# )
# gen_code_matrix_regrP<- cbind(
#   aa_logreg_matrix[gen_code_matrix[,1],"p"],
#   aa_logreg_matrix[gen_code_matrix[,2],"p"],
#   c(aa_logreg_matrix[gen_code_matrix[1:2,3],"p"],NA,NA,aa_logreg_matrix[gen_code_matrix[5:16,3],"p"]),
#   c(aa_logreg_matrix[gen_code_matrix[1:2,4],"p"],NA,aa_logreg_matrix[gen_code_matrix[4:16,4],"p"])
# )
# 
# for(i in 1:2){
#   if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_logreg_aa_HLABinder_regrCoeff.svg"))
#   else pdf(paste0("results/figs/pdf/fig",fig_nr,"_logreg_aa_HLABinder_regrCoeff.pdf"))
#   corrplot(cbind(gen_code_matrix_regrCoeff,0,scale=seq(-1,1,by=0.25)),p.mat = cbind(gen_code_matrix_regrP,1,0),sig.level=0.05,insig = "blank", col=brewer.pal(n = 9, name = "PuOr"),cl.pos = "b")
#   dev.off()
# }
# 
# # Visualize regression lines for ind aa
# for(i in 1:2){
#   if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_logreg_aa_HLABinder_regrLines.svg"))
#   else pdf(paste0("results/figs/pdf/fig",fig_nr,"_logreg_aa_HLABinder_regrLines.pdf"))
#   x<-seq(0,9,0.1)
#   for(j in 1:length(aa_all)){
#     aa<- aa_all[j]
#     if(aa%in%aa_hp) col_aa<- "green"
#     else if(aa%in%aa_polar) col_aa<- "blue"
#     else col_aa<- "red"
#     
#     n_aa<- as.numeric(pep_matrix[,aa])
#     glm_model<- glm(pep_isHA~n_aa,family=binomial(link = "logit"))
#     pred1<-predict(glm_model,data.frame(n_aa=x),type="response",se.fit=TRUE)
#     
#     if(j==1){
#       plot(x,100*pred1$fit,type="l",frame.plot=F,xlab="Number of amino acids in nonapeptide",ylab="% HLA binding nonapeptides",ylim=c(0,100),col=col_aa)
#       abline(h=100*mean(pep_isHA,na.rm=T),lty=2)
#     }
#     else points(x,100*pred1$fit,type="l",col=col_aa)
#     text(x = 9.1,y = 100*pred1$fit[length(x)],labels = aa,xpd=NA)
#   }
#   dev.off()
# }
# 
# # Multiple regression
# pep_matrix_fused<- as.data.frame(cbind(pep_isHA,pep_matrix))
# glm_model<- glm(pep_isHA~A+F+G+I+L+M+P+V+W+C+N+Q+S+T+Y+D+E+H+K+R,data = pep_matrix_fused,family=binomial(link = "logit"))
# glm_multiple_model<- glm_model
# glm_coeff<- summary(glm_model)$coefficients
# glm_coeff<- as.data.frame(glm_coeff[order(glm_coeff[,"Estimate"]),])
# cat("Multiple logistic regression results:", " \n")
# print.data.frame(format(glm_coeff,digits=3))
# 
# # # 1 HLA type? Checked, no difference

# Repeat for main classes only
aa_classes<- list(aa_hp,aa_polar,aa_charged)
aa_logreg_matrix_classes<- matrix(NA,length(aa_classes),2,dimnames = list(c("hp","polar","charged"),c("estimate","p")))
for(i in 1:length(aa_classes)){
  aa<- aa_classes[[i]]
  n_aa<- rowSums(pep_matrix[,aa])
  glm_model<- glm(pep_isHA~n_aa,family=binomial(link = "logit"))
  aa_logreg_matrix_classes[i,"estimate"]<- glm_model$coefficients["n_aa"]
  aa_logreg_matrix_classes[i,"p"]<- summary(glm_model)$coefficients["n_aa",4]
}
cat("Logistic regression results for main aa classes:", " \n")
print.data.frame(format(as.data.frame.matrix(aa_logreg_matrix_classes),digits=3))  

# Visualize regression lines for aa classes
for(i in 1:2){
  if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_logreg_aa_class_HLABinder.svg"))
  else pdf(paste0("results/figs/pdf/fig",fig_nr,"_logreg_aa_class_HLABinder.pdf"))
  x<-seq(0,9,0.1)
  for(j in 1:length(aa_classes)){
    aa<- aa_classes[[j]]
    if(j==1) col_aa<- "green"
    else if(j==2) col_aa<- "blue"
    else col_aa<- "red"
    
    n_aa<- rowSums(pep_matrix[,aa])
    glm_model<- glm(pep_isHA~n_aa,family=binomial(link = "logit"))
    pred1<-predict(glm_model,data.frame(n_aa=x),type="response",se.fit=TRUE)
    
    if(j==1){
      plot(x,100*pred1$fit,type="l",frame.plot=F,xlab="Number of amino acids in nonapeptide",ylab="% HLA binding nonapeptides",col=col_aa)
      abline(h=100*mean(pep_isHA,na.rm=T),lty=2)
    }
    else points(x,100*pred1$fit,type="l",col=col_aa)
    text(x = 9.1,y = 100*pred1$fit[length(x)],labels = c("hp","polar","charged")[j],xpd=NA)
  }
  dev.off()
}

# Save all data
################
{
  FULL_FILENAME <- parent.frame(2)$ofile # Only works when sourced, otherwise is.null
  FULL_FILENAME<- gsub("scripts/","",FULL_FILENAME)
  FULL_FILENAME<- paste0("fig",fig_nr,"_",FULL_FILENAME)
  # save(glm_multiple_model,aa_logreg_matrix,aa_logreg_matrix_classes,gen_code_matrix,file=paste0("results/data/",FULL_FILENAME,"Data"))
  save(aa_logreg_matrix_classes,file=paste0("results/data/",FULL_FILENAME,"Data"))
}

