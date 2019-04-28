#########################################
# manuscript_cor_triNT_HBMR_aaClass.R
#########################################

# Load
#######
{
  GPPM<- readRDS(file="data/GPPM_inclHLAAlleles.rds") # Takes some time to load
}

# Libraries
############
{
  library(svglite)
}

# Preprocess and filtering
###########################
{
  length(GPPM) # 93739662
  
  # Filter 1) Restrict analysis to synonymous and non-synonymous mutations (remove nonsense)
  idx_NS<- which(GPPM$variant%in%c("nonsynonymous SNV","synonymous SNV"))
  GPPM<- GPPM[idx_NS]

  # Define high affinity HLA binders (cut-off: Kd<500nM)
  GPPM$isHA<- GPPM$HLA_aff_mean<500
}

# Some numbers about GPPM object
################################
{
  # General numbers
  cat("n genes:",length(unique(GPPM$gene)),"  \n") # 17,992
  cat("n mutations:",length(GPPM),"  \n") # n mutations: 88,970,246   
  var_t<- table(GPPM$variant)
  cat("n syn mutations:",var_t["synonymous SNV"],"  \n") # n syn mutations: 21,203,704   
  cat("n nonsyn mutations:",var_t["nonsynonymous SNV"],"  \n") # n nonsyn mutations: 67,766,542
  cat("prop HLA binding:",format(mean(GPPM$isHA),digits=3),"  \n","  \n") # 0.222
}

# 1) n/s analysis
#############################
{
  # Nonsynonymous to synonymous? (ratio N/S) 
  sim_all_isHA_variant_t<- table(variant=GPPM$variant,isHA=GPPM$isHA)
  sim_all_isHA_variant_t<- sim_all_isHA_variant_t[c("synonymous SNV","nonsynonymous SNV"),]
  colnames(sim_all_isHA_variant_t)<- c("non-HLA binding","HLA binding")
  print.data.frame(as.data.frame.matrix(sim_all_isHA_variant_t))  
  # non-HLA binding HLA binding
  # synonymous SNV           16194676     5009028
  # nonsynonymous SNV        53036068    14730474
  
  sim_all_ns<- sim_all_isHA_variant_t["nonsynonymous SNV",]/sim_all_isHA_variant_t["synonymous SNV",]
  sim_all_fisher_test<- fisher.test(sim_all_isHA_variant_t)
  sim_all_p<- format(sim_all_fisher_test$p.value,digits = 3)
  sim_all_HBMR<- format(sim_all_fisher_test$estimate,digits = 3)
  sim_all_prop_HLA<- prop.table(sim_all_isHA_variant_t,1)[,"HLA binding"] 
  
  cat("n/s HLA binding:",sim_all_ns["HLA binding"],"  \n") # 2.940785   
  cat("n/s non-HLA binding:",sim_all_ns["non-HLA binding"],"  \n") # 3.274908 
  # Note these ratios are not usefull as biased themselves by mut signatures!
  cat("p=",sim_all_p,"  \n") # p=0
  cat("N/S HBMR=",sim_all_HBMR,"  \n") # N/S HBMR=0.898
  cat("prop HLA synonymous:",format(sim_all_prop_HLA["synonymous SNV"],digits=3),"  \n") # 0.236
  cat("prop HLA non-synonymous:",format(sim_all_prop_HLA["nonsynonymous SNV"],digits=3),"  \n") #  0.217
  
  # Put in pie chart
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"_NS_pie.pdf"))
    else svglite(paste0("results/figs/svg/fig",fig_nr,"_NS_pie.svg"))
    par(mfrow=c(2,1))
    pie(sim_all_isHA_variant_t["synonymous SNV",],col = c("white","grey"),labels = c("","HLA binding"),main="Synonymous")
    pie(sim_all_isHA_variant_t["nonsynonymous SNV",],col = c("white","grey"),labels = c("","HLA binding"),main="Non-synonymous")
    dev.off()
  }
}

# Effect of triNT signature for silent/non-silent
#################################################
{
  # Color for plots
  substTypes<-c("C>A","C>G","C>T","T>A","T>C","T>G")
  substTypes_cols<-c("blue","black","red","grey","green","pink")
  names(substTypes_cols)<-substTypes
  
  # subst type effect
  sim_all_isHA_substType_variant_t<- table(subst_type=GPPM$subst_type,variant=factor(GPPM$variant,levels=c("synonymous SNV","nonsynonymous SNV")),isHA=GPPM$isHA)
  substType_matrix<- matrix(NA,length(substTypes),6,dimnames = list(substTypes,c("ns-","ns+","s+","n+","HBMR","p")))
  for(substType in substTypes){
    sim_all_isHA_substType_variant_t_temp<- sim_all_isHA_substType_variant_t[substType,,]
    colnames(sim_all_isHA_substType_variant_t_temp)<- c("non-HLA binding","HLA binding")
    substType_matrix[substType,c("ns-","ns+")]<- as.numeric(format(sim_all_isHA_substType_variant_t_temp["nonsynonymous SNV",]/sim_all_isHA_substType_variant_t_temp["synonymous SNV",],digits=3))
    sim_all_fisher_test_temp<- fisher.test(sim_all_isHA_substType_variant_t_temp)
    substType_matrix[substType,"p"]<- as.numeric(format(sim_all_fisher_test_temp$p.value,digits = 3))
    substType_matrix[substType,"HBMR"]<- as.numeric(format(sim_all_fisher_test_temp$estimate,digits = 3))
    substType_matrix[substType,c("s+","n+")]<- as.numeric(format(prop.table(sim_all_isHA_substType_variant_t_temp,1)[,"HLA binding"],digits=2)) 
  }
  print.data.frame(as.data.frame.matrix(substType_matrix))  
  #     ns-  ns+   s+   n+  HBMR       p
  # C>A 3.95 3.00 0.27 0.22 0.759 0.00e+00
  # C>G 5.11 3.75 0.27 0.21 0.734 0.00e+00
  # C>T 1.59 1.31 0.24 0.21 0.823 0.00e+00
  # T>A 4.82 4.86 0.22 0.22 1.010 2.69e-05
  # T>C 2.25 2.57 0.21 0.23 1.140 0.00e+00
  # T>G 5.03 5.06 0.22 0.22 1.000 3.76e-03
  
  # triNT effect
  sim_all_isHA_triNT_variant_t<- table(subst_type=GPPM$subst_type3,variant=factor(GPPM$variant,levels=c("synonymous SNV","nonsynonymous SNV")),isHA=GPPM$isHA)
  substTypes3<- unique(GPPM$subst_type3)
  triNT_matrix<- matrix(NA,length(substTypes3),6,dimnames = list(substTypes3,c("ns-","ns+","s+","n+","HBMR","p")))
  for(substType3 in substTypes3){
    cat(substType3," \n")
    sim_all_isHA_triNT_variant_t_temp<- sim_all_isHA_triNT_variant_t[substType3,,]
    colnames(sim_all_isHA_triNT_variant_t_temp)<- c("non-HLA binding","HLA binding")
    triNT_matrix[substType3,c("ns-","ns+")]<- as.numeric(format(sim_all_isHA_triNT_variant_t_temp["nonsynonymous SNV",]/sim_all_isHA_triNT_variant_t_temp["synonymous SNV",],digits=3))
    sim_all_fisher_test_temp<- fisher.test(sim_all_isHA_triNT_variant_t_temp)
    triNT_matrix[substType3,"p"]<- as.numeric(format(sim_all_fisher_test_temp$p.value,digits = 3))
    triNT_matrix[substType3,"HBMR"]<- as.numeric(format(sim_all_fisher_test_temp$estimate,digits = 3))
    triNT_matrix[substType3,c("s+","n+")]<- as.numeric(format(prop.table(sim_all_isHA_triNT_variant_t_temp,1)[,"HLA binding"],digits=2)) 
  }
  print.data.frame(as.data.frame.matrix(triNT_matrix[order(triNT_matrix[,"HBMR"]),]))  

  # Exclude 4 subs types without synonymous mut 
  cat("No synonymous: ",names(which(is.na(triNT_matrix[,"s+"]))),"  \n") # No synonymous:  ATT>AAT ATT>AGT ACT>AGT ACT>AAT
  triNT_matrix<- triNT_matrix[!is.na(triNT_matrix[,"s+"]),]

  # Plot HBMR
  for(i in 1:2){
    if(i==1) pdf(paste0("results/figs/pdf/fig",fig_nr,"_HBMR_triNT.pdf"),width = 14)
    else svglite(paste0("results/figs/svg/fig",fig_nr,"_HBMR_triNT.svg"),width = 20)
    par(oma=c(2,2,2,2))
    barplot(sort(triNT_matrix[,"HBMR"]),beside = T,las=2,col=substTypes_cols[paste0(substr(names(sort(triNT_matrix[,"HBMR"])),2,2),">",substr(names(sort(triNT_matrix[,"HBMR"])),6,6))],ylab="HBMR",border = NA)
    legend("topleft",legend = substTypes,text.col = substTypes_cols,bty="n")
    dev.off()
  }
}

# Save all data
################
{
  FULL_FILENAME <- parent.frame(2)$ofile # Only works when sourced, otherwise is.null
  FULL_FILENAME<- gsub("scripts/","",FULL_FILENAME)
  FULL_FILENAME<- paste0("fig",fig_nr,"_",FULL_FILENAME)
  save(triNT_matrix,substType_matrix,file = paste0("results/data/",FULL_FILENAME,"Data"))
}
