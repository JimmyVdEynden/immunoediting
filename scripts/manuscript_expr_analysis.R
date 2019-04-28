###################################
# manuscript_expr_analysis.R
###################################

# Load data
###########
{
  TCGA_maf<- readRDS("data/TCGA_maf_mainHLAType.rds")
  load("data/GO_v62.RData")
}

# Libraries
############
{
  library(svglite)
  library(plotrix)
  source("scripts/functions/do_GSEA2.R")
  library(WriteXLS)
}

# 0) Preprocess & Filtering
#############################
{
  # Filter 1) Restrict analysis to synonymous and non-synonymous mutations (indels were removed previously, also remove nonsense)
  TCGA_maf<- TCGA_maf[!is.na(TCGA_maf[,"Variant_Classification"])&TCGA_maf[,"Variant_Classification"]%in%c("synonymous SNV","nonsynonymous SNV"),]

  # Turn barcodes in factors
  TCGA_maf[,"Tumor_Sample_Barcode"]<- factor(TCGA_maf[,"Tumor_Sample_Barcode"])

  # Define high affinity HLA binders (cut-off: Kd<500nM)
  isHA_wt<- TCGA_maf[,"wt_HLA_mean_aff"]<500
  isHA_mut<- TCGA_maf[,"mut_HLA_mean_aff"]<500
  TCGA_maf<- cbind(TCGA_maf,isHA_wt=isHA_wt,isHA_mut=isHA_mut)
}

# 1) GSEA on non-expressed genes
#################################
{
  # Focus on cancers from rest of study
  cancers<- sort(names(table(TCGA_maf[,"Cancer"]))[table(TCGA_maf[,"Cancer"])>10000])
  
  # Get expressed & non-expressed genes 
  TCGA_maf_nonExpr<- TCGA_maf[!is.na(TCGA_maf[,"mRNA"])&TCGA_maf[,"mRNA"]<=expr_cu,]
  TCGA_maf_expr<- TCGA_maf[!is.na(TCGA_maf[,"mRNA"])&TCGA_maf[,"mRNA"]>expr_cu,]
  prop.table(table(TCGA_maf_nonExpr[,"isHA_wt"])) 
  prop.table(table(TCGA_maf_expr[,"isHA_wt"])) 
  TCGA_maf_nonExpr_genes<- unique(TCGA_maf_nonExpr[,1]) 
  TCGA_maf_expr_genes<- unique(TCGA_maf_expr[,1])
  cat("nr of expressed genes: ", length(TCGA_maf_expr_genes)," (HLA binding:",format(100*prop.table(table(TCGA_maf_expr[,"isHA_wt"]))["TRUE"],digits=3),"%)","  \n")
  cat("nr of non-expressed genes: ", length(TCGA_maf_nonExpr_genes)," (HLA binding:",format(100*prop.table(table(TCGA_maf_nonExpr[,"isHA_wt"]))["TRUE"],digits=3),"%)","  \n")
  
  # Pan cancer
  expr_can_t<- table(Cancer=TCGA_maf$Cancer, isExpr=TCGA_maf$mRNA>expr_cu, isHA=TCGA_maf$isHA_wt)
  expr_Pancan_t<- colSums(expr_can_t)
  #         isHA
  # isExpr  FALSE   TRUE
  # FALSE  152093  60144
  # TRUE  1154891 319766
  cat("Pan: P=",fisher.test(expr_Pancan_t)$p.value," \n")
  print.data.frame(as.data.frame.matrix(expr_Pancan_t))  

  # Barplots per cancer type
  expr_can_prop_HA<- prop.table(expr_can_t[cancers,"TRUE",],1)[,"TRUE"]
  expr_can_prop_nonHA<- prop.table(expr_can_t[cancers,"FALSE",],1)[,"TRUE"]
  p<- NULL
  for(cancer in cancers){
    p<- c(p,fisher.test(rbind(expr_can_t[cancer,"TRUE",],expr_can_t[cancer,"FALSE",]))$p.value)
  }
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_expr.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_expr.pdf"))
    bp<- barplot(rbind(expr_can_prop_HA,expr_can_prop_nonHA),beside=T,las=2,ylab="% HLA binding",col=c("black","grey"))
    text(bp[1,],0.3,xpd=NA,labels = signif(p,2),srt=90)
    dev.off()
  }
  
  # Split by synonymous/non-synonymous
  
  # Pan cancer
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_expr_ns_pan.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_expr_ns_pan.pdf"))
    # par(mfrow=c(1,2))
    expr_can_t<- table(Cancer=TCGA_maf$Cancer, isExpr=TCGA_maf$mRNA>expr_cu, isHA=TCGA_maf$isHA_wt,TCGA_maf$Variant_Classification)
    expr_Pancan_t<- colSums(expr_can_t)
    expr_pan_prop_s<- prop.table(expr_Pancan_t[,,"synonymous SNV"],1)[,"TRUE"]
    expr_pan_prop_n<- prop.table(expr_Pancan_t[,,"nonsynonymous SNV"],1)[,"TRUE"]
    p<- c(fisher.test(expr_Pancan_t[,,"synonymous SNV"])$p.value,fisher.test(expr_Pancan_t[,,"nonsynonymous SNV"])$p.value)
    bp<- barplot(100*cbind(s=expr_pan_prop_s,n=expr_pan_prop_n),beside=T,ylab="% HLA binding",col=c("black","grey"))
    text(bp,-2,xpd=NA,labels = paste0(100*cbind(signif(expr_pan_prop_s,2),signif(expr_pan_prop_n,2)),"%"))
    text(bp[1,],35,xpd=NA,labels = paste0("P=",signif(p,2)))
    dev.off()
  }

  # Per cancer
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_expr_ns.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_expr_ns.pdf"))
    par(mfrow=c(2,1))
    for(variant in c("synonymous SNV","nonsynonymous SNV")){
      p<- NULL
      for(cancer in cancers){
        p<- c(p,fisher.test(rbind(expr_can_t[cancer,"TRUE",,variant],expr_can_t[cancer,"FALSE",,variant]))$p.value)
      }
      expr_can_prop_HA_s<- prop.table(expr_can_t[cancers,"TRUE",,variant],1)[,"TRUE"]
      expr_can_prop_nonHA_s<- prop.table(expr_can_t[cancers,"FALSE",,variant],1)[,"TRUE"]
      bp<- barplot(rbind(expr_can_prop_HA_s,expr_can_prop_nonHA_s),beside=T,las=2,ylab="% HLA binding",col=c("black","grey"),main=variant)
      text(bp[1,],0.3,xpd=NA,labels = signif(p,2),srt=90)
    }
    dev.off()
  }

  # GSEA
  GSEA_nonExpr_GO<- do_GSEA2(TCGA_maf_nonExpr_genes,union(TCGA_maf_nonExpr_genes,TCGA_maf_expr_genes),GO_ls,isList = T,min_genes = 2)
  GSEA_expr_GO<- do_GSEA2(TCGA_maf_expr_genes,union(TCGA_maf_nonExpr_genes,TCGA_maf_expr_genes),GO_ls,isList = T,min_genes = 2)
  WriteXLS(c("GSEA_nonExpr_GO","GSEA_expr_GO"),ExcelFileName = paste0("results/tables/table_",fig_nr,".xls"),SheetNames = c("GO not expressed","GO expressed"),row.names = TRUE)
  
  # Plot GSEA
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_expr_GSEA.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_expr_GSEA.pdf"))
    par(mfrow=c(2,1))
    barplot(-log10(as.numeric(as.character(GSEA_nonExpr_GO[10:1,"q"]))),horiz = T,xlim=c(0,200),xlab="GSEA -log10(q) value")
    text(rep(0,10),seq(11.5,0.7,by = -1.2),labels = rownames(GSEA_nonExpr_GO)[1:10],adj = 0)  
    barplot(-log10(as.numeric(as.character(GSEA_expr_GO[10:1,"q"]))),horiz = T,xlim=c(0,200),xlab="GSEA -log10(q) value")
    text(rep(0,10),seq(11.5,0.7,by = -1.2),labels = rownames(GSEA_expr_GO)[1:10],adj = 0)  
    dev.off()
  }
}

# 2) nHP AA in each group?
###########################
{
  GPPM_aa_t<- readRDS("data/GPPM_gene_aa_t.rds")
  load("data/aa_classes.RData")
  GPPM_aa_t_propHP<- rowSums(GPPM_aa_t[,aa_hp])/rowSums(GPPM_aa_t)
  GPPM_aa_t_HP<- rowSums(GPPM_aa_t[,aa_hp])
  GPPM_aa_t_all<- rowSums(GPPM_aa_t) 
  
  # Average per enriched pathway?
  propHP_nonExpr_ls<- list()
  for(i in 10:1){
    propHP_nonExpr_ls<- c(propHP_nonExpr_ls,list(GPPM_aa_t_propHP[GO_ls[[rownames(GSEA_nonExpr_GO)[i]]]]))
  }
  
  propHP_expr_ls<- list()
  for(i in 10:1){
    propHP_expr_ls<- c(propHP_expr_ls,list(GPPM_aa_t_propHP[GO_ls[[rownames(GSEA_expr_GO)[i]]]]))
  }
  
  # Plot
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_expr_GSEA_isHP.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_expr_GSEA_isHP.pdf"))
    par(mfcol=c(2,2))
    boxplot(propHP_nonExpr_ls,horizontal=T,outline=F,ylim=c(0.3,0.7),axes=F)  
    axis(1)
    boxplot(propHP_expr_ls,horizontal=T,outline=F,ylim=c(0.3,0.7),axes=F)  
    axis(1)
    dev.off()
  }
}

# 3) Drivers genes?
####################
{
  gene_var_HA_t<- readRDS(file = "data/GPPM_gene_var_HA_t.rds")
  TCGA_t<- table(TCGA_maf$Hugo_Symbol,TCGA_maf$Variant_Classification,TCGA_maf$wt_HLA_mean_aff<500)
  CGC_v83 <- readRDS("data/CGC_v83.rds")
  
  # %CGC?
  TCGA_CGC_t<- table(TCGA_maf$Hugo_Symbol%in%CGC_v83,TCGA_maf$wt_HLA_mean_aff<500)
  cat("nr of driver gene muts: ", sum(TCGA_CGC_t["TRUE",])," (HLA binding:",format(100*prop.table(TCGA_CGC_t["TRUE",])["TRUE"],digits=3),"%)","  \n")
  cat("nr of non-driver gene muts: ", sum(TCGA_CGC_t["FALSE",])," (HLA binding:",format(100*prop.table(TCGA_CGC_t["FALSE",])["TRUE"],digits=3),"%)","  \n")
  cat("Pan: P=",fisher.test(TCGA_CGC_t)$p.value," \n")
  print.data.frame(as.data.frame.matrix(TCGA_CGC_t))  
  
  for(i in 1:2){
    if(i==1) svglite(paste0("results/figs/svg/fig",fig_nr,"_driver_pan.svg"))
    else pdf(paste0("results/figs/pdf/fig",fig_nr,"_driver_pan.pdf"))
    par(mfrow=c(1,2))
    # TCGA data

    variant<- "synonymous SNV"
    propHLA<- 100*TCGA_t[,variant,"TRUE"]/(TCGA_t[,variant,"FALSE"]+TCGA_t[,variant,"TRUE"])
    isDriver<- names(propHLA)%in%CGC_v83
    boxplot(list(propHLA[!isDriver],propHLA[isDriver],NA,NA,NA),outline=F,frame.plot=F,axes=F,ylab="% HLA binding",main="TCGA",ylim=c(0,65))
    p_temp<- signif(wilcox.test(propHLA[!isDriver],propHLA[isDriver])$p.value,3)
    text(1.5,60,labels = paste0("P=",p_temp),xpd=NA)
    text(c(1,2),c(-2,-2),labels = c(paste0(signif(median(propHLA[!isDriver],na.rm=T),3),"%"),paste0(signif(median(propHLA[isDriver],na.rm=T),3),"%")),xpd=NA)
    
    variant<- "nonsynonymous SNV"
    propHLA<- 100*TCGA_t[,variant,"TRUE"]/(TCGA_t[,variant,"FALSE"]+TCGA_t[,variant,"TRUE"])
    isDriver<- names(propHLA)%in%CGC_v83
    boxplot(list(NA,NA,NA,propHLA[!isDriver],propHLA[isDriver]),frame.plot=F,outline=F,axes=F,add=T)
    p_temp<- signif(wilcox.test(propHLA[!isDriver],propHLA[isDriver])$p.value,3)
    text(4.5,60,labels = paste0("P=",p_temp),xpd=NA)
    text(c(4,5),c(-2,-2),labels = c(paste0(signif(median(propHLA[!isDriver],na.rm=T),3),"%"),paste0(signif(median(propHLA[isDriver],na.rm=T),3),"%")),xpd=NA)
    
    axis(1,at=c(1,2,4,5),labels = c("s-","s+","n-","n+")) 
    axis(2)      
    
    # GPPM data
    variant<- "synonymous SNV"
    propHLA<- 100*gene_var_HA_t[,variant,"TRUE"]/(gene_var_HA_t[,variant,"FALSE"]+gene_var_HA_t[,variant,"TRUE"])
    isDriver<- names(propHLA)%in%CGC_v83
    boxplot(list(propHLA[!isDriver],propHLA[isDriver],NA,NA,NA),outline=F,frame.plot=F,axes=F,ylab="% HLA binding",main="Whole exome",ylim=c(0,65))
    p_temp<- signif(wilcox.test(propHLA[!isDriver],propHLA[isDriver])$p.value,3)
    text(1.5,60,labels = paste0("P=",p_temp),xpd=NA)
    text(c(1,2),c(-2,-2),labels = c(paste0(signif(median(propHLA[!isDriver],na.rm=T),3),"%"),paste0(signif(median(propHLA[isDriver],na.rm=T),3),"%")),xpd=NA)
    
    variant<- "nonsynonymous SNV"
    propHLA<- 100*gene_var_HA_t[,variant,"TRUE"]/(gene_var_HA_t[,variant,"FALSE"]+gene_var_HA_t[,variant,"TRUE"])
    isDriver<- names(propHLA)%in%CGC_v83
    boxplot(list(NA,NA,NA,propHLA[!isDriver],propHLA[isDriver]),frame.plot=F,outline=F,axes=F,add=T)
    p_temp<- signif(wilcox.test(propHLA[!isDriver],propHLA[isDriver])$p.value,3)
    text(4.5,60,labels = paste0("P=",p_temp),xpd=NA)
    text(c(4,5),c(-2,-2),labels = c(paste0(signif(median(propHLA[!isDriver],na.rm=T),3),"%"),paste0(signif(median(propHLA[isDriver],na.rm=T),3),"%")),xpd=NA)
    
    axis(1,at=c(1,2,4,5),labels = c("s-","s+","n-","n+")) 
    axis(2)      
    
    dev.off()
  }
}

# Save all data
################
{
  FULL_FILENAME <- parent.frame(2)$ofile # Only works when sourced, otherwise is.null
  FULL_FILENAME<- gsub("scripts/","",FULL_FILENAME)
  FULL_FILENAME<- paste0("fig",fig_nr,"_",FULL_FILENAME)
  save(GSEA_nonExpr_GO,GSEA_expr_GO,file = paste0("results/data/",FULL_FILENAME,"Data"))
}

