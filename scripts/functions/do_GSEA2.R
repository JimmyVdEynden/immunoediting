do_GSEA2<- function(genes_retrieved,genes_all,GSEA_db,min_genes,isList=FALSE){
  if(isList) GSEA_names<-names(GSEA_db)
  else GSEA_names<-rownames(GSEA_db)
  GSEA_table<-matrix(NA, length(GSEA_names), 8, dimnames=list(GSEA_names,c("n_genes_pw","n_genes_pw_neg","n_genes_pw_pos","prop_neg","prop_pos","p","q","genes")))
  for(i in 1:nrow(GSEA_table)){
    # cat(i," ")
    if(isList) genes_temp<-unlist(GSEA_db[[i]])
    else{
      genes_temp<-unique(as.character(GSEA_db[i,-1]))
      genes_temp<-genes_temp[genes_temp!=""]
    }
    genes_overlap<-genes_temp[genes_temp%in%genes_all]
    GSEA_table[i,"n_genes_pw"]<-length(genes_overlap)
    if(length(genes_overlap)<min_genes) next #No use analysing 
    else{
      isRetrieved<- genes_all%in%genes_retrieved  
      inDB=factor(genes_all%in%genes_overlap,levels=c(FALSE,TRUE))
      compare_pw_temp_t<-table(isRetrieved,inDB)
      GSEA_table[i,c("n_genes_pw_neg","n_genes_pw_pos")]<-compare_pw_temp_t[,"TRUE"]
      GSEA_table[i,c("prop_neg","prop_pos")]<-round(100*prop.table(compare_pw_temp_t,1)[,"TRUE"],1)
      GSEA_table[i,"genes"]<-paste(genes_all[isRetrieved&inDB==TRUE],collapse=",")
      GSEA_table[i,"p"]<-fisher.test(compare_pw_temp_t,alternative = "greater")$p.value   
    }
  }
  GSEA_table<-GSEA_table[!is.na(GSEA_table[,"p"])&!duplicated(rownames(GSEA_table)),]
  # GSEA_table[,"q"]<- qvalue(as.numeric(GSEA_table[,"p"]))$qvalues # To use Storey method (library qvalue)
  GSEA_table[,"q"]<-p.adjust(as.numeric(GSEA_table[,"p"]),"fdr")
  GSEA_table<-as.data.frame(GSEA_table[order(as.numeric(GSEA_table[,"p"])),])
  return(GSEA_table)
}