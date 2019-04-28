get_peptides_from_GPP<- function(GPP,peptide_length,return_ls=TRUE){
  GPP_aaPos<- GPP$aaPos 
  GPP_aaSeq<- GPP$aa
  GPP_pept<- NULL
  for(i in 1:ncol(GPP_aaPos)){
    
    # Get protein information
    GPP_aaPos_temp<- GPP_aaPos[,i]
    GPP_aaSeq_temp<- GPP_aaSeq[,i]
    names(GPP_aaSeq_temp)<- GPP_aaPos_temp
    GPP_aaSeq_temp<- GPP_aaSeq_temp[GPP_aaSeq_temp!=""]
    GPP_aaSeq_temp<- GPP_aaSeq_temp[!duplicated(names(GPP_aaSeq_temp))]
    GPP_aaSeq_temp<- GPP_aaSeq_temp[order(as.numeric(names(GPP_aaSeq_temp)),decreasing = F)]
    
    # Get peptide information
    if(return_ls){
      GPP_pept_temp_ls<- list() 
      for(a in 1:length(GPP_aaSeq_temp)){
        GPP_pept_temp<- rep(NA,peptide_length)
        for(b in 1:peptide_length){
          if((a-b+1)<=0) next # Not before first aa
          if((a-b+peptide_length)>length(GPP_aaSeq_temp)) next # Not after last aa
          GPP_pept_temp[b]<- paste(GPP_aaSeq_temp[(a-b+1):((a-b+1)+peptide_length-1)],collapse="")
        }
        GPP_pept_temp_ls<- c(GPP_pept_temp_ls,list(GPP_pept_temp))
      }
      names(GPP_pept_temp_ls)<- names(GPP_aaSeq_temp)
      GPP_pept_temp_ls<- c(GPP_pept_temp_ls,"0"=list(rep(NA,9)))
      GPP_pept_temp_ls<- GPP_pept_temp_ls[as.character(GPP_aaPos_temp)]
      GPP_pept<- cbind(GPP_pept,GPP_pept_temp_ls)
    }
    else{
      GPP_pept_temp<- GPP_aaSeq_temp
      GPP_pept_temp[]<- NA 
      for(a in 1:(length(GPP_pept_temp)-peptide_length+1)){
        GPP_pept_temp[a]<- paste(GPP_aaSeq_temp[a:(a+peptide_length-1)],collapse="")
      }
      GPP_pept_temp<- c(GPP_pept_temp,"0"=NA)
      GPP_pept_temp<- GPP_pept_temp[as.character(GPP_aaPos_temp)]
      GPP_pept<- cbind(GPP_pept,GPP_pept_temp)
    }
  }
  colnames(GPP_pept)<- colnames(GPP_aaSeq)
  
  return(GPP_pept)
}
