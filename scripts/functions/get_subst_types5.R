get_subst_types5<- function(chr,pos,alt_allele){
  
  # Libraries
  require(BSgenome.Hsapiens.UCSC.hg19)
  
  # Get PentaNT ref information
  chr<- gsub("chr","",chr) # Remove in case format is alread "chr"
  chr<- paste0("chr",chr)
  mutPentaNT_Nl<- getSeq(Hsapiens,chr,as.numeric(pos)-2,as.numeric(pos)+2)
  
  # Get PentaNT var information
  mutPentaNT_T<- mutPentaNT_Nl
  mutNT_T<- alt_allele
  mutNT_T[mutNT_T=="0"]<- "" #Zero not recognized by DNAString
  mutNT_T[grepl(",",mutNT_T)]<- "" # Few are undefined, remove
  mutNT_T<- DNAStringSet(mutNT_T)
  mutNT_T[width(mutNT_T)!=1]<- ""
  subseq(mutPentaNT_T,3,3)<- mutNT_T
  
  # Take pyrimidine as ref
  idx_compl<-which(substr(mutPentaNT_Nl,3,3)%in%c("A","G"))
  mutPentaNT_Nl[idx_compl]<-reverseComplement(mutPentaNT_Nl[idx_compl])
  mutPentaNT_T[idx_compl]<-reverseComplement(mutPentaNT_T[idx_compl])
  mutPentaNT_Nl[width(mutPentaNT_Nl)!=5]<- DNAString("N")
  mutPentaNT_T[width(mutPentaNT_T)!=5]<- DNAString("N")
  
  # Get PentaNT substitution type
  subst_type5<- paste0(mutPentaNT_Nl,">",mutPentaNT_T)

  # Add information to maf
  return(subst_type5)

    
}
