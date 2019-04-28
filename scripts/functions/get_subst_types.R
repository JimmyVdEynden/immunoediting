get_subst_types<- function(chr,pos,alt_allele){
  
  # Libraries
  require(BSgenome.Hsapiens.UCSC.hg19)
  
  # Get triNT ref information
  chr<- gsub("chr","",chr) # Remove in case format is alread "chr"
  chr<- paste0("chr",chr)
  mutTriNT_Nl<- getSeq(Hsapiens,chr,as.numeric(pos)-1,as.numeric(pos)+1)
  
  # Get triNT var information
  mutTriNT_T<- mutTriNT_Nl
  mutNT_T<- alt_allele
  mutNT_T[mutNT_T=="0"]<- "" #Zero not recognized by DNAString
  mutNT_T[grepl(",",mutNT_T)]<- "" # Few are undefined, remove
  mutNT_T<- DNAStringSet(mutNT_T)
  mutNT_T[width(mutNT_T)!=1]<- ""
  subseq(mutTriNT_T,2,2)<- mutNT_T
  
  # Take pyrimidine as ref
  idx_compl<-which(substr(mutTriNT_Nl,2,2)%in%c("A","G"))
  mutTriNT_Nl[idx_compl]<-reverseComplement(mutTriNT_Nl[idx_compl])
  mutTriNT_T[idx_compl]<-reverseComplement(mutTriNT_T[idx_compl])
  mutTriNT_Nl[width(mutTriNT_Nl)!=3]<- DNAString("N")
  mutTriNT_T[width(mutTriNT_T)!=3]<- DNAString("N")
  
  # Get main substitution type
  subst_type<- paste0(substr(mutTriNT_Nl,2,2),">",substr(mutTriNT_T,2,2))
  
  # Get triNT substitution type
  subst_type3<- paste0(mutTriNT_Nl,">",mutTriNT_T)
  
  # Add information to maf
  return(list(subst_type=subst_type,subst_type3=subst_type3))
  
}
