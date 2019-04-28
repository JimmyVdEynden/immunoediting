create_GPosProt_from_gene<- function(gene,edb){

  # Libraries
#   library(ensembldb)
#   library(BSgenome.Hsapiens.UCSC.hg19)

  # Get all uniprot ids
  prts_id<- genes(edb,filter=GenenameFilter(gene),columns = listColumns(edb,"uniprot"))
  prts_id<- prts_id[!is.na(prts_id$uniprot_id)]
  prts_id<- sort(unique(prts_id$protein_id))
  
  # Get aa sequences + width
  prts_aa <- proteins(edb, filter = ProteinIdFilter(prts_id),
                      return.type = "AAStringSet")
  prts_width<- width(prts_aa)
  names(prts_width)<- prts_id
  
  # Get CDS gene coordinates
  cds_gene <- cdsBy(edb, by="gene",filter = ProteinIdFilter(prts_id),
                    columns = "protein_id")
  cds_gene_gpos <- GPos(cds_gene[[1]]) # Turn into GPos object
  seqlevels(cds_gene_gpos)<- paste("chr",seqlevels(cds_gene_gpos),sep="")
  
  # Get CDS transcript coordinates
  cds <- cdsBy(edb, by="tx",filter = ProteinIdFilter(prts_id),
               columns = "protein_id")
  seqlevels(cds)<- paste("chr",seqlevels(cds),sep="")
  
  # Get NT sequences
  seq_temp<- getSeq(Hsapiens,seqnames(cds_gene_gpos),pos(cds_gene_gpos),pos(cds_gene_gpos))
  cds_gene_gpos$nt<- seq_temp
  seq_temp<- getSeq(Hsapiens,seqnames(cds_gene_gpos),pos(cds_gene_gpos)-1,pos(cds_gene_gpos)+1)
  cds_gene_gpos$triNt<- seq_temp
  
  # Calculate aa positions
  aaPos<- NULL
  aaSeq<- NULL
  isNegStrand<- as.character(unique(strand(cds_gene_gpos)))=="-"
  for(tx in names(cds)){
    
    cds_temp<- cds[[tx]]
    prt_id<- unique(cds_temp$protein_id)
    # cat(tx, " ",prt_id,"\n")
    cds_temp_gpos <- GPos(cds_temp)
    
    # Determine AA positions in prot
    isOverlap_temp<- countOverlaps(cds_gene_gpos,cds_temp_gpos)
    aaPos_temp<- isOverlap_temp
    aaCo<- ceiling(seq(1,sum(aaPos_temp))/3)
    if(isNegStrand) aaCo<- rev(aaCo)  
    aaPos_temp[aaPos_temp==1]<- aaCo
    
    # Remove stopcodon: 
    # Sometimes not present in coodinates for some reason
    # e.g. TP53 ENSP00000425104 & ENSP00000473895
    # In this case don't remove
    if(max(aaPos_temp)!=prts_width[prt_id]) aaPos_temp[aaPos_temp==max(aaPos_temp)]<- 0 # Remove stopcodon
    
    # Check 1: similar protein lengths?
    if(max(aaPos_temp)!=prts_width[prt_id]){
      warning(paste(prt_id," ",tx," Check proteins: seqlengths ",max(aaPos_temp)," vs. ",prts_width[prt_id]," REMOVING PROTEIN"))
      next
    }
    # Determine precise AA in protein
    seq_prt<- as.character(prts_aa[prt_id])
    get_aaSeq<- function(aa_idx){
      if(aa_idx==0) aa_temp<- ""
      else aa_temp<- substr(seq_prt,aa_idx,aa_idx)
      return(aa_temp)
    }
    aa_seq_temp<- sapply(aaPos_temp,get_aaSeq)
    
    # Check 2: compare to aaSeq to translated sequence
    # If different don't use protein information (e.g. ENST00000436599 starts with X, gives wrong results)
    nt_seq_prot<- getSeq(Hsapiens,seqnames(cds_gene_gpos),pos(cds_gene_gpos),pos(cds_gene_gpos),strand=strand(cds_gene_gpos))
    nt_seq_prot<- nt_seq_prot[aaPos_temp!=0]
    nt_seq_prot<- DNAString(paste(nt_seq_prot,collapse=""))
    if(isNegStrand) nt_seq_prot<- rev(nt_seq_prot)  
    nt_seq_prot_translated<- as.character(translate(nt_seq_prot))
    aa_seq_temp2<- aa_seq_temp[aa_seq_temp!=""]
    aa_seq_temp2<- aa_seq_temp2[seq(1,length(aa_seq_temp2),3)]
    if(isNegStrand) aa_seq_temp2<- rev(aa_seq_temp2)  
    aa_seq_temp2<- paste(aa_seq_temp2,collapse = "")
    if(aa_seq_temp2!=nt_seq_prot_translated){
      warning(paste(prt_id," ",tx," Check proteins: seq ",nt_seq_prot_translated," vs. ", aa_seq_temp2, "REMOVING PROTEIN"))
      next
    }
    
    # Add protein information
    aaPos<- cbind(aaPos,aaPos_temp)
    colnames(aaPos)[ncol(aaPos)]<- prt_id
    aaSeq<- cbind(aaSeq,aa_seq_temp)
    colnames(aaSeq)[ncol(aaSeq)]<- prt_id
    
  }
  rownames(aaSeq)<- NULL
  cds_gene_gpos$aaPos<- aaPos
  cds_gene_gpos$aa<- aaSeq
  
  #   # Check 3: aa same for each protein?
  #   for(i in 1:nrow(aaSeq)){
  #     if(length(setdiff(unique(aaSeq[i,]),""))>1) cat(i," ")
  #   }
  #   # Not always, check IGV for those
  #   cds_gene_gpos[597:598]
  
  return(cds_gene_gpos)
}

#   create_GPosProt_from_gene<- function(gene,edb){
#     
#     # Get all uniprot ids
#     prts_id<- genes(edb,filter=GenenameFilter(gene),columns = listColumns(edb,"uniprot"))
#     prts_id<- prts_id[!is.na(prts_id$uniprot_id)]
#     prts_id<- sort(unique(prts_id$protein_id))
#     
#     # Get aa sequences + width
#     prts_aa <- proteins(edb, filter = ProteinIdFilter(prts_id),
#                         return.type = "AAStringSet")
#     prts_width<- width(prts_aa)
#     names(prts_width)<- prts_id
#     
#     # Get CDS gene coordinates
#     cds_gene <- cdsBy(edb, by="gene",filter = ProteinIdFilter(prts_id),
#                       columns = "protein_id")
#     cds_gene_gpos <- GPos(cds_gene[[1]]) # Turn into GPos object
#     seqlevels(cds_gene_gpos)<- paste("chr",seqlevels(cds_gene_gpos),sep="")
#     
#     # Get CDS transcript coordinates
#     cds <- cdsBy(edb, by="tx",filter = ProteinIdFilter(prts_id),
#                  columns = "protein_id")
#     seqlevels(cds)<- paste("chr",seqlevels(cds),sep="")
#     
#     # Calculate aa positions
#     aaPos<- NULL
#     aaSeq<- NULL
#     isNegStrand<- as.character(unique(strand(cds_gene_gpos)))=="-"
#     for(tx in names(cds)){
#       
#       cds_temp<- cds[[tx]]
#       prt_id<- unique(cds_temp$protein_id)
#       cat(tx, " ",prt_id,"\n")
#       cds_temp_gpos <- GPos(cds_temp)
#       
#       # Determine AA positions in prot
#       isOverlap_temp<- countOverlaps(cds_gene_gpos,cds_temp_gpos)
#       aaPos_temp<- isOverlap_temp
#       aaCo<- ceiling(seq(1,sum(aaPos_temp))/3)
#       if(isNegStrand) aaCo<- rev(aaCo)  
#       aaPos_temp[aaPos_temp==1]<- aaCo
#       aaPos<- cbind(aaPos,aaPos_temp)
#       colnames(aaPos)[ncol(aaPos)]<- prt_id
#       
#       # Remove stopcodon: 
#       # Sometimes not present in coodinates for some reason
#       # e.g. TP53 ENSP00000425104 & ENSP00000473895
#       # In this case don't remove
#       if(max(aaPos_temp)!=prts_width[prt_id]) aaPos_temp[aaPos_temp==max(aaPos_temp)]<- 0 # Remove stopcodon
#       
#       # Check 1: similar protein lengths?
#       if(max(aaPos_temp)!=prts_width[prt_id]) warning(paste(prt_id," ",tx," Check proteins: seqlengths ",max(aaPos_temp)," vs. ",prts_width[prt_id]))
#       
#       # Determine precise AA in protein
#       seq_prt<- as.character(prts_aa[prt_id])
#       #     aa_seq_temp<- NULL
#       #     for(aa_idx in aaPos_temp){
#       #       cat(aa_idx," ")
#       #       if(aa_idx==0) aa_temp<- ""
#       #       else aa_temp<- substr(seq_prt,aa_idx,aa_idx)
#       #       aa_seq_temp<- c(aa_seq_temp,aa_temp)
#       #     }
#       #     aaSeq<- cbind(aaSeq,aa_seq_temp)
#       #     colnames(aaSeq)[ncol(aaSeq)]<- prt_id
#       
#       get_aaSeq<- function(aa_idx){
#         if(aa_idx==0) aa_temp<- ""
#         else aa_temp<- substr(seq_prt,aa_idx,aa_idx)
#         return(aa_temp)
#       }
#       aa_seq_temp<- sapply(aaPos_temp,get_aaSeq)
#       aaSeq<- cbind(aaSeq,aa_seq_temp)
#       colnames(aaSeq)[ncol(aaSeq)]<- prt_id
#       
#       # Check 2: compare to aaSeq to translated sequence
#       nt_seq_prot<- getSeq(Hsapiens,seqnames(cds_gene_gpos),pos(cds_gene_gpos),pos(cds_gene_gpos),strand=strand(cds_gene_gpos))
#       nt_seq_prot<- nt_seq_prot[aaPos_temp!=0]
#       nt_seq_prot<- DNAString(paste(nt_seq_prot,collapse=""))
#       if(isNegStrand) nt_seq_prot<- rev(nt_seq_prot)  
#       nt_seq_prot_translated<- as.character(translate(nt_seq_prot))
#       aa_seq_temp<- aa_seq_temp[aa_seq_temp!=""]
#       aa_seq_temp<- aa_seq_temp[seq(1,length(aa_seq_temp),3)]
#       if(isNegStrand) aa_seq_temp<- rev(aa_seq_temp)  
#       aa_seq_temp<- paste(aa_seq_temp,collapse = "")
#       cat(substr(aa_seq_temp,1,1)," ")
#       if(aa_seq_temp!=nt_seq_prot_translated) warning(paste(prt_id," ",tx," Check proteins: seq ",nt_seq_prot_translated," vs. ", aa_seq_temp))
#     }
#     rownames(aaSeq)<- NULL
#     cds_gene_gpos$aaPos<- aaPos
#     cds_gene_gpos$aaSeq<- aaSeq
#     
#     #   # Check 3: aa same for each protein?
#     #   for(i in 1:nrow(aaSeq)){
#     #     if(length(setdiff(unique(aaSeq[i,]),""))>1) cat(i," ")
#     #   }
#     #   # Not always, check IGV for those
#     #   cds_gene_gpos[597:598]
#     
#     # return(cds_gene_gpos)
#   }
#   
