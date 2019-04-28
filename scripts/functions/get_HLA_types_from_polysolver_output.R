get_HLA_types_from_polysolver_output<- function(polysolver_output_folder,convert_HLAtypes_to_netMHC_format=FALSE){
  filenames<- list.files(polysolver_output_folder)
  TCGA_HLA_types<- matrix(NA,length(filenames),6,dimnames = list(substr(filenames,1,12),c("HLA-A1","HLA-A2","HLA-B1","HLA-B2","HLA-C1","HLA-C2")))
  for(filename in filenames){
    if(length(readLines(paste(polysolver_output_folder,filename,sep="")))==0) next
    barcode<- substr(filename,1,12)
    TCGA_HLA_type<- read.table(paste(polysolver_output_folder,filename,sep=""),stringsAsFactors = F)
    TCGA_HLA_types[barcode,"HLA-A1"]<- TCGA_HLA_type[1,2]
    TCGA_HLA_types[barcode,"HLA-A2"]<- TCGA_HLA_type[1,3]
    TCGA_HLA_types[barcode,"HLA-B1"]<- TCGA_HLA_type[2,2]
    TCGA_HLA_types[barcode,"HLA-B2"]<- TCGA_HLA_type[2,3]
    TCGA_HLA_types[barcode,"HLA-C1"]<- TCGA_HLA_type[3,2]
    TCGA_HLA_types[barcode,"HLA-C2"]<- TCGA_HLA_type[3,3]
  }
  # Convert HLA types to netMHCPan format (e.g. HLA-A33:03)
  if(convert_HLAtypes_to_netMHC_format){
    TCGA_HLA_types<- toupper(substr(TCGA_HLA_types,1,11))
    TCGA_HLA_types<- gsub('^(.{8}).', '\\1:', TCGA_HLA_types)
    TCGA_HLA_types<- gsub("\\_","",TCGA_HLA_types)
    TCGA_HLA_types<- gsub("HLA", "HLA-", TCGA_HLA_types)
  }
  return(TCGA_HLA_types)
}

