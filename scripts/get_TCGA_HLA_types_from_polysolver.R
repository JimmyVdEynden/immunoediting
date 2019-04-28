########################################################
# get_TCGA_HLA_types_from_polysolver.R
########################################################

# See run_polysolver_TCGA.sh

# Load functions
#################
source("scripts/functions/get_HLA_types_from_polysolver_output.R")

# Get HLA types
###############
TCGA_HLA_types<- get_HLA_types_from_polysolver_output("temp/polysolver/TCGA_HLA_types/")
TCGA_HLA_types_netMHC<- get_HLA_types_from_polysolver_output("temp/polysolver/TCGA_HLA_types/",convert_HLAtypes_to_netMHC_format = T)

# Some HLA types were default, because of abscent reads! Exclude those (19 in total)
polysolver_folders<- list.files("temp/polysolver/output",full.names = T)
sample_to_exclude<- NULL
sample_to_exclude2<- NULL
for(i in 1:length(polysolver_folders)){
  cat(i, " ")
  if(i%in%c(4787,7595,2289)) next # Some problem, check manually, ok!
  folder_temp<- polysolver_folders[i]
  if(sum(read.table(paste0(folder_temp,"/counts1.R0k6"),sep = "\t")[,2],na.rm=T)==0) sample_to_exclude<- append(sample_to_exclude,gsub(".*\\/","",folder_temp))
  if(sum(read.table(paste0(folder_temp,"/counts2.R0k6"),sep = "\t")[,2],na.rm=T)==0) sample_to_exclude2<- append(sample_to_exclude2,gsub(".*\\/","",folder_temp))
}
identical(sample_to_exclude,sample_to_exclude2) # TRUE
TCGA_HLA_types<- TCGA_HLA_types[!rownames(TCGA_HLA_types)%in%substr(sample_to_exclude,1,12),]
TCGA_HLA_types_netMHC<- TCGA_HLA_types_netMHC[!rownames(TCGA_HLA_types_netMHC)%in%substr(sample_to_exclude,1,12),]

# Save
save(TCGA_HLA_types,TCGA_HLA_types_netMHC,file="data/TCGA_HLA_types.RData")

# Put different HLA alleles in 1 txt file to be used by netMHCPan
##################################################################

# Spread over multiple lines to avoid netMHCPan problems and to parallelize scripts
HLA_types<- sort(unique(as.character(TCGA_HLA_types_netMHC)))
HLA_types<- setdiff(HLA_types,"HLA-B13:07") # HLA-B13:07 not in netMHC database!
max_per_line<- round(length(HLA_types)/10)
write.table(
  c(
  paste(HLA_types[1:max_per_line],collapse=","),
  paste(HLA_types[(max_per_line+1):(2*max_per_line)],collapse=","),
  paste(HLA_types[(2*max_per_line+1):(3*max_per_line)],collapse=","),
  paste(HLA_types[(3*max_per_line+1):(4*max_per_line)],collapse=","),
  paste(HLA_types[(4*max_per_line+1):(5*max_per_line)],collapse=","),
  paste(HLA_types[(5*max_per_line+1):(6*max_per_line)],collapse=","),
  paste(HLA_types[(6*max_per_line+1):(7*max_per_line)],collapse=","),
  paste(HLA_types[(7*max_per_line+1):(8*max_per_line)],collapse=","),
  paste(HLA_types[(8*max_per_line+1):(9*max_per_line)],collapse=","),
  paste(HLA_types[(9*max_per_line+1):length(HLA_types)],collapse=",")
  ),
  "temp/HLA_types_all.txt",sep=",",quote=FALSE,row.names = F,col.names = F
  )
