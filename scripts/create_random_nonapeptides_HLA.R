###################################
# create_random_nonapeptides_HLA.R
###################################

# Create set of 10e6 random nonapeptides and calculate HLA affinities for canonical HLA type

# Define aa classes
aa_hp<- sort(c("G","A","P","V","L","I","M","W","F"))
aa_polar<- sort(c("S","T","Y","N","Q","C"))
aa_charged<- sort(c("K","R","H","D","E"))
aa_all<- c(aa_hp,aa_polar,aa_charged)

# Create random nonapeptides
nonapep_random<- rep(NA,1000000)
for(i in 1:1000000){
  cat(i," ")
  nonapep_random[i]<- paste0(sample(aa_all,9,replace = T),collapse = "")
}
write.table(nonapep_random,file = "temp/nonapep_random.pep",quote=F,row.names = F,col.names = F)

# Run NetMHC Pan (10x)
system("scripts/other/runNetMHCPan_nonapep_random.sh") # Ignore warnings,works!

# Create matrix of amino acids in peptides
pep_matrix<- matrix(0,length(nonapep_random),length(aa_all),dimnames = list(nonapep_random,aa_all))
for(i in 1:nrow(pep_matrix)){
  cat(i," ")
  pep_matrix[i,]<- table(factor(unlist(strsplit(nonapep_random[i],"")),levels=aa_all))
}

# Create matrix of pep x HLA types from netMHCPan results
for(i in 1:6){
  cat(i, " ")
  gene_mhc_colnames<- setdiff(unlist(strsplit(readLines(paste0("temp/nonapep_random_",i,".xls"),n=1),"\t")),"")
  gene_mhc_nM_temp<- read.table(paste0("temp/nonapep_random_",i,".xls"),header=T,sep="\t",skip=1,colClasses = c(c("NULL","NULL","NULL"),rep(c("NULL","NULL","numeric","NULL"),length(gene_mhc_colnames)),c("NULL","NULL")))
  colnames(gene_mhc_nM_temp)<- gene_mhc_colnames
  if(i==1) gene_mhc_nM<- gene_mhc_nM_temp
  else gene_mhc_nM<- cbind(gene_mhc_nM,gene_mhc_nM_temp)
}
saveRDS(gene_mhc_nM,file = "temp/nonapep_random.rds")

# Calculate harmonic mean
source("scripts/functions/harmonic_mean.R")
nonapep_random_HLA_mean_aff<- apply(gene_mhc_nM,1,harmonic_mean)
saveRDS(nonapep_random_HLA_mean_aff,file = "temp/nonapep_random_HLA_aff.rds")

# Save data
nonapep_random_HLA_aff<- gene_mhc_nM
save(nonapep_random_HLA_aff,nonapep_random_HLA_mean_aff,pep_matrix,file="data/nonapep_random.RData")
