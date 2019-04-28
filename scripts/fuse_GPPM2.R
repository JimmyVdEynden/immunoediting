# Called from other/fuse_GPPM.sh

# Select files to fuse
GPPM_files<- list.files("temp/",full.names = T)
GPPM_files<- GPPM_files[grepl("GPPM_all",GPPM_files)]

# Fuse
pb <- txtProgressBar(min = 0, max = length(GPPM_files), style = 3)
for(i in 1:length(GPPM_files)){
  # cat(i, " ")
  setTxtProgressBar(pb, i)
  GPPM<- readRDS(GPPM_files[i])
  if(i==1) GPPM_all<- GPPM
  else GPPM_all<- c(GPPM_all,GPPM)
}
saveRDS(GPPM_all,file = "temp/GPPM_all.rds")
