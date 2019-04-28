# Called from other/fuse_GPPM.sh

# Adding one by one slows down too much: solution: 1) max 1000 at a time and 2) parallelize shell
# Therefore need to pass command args
args<- commandArgs(trailingOnly = T) # TRailing only to only use args supplied
j<- as.numeric(args[1])

# Select files to fuse
GPPM_files<- list.files("temp/GPPM/",full.names = T)
max_per_run<- round(length(GPPM_files)/20)
if(j!=20) GPPM_files_sel<- GPPM_files[((j-1)*max_per_run+1):(j*max_per_run)]
if(j==20) GPPM_files_sel<- GPPM_files[((j-1)*max_per_run+1):length(GPPM_files)]

# Fuse
pb <- txtProgressBar(min = 0, max = length(GPPM_files_sel), style = 3)
for(i in 1:length(GPPM_files_sel)){
  # cat(i, " ")
  setTxtProgressBar(pb, i)
  GPPM<- readRDS(GPPM_files_sel[i])
  if(i==1) GPPM_all<- GPPM
  else GPPM_all<- c(GPPM_all,GPPM)
}
saveRDS(GPPM_all,file = paste0("temp/GPPM_all_",j,".rds"))
