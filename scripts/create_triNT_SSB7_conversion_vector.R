#############################################
# create_triNT_SSB7_conversion_vector.R
#############################################

# Create convert that allows converting triNT to SSB7 class as described by Zapata et al.

# Download SSB7 information from synapse (https://www.synapse.org/; id syn11935058)
############################################################################################################

CpG_id<- readRDS("data/SSB7_id.rds")

# Add 96-class and 6-class substititon types in our format
triNT_from<- gsub("_.*","",CpG_id[,1])
subst_type<- paste0(substr(CpG_id[,1],2,2),">",substr(CpG_id[,1],5,5))
subst_type_convertIdx<- which(subst_type%in%c("A>C","A>G","A>T","G>A","G>C","G>T"))
subst_type[subst_type=="A>C"]<- "T>G"
subst_type[subst_type=="A>G"]<- "T>C"
subst_type[subst_type=="A>T"]<- "T>A"
subst_type[subst_type=="G>A"]<- "C>T"
subst_type[subst_type=="G>C"]<- "C>G"
subst_type[subst_type=="G>T"]<- "C>A"
triNT_from[subst_type_convertIdx]<- reverseComplement(DNAStringSet(triNT_from[subst_type_convertIdx]))
subst_type3<- paste0(triNT_from,">",substr(triNT_from,1,1),substr(subst_type,3,3),substr(triNT_from,3,3))
CpG_id<- cbind(CpG_id,subst_type,subst_type3)
SSB7<- subst_type
SSB7[CpG_id[,2]=="cpg"]<- "cpg"
CpG_id<- cbind(CpG_id,SSB7)
CpG_conversion<- CpG_id[,c(3,4,5)]
CpG_conversion<- CpG_conversion[which(!duplicated(CpG_conversion[,"subst_type3"])),]
CpG_conversion_vector<- as.character(CpG_conversion[,"SSB7"])
names(CpG_conversion_vector)<- CpG_conversion[,"subst_type3"]
saveRDS(CpG_conversion_vector,file="data/SSB7_conversion_vector.rds")

