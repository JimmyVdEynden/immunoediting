######################################
# get_gene_NS_sites.R
######################################

# Aim: Get N & S sites for all genes & subst types

# Load
GPPM<- readRDS("data/GPPM_inclHLAAlleles.rds")
load("temp/CpG_conversion.RData")
GPPM$SSB7<- CpG_conversion_vector[GPPM$subst_type3]

# Get sites
NS_t<- list()
NS_t$main<- table(GPPM$gene,GPPM$subst_type,GPPM$variant)
NS_t$SSB7<- table(GPPM$gene,GPPM$SSB7,GPPM$variant)
NS_t$triNT<- table(GPPM$gene,GPPM$subst_type3,GPPM$variant)
saveRDS(NS_t,file = "data/gene_NS_sites.rds")


