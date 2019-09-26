#####################################
# manuscript_Rsim_concept.R
#####################################

# Load data
###########

# Expected numbers for prototypical type for 10000 random mut per subst type
GPPM_subset_subst_type_matrix_ls <- readRDS("data/GPPM_subset_subst_type_matrix_ls.rds")

# Mutation data
TCGA_maf <- readRDS("data/TCGA_maf.rds")

# Functions
############
source("scripts/functions/context3d.R")
library(svglite)

# Plot aggregated mutational signature
#######################################

# Select cancer as example
cancer<- "SKCM"
TCGA_maf_tmp<- TCGA_maf[TCGA_maf$Cancer==cancer,] 

# Get ms
ms_tmp<- prop.table(table(TCGA_maf_tmp$subst_type3))

# Plot
plot_tmp<- context3d(ms_tmp,scalexy = 1000,scalez = 50000,alpha=0.4)

# Save figure
rgl.postscript("results/figs/pdf/fig5_concept_ms.pdf",fmt="pdf")
rgl.postscript("results/figs/svg/fig5_concept_ms.svg",fmt="svg")

# Plot n+/n- per ms
###################

# Take prototypical type as an example
exp_per_subst<- GPPM_subset_subst_type_matrix_ls$triNT_mut[,"nonsynonymous SNV","TRUE"]/GPPM_subset_subst_type_matrix_ls$triNT_mut[,"nonsynonymous SNV","FALSE"]
plot_tmp<- context3d(exp_per_subst,scalexy = 1000,plot_heatmap = T)

# Save figure
rgl.postscript("results/figs/pdf/fig5_concept_exp.pdf",fmt="pdf")
rgl.postscript("results/figs/svg/fig5_concept_exp.svg",fmt="svg")

# Save scale
svglite("results/figs/svg/fig5_concept_exp_scale.svg")
plot(1:100,rep(1,100),pch="|",col=colorRampPalette(colors = c("white","yellow","red"))(100)[1:100],cex=5,axes=F,xlab=NA,ylab=NA)
points(50,1,pch="|")
text(50,0.5,round(mean(exp_per_subst),2),xpd=NA)
points(100,1,pch="|")
text(100,0.5,2*round(mean(exp_per_subst),2),xpd=NA)
points(0,1,pch="|")
text(0,0.5,0.5*round(mean(exp_per_subst),2),xpd=NA)
dev.off()
