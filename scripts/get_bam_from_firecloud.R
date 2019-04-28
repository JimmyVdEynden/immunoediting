###################################################################################################
# get_bam_from_firecloud: script to extract and download specific region from bam file in Firecloud
###################################################################################################

# 1.1 Download sample metadata for specific cancer from firecloud & put file in downloads/firecloud/metadata_downloads/

# 1.2 Create own metadata based on the downloaded ones
{
  sample_files<- list.files("downloads/firecloud/metadata_downloads/",full.names = T)
  
  metadata_participant<- NULL
  metadata_sample<- NULL
  metadata_sampleSet<- NULL
  
  for(sample_file in sample_files){
    
    # Read file
    metadata_temp<- read.table(sample_file,header=TRUE,sep="\t",stringsAsFactors = F)
    cancer<- gsub(".*sample_|\\.txt","",sample_file)
    
    # Only select metadata with gs links: 335, exactly the same as using TCGAbiolinks selection (see download_gdc_bam_10.R in arch/)
    if(cancer=="LAML") metadata_selected_temp<- metadata_temp[metadata_temp[,"sample_type"]=="NT",] # NT form LAML
    else metadata_selected_temp<- metadata_temp[metadata_temp[,"sample_type"]=="NB",] # NT form LAML
    metadata_selected_temp<- metadata_selected_temp[!is.na(metadata_selected_temp[,"WXS_bai_path"])&metadata_selected_temp[,"WXS_bai_path"]!="",]
    
    # Extract columns to keep
    sample_id<- metadata_selected_temp[,1]
    participant<- substr(metadata_selected_temp[,"tcga_sample_id"],1,12)
    BAM<- metadata_selected_temp[,"WXS_bam_path"]
    BAM_IDX<- metadata_selected_temp[,"WXS_bai_path"]
    tcga_sample_id<- metadata_selected_temp[,"tcga_sample_id"]
    metadata_selected_temp<- cbind("entity:sample_id"=sample_id,BAM,BAM_IDX,participant,tcga_sample_id)
    
    # Create metadata files
    metadata_participant_temp<- metadata_selected_temp[,"participant",drop=FALSE]
    colnames(metadata_participant_temp)<- "entity:participant_id"
    metadata_sample_temp<- metadata_selected_temp
    metadata_sampleSet_temp<- cbind("membership:sample_set_id"=cancer,sample=metadata_selected_temp[,"entity:sample_id"])
    
    # Check metadata
    if(unique(metadata_sampleSet_temp[,1])!=unique(gsub("-.*","",metadata_sampleSet_temp[,2]))) stop("Check cancer ",cancer)
    
    # Add to other metadata
    metadata_participant<- rbind(metadata_participant,metadata_participant_temp) 
    metadata_sample<- rbind(metadata_sample,metadata_sample_temp) 
    metadata_sampleSet<- rbind(metadata_sampleSet,metadata_sampleSet_temp) 
  }
  # Save
  write.table(metadata_participant,"downloads/firecloud/metadata/metadata_participant.txt",sep="\t",row.names = F,col.names = T,quote=F)
  write.table(metadata_sample,"downloads/firecloud/metadata/metadata_sample.txt",sep="\t",row.names = F,col.names = T,quote=F)
  write.table(metadata_sampleSet,"downloads/firecloud/metadata/metadata_sampleSet.txt",sep="\t",row.names = F,col.names = T,quote=F)
}

# 1.3 Reduce selected bams to HLA regions on chr6 on firecloud 
{
  # 1.3.1 Upload metadata: 
  # 1) metadata_participant.txt
  # 2) metadata_sample.txt
  # 3) metadata_sampleSet.txt

  # 1.3.2 Upload script bamSplit_HLA.wdl in method repository

  # 1.3.3 Method configuration
  # bamSplit_HLA.wdl
  # root entity: sample
  # bamSplit_HLA.splitBam_HLA.sampleName: (String): this.tcga_sample_id
  # bamSplit_HLA.splitBam_HLA.inputBAM: (File): this.BAM
  # bamSplit_HLA.splitBam_HLA.inputBAI: (File): this.BAM_IDX
  # bamSplit_HLA.splitBam_HLA.splittedBAM: (File):this.splittedBAM
  # bamSplit_HLA.splitBam_HLA.splittedBAM_IDX: (File):this.splittedBAM_IDX

  # 1.3.4 Run all cancer samples: sampleset: e.g. COAD & this.samples
}

# 1.4 Download bucket & put bamfiles in downloads/firecloud/bamSplit_HLA/bam/
{
  # run scripts/other/download_bam_from_firecloud.sh
}
