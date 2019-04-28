Van den Eynden J, Jiménez-Sánchez A, Miller ML and Larsson E. Lack of detectable neoantigen depletion signals in the untreated cancer genome. 2019.
############################################################################################################

Source code used to produce the results as reported in the manuscript.

The entire pipeline is described in the R markdown notebook "immunoediting.Rmd" (and html "immunoediting.nb.html", which also contains some core results). 

After downloading data available at ... data/ folder, running section 3 ("Analysis") will produce all figures, tables, ... as reported in the manuscripts. We recommend not to run section 2 ("Data processing") as such due to long processing times.    

Note that we are not allowed to make protected TCGA data or any downstream data available. Concretely this implies that the files data/TCGA_maf.rds & data/TCGA_HLA_types.RData are not provided and any script depending in this file will return an error.

To use the HLA-binding annotation and/or metrics reported in the manuscript on other mutation data, follow instructions and example file provided in folder "running_analysis_on_other_mutation_data/" 

sessionInfo()

R version 3.4.4 (2018-03-15)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS release 6.9 (Final)

Matrix products: default
BLAS: /usr/lib64/R/lib/libRblas.so
LAPACK: /usr/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grDevices datasets  stats4    parallel  stats     graphics  utils     methods   base     

other attached packages:
 [1] exact2x2_1.6.3                          exactci_1.3-3                          
 [3] ssanv_1.1                               WriteXLS_4.0.0                         
 [5] RColorBrewer_1.1-2                      corrplot_0.84                          
 [7] plotrix_3.7-2                           gdtools_0.1.7                          
 [9] svglite_1.2.1                           TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
[11] VariantAnnotation_1.22.3                Rsamtools_1.28.0                       
[13] SummarizedExperiment_1.6.4              DelayedArray_0.2.7                     
[15] matrixStats_0.52.2                      BSgenome.Hsapiens.UCSC.hg19_1.4.0      
[17] BSgenome_1.44.1                         rtracklayer_1.36.4                     
[19] Biostrings_2.44.2                       XVector_0.16.0                         
[21] EnsDb.Hsapiens.v75_2.1.0                ensembldb_2.0.4                        
[23] AnnotationFilter_1.0.0                  GenomicFeatures_1.28.5                 
[25] AnnotationDbi_1.38.2                    Biobase_2.36.2                         
[27] GenomicRanges_1.28.5                    GenomeInfoDb_1.12.2                    
[29] IRanges_2.10.3                          S4Vectors_0.14.4                       
[31] BiocGenerics_0.22.1                    

loaded via a namespace (and not attached):
 [1] lattice_0.20-35               htmltools_0.3.6               yaml_2.1.18                  
 [4] interactiveDisplayBase_1.14.0 blob_1.1.0                    XML_3.98-1.9                 
 [7] rlang_0.1.2                   DBI_0.7                       BiocParallel_1.10.1          
[10] bit64_0.9-7                   GenomeInfoDbData_0.99.0       zlibbioc_1.22.0              
[13] ProtGenerics_1.8.0            memoise_1.1.0                 biomaRt_2.32.1               
[16] httpuv_1.3.5                  BiocInstaller_1.26.1          curl_2.8.1                   
[19] Rcpp_0.12.16                  xtable_1.8-2                  mime_0.5                     
[22] bit_1.1-12                    AnnotationHub_2.8.2           packrat_0.4.9-2              
[25] digest_0.6.15                 shiny_1.0.5                   grid_3.4.4                   
[28] tools_3.4.4                   bitops_1.0-6                  lazyeval_0.2.0               
[31] RCurl_1.95-4.8                tibble_1.3.4                  RSQLite_2.0                  
[34] pkgconfig_2.0.1               Matrix_1.2-11                 httr_1.3.1                   
[37] R6_2.2.2                      GenomicAlignments_1.12.2      compiler_3.4.4  