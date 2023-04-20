################################################################
# This script runs an over representation analysis
# for deferentially methylated gene sets against several
# gene set databases. The results and network figures are
# automatically saved to a folder titled geneset_overrep_analysis.
# This script must be run after the main script.

# Author: Parker Grosjean
# March 2021
################################################################

library(WGCNA)
library(minfi)
library(tidyverse)
library(missMethyl)
library(org.Hs.eg.db)
library(data.table)
library(ggplot2)
library(ggridges)
library(methylGSA)
library(clusterProfiler)
library(MeSH.Hsa.eg.db)
library(ReactomePA)
library(meshes)
library(msigdbr)
library(DOSE)
library(enrichplot)
library(ggupset)
library(ggnewscale)

####
####
####
#### Defining functions for testing over representation
####
####
####

go_over_representation <- function(probe_list, total_probe_list, genomic_features, annotation){
  # Generating map between refseq IDs and Entrez IDs
  map <- org.Hs.egREFSEQ2EG
  mapped_genes <- mappedkeys(map)
  gene_map <- as.list(map[mapped_genes])
  dname.info <- annotation[annotation$Name %in% probe_list, ]
  # filtering for genomic features of interest
  if (genomic_features != "ALL"){
    dname.info$UCSC_RefGene_Group <- Map(strsplit, dname.info$UCSC_RefGene_Group, ';') %>% Map(unlist, .)
    dname.info$to_keep <- Map(intersect, dname.info$UCSC_RefGene_Group, genomic_features) %>% Map(length, .)
    dname.info <- dname.info[dname.info$to_keep > 0, ]
  }
  # finding refseq ids of interest
  refseq.list <- dname.info$UCSC_RefGene_Accession[dname.info$UCSC_RefGene_Accession != ""]
  refseq.list <- strsplit(refseq.list, ';')
  refseq.list <- unlist(unique(unlist(refseq.list)))
  # converting RefSeq IDs to EntrezIDs
  gene.list <- gene_map[refseq.list]
  # Getting full universe of genes from annotation
  all.cpg.info <-  annotation[annotation$Name %in% total_probe_list, ]
  gene.universe <- all.cpg.info$UCSC_RefGene_Accession[all.cpg.info$UCSC_RefGene_Accession != ""]
  gene.universe <- strsplit(gene.universe, ';')
  gene.universe <- unlist(unique(unlist(gene.universe)))
  # converting RefSeq IDs to EntrezIDs
  gene.universe <- gene_map[gene.universe]
  # Looking for enriched GO terms
  ego <- enrichGO(gene          = as.character(gene.list),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.1,
                  qvalueCutoff  = 0.1,
                  universe      = as.character(gene.universe),
                  readable      = TRUE)
  # saving result
  item_plot = NULL
  try(item_plot <- enrichplot::cnetplot(ego))
  pre_result <- ego@result
  pre_result <- pre_result[pre_result$p.adjust < 0.1, ]
  result <- pre_result[order(pre_result$p.adjust), ]
  return(list(item_plot, result))
}


kegg_over_representation <- function(probe_list, total_probe_list, genomic_features, annotation){
  # Generating map between refseq IDs and Entrez IDs
  map <- org.Hs.egREFSEQ2EG
  mapped_genes <- mappedkeys(map)
  gene_map <- as.list(map[mapped_genes])
  dname.info <- annotation[annotation$Name %in% probe_list, ]
  # filtering for genomic features of interest
  if (genomic_features != "ALL"){
    dname.info$UCSC_RefGene_Group <- Map(strsplit, dname.info$UCSC_RefGene_Group, ';') %>% Map(unlist, .)
    dname.info$to_keep <- Map(intersect, dname.info$UCSC_RefGene_Group, genomic_features) %>% Map(length, .)
    dname.info <- dname.info[dname.info$to_keep > 0, ]
  }
  # finding refseq ids of interest
  refseq.list <- dname.info$UCSC_RefGene_Accession[dname.info$UCSC_RefGene_Accession != ""]
  refseq.list <- strsplit(refseq.list, ';')
  refseq.list <- unlist(unique(unlist(refseq.list)))
  # converting RefSeq IDs to EntrezIDs
  gene.list <- gene_map[refseq.list]
  # Getting full universe of genes from annotation
  all.cpg.info <-  annotation[annotation$Name %in% total_probe_list, ]
  gene.universe <- all.cpg.info$UCSC_RefGene_Accession[all.cpg.info$UCSC_RefGene_Accession != ""]
  gene.universe <- strsplit(gene.universe, ';')
  gene.universe <- unlist(unique(unlist(gene.universe)))
  # converting RefSeq IDs to EntrezIDs
  gene.universe <- gene_map[gene.universe]
  # Looking for enriched GO terms
  ekegg <- enrichKEGG(gene          = as.character(gene.list),
                      organism      = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.1,
                      qvalueCutoff  = 0.1,
                      universe      = as.character(gene.universe))
  # saving result
  ekegg <- setReadable(ekegg, 'org.Hs.eg.db', 'ENTREZID')
  item_plot <- NULL
  try(item_plot <- enrichplot::cnetplot(ekegg, colorEdge = TRUE, cex_label_category = 1.5))
  pre_result <- ekegg@result
  pre_result <- pre_result[pre_result$p.adjust < 0.1, ]
  result <- pre_result[order(pre_result$p.adjust), ]
  return(list(item_plot, result))
}


reactome_over_representation <- function(probe_list, total_probe_list, genomic_features, annotation){
  # Generating map between refseq IDs and Entrez IDs
  map <- org.Hs.egREFSEQ2EG
  mapped_genes <- mappedkeys(map)
  gene_map <- as.list(map[mapped_genes])
  dname.info <- annotation[annotation$Name %in% probe_list, ]
  # filtering for genomic features of interest
  if (genomic_features != "ALL"){
    dname.info$UCSC_RefGene_Group <- Map(strsplit, dname.info$UCSC_RefGene_Group, ';') %>% Map(unlist, .)
    dname.info$to_keep <- Map(intersect, dname.info$UCSC_RefGene_Group, genomic_features) %>% Map(length, .)
    dname.info <- dname.info[dname.info$to_keep > 0, ]
  }
  # finding refseq ids of interest
  refseq.list <- dname.info$UCSC_RefGene_Accession[dname.info$UCSC_RefGene_Accession != ""]
  refseq.list <- strsplit(refseq.list, ';')
  refseq.list <- unlist(unique(unlist(refseq.list)))
  # converting RefSeq IDs to EntrezIDs
  gene.list <- gene_map[refseq.list]
  # Getting full universe of genes from annotation
  all.cpg.info <-  annotation[annotation$Name %in% total_probe_list, ]
  gene.universe <- all.cpg.info$UCSC_RefGene_Accession[all.cpg.info$UCSC_RefGene_Accession != ""]
  gene.universe <- strsplit(gene.universe, ';')
  gene.universe <- unlist(unique(unlist(gene.universe)))
  # converting RefSeq IDs to EntrezIDs
  gene.universe <- gene_map[gene.universe]
  # Looking for enriched GO terms
  ereact <- enrichPathway(gene          = as.character(gene.list),
                      organism      = "human",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.1,
                      qvalueCutoff  = 0.1,
                      universe      = as.character(gene.universe),
                      readable = TRUE)
  # saving result
  ereact <- setReadable(ereact, 'org.Hs.eg.db', 'ENTREZID')
  item_plot = NULL
  try(item_plot <- enrichplot::cnetplot(ereact, colorEdge = TRUE, cex_label_category = 1.5))
  pre_result <- ereact@result
  pre_result <- pre_result[pre_result$p.adjust < 0.1, ]
  result <- pre_result[order(pre_result$p.adjust), ]
  return(list(item_plot, result))
}


mesh_over_representation <- function(probe_list, total_probe_list, genomic_features, annotation){
  # Generating map between refseq IDs and Entrez IDs
  map <- org.Hs.egREFSEQ2EG
  mapped_genes <- mappedkeys(map)
  gene_map <- as.list(map[mapped_genes])
  dname.info <- annotation[annotation$Name %in% probe_list, ]
  # filtering for genomic features of interest
  if (genomic_features != "ALL"){
    dname.info$UCSC_RefGene_Group <- Map(strsplit, dname.info$UCSC_RefGene_Group, ';') %>% Map(unlist, .)
    dname.info$to_keep <- Map(intersect, dname.info$UCSC_RefGene_Group, genomic_features) %>% Map(length, .)
    dname.info <- dname.info[dname.info$to_keep > 0, ]
  }
  # finding refseq ids of interest
  refseq.list <- dname.info$UCSC_RefGene_Accession[dname.info$UCSC_RefGene_Accession != ""]
  refseq.list <- strsplit(refseq.list, ';')
  refseq.list <- unlist(unique(unlist(refseq.list)))
  # converting RefSeq IDs to EntrezIDs
  gene.list <- gene_map[refseq.list]
  # Getting full universe of genes from annotation
  all.cpg.info <-  annotation[annotation$Name %in% total_probe_list, ]
  gene.universe <- all.cpg.info$UCSC_RefGene_Accession[all.cpg.info$UCSC_RefGene_Accession != ""]
  gene.universe <- strsplit(gene.universe, ';')
  gene.universe <- unlist(unique(unlist(gene.universe)))
  # converting RefSeq IDs to EntrezIDs
  gene.universe <- gene_map[gene.universe]
  # Looking for enriched GO terms
  emesh <- enrichMeSH(gene          = as.character(gene.list),
                      MeSHDb      = "MeSH.Hsa.eg.db",
                      database   = "gendoo", 
                      category    = "C",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.1,
                      qvalueCutoff  = 0.1,
                      universe      = as.character(gene.universe))
  # saving result
  item_plot = NULL
  try(item_plot <- enrichplot::cnetplot(emesh))
  pre_result <- emesh@result
  pre_result <- pre_result[pre_result$p.adjust < 0.1, ]
  result <- pre_result[order(pre_result$p.adjust), ]
  return(list(item_plot, result))
}


msigdb_over_representation <- function(probe_list, total_probe_list, genomic_features, annotation, all_gene_sets, c_number){
  # Generating map between refseq IDs and Entrez IDs
  map <- org.Hs.egREFSEQ2EG
  mapped_genes <- mappedkeys(map)
  gene_map <- as.list(map[mapped_genes])
  dname.info <- annotation[annotation$Name %in% probe_list, ]
  # filtering for genomic features of interest
  if (genomic_features != "ALL"){
    dname.info$UCSC_RefGene_Group <- Map(strsplit, dname.info$UCSC_RefGene_Group, ';') %>% Map(unlist, .)
    dname.info$to_keep <- Map(intersect, dname.info$UCSC_RefGene_Group, genomic_features) %>% Map(length, .)
    dname.info <- dname.info[dname.info$to_keep > 0, ]
  }
  # finding refseq ids of interest
  refseq.list <- dname.info$UCSC_RefGene_Accession[dname.info$UCSC_RefGene_Accession != ""]
  refseq.list <- strsplit(refseq.list, ';')
  refseq.list <- unlist(unique(unlist(refseq.list)))
  # converting RefSeq IDs to EntrezIDs
  gene.list <- gene_map[refseq.list]
  # Getting full universe of genes from annotation
  all.cpg.info <-  annotation[annotation$Name %in% total_probe_list, ]
  gene.universe <- all.cpg.info$UCSC_RefGene_Accession[all.cpg.info$UCSC_RefGene_Accession != ""]
  gene.universe <- strsplit(gene.universe, ';')
  gene.universe <- unlist(unique(unlist(gene.universe)))
  # converting RefSeq IDs to EntrezIDs
  gene.universe <- gene_map[gene.universe]
  # Looking for enriched GO terms
  msigdbr_df <- all_gene_sets[all_gene_sets$gs_cat == c_number, ]
  msigdbr_t2g = msigdbr_df %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
  emsigdb <- enricher(gene = as.character(gene.list), 
                      TERM2GENE = msigdbr_t2g,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.1,
                      qvalueCutoff  = 0.1,
                      universe      = as.character(gene.universe))
  
  # saving result
  item_plot = NULL
  try(item_plot <- enrichplot::cnetplot(emsigdb))
  pre_result <- emsigdb@result
  pre_result <- pre_result[pre_result$p.adjust < 0.1, ]
  result <- pre_result[order(pre_result$p.adjust), ]
  return(list(item_plot, result))
}



####
####
####
#### Testing Nominally differentially methylated probes for pathway over representation
####
####
####

# Get Epic Annotations
annEpic = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# Defining files to test
files.to.test <- c(list.files(path = "./unfiltered_results", pattern="*.csv"))
files.to.test <- paste(getwd(),"/unfiltered_results/", files.to.test, sep="")
# Defining nominal p-value significance threshold
sig_thresh <- 0.05
# Defining Genomic Locations of Interest
genomic_features <- c("TSS200" , "TSS1500", "5'UTR") # upstream
# genomic_features <- c("1stExon" , "Body") # gene body
# genomic_features <- "ALL"
# Collapsing Genomic Features for use in csv filename
gfs <- paste(genomic_features, collapse='_') # genomic features string
# MSigDb gene sets
all_gene_sets = msigdbr(species = "Homo sapiens")
# setting up directories for saving files
if (dir.exists(paste(getwd(), "/geneset_overrep_analysis", sep="")) == FALSE){
  dir.create(paste(getwd(), "/geneset_overrep_analysis", sep=""))
}
# Iterating through files for testing
for (genomic_features in  list(c("ALL"), c("TSS200" , "TSS1500", "5'UTR", "1stExon"), c("Body"))){
  gfs <- paste(genomic_features, collapse='_') # genomic features string
  for (file in files.to.test){
    print(paste("analyzing file:", file))
    readin_filepath <- unlist(strsplit(file, "/"))
    readin_filename <- readin_filepath[length(readin_filepath)]
    readin_filename <- substr(readin_filename[1], 1, nchar(readin_filename)-4)
    # making folder for saving results
    if (dir.exists(paste(getwd(), "/geneset_overrep_analysis/", readin_filename, sep="")) == FALSE){
      dir.create(paste(getwd(), "/geneset_overrep_analysis/", readin_filename, sep=""))
    }
    if (dir.exists(paste(getwd(), "/geneset_overrep_analysis/", readin_filename, "/", gfs, sep="")) == FALSE){
      dir.create(paste(getwd(), "/geneset_overrep_analysis/", readin_filename, "/", gfs, sep=""))
    }
    # running all over representation analyses
    df <- fread(file)
    # finding nominal p-values less than 0.05
    sig_df <- df[df$P.Value < sig_thresh]
    sig_probe_list <- sig_df$Name
    all_probe_list <- df$Name
    # pdf(file=paste("./geneset_overrep_analysis/", readin_filename, "/", gfs, "/go_ora.pdf", sep=""), width=10, height=8)
    go_res <- go_over_representation(sig_probe_list, all_probe_list, genomic_features, annEpic)
    # print(go_res[1])
    # dev.off()
    go_res <- as.data.frame(go_res[2])
    if (dim(go_res)[1] > 0){
      csv_name <- paste("./geneset_overrep_analysis/", readin_filename, "/", gfs, "/go_ora.csv", sep="")
      write.csv(go_res, csv_name)
    }
    pdf(file=paste("./geneset_overrep_analysis/", readin_filename, "/", gfs, "/kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(sig_probe_list, all_probe_list, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("./geneset_overrep_analysis/", readin_filename, "/", gfs, "/kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    pdf(file=paste("./geneset_overrep_analysis/", readin_filename, "/", gfs, "/reactome_ora.pdf", sep=""), width=10, height=8)
    reactome_res <- reactome_over_representation(sig_probe_list, all_probe_list, genomic_features, annEpic)
    print(reactome_res[1])
    reactome_res <- as.data.frame(reactome_res[2])
    if (dim(reactome_res)[1] > 0){
      csv_name <- paste("./geneset_overrep_analysis/", readin_filename, "/", gfs, "/reactome_ora.csv", sep="")
      write.csv(reactome_res, csv_name)
    }
    mesh_res <- mesh_over_representation(sig_probe_list, all_probe_list, genomic_features, annEpic)
    # mesh_res[1]
    mesh_res <- as.data.frame(mesh_res[2])
    if (dim(mesh_res)[1] > 0){
      csv_name <- paste("./geneset_overrep_analysis/", readin_filename, "/", gfs, "/mesh_ora.csv", sep="")
      write.csv(mesh_res, csv_name)
    }
    msigdbC7_res <- msigdb_over_representation(sig_probe_list, all_probe_list, genomic_features, annEpic, all_gene_sets, "C7")
    # msigdbC7_res[1]
    msigdbC7_res <- as.data.frame(msigdbC7_res[2])
    if (dim( msigdbC7_res)[1] > 0){
      csv_name <- paste("./geneset_overrep_analysis/", readin_filename, "/", gfs, "/ msigdbC7_ora.csv", sep="")
      write.csv(msigdbC7_res, csv_name)
    }
    msigdbC8_res <- msigdb_over_representation(sig_probe_list, all_probe_list, genomic_features, annEpic, all_gene_sets, "C8")
    # msigdbC8_res[1]
    msigdbC8_res <- as.data.frame(msigdbC8_res[2])
    if (dim( msigdbC8_res)[1] > 0){
      csv_name <- paste("./geneset_overrep_analysis/", readin_filename, "/", gfs, "/ msigdbC8_ora.csv", sep="")
      write.csv(msigdbC8_res, csv_name)
    }
  }
}
