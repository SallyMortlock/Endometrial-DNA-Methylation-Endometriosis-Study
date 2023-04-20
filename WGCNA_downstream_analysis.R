# Final script for WGCNA linear models with limma package for finding significant p-values
################################################################
# This script runs the downstream analysis associated with the
# WGCNA. The downstream analysis consists of building a set of
# linear models with Limma where each module (1st PC of WGCNA eigen module)
# becomes a function of the explanatory variables of cycle phase,
# case control status or case stage, surrogate variables, and
# covariates. The nominally significant WGCNA modules are then
# used for pathway over representation analysis using the KEGG
# and Reactome databases, where the results are automatically saved 
# to a folder titled wgcna_overrep_analysis. This script must be
# run after SVA_for_WGCNA.R

# Author: Parker Grosjean
# March 2021
################################################################

library(limma)
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
library(MeSH.Hsa.eg.db)#
library(ReactomePA)
library(meshes)
library(msigdbr)
library(DOSE)
library(enrichplot)
library(ggupset)
library(ggnewscale)

# Function to label significance of p-values when plotting heatmap
signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
         symbols = c("***", "**", "*", " "))
}

# defining functions for over represenation analysis
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
  item_plot <- enrichplot::cnetplot(ego)
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
  if (is.null(ekegg)){
    return(list(NULL, NULL))
  }
  ekegg <- setReadable(ekegg, 'org.Hs.eg.db', 'ENTREZID')
  print(dim(ekegg@result))
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
  item_plot <- NULL
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
  item_plot <- enrichplot::cnetplot(emesh)
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
  item_plot <- enrichplot::cnetplot(emsigdb)
  pre_result <- emsigdb@result
  pre_result <- pre_result[pre_result$p.adjust < 0.1, ]
  result <- pre_result[order(pre_result$p.adjust), ]
  return(list(item_plot, result))
}

############ Loading in the WGCNA Files ############
bwnet_raw_50_var <- readRDS("WGCNA_bwnet.rds")

############ Get Epic chip annotations ############
annEpic = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Change this below to run enrichment analysis for different WGCNA network instances
net <- bwnet_raw_50_var 

############ Loading in SVA data ############
svaout_50_var_cc_sv <- readRDS("svaout_50_var_cc_sv.rds")
svaout_50_var_stage_sv <- readRDS("svaout_50_var_stage_sv.rds")

############ Preprocessing for Enrichment Analysis ############
# Extracting Module EigenProbes and Probe Module Links
colors_raw = net$colors
# converting ME names for better figures
mod_me_names <- function(x){
  if (nchar(x) < 4){
    return(paste(substr(x, 1, 2), '0', substr(x, 3, nchar(x)), sep=''))
  }
  else{
    return(x)
  }
}

colors_raw <- paste(rep(c("ME"), length(colors_raw)), colors_raw, sep="")
colors_raw <- lapply(colors_raw, mod_me_names)
names(colors_raw) <- names(net$colors)

# setting List of all tested cpgs
probe_list = names(colors_raw)
names(probe_list) = NULL
ME_df <- net$MEs
colnames(ME_df) <- lapply(colnames(ME_df), mod_me_names)

# Merging MEs from WGCNA with Clinical Variables
load("SVA/for_sv_calc_WGCNA.RData")
head(cov_select_no_ffpe)
cov_select_no_ffpe$ID <- rownames(cov_select_no_ffpe)
phenotype_mappings = cov_select_no_ffpe[, c("ID", "Batch", "Institute_for_Analysis", "Endometriosis_Yes_No_", "Cycle_phase_for_Analysis", "Endometriosis_stage_grouped_I_II__III_IV_")]
ME_clinical_df <- merge(ME_df, phenotype_mappings, by.x="row.names", by.y="ID", how="inner")
ME_clinical_df <- ME_clinical_df %>% rename(mapping_ids=Row.names)

###### Preprocessing Pipeline for Case Control and Cycle Phase #######
# Merging clinical and WGCNA with SVs
sv_vec <- rep(c("SV"), dim(svaout_50_var_cc_sv)[2])
range_vec <- as.character(c(1:dim(svaout_50_var_cc_sv)[2]))
colnames(svaout_50_var_cc_sv) <- paste(sv_vec, range_vec, sep="")
ME_clinical_df <- merge(ME_clinical_df, svaout_50_var_cc_sv, by.x="mapping_ids", by.y="row.names", how="inner")

# Filtering Clinical and ME DataFrame
ME_names <- sort(colnames(ME_clinical_df)[grepl("ME*", colnames(ME_clinical_df))])
SV_names <- sort(colnames(ME_clinical_df)[grepl("SV*", colnames(ME_clinical_df))])
covariates <- c("Batch", "Institute_for_Analysis")
ME_clinical_df_mod <- ME_clinical_df[,c(ME_names, SV_names, "Endometriosis_Yes_No_", "Cycle_phase_for_Analysis", covariates)]
ME_clinical_df_mod <- ME_clinical_df_mod[ME_clinical_df_mod[,"Cycle_phase_for_Analysis"] != "SE",]

# Turning Clinical Covariates and Variables to Factors
ME_clinical_df_mod$Batch <- as.factor(ME_clinical_df_mod$Batch)
ME_clinical_df_mod$Institute_for_Analysis <- as.factor(ME_clinical_df_mod$Institute_for_Analysis)
ME_clinical_df_mod$Endometriosis_Yes_No_ <- as.factor(ME_clinical_df_mod$Endometriosis_Yes_No_)
ME_clinical_df_mod$Cycle_phase_for_Analysis <- factor(ME_clinical_df_mod$Cycle_phase_for_Analysis, levels = c("Menstrual", "PE", "ESE", "MSE", "LSE"))


###### Preprocessing Pipeline for Endometriosis Stage #######
# Merging clinical and WGCNA with SVs

ME_clinical_df <- merge(ME_df, phenotype_mappings, by.x="row.names", by.y="ID", how="inner")
ME_clinical_df <- ME_clinical_df %>% rename(mapping_ids=Row.names)

sv_vec <- rep(c("SV"), dim(svaout_50_var_stage_sv)[2])
range_vec <- as.character(c(1:dim(svaout_50_var_stage_sv)[2]))
colnames(svaout_50_var_stage_sv) <- paste(sv_vec, range_vec, sep="")
ME_clinical_df_stage <- merge(ME_clinical_df, svaout_50_var_stage_sv, by.x="mapping_ids", by.y="row.names", how="inner")

# Filtering Clinical and ME DataFrame
ME_names <- sort(colnames(ME_clinical_df_stage)[grepl("ME*", colnames(ME_clinical_df_stage))])
SV_names_stage <- sort(colnames(ME_clinical_df_stage)[grepl("SV*", colnames(ME_clinical_df_stage))])
covariates <- c("Batch", "Institute_for_Analysis")
ME_clinical_df_stage_mod <- ME_clinical_df_stage[,c(ME_names, SV_names_stage, "Endometriosis_stage_grouped_I_II__III_IV_", "Cycle_phase_for_Analysis", covariates)]
ME_clinical_df_stage_mod <- ME_clinical_df_stage_mod[ME_clinical_df_stage_mod[,"Endometriosis_stage_grouped_I_II__III_IV_"] != "UNKNOWN",]
ME_clinical_df_stage_mod <- ME_clinical_df_stage_mod[ME_clinical_df_stage_mod[,"Cycle_phase_for_Analysis"] != "SE",]

# Turning Clinical Covariates and Variables to Factors
ME_clinical_df_stage_mod$Batch <- as.factor(ME_clinical_df_stage_mod$Batch)
ME_clinical_df_stage_mod$Institute_for_Analysis <- as.factor(ME_clinical_df_stage_mod$Institute_for_Analysis)
ME_clinical_df_stage_mod$Endometriosis_stage_grouped_I_II__III_IV_ <- as.factor(ME_clinical_df_stage_mod$Endometriosis_stage_grouped_I_II__III_IV_)
ME_clinical_df_stage_mod$Cycle_phase_for_Analysis <- factor(ME_clinical_df_stage_mod$Cycle_phase_for_Analysis, levels = c("Menstrual", "PE", "ESE", "MSE", "LSE"))

# Limma Linear Modeling Case Control
colnames(ME_clinical_df_mod)
design_cc <- ME_clinical_df_mod[,c(SV_names, "Cycle_phase_for_Analysis", "Endometriosis_Yes_No_")]
head(design_cc)
mymod_no_ffpe_cc <- model.matrix(~Cycle_phase_for_Analysis+Endometriosis_Yes_No_+., design_cc)
head(mymod_no_ffpe_cc)
ME_vals_cc <- ME_clinical_df_mod[,ME_names]
fit_no_ffpe_cc <- lmFit(t(ME_vals_cc), mymod_no_ffpe_cc)
for_contrasting <-  c("Intercept", "Cycle_phase_for_AnalysisPE", "Cycle_phase_for_AnalysisESE", "Cycle_phase_for_AnalysisMSE", "Cycle_phase_for_AnalysisLSE", "Endometriosis_Yes_No_Yes") 
contrasts_case_control_no_ffpe <- makeContrasts("Endometriosis_Yes_No_Yes", levels = mymod_no_ffpe_cc)
fit_case_control_no_ffpe <- contrasts.fit(fit_no_ffpe_cc, contrasts_case_control_no_ffpe)
fit_case_control_no_ffpe <- eBayes(fit_case_control_no_ffpe)
summary(decideTests(fit_case_control_no_ffpe),adjust.method="BH",p.value=0.05)
cc_up <- rownames(decideTests(fit_case_control_no_ffpe)[decideTests(fit_case_control_no_ffpe) > 0, ])
cc_down <- rownames(decideTests(fit_case_control_no_ffpe)[decideTests(fit_case_control_no_ffpe) < 0, ])

# Limma Linear Modeling Stage and Cycle Stage
design_cc <- ME_clinical_df_stage_mod[,c(SV_names, "Cycle_phase_for_Analysis", "Endometriosis_stage_grouped_I_II__III_IV_")]
mymod_no_ffpe_stage <- model.matrix(~Cycle_phase_for_Analysis+Endometriosis_stage_grouped_I_II__III_IV_+., design_cc)
colnames(mymod_no_ffpe_stage) <- c("Intercept", "PE", "ESE", "MSE", "LSE", "I_II", "III_IV", "unknown", SV_names)
ME_vals_stage <- ME_clinical_df_stage_mod[,ME_names]
fit_no_ffpe_stage <- lmFit(t(ME_vals_stage), mymod_no_ffpe_stage)
for_contrasting <- c("III_IV - I_II", "III_IV", "I_II", "PE - (ESE + LSE + MSE)/3", "PE", "(ESE + LSE + MSE)/3")
# Testing for Stage III-IV vs Control
contrasts_stage_no_ffpe <- makeContrasts(for_contrasting[2], levels = mymod_no_ffpe_stage)
fit_stage_no_ffpe <- contrasts.fit(fit_no_ffpe_stage, contrasts_stage_no_ffpe)
fit_stage_no_ffpe <- eBayes(fit_stage_no_ffpe)
summary(decideTests(fit_stage_no_ffpe),adjust.method="BH",p.value=0.05)
stage_up <- rownames(decideTests(fit_stage_no_ffpe)[decideTests(fit_stage_no_ffpe) > 0, ])
stage_down <- rownames(decideTests(fit_stage_no_ffpe)[decideTests(fit_stage_no_ffpe) < 0, ])

# Testing for Stage I-II vs Control
contrasts_stagea_no_ffpe <- makeContrasts(for_contrasting[3], levels = mymod_no_ffpe_stage)
fit_stagea_no_ffpe <- contrasts.fit(fit_no_ffpe_stage, contrasts_stagea_no_ffpe)
fit_stagea_no_ffpe <- eBayes(fit_stagea_no_ffpe)
summary(decideTests(fit_stagea_no_ffpe),adjust.method="BH",p.value=0.05)
stagea_up <- rownames(decideTests(fit_stagea_no_ffpe)[decideTests(fit_stagea_no_ffpe) > 0, ])
stagea_down <- rownames(decideTests(fit_stagea_no_ffpe)[decideTests(fit_stagea_no_ffpe) < 0, ])

# Testing for Stage PE vs SE
contrasts_stage_no_ffpe_cp <- makeContrasts(for_contrasting[4], levels = mymod_no_ffpe_stage)
fit_stage_no_ffpe_cp <- contrasts.fit(fit_no_ffpe_stage, contrasts_stage_no_ffpe_cp)
fit_stage_no_ffpe_cp <- eBayes(fit_stage_no_ffpe_cp)
summary(decideTests(fit_stage_no_ffpe_cp),adjust.method="BH",p.value=0.05)
cp_up <- rownames(decideTests(fit_stage_no_ffpe_cp)[decideTests(fit_stage_no_ffpe_cp) > 0, ])
cp_down <- rownames(decideTests(fit_stage_no_ffpe_cp)[decideTests(fit_stage_no_ffpe_cp) < 0, ])

# Testing for Menstrual vs PE
contrasts_stage_no_ffpe_cp1 <- makeContrasts(for_contrasting[5], levels = mymod_no_ffpe_stage)
fit_stage_no_ffpe_cp1 <- contrasts.fit(fit_no_ffpe_stage, contrasts_stage_no_ffpe_cp1)
fit_stage_no_ffpe_cp1 <- eBayes(fit_stage_no_ffpe_cp1)
summary(decideTests(fit_stage_no_ffpe_cp1),adjust.method="BH",p.value=0.05)
cp1_up <- rownames(decideTests(fit_stage_no_ffpe_cp1)[decideTests(fit_stage_no_ffpe_cp1) > 0, ])
cp1_down <- rownames(decideTests(fit_stage_no_ffpe_cp1)[decideTests(fit_stage_no_ffpe_cp1) < 0, ])

# Testing for Menstrual vs SE
contrasts_stage_no_ffpe_cp2 <- makeContrasts(for_contrasting[6], levels = mymod_no_ffpe_stage)
fit_stage_no_ffpe_cp2 <- contrasts.fit(fit_no_ffpe_stage, contrasts_stage_no_ffpe_cp2)
fit_stage_no_ffpe_cp2 <- eBayes(fit_stage_no_ffpe_cp2)
summary(decideTests(fit_stage_no_ffpe_cp2),adjust.method="BH",p.value=0.05)
cp2_up <- rownames(decideTests(fit_stage_no_ffpe_cp2)[decideTests(fit_stage_no_ffpe_cp2) > 0, ])
cp2_down <- rownames(decideTests(fit_stage_no_ffpe_cp2)[decideTests(fit_stage_no_ffpe_cp2) < 0, ])

### Generating WGCNA Heatmap
cc_summary_MEs <- fit_case_control_no_ffpe$t
cc_pvals_MEs <- p.adjust(fit_case_control_no_ffpe$p.value, 'BH', length(cc_summary_MEs))
stage_summary_MEs <- fit_stage_no_ffpe$t
stage_pvals_MEs <- p.adjust(fit_stage_no_ffpe$p.value, 'BH', length(stage_summary_MEs))
stagea_summary_MEs <- fit_stagea_no_ffpe$t
stagea_pvals_MEs <- p.adjust(fit_stagea_no_ffpe$p.value, 'BH', length(stagea_summary_MEs))
cp_summary_MEs <- fit_stage_no_ffpe_cp$t
cp_pvals_MEs <- p.adjust(fit_stage_no_ffpe_cp$p.value, 'BH', length(cp_summary_MEs))
cp1_summary_MEs <- fit_stage_no_ffpe_cp1$t
cp1_pvals_MEs <- p.adjust(fit_stage_no_ffpe_cp1$p.value, 'BH', length(cp1_summary_MEs))
cp2_summary_MEs <- fit_stage_no_ffpe_cp2$t
cp2_pvals_MEs <- p.adjust(fit_stage_no_ffpe_cp2$p.value, 'BH', length(cp2_summary_MEs))
module_es <- matrix(data = c(cc_summary_MEs[,1], stagea_summary_MEs[,1], stage_summary_MEs[,1], cp_summary_MEs[,1], cp1_summary_MEs[,1], cp2_summary_MEs[,1]), nrow=length(cp_summary_MEs[,1]), ncol=6)
rownames(module_es) <- rownames(cc_summary_MEs)
text_matrix <- matrix(rep(c(""), length(module_es)*3), nrow=dim(module_es)[1], ncol=6)
# Modifying for Case Control
for (ind in seq(1, length(rownames(module_es)))){
  text_matrix[ind, 1] = signif.num(cc_pvals_MEs[ind])
}
# Modifying for Stage
for (ind in seq(1, length(rownames(module_es)))){
  text_matrix[ind, 2] = signif.num(stagea_pvals_MEs[ind])
}
# Modifying for Stage
for (ind in seq(1, length(rownames(module_es)))){
  text_matrix[ind, 3] = signif.num(stage_pvals_MEs[ind])
}
# Modifying for Cycle Phase
for (ind in seq(1, length(rownames(module_es)))){
  text_matrix[ind,4] = signif.num(cp_pvals_MEs[ind])
}
for (ind in seq(1, length(rownames(module_es)))){
  text_matrix[ind, 5] = signif.num(cp1_pvals_MEs[ind])
}
for (ind in seq(1, length(rownames(module_es)))){
  text_matrix[ind, 6] = signif.num(cp2_pvals_MEs[ind])
}
png("WGCNA/Module_Heatmap.png", height = 3000, width = 2000, res = 300)
labeledHeatmap(Matrix = module_es,
               xLabels = c("Case - Control", "St. I-II - Control", "St. III-IV - Control", "PE - SE", "PE - Menstrual", "SE - Menstrual"),
               yLabels = rownames(module_es),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               setStdMargins = TRUE,
               cex.text = 2,
               textAdj = c(0.5, 0.75),
               zlim = c(min(module_es),max(module_es)),
               textMatrix = text_matrix,
               verticalSeparator.interval = 1,
               horizontalSeparator.interval = 1)
dev.off()

### Testing for Pathway Over Representation
setwd("WGCNA")
if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis", sep="")) == FALSE){
  dir.create(paste(getwd(), "/wgcna_overrep_analysis", sep=""))
}

# defining all probes tested
all_probes <- names(colors_raw)
# defining genomic_features of interest
for (genomic_features in  list(c("ALL"), c("TSS200" , "TSS1500", "5'UTR", "1stExon"), c("Body"))){
  print(paste("genomic_features", genomic_features))
  # Collapsing Genomic Features for use in csv filename
  gfs <- paste(genomic_features, collapse='_') # genomic features string
  
  ### Case Control
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cc_up", sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cc_up", sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cc_down", sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cc_down", sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cc_up/", gfs, sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cc_up/", gfs,  sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cc_down/", gfs, sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cc_down/", gfs, sep=""))
  }
  ## Testing increased methylation
  cc_up_probes <- names(colors_raw[colors_raw %in% cc_up])
  if (length(cc_up_probes) > 0){
    # KEGG
    pdf(file=paste("./wgcna_overrep_analysis/cc_up/", gfs, "/kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(cc_up_probes, all_probes, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cc_up/", gfs, "/kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    # Reactome
    pdf(file=paste("./wgcna_overrep_analysis/cc_up/", gfs, "/reactome_ora.pdf", sep=""), width=10, height=8)
    react_res <- reactome_over_representation(cc_up_probes, all_probes, genomic_features, annEpic)
    print(react_res[1])
    dev.off()
    react_res <- as.data.frame(react_res[2])
    if (dim(react_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cc_up/", gfs, "/reactome_ora.csv", sep="")
      write.csv(react_res, csv_name)
    }
  }
  ## Testing decreased methylation
  cc_down_probes <- names(colors_raw[colors_raw %in% cc_down])
  if (length(cc_down_probes) > 0){
    # KEGG
    pdf(file=paste("./wgcna_overrep_analysis/cc_down/", gfs, "/kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(cc_down_probes, all_probes, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cc_down/", gfs, "/kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    # Reactome
    pdf(file=paste("./wgcna_overrep_analysis/cc_down/", gfs, "/reactome_ora.pdf", sep=""), width=10, height=8)
    react_res <- reactome_over_representation(cc_down_probes, all_probes, genomic_features, annEpic)
    print(react_res[1])
    dev.off()
    react_res <- as.data.frame(react_res[2])
    if (dim(react_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cc_down/", gfs, "/reactome_ora.csv", sep="")
      write.csv(react_res, csv_name)
    }
  }
  
  ## Stage
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/stage_up", sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/stage_up", sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/stage_down", sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/stage_down", sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/stage_up/", gfs, sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/stage_up/", gfs, sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/stage_down/", gfs, sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/stage_down/", gfs, sep=""))
  }
  # Testing increased methylation
  stage_up_probes <- names(colors_raw[colors_raw %in% stage_up])
  if (length(stage_up_probes) > 0){
    # KEGG
    pdf(file=paste("./wgcna_overrep_analysis/stage_up/", gfs, "/kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(stage_up_probes, all_probes, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/stage_up/", gfs, "/kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    # Reactome
    pdf(file=paste("./wgcna_overrep_analysis/stage_up/", gfs, "/reactome_ora.pdf", sep=""), width=10, height=8)
    react_res <- reactome_over_representation(stage_up_probes, all_probes, genomic_features, annEpic)
    print(react_res[1])
    dev.off()
    react_res <- as.data.frame(react_res[2])
    if (dim(react_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/stage_up/", gfs, "/reactome_ora.csv", sep="")
      write.csv(react_res, csv_name)
    }
  }
  # Testing decreased methylation
  stage_down_probes <- names(colors_raw[colors_raw %in% stage_down])
  if (length(stage_down_probes) > 0){
    # KEGG
    pdf(file=paste("./wgcna_overrep_analysis/stage_down/", gfs, "/kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(stage_down_probes, all_probes, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/stage_down/", gfs, "/kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    # Reactome
    pdf(file=paste("./wgcna_overrep_analysis/stage_down/", gfs, "/reactome_ora.pdf", sep=""), width=10, height=8)
    react_res <- reactome_over_representation(stage_down_probes, all_probes, genomic_features, annEpic)
    print(react_res[1])
    dev.off()
    react_res <- as.data.frame(react_res[2])
    if (dim(react_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/stage_down/", gfs, "/reactome_ora.csv", sep="")
      write.csv(react_res, csv_name)
    }
  }
  
  ## Cycle Phase PE - SE
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_se_up", sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_se_up", sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_se_down", sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_se_down", sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_se_up/", gfs, sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_se_up/", gfs, sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_se_down/", gfs, sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_se_down/", gfs, sep=""))
  }
  # Testing increased methylation
  cp_up_probes <- names(colors_raw[colors_raw %in% cp_up])
  if (length(cp_up_probes) > 0){
    # KEGG
    pdf(file=paste("./wgcna_overrep_analysis/cp_pe_se_up/", gfs, "/kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(cp_up_probes, all_probes, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_pe_se_up/", gfs, "/kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    # Reactome
    pdf(file=paste("./wgcna_overrep_analysis/cp_pe_se_up/", gfs, "/reactome_ora.pdf", sep=""), width=10, height=8)
    react_res <- reactome_over_representation(cp_up_probes, all_probes, genomic_features, annEpic)
    print(react_res[1])
    dev.off()
    react_res <- as.data.frame(react_res[2])
    if (dim(react_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_pe_se_up/", gfs, "/reactome_ora.csv", sep="")
      write.csv(react_res, csv_name)
    }
  }
  # Testing decreased methylation
  cp_down_probes <- names(colors_raw[colors_raw %in% cp_down])
  if (length(cp_down_probes) > 0){
    # KEGG
    pdf(file=paste("./wgcna_overrep_analysis/cp_pe_se_down/", gfs, "/kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(cp_down_probes, all_probes, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_pe_se_down/", gfs, "/kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    # Reactome
    pdf(file=paste("./wgcna_overrep_analysis/cp_pe_se_down/", gfs, "/reactome_ora.pdf", sep=""), width=10, height=8)
    react_res <- reactome_over_representation(cp_down_probes, all_probes, genomic_features, annEpic)
    print(react_res[1])
    dev.off()
    react_res <- as.data.frame(react_res[2])
    if (dim(react_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_pe_se_down/", gfs, "/reactome_ora.csv", sep="")
      write.csv(react_res, csv_name)
    }
  }
  
  ## Cycle Phase PE - Menstrual
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_menstrual_up", sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_menstrual_up", sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_menstrual_down", sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_menstrual_down", sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_menstrual_up/", gfs, sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_menstrual_up/", gfs,  sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_menstrual_down/", gfs,  sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_pe_menstrual_down/", gfs,  sep=""))
  }
  # Testing increased methylation
  cp1_up_probes <- names(colors_raw[colors_raw %in% cp1_up])
  if (length(cp1_up_probes) > 0){
    # KEGG
    pdf(file=paste("./wgcna_overrep_analysis/cp_pe_menstrual_up/", gfs, "/kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(cp1_up_probes, all_probes, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_pe_menstrual_up/", gfs, "/kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    # Reactome
    pdf(file=paste("./wgcna_overrep_analysis/cp_pe_menstrual_up/", gfs, "/reactome_ora.pdf", sep=""), width=10, height=8)
    react_res <- reactome_over_representation(cp1_up_probes, all_probes, genomic_features, annEpic)
    print(react_res[1])
    dev.off()
    react_res <- as.data.frame(react_res[2])
    if (dim(react_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_pe_menstrual_up/", gfs, "/reactome_ora.csv", sep="")
      write.csv(react_res, csv_name)
    }
  }
  # Testing decreased methylation
  cp1_down_probes <- names(colors_raw[colors_raw %in% cp1_down])
  if (length(cp1_down_probes) > 0){
    # KEGG
    pdf(file=paste("./wgcna_overrep_analysis/cp_pe_menstrual_down/", gfs, "/kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(cp1_down_probes, all_probes, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_pe_menstrual_down/", gfs, "/kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    # Reactome
    pdf(file=paste("./wgcna_overrep_analysis/cp_pe_menstrual_down/", gfs, "/reactome_ora.pdf", sep=""), width=10, height=8)
    react_res <- reactome_over_representation(cp1_down_probes, all_probes, genomic_features, annEpic)
    print(react_res[1])
    dev.off()
    react_res <- as.data.frame(react_res[2])
    if (dim(react_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_pe_menstrual_down/", gfs, "/reactome_ora.csv", sep="")
      write.csv(react_res, csv_name)
    }
  }
  
  ## Cycle Phase SE - Menstrual
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_se_menstrual_up", sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_se_menstrual_up", sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_se_menstrual_down", sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_se_menstrual_down", sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_se_menstrual_up/", gfs, sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_se_menstrual_up/", gfs, sep=""))
  }
  if (dir.exists(paste(getwd(), "/wgcna_overrep_analysis/cp_se_menstrual_down/", gfs, sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_overrep_analysis/cp_se_menstrual_down/", gfs, sep=""))
  }
  # Testing increased methylation
  cp2_up_probes <- names(colors_raw[colors_raw %in% cp2_up])
  if (length(cp2_up_probes) > 0){
    # KEGG
    pdf(file=paste("./wgcna_overrep_analysis/cp_se_menstrual_up/", gfs, "/kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(cp2_up_probes, all_probes, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_se_menstrual_up/", gfs, "/kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    # Reactome
    pdf(file=paste("./wgcna_overrep_analysis/cp_se_menstrual_up/", gfs, "/reactome_ora.pdf", sep=""), width=10, height=8)
    react_res <- reactome_over_representation(cp2_up_probes, all_probes, genomic_features, annEpic)
    print(react_res[1])
    dev.off()
    react_res <- as.data.frame(react_res[2])
    if (dim(react_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_se_menstrual_up/", gfs, "/reactome_ora.csv", sep="")
      write.csv(react_res, csv_name)
    }
  }
  # Testing decreased methylation
  cp2_down_probes <- names(colors_raw[colors_raw %in% cp2_down])
  if (length(cp2_down_probes) > 0){
    # KEGG
    pdf(file=paste("./wgcna_overrep_analysis/cp_se_menstrual_down/", gfs, "/kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(cp2_down_probes, all_probes, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_se_menstrual_down/", gfs, "/kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    # Reactome
    pdf(file=paste("./wgcna_overrep_analysis/cp_se_menstrual_down/", gfs, "/reactome_ora.pdf", sep=""), width=10, height=8)
    react_res <- reactome_over_representation(cp2_down_probes, all_probes, genomic_features, annEpic)
    print(react_res[1])
    dev.off()
    react_res <- as.data.frame(react_res[2])
    if (dim(react_res)[1] > 0){
      csv_name <- paste("./wgcna_overrep_analysis/cp_se_menstrual_down/", gfs, "/reactome_ora.csv", sep="")
      write.csv(react_res, csv_name)
    }
  }
}


############## Per Modules Over Representation Analysis ##############
if (dir.exists(paste(getwd(), "/wgcna_module_overrep_analysis", sep="")) == FALSE){
  dir.create(paste(getwd(), "/wgcna_module_overrep_analysis", sep=""))
}

# iterating through all modules
for (genomic_features in  list(c("ALL"), c("TSS200" , "TSS1500", "5'UTR", "1stExon"), c("Body"))){
  print(paste("genomic_features", genomic_features))
  # Collapsing Genomic Features for use in csv filename
  gfs <- paste(genomic_features, collapse='_') # genomic features string
  if (dir.exists(paste(getwd(), "/wgcna_module_overrep_analysis/", gfs, "/", sep="")) == FALSE){
    dir.create(paste(getwd(), "/wgcna_module_overrep_analysis/", gfs, "/", sep=""))
  }
  for (module in unique(colors_raw)){
    probes <- names(colors_raw[colors_raw==module])
    # KEGG
    pdf(file=paste("./wgcna_module_overrep_analysis/", gfs, "/", module, "_kegg_ora.pdf", sep=""), width=10, height=8)
    kegg_res <- kegg_over_representation(probes, all_probes, genomic_features, annEpic)
    print(kegg_res[1])
    dev.off()
    kegg_res <- as.data.frame(kegg_res[2])
    if (dim(kegg_res)[1] > 0){
      csv_name <- paste("wgcna_module_overrep_analysis/", gfs, "/", module, "_kegg_ora.csv", sep="")
      write.csv(kegg_res, csv_name)
    }
    # Reactome
    pdf(file=paste("./wgcna_module_overrep_analysis/", gfs, "/", module, "_reactome_ora.pdf", sep=""), width=10, height=8)
    react_res <- kegg_over_representation(probes, all_probes, genomic_features, annEpic)
    print(react_res[1])
    dev.off()
    react_res <- as.data.frame(react_res[2])
    if (dim(react_res)[1] > 0){
      csv_name <- paste("wgcna_module_overrep_analysis/", gfs, "/", module, "_reactome_ora.csv", sep="")
      write.csv(react_res, csv_name)
    }
  }
}
