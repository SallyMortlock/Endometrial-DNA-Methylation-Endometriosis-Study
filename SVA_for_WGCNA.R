################################################################
# This script runs a surrogate variable analysis (SVA) for use
# in linear models used to find statistically significant
# deferentially methylated WGCNA modules.
# The output of this script is used in the WGCNA_downstream_analysis.R
# script and thus this scripts must be run prior to running
# the downstream WGCNA analysis.

# Author: Parker Grosjean
# March 2021
################################################################

library(isva)
library(SmartSVA)
library(openxlsx)
library(data.table)
library(tidyverse)

############ Calculating SVs to include in Linear Models ############
annotation_file <- data.frame(read.xlsx("Data_Annotation.xlsx"),
                              stringsAsFactors = FALSE)
annotation_file$Sample.type<-gsub("\\ ","_",annotation_file$Sample.type)
annotation_file<-annotation_file[!is.na(annotation_file$Epic_Complete.Bar.code),]
rownames(annotation_file)<-annotation_file$Epic_Complete.Bar.code

# Reading in M-values
M_batch1_batch2 <- data.frame(fread("M_batch1_batch2.csv"))
rownames(M_batch1_batch2)<-M_batch1_batch2[,1]
M_batch1_batch2<-M_batch1_batch2[,-c(1)]
colnames(M_batch1_batch2)<-gsub("X","",colnames(M_batch1_batch2))

cov<-annotation_file[colnames(M_batch1_batch2),c("Batch", "Institute.for.Analysis", "Endometriosis..Yes.No.", "Cycle.phase.for.Analysis", "Endometriosis.stage.grouped..I.II....III.IV.", "Sample.type")]
levs <- c(" Menstrual", "PE", "ESE", "MSE", "LSE")
cov$Cycle.phase.for.Analysis<- factor(cov$Cycle.phase.for.Analysis, levs)
cov$Endometriosis.stage.grouped..I.II....III.IV. <- as.factor(cov$Endometriosis.stage.grouped..I.II....III.IV.)

# Taking the top 50% most variable CpG probes for running SVA and WGCNA
var_probes <- apply(M_batch1_batch2, FUN = var, MARGIN = 1)
var_quantile <- quantile(var_probes, na.rm = T)
M_batch1_batch2_50_var <- M_batch1_batch2[which(var_probes>var_quantile[3]),]  ## remove bottom 50% variance probes

# saving checkpoint
save(cov, M_batch1_batch2_50_var, file="for_sv_calc_unfiltered.RData")
load("for_sv_calc_unfiltered.RData")

## SV Calc for 50 var (Case Control and Cycle Phase)
Y.r_no_ffpe_50_var <- t(resid(lm(t(M_batch1_batch2_50_var) ~ Endometriosis..Yes.No. + Cycle.phase.for.Analysis, data=cov)))
mod_no_ffpe_50_var <- model.matrix(~Endometriosis..Yes.No. + Cycle.phase.for.Analysis, data=cov)
mod_0_no_ffpe_50_var <- model.matrix(~1, data=cov)
num_sv_no_ffpe_50_var <- EstDimRMT(Y.r_no_ffpe_50_var, FALSE)$dim + 1 
svaout_50_var_cc <- smartsva.cpp(as.matrix(M_batch1_batch2_50_var), n.sv=49, mod_no_ffpe_50_var, mod_0_no_ffpe_50_var, VERBOSE = T, epsilon=1e-2)
svaout_50_var_cc_sv <- svaout_50_var_cc$sv 
rownames(svaout_50_var_cc_sv)<-colnames(M_batch1_batch2_50_var)

## SV Calc for 50 var (Stage)
Y.r_no_ffpe_50_var <- t(resid(lm(t(M_batch1_batch2_50_var) ~ Cycle.phase.for.Analysis + Endometriosis.stage.grouped..I.II....III.IV., data=cov)))
mod_no_ffpe_50_var <- model.matrix(~Cycle.phase.for.Analysis + Endometriosis.stage.grouped..I.II....III.IV., data=cov)
mod_0_no_ffpe_50_var <- model.matrix(~1, data=cov)
num_sv_no_ffpe_50_var <- EstDimRMT(Y.r_no_ffpe_50_var, FALSE)$dim + 1
svaout_no_ffpe_50_var <- smartsva.cpp(as.matrix(M_batch1_batch2_50_var),n.sv=num_sv_no_ffpe_50_var, mod_no_ffpe_50_var, mod_0_no_ffpe_50_var, VERBOSE = T, epsilon=1e-2)
svaout_50_var_stage <- svaout_no_ffpe_50_var
svaout_50_var_stage_sv <- svaout_50_var_stage$sv 
rownames(svaout_50_var_stage_sv)<-colnames(M_batch1_batch2_50_var)

# Saving Data for 50 var
save(svaout_50_var_stage_sv, svaout_50_var_cc_sv, file="SVs_50_var_unfiltered.RData")

