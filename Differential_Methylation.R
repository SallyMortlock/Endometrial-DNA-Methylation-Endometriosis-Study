
###################################################
#
# analysis_model_SVA_final.R
# Model fitting with SVA
#
# Idit Kosti, PhD,  Marina Sirota, PhD, Parker Grosjean
# Mar 2021
# Rerun and edited Jul 2022 by Sally Mortlock, PhD
#
###################################################

install.packages('qqman')
install.packages("readxl")
install.packages("ggcorrplot")
install.packages("UpSetR")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap") #sudo apt-get install libcairo2-dev
BiocManager::install("DMRcate")
BiocManager::install( "impute")
install.packages("WGCNA")
install.packages("SmartSVA")
install.packages("tidyr")

library(tidyr)
library(WGCNA)
library(dplyr)
library(DMRcate)
library("UpSetR")
library(ggcorrplot)
library("readxl")
library(limma)
library(ggplot2)
library(EnhancedVolcano)
library(Matrix)
library(corrplot)
library(ComplexHeatmap)
library(qqman)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(lattice)
library(vioplot)
library(bumphunter)
library(bacon)
library(SmartSVA)
library(caret)
library(stringr)
library(missMethyl)
library(UpSetR)
library(sva)


### FUNCTIONS ####

plot_qq_manhattan <- function(name, fit) {
  #res1<-data.frame(topTable(fit,number = 759345))
  
  png(paste(name, "_qq_t.png", sep = ""))
  main<-paste(gsub("_", " ", name)," Lambda = ",median(qchisq(fit$p.value, df=1, lower.tail=FALSE)) / qchisq(0.5, 1),sep="")
  qqt(fit$t,df=fit$df.residual+fit$df.prior,main = main,cex.main=0.7)
  abline(0,1)
  dev.off()
  
  
  #png(paste(name, "_manhattan.png", sep = ""), w = 12, h = 8, units = "in", res = 200)
  #res1[,"CHR"] <- genomic_locations[rownames(res1),"chr"]
  #res1$CHR<-gsub("chr","",res1$CHR)
  #res1$CHR<-as.numeric(res1$CHR)
  #res1[,"pos"] <- as.numeric(genomic_locations[rownames(res1),"pos"])
  #res1<-res1[complete.cases(res1),]
  
  #manhattan(res1, chr="CHR", bp="pos",  p="adj.P.Val", suggestiveline = T, genomewideline = T, main = gsub("_", " ", name))
  
  #dev.off()
}

plot_heatmap_agg <- function(name, data, cov, probes, samples) {
  png(paste(name, ".png", sep = ""), width = 4000, height = 3000, res = 300)
  cov_select2 <- data.frame(cov[samples,])
  col_cyclephase = RColorBrewer::brewer.pal(name = "Set3",n=6)
  names(col_cyclephase) <- c("PE","ESE","MSE","LSE","Menstrual","SE")
  col_institute = RColorBrewer::brewer.pal(name = "Set1",n=4)
  names(col_institute)<-c("UM","UCSF","EndOx","CIR_University_of_Edinburgh")
  col_Batch = c("1" = "#9E3CC3", "2" = "#7BC213")
  names(col_Batch)<-c("1","2")
  col_plate = c("#006A40FF", "#F08892FF", "#75B41EFF", "#95828DFF", "#708C98FF", "#8AB8CFFF", "#007E7FFF", "#358359FF", "#8BA1BCFF", "#5A5895FF", "#F2990CFF", "#5A5895FF")
  names(col_plate)<-c("1473","1474","1475","1476","1477","1478","1479","1480","1547","1548","1549","1550")
  col_status = c("#FFA500","#4169E1")
  names(col_status)<-c("Yes","No")
  col_ethnicity=RColorBrewer::brewer.pal(name = "Set2",n=6)
  names(col_ethnicity)<-c("EUR","ADMIX","AMR","AFR","SAS","EAS")
  colors<-list(CyclePhase =col_cyclephase ,Institute = col_institute, Batch =col_Batch, Plate =col_plate, Status =col_status, Ethnicity=col_ethnicity )
  column_ha = HeatmapAnnotation(CyclePhase = cov_select2$Cycle_phase_for_Analysis_agg,Institute = cov_select2$Institute_for_Analysis,Batch = cov_select2$Batch, Plate = cov_select2$BS_Plate_No_, Status = cov_select2$Endometriosis_Yes_No_,Ethnicity=cov_select2$Genetic_Ancestry,col=colors)
  h<-Heatmap(as.matrix(data[probes,samples]), name = "mat", top_annotation = column_ha, cluster_rows = TRUE, show_column_names = FALSE, show_row_names = FALSE,column_title = gsub("_", " ", name))
  print(h)
  dev.off()
}

plot_heatmap_agg2 <- function(name, data, cov, probes, samples) {
  png(paste(name, ".png", sep = ""), width = 4000, height = 3000, res = 300)
  cov_select2 <- data.frame(cov[samples,])
  col_cyclephase = RColorBrewer::brewer.pal(name = "Set3",n=6)
  names(col_cyclephase) <- c("PE","ESE","MSE","LSE","Menstrual","SE")
  colors<-list(CyclePhase =col_cyclephase)
  column_ha = HeatmapAnnotation(CyclePhase = cov_select2$Cycle_phase_for_Analysis_agg,col=colors, show_legend = FALSE)
  h<-Heatmap(as.matrix(data[probes,samples]), name = "mat", top_annotation = column_ha, cluster_rows = TRUE, show_heatmap_legend = FALSE, show_column_names = FALSE, show_row_names = FALSE,column_title = gsub("_", " ", name))
  print(h)
  dev.off()
}

plot_heatmap <- function(name, data, cov, probes, samples) {
  png(paste(name, ".png", sep = ""), width = 4000, height = 3000, res = 300)
  cov_select2 <- data.frame(cov[samples,])
  col_cyclephase = RColorBrewer::brewer.pal(name = "Set3",n=6)
  names(col_cyclephase) <- c("PE","ESE","MSE","LSE","Menstrual","SE")
  col_institute = RColorBrewer::brewer.pal(name = "Set1",n=4)
  names(col_institute)<-c("UM","UCSF","EndOx","CIR_University_of_Edinburgh")
  col_Batch = c("1" = "#9E3CC3", "2" = "#7BC213")
  names(col_Batch)<-c("1","2")
  col_plate = c("#006A40FF", "#F08892FF", "#75B41EFF", "#95828DFF", "#708C98FF", "#8AB8CFFF", "#007E7FFF", "#358359FF", "#8BA1BCFF", "#5A5895FF", "#F2990CFF", "#5A5895FF")
  names(col_plate)<-c("1473","1474","1475","1476","1477","1478","1479","1480","1547","1548","1549","1550")
  col_status = c("#FFA500","#4169E1")
  names(col_status)<-c("Yes","No")
  col_ethnicity=RColorBrewer::brewer.pal(name = "Set2",n=6)
  names(col_ethnicity)<-c("EUR","ADMIX","AMR","AFR","SAS","EAS")
  colors<-list(CyclePhase =col_cyclephase ,Institute = col_institute, Batch =col_Batch, Plate =col_plate, Status =col_status, Ethnicity=col_ethnicity )
  column_ha = HeatmapAnnotation(CyclePhase = cov_select2$Cycle_phase_for_Analysis,Institute = cov_select2$Institute_for_Analysis,Batch = cov_select2$Batch, Plate = cov_select2$BS_Plate_No_, Status = cov_select2$Endometriosis_Yes_No_,Ethnicity=cov_select2$Genetic_Ancestry,col=colors)
  h<-Heatmap(as.matrix(data[probes,samples]), name = "mat", top_annotation = column_ha, cluster_rows = TRUE, show_column_names = FALSE, show_row_names = FALSE,column_title = gsub("_", " ", name))
  print(h)
  dev.off()
}


plot_pca <- function(name, PCs, cov_select) {
  mynames<-c("BS_Plate_No_","Institute_for_Analysis","Endometriosis_Yes_No_","Cycle_phase_for_Analysis","Genetic_Ancestry","Batch", "Sample_type")
  
  col_cyclephase = RColorBrewer::brewer.pal(name = "Set3",n=6)
  names(col_cyclephase) <- c("PE","ESE","MSE","LSE","Menstrual","SE")
  col_institute = RColorBrewer::brewer.pal(name = "Set1",n=4)
  names(col_institute)<-c("UM","UCSF","EndOx","CIR_University_of_Edinburgh")
  col_batch = c("1" = "#9E3CC3", "2" = "#7BC213")
  names(col_batch)<-c("1","2")
  col_plate = c("#006A40FF", "#F08892FF", "#75B41EFF", "#95828DFF", "#708C98FF", "#8AB8CFFF", "#007E7FFF", "#358359FF", "#8BA1BCFF", "#5A5895FF", "#F2990CFF", "#5A5895FF")
  names(col_plate)<-c("1473","1474","1475","1476","1477","1478","1479","1480","1547","1548","1549","1550")
  col_status = c("#FFA500","#4169E1")
  names(col_status)<-c("Yes","No")
  col_ethnicity=RColorBrewer::brewer.pal(name = "Set2",n=6)
  names(col_ethnicity)<-c("EUR","ADMIX","AMR","AFR","SAS","EAS")
  col_sample_type <- c('black','chartreuse3')
  names(col_sample_type) <- c("FFPE", "Fresh Frozen")
  colors<-list(Plate =col_plate, Institute = col_institute,  Status =col_status, CyclePhase =col_cyclephase, Ethnicity=col_ethnicity, Batch =col_batch, Sample_Type = col_sample_type)
  
  pdf(name)
  PropEx <- round((PCs$sdev^2/sum(PCs$sdev^2))*100,digits=2)
  head(PropEx)
  barplot(PropEx[1:8], las=2, xlab='', ylab='% Variance Explained')
  for (i in 1:length(mynames)){
    print(i)
    myvar<-cov_select[,mynames[i]]
    
    if(class(myvar)=="factor"){levels(myvar)<-c(levels(myvar),"Missing")}
    
    myvar[which(is.na(myvar))]<-"Missing"
    myvar<-factor(myvar,levels=unique(myvar))
    colorlist<-colors[[i]]
    
    mycolors<-unlist(lapply(myvar,function(x){return(colorlist[which(levels(myvar)==x)])}))
    plot.new()
    par(mar=c(0, 0, 0, 0))
    legend("center",levels(myvar), title=mynames[i],fill=colorlist[1:length(levels(myvar))])
    pairs(PCs$x[,1:3],labels=paste0("PC",1:3),col=mycolors,pch=20,cex=2)
  }
  dev.off()
}

plot_corr <- function(PCs, cov_select) {
  x<-cov_select
  x[,1] <- as.numeric(x[,1])
  for(i in 2:dim(x)[2]) {
    x[,i] <- as.factor(x[,i])
  }
  x$PC1<-as.numeric(PCs$x[,1])
  x$PC2<-as.numeric(PCs$x[,2])
  x$PC3<-as.numeric(PCs$x[,3])
  x$PC4<-as.numeric(PCs$x[,4])
  x$PC5<-as.numeric(PCs$x[,5])
  x$PC6<-as.numeric(PCs$x[,6])
  x$PC7<-as.numeric(PCs$x[,7])
  x$PC8<-as.numeric(PCs$x[,8])
  x$PC9<-as.numeric(PCs$x[,9])
  x$PC10<-as.numeric(PCs$x[,10])
  
  model.matrix(~0+., data=x) %>% 
    cor(use="pairwise.complete.obs") %>% 
    ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size=2)
}

get_corrected_data <- function(y, mod, svaobj) {
  X=cbind(mod,svaobj$sv) 
  Hat=solve(t(X)%*%X)%*%t(X) 
  beta=(Hat%*%t(y)) 
  P=ncol(mod) 
  cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),]) 
  return(cleany) 
}


# Get Epic chip annotation
annEpic = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Pull genomic locations
genomic_locations<-annEpic[rownames(M_batch1_batch2),c("chr","pos")]

# Load Data
M_batch1_batch2 <- readRDS("M_batch1_batch2_July2022.rds")
cov <- read.csv("cov_July2022.csv", header=T)
rownames(cov) <- cov$X
cov$X <- NULL

# Raw Data PCA and Correlation Plot
PCs <- prcomp(t(M_batch1_batch2),center=F, scale.=F)
saveRDS(PCs,"PCs.rds")

plot_pca("PCA_presva.pdf", PCs, cov)

pdf("corrplot_presva.pdf", w = 10, h = 10)
plot_corr(PCs, cov[,-6])
dev.off()

# Run SVA 
mod <- model.matrix(~Endometriosis_Yes_No_ + Cycle_phase_for_Analysis, data=cov)
mod_0 <- model.matrix(~1, data=cov)

Y.r <- t(resid(lm(t(M_batch1_batch2) ~ Endometriosis_Yes_No_ + Cycle_phase_for_Analysis, data=cov)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1 

svaout <- smartsva.cpp(as.matrix(M_batch1_batch2),n.sv=n.sv, mod, mod_0)
svaout_sv <- svaout$sv
saveRDS(svaout_sv,"svaout_sv.rds")

rownames(svaout_sv)<-colnames(M_batch1_batch2)

# Corrected Data PCA and Correlation Plot
M_batch1_batch2_corrected <- get_corrected_data(as.matrix(M_batch1_batch2),mod,svaout)

PCs_corrected <- prcomp(t(M_batch1_batch2_corrected),center=F, scale.=F)

plot_pca("PCA_corrected_postsva.pdf", PCs_corrected, cov)

pdf("corrplot_corrected_postsva.pdf", w = 10, h = 10)
plot_corr(PCs_corrected, cov[,-6])
dev.off()

# Single Site Analysis

# Raw Data + SVs
regressme <- data.frame(cov[,4], cov[,3],svaout_sv)
colnames(regressme)[1] <- "Cycle_phase_for_Analysis"
colnames(regressme)[2] <- "Endometriosis_Yes_No_"

mymod_cycle_phase <- model.matrix(~0+Cycle_phase_for_Analysis+Endometriosis_Yes_No_+., regressme)

fit_cycle_phase <- lmFit(M_batch1_batch2, mymod_cycle_phase)

contrasts_cyclephase<- makeContrasts(Cycle_phase_for_AnalysisESE-Cycle_phase_for_AnalysisLSE, Cycle_phase_for_AnalysisESE-Cycle_phase_for_AnalysisMSE, Cycle_phase_for_AnalysisESE-Cycle_phase_for_AnalysisPE, Cycle_phase_for_AnalysisLSE-Cycle_phase_for_AnalysisMSE,Cycle_phase_for_AnalysisLSE-Cycle_phase_for_AnalysisPE, Cycle_phase_for_AnalysisMSE-Cycle_phase_for_AnalysisPE, Cycle_phase_for_AnalysisMSE-Cycle_phase_for_AnalysisMenstrual, Cycle_phase_for_AnalysisESE-Cycle_phase_for_AnalysisMenstrual, Cycle_phase_for_AnalysisLSE-Cycle_phase_for_AnalysisMenstrual, Cycle_phase_for_AnalysisPE-Cycle_phase_for_AnalysisMenstrual, levels = mymod_cycle_phase)
 
fit_cycle_phase <-contrasts.fit(fit_cycle_phase, contrasts_cyclephase)
fit_cycle_phase <- eBayes(fit_cycle_phase)

summary(decideTests(fit_cycle_phase),adjust.method="BH",p.value=0.05)
saveRDS(fit_cycle_phase,"fit_cycle_phase.rds")

mymod_cc <- model.matrix(~Cycle_phase_for_Analysis+Endometriosis_Yes_No_+., regressme)
fit_cc <- lmFit(M_batch1_batch2,mymod_cc)

contrasts_case_control <- makeContrasts(Endometriosis_Yes_No_Yes,levels = mymod_cc)

fit_case_control<-contrasts.fit(fit_cc, contrasts_case_control)
fit_case_control<- eBayes(fit_case_control)

summary(decideTests(fit_case_control),adjust.method="BH",p.value=0.05)
#0

#Cycle Phase Analysis
cycle_phase_names <- c("ESE_LSE", "ESE_MSE", "ESE_PE", "MSE_LSE", "LSE_PE", "MSE_PE", "MSE_Menstrual", "ESE_Menstrual", "LSE_Menstrual", "PE_Menstrual")
fit_cycle_phase <- readRDS("fit_cycle_phase.rds")


# Single Site
for (i in 1:ncol(fit_cycle_phase))
{
  print(paste(colnames(fit_cycle_phase)[i], cycle_phase_names[i]))
  
  ## Plot qq plots and manhattam plots
  plot_qq_manhattan(paste("CyclePhase_M_values",cycle_phase_names[i], sep=""), fit_cycle_phase[,colnames(fit_cycle_phase)[i]])
}

for (i in 1:ncol(fit_cycle_phase))
{
  print(paste(colnames(fit_cycle_phase)[i], cycle_phase_names[i]))
  
  ## Get and save pvals and delta beta
  sig_sites <- topTable(fit_cycle_phase[,colnames(fit_cycle_phase)[i]],number = nrow(M_batch1_batch2))
  split<-strsplit(cycle_phase_names[i],"\\_")
  sig_sites <- sig_sites[sig_sites$adj.P.Val<0.05,]
  sig_sites <- sig_sites[order(sig_sites$P.Value),]
  sig_sites2 <- sig_sites[1:50,]
  cpgs <- which(rownames(M_batch1_batch2) %in% rownames(sig_sites2))
  
  ## Plot heatmap
  plot_heatmap(name = paste(cycle_phase_names[i], "_heatmap", sep = ""),data = M_batch1_batch2,cov = cov,probes = cpgs,samples = rownames(cov[cov$Cycle_phase_for_Analysis %in% c(split[[1]][1],split[[1]][2]),]))
  
  write.csv(sig_sites, paste(cycle_phase_names[i],".csv", sep = ""))
  
} 

# DMR
for (i in 1:ncol(fit_cycle_phase)) {
  print(cycle_phase_names[i])
  
  myAnnotation <- DMRcate::cpg.annotate(object = as.matrix(M_batch1_batch2), datatype = "array", what = "M",
                                        analysis.type = "differential", design = as.matrix(mymod_cycle_phase),
                                        contrasts = TRUE, cont.matrix = contrasts_cyclephase,
                                        coef = colnames(fit_cycle_phase)[i], arraytype = "EPIC", fdr = 0.1)
  DMRs <- dmrcate(myAnnotation, lambda=1000, C=2,pcutoff = 0.05)
  
  results.ranges <- extractRanges(DMRs, genome = "hg19")
  results.ranges_select<-results.ranges[results.ranges$Fisher<0.1,]
  
  print(paste("DMRs", length(results.ranges_select$no.cpgs)))
  
  write.csv(results.ranges_select,paste("DMRcate_results_select", cycle_phase_names[i], ".csv", sep = ""))

 
### CASE / CONTROL ###
plot_qq_manhattan(paste("CaseControl_M_values",cycle_phase_names[1], sep=""), fit_case_control[,colnames(fit_case_control)[1]])

# Single Sites
CC <- topTable(fit_case_control[,colnames(fit_case_control)[1]],number = nrow(M_batch1_batch2))
write.csv(CC, "CC.csv")

# DMR
myAnnotation_CC <- DMRcate::cpg.annotate(object = as.matrix(M_batch1_batch2), datatype = "array", what = "M",
                                         analysis.type = "differential", design = as.matrix(mymod_cc),
                                         contrasts = TRUE, cont.matrix = contrasts_case_control,
                                         coef = "Endometriosis_Yes_No_Yes", arraytype = "EPIC", fdr = 0.1)

DMRs_CC <- dmrcate(myAnnotation_CC, lambda=1000, C=2,pcutoff = 0.05)
results.ranges_CC <- extractRanges(DMRs_CC, genome = "hg19")
results.ranges_CC_select<-results.ranges_CC[results.ranges_CC$Fisher<0.1,]
write.csv(results.ranges_CC_select,"DMRcate_results_0_5_CC.csv")

### AGGREGATED ANALYSIS ####

cov[,"Cycle_phase_for_Analysis_agg"]<-cov$Cycle_phase_for_Analysis
cov$Cycle_phase_for_Analysis_agg<-gsub("ESE","SE",cov$Cycle_phase_for_Analysis_agg)
cov$Cycle_phase_for_Analysis_agg<-gsub("MSE","SE",cov$Cycle_phase_for_Analysis_agg)
cov$Cycle_phase_for_Analysis_agg<-gsub("LSE","SE",cov$Cycle_phase_for_Analysis_agg)

regressme_agg <- data.frame(cov$Cycle_phase_for_Analysis_agg,cov$Endometriosis_Yes_No_,svaout_sv)
colnames(regressme_agg)[1] <- "Cycle_phase_for_Analysis_agg"
colnames(regressme_agg)[2] <- "Endometriosis_Yes_No_"

mymod_agg_cycle_phase <- data.frame(model.matrix(~0+Cycle_phase_for_Analysis_agg+Endometriosis_Yes_No_+., regressme_agg))
fit_agg_cycle_phase <- lmFit(M_batch1_batch2,mymod_agg_cycle_phase)

contrasts_cyclephase_agg<- makeContrasts(Cycle_phase_for_Analysis_aggSE-Cycle_phase_for_Analysis_aggPE, Cycle_phase_for_Analysis_aggSE-Cycle_phase_for_Analysis_aggMenstrual, Cycle_phase_for_Analysis_aggPE-Cycle_phase_for_Analysis_aggMenstrual, levels = mymod_agg_cycle_phase)

fit_cycle_phase_agg <-contrasts.fit(fit_agg_cycle_phase, contrasts_cyclephase_agg)
fit_cycle_phase_agg <- eBayes(fit_cycle_phase_agg)

summary(decideTests(fit_cycle_phase_agg),adjust.method="BH",p.value=0.05)

mymod_agg_cc <- data.frame(model.matrix(~Cycle_phase_for_Analysis_agg+Endometriosis_Yes_No_+., regressme_agg))
fit_agg_cc <- lmFit(M_batch1_batch2, mymod_agg_cc)

contrasts_case_control_agg <- makeContrasts(Endometriosis_Yes_No_Yes,levels = mymod_agg_cc)

fit_case_control_agg <- contrasts.fit(fit_agg_cc, contrasts_case_control_agg)
fit_case_control_agg <- eBayes(fit_case_control_agg)

summary(decideTests(fit_case_control_agg),adjust.method="BH",p.value=0.05)

cycle_phase_names_agg <- c("SE_PE", "SE_Menstrual", "PE_Menstrual")

# Single Site
for (i in 1:ncol(fit_cycle_phase_agg))
{
  print(paste(colnames(fit_cycle_phase_agg)[i], cycle_phase_names_agg[i]))
  
  plot_qq_manhattan(paste("CyclePhase_M_values_agg",cycle_phase_names_agg[i], sep=""), fit_cycle_phase_agg[,colnames(fit_cycle_phase_agg)[i]])
  
  ## Get and save pvals and delta beta
  sig_sites <- topTable(fit_cycle_phase_agg[,colnames(fit_cycle_phase_agg)[i]],number = nrow(M_batch1_batch2))
  split<-strsplit(cycle_phase_names_agg[i],"\\_")
  sig_sites <- sig_sites[sig_sites$adj.P.Val<0.05,]
  sig_sites <- sig_sites[order(sig_sites$P.Value),]
  sig_sites2 <- sig_sites[1:50,]
  cpgs <- which(rownames(M_batch1_batch2) %in% rownames(sig_sites2))
  
  ## Plot heatmap
  plot_heatmap_agg(name = paste(cycle_phase_names_agg[i], "_heatmap_agg", sep = ""),data = M_batch1_batch2,cov = cov,probes = cpgs,samples = rownames(cov[cov$Cycle_phase_for_Analysis_agg %in% c(split[[1]][1],split[[1]][2]),]))
  
  write.csv(sig_sites, paste(cycle_phase_names_agg[i],"_agg.csv", sep = ""))
  
} 

# DMR
for (i in 1:ncol(fit_cycle_phase_agg)) {
  print(cycle_phase_names_agg[i])
  
  myAnnotation <- DMRcate::cpg.annotate(object = as.matrix(M_batch1_batch2), datatype = "array", what = "M",
                                        analysis.type = "differential", design = as.matrix(mymod_agg_cycle_phase),
                                        contrasts = TRUE, cont.matrix = contrasts_cyclephase_agg,
                                        coef = colnames(fit_cycle_phase_agg)[i], arraytype = "EPIC", fdr = 0.1)
  DMRs <- dmrcate(myAnnotation, lambda=1000, C=2,pcutoff = 0.05)
  
  results.ranges <- extractRanges(DMRs, genome = "hg19")
  results.ranges_select<-results.ranges[results.ranges$Fisher<0.1,]
  
  print(paste("DMRs", length(results.ranges_select$no.cpgs)))
  
  write.csv(results.ranges_select,paste("DMRcate_results_select", cycle_phase_names_agg[i], "_agg.csv", sep = ""))
  }

## CASE CONTROL ##
plot_qq_manhattan(paste("CaseControl_M_values_agg",colnames(fit_case_control_agg)[1], sep=""), fit_case_control_agg[,colnames(fit_case_control_agg)[1]])

# DMR
myAnnotation_CC <- DMRcate::cpg.annotate(object = as.matrix(M_batch1_batch2), datatype = "array", what = "M",
                                            analysis.type = "differential", design = as.matrix(mymod_agg_cc),
                                            contrasts = TRUE, cont.matrix = contrasts_case_control_agg,
                                            coef = "Endometriosis_Yes_No_Yes", arraytype = "EPIC", fdr = 0.1)
DMRs_CC <- dmrcate(myAnnotation_CC, lambda=1000, C=2,pcutoff = 0.05)
results.ranges_CC <- extractRanges(DMRs_CC)
results.ranges_CC_select<-results.ranges_CC[results.ranges_CC$Fisher<0.1,]
write.csv(results.ranges_CC_select,"DMRcate_results_select_CC_agg.csv")


### STAGE ANALYSIS ###

cov[,"Stage_grouped"]<-annotation_file[rownames(cov),colnames(annotation_file) %in% c("Endometriosis.stage.grouped..I.II....III.IV.")]
cov$Stage_grouped<-gsub("-","_",cov$Stage_grouped)

mod_stage <- model.matrix(~0+Stage_grouped + Cycle_phase_for_Analysis, data=cov)
mod_0_stage <- model.matrix(~1, data=cov)

Y.r <- t(resid(lm(t(M_batch1_batch2) ~ Stage_grouped + Cycle_phase_for_Analysis, data=cov)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1

svaout_stage <- smartsva.cpp(as.matrix(M_batch1_batch2),n.sv=n.sv, mod_stage, mod_0_stage)
svaout_sv_stage <- svaout_stage$sv
svaout_sv_stage <- readRDS("svaout_sv_stage.rds")

## I_II vs Control
index <- which(cov$Stage_grouped=="CONTROL"|cov$Stage_grouped=="I_II")
cov_C_I_II <- cov[index,]
M_batch1_batch2_C_I_II <- M_batch1_batch2[,index]

mod_stage <- model.matrix(~0+Stage_grouped + Cycle_phase_for_Analysis, data=cov_C_I_II)
mod_0_stage <- model.matrix(~1, data=cov_C_I_II)

Y.r <- t(resid(lm(t(M_batch1_batch2_C_I_II) ~ Stage_grouped + Cycle_phase_for_Analysis, data=cov_C_I_II)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1

svaout_stage_C_I_II <- smartsva.cpp(as.matrix(M_batch1_batch2_C_I_II),n.sv=n.sv, mod_stage, mod_0_stage)
svaout_sv_stage_C_I_II <- svaout_stage_C_I_II$sv
rownames(svaout_sv_stage_C_I_II) <- rownames(cov_C_I_II)
saveRDS(svaout_sv_stage_C_I_II, "svaout_sv_stage_C_I_II.rds")
write.table(svaout_sv_stage_C_I_II,"svaout_sv_stage_C_I_II.txt", quote=F, row.names=T, sep="\t")

svaout_sv_stage <- readRDS("svaout_sv_stage_C_I_II.rds")
index <- which(cov$Stage_grouped=="CONTROL"|cov$Stage_grouped=="I_II")
cov_stage <- cov[index,]
M_batch1_batch2_stage <- M_batch1_batch2[,index]

regressme_stage <- data.frame(cov_stage$Cycle_phase_for_Analysis,cov_stage$Endometriosis_Yes_No_,svaout_sv_stage)
colnames(regressme_stage)[1] <- "Cycle_phase_for_Analysis"
colnames(regressme_stage)[2] <- "Endometriosis_Yes_No_"

mymod_stage_cc <- data.frame(model.matrix(~Cycle_phase_for_Analysis+Endometriosis_Yes_No_+., regressme_stage))
fit_stage_cc <- lmFit(M_batch1_batch2_stage, mymod_stage_cc)

contrasts_case_control_stage <- makeContrasts(Endometriosis_Yes_No_Yes,levels = mymod_stage_cc)

fit_case_control_stage <- contrasts.fit(fit_stage_cc, contrasts_case_control_stage)
fit_case_control_stage <- eBayes(fit_case_control_stage)

summary(decideTests(fit_case_control_stage),adjust.method="BH",p.value=0.05)
saveRDS(fit_case_control_stage, "fit_C_I_II.rds")

plot_qq_manhattan(paste("CaseControl_M_values_C_I_II",colnames(fit_case_control_stage)[1], sep=""), fit_case_control_stage[,colnames(fit_case_control_stage)[1]])

fit_case_control_stage <- readRDS("fit_C_I_II.rds")
res <- topTable(fit_case_control_stage[,colnames(fit_case_control_stage)[1]],number = nrow(M_batch1_batch2))
write.csv(res,"C_I_II.csv", row.names=T, quote = F)

# DMR
myAnnotation_CC <- DMRcate::cpg.annotate(object = as.matrix(M_batch1_batch2_stage), datatype = "array", what = "M",
                                         analysis.type = "differential", design = as.matrix(mymod_stage_cc),
                                         contrasts = TRUE, cont.matrix = contrasts_case_control_stage,
                                         coef = "Endometriosis_Yes_No_Yes", arraytype = "EPIC", fdr = 0.1)
DMRs_CC <- dmrcate(myAnnotation_CC, lambda=1000, C=2,pcutoff = 0.05)
results.ranges_CC <- extractRanges(DMRs_CC)
results.ranges_CC_select<-results.ranges_CC[results.ranges_CC$Fisher<0.1,]
write.csv(results.ranges_CC_select,"DMRcate_results_select_C_I_II.csv")

## III_IV vs Control
index <- which(cov$Stage_grouped=="CONTROL"|cov$Stage_grouped=="III_IV")
cov_C_III_IV <- cov[index,]
M_batch1_batch2_C_III_IV <- M_batch1_batch2[,index]

mod_stage <- model.matrix(~0+Stage_grouped + Cycle_phase_for_Analysis, data=cov_C_III_IV)
mod_0_stage <- model.matrix(~1, data=cov_C_III_IV)

Y.r <- t(resid(lm(t(M_batch1_batch2_C_III_IV) ~ Stage_grouped + Cycle_phase_for_Analysis, data=cov_C_III_IV)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1

svaout_stage_C_III_IV <- smartsva.cpp(as.matrix(M_batch1_batch2_C_III_IV),n.sv=n.sv, mod_stage, mod_0_stage)
svaout_sv_stage_C_III_IV <- svaout_stage_C_III_IV$sv
rownames(svaout_sv_stage_C_III_IV) <- rownames(cov_C_III_IV)
saveRDS(svaout_sv_stage_C_III_IV, "svaout_sv_stage_C_III_IV.rds")
write.table(svaout_sv_stage_C_III_IV,"svaout_sv_stage_C_III_IV.txt", quote=F, row.names=T, sep="\t")

svaout_sv_stage <- readRDS("svaout_sv_stage_C_III_IV.rds")

index <- which(cov$Stage_grouped=="CONTROL"|cov$Stage_grouped=="III_IV")
cov_stage <- cov[index,]
M_batch1_batch2_stage <- M_batch1_batch2[,index]

regressme_stage <- data.frame(cov_stage$Cycle_phase_for_Analysis,cov_stage$Endometriosis_Yes_No_,svaout_sv_stage)
colnames(regressme_stage)[1] <- "Cycle_phase_for_Analysis"
colnames(regressme_stage)[2] <- "Endometriosis_Yes_No_"

mymod_stage_cc <- data.frame(model.matrix(~Cycle_phase_for_Analysis+Endometriosis_Yes_No_+., regressme_stage))
fit_stage_cc <- lmFit(M_batch1_batch2_stage, mymod_stage_cc)

contrasts_case_control_stage <- makeContrasts(Endometriosis_Yes_No_Yes,levels = mymod_stage_cc)

fit_case_control_stage <- contrasts.fit(fit_stage_cc, contrasts_case_control_stage)
fit_case_control_stage <- eBayes(fit_case_control_stage)

summary(decideTests(fit_case_control_stage),adjust.method="BH",p.value=0.05)
saveRDS(fit_case_control_stage, "fit_C_III_IV.rds")

plot_qq_manhattan(paste("CaseControl_M_values_C_III_IV",colnames(fit_case_control_stage)[1], sep=""), fit_case_control_stage[,colnames(fit_case_control_stage)[1]])

fit_case_control_stage <- readRDS("fit_C_III_IV.rds")
res <- topTable(fit_case_control_stage[,colnames(fit_case_control_stage)[1]],number = nrow(M_batch1_batch2))
write.csv(res,"C_III_IV.csv", row.names=T, quote = F)
index <- which(rownames(annEpic) %in% rownames(res))
head(annEpic)
res$CpG <- rownames(res)
res2 <- as.data.frame(merge(res,annEpic,by.x="CpG",by.y="Name"))
res2 <- res2[order(res2$P.Value),]
write.csv(res2[res2$adj.P.Val<0.05,],"C_III_IV_FDR.csv", row.names=T, quote = F)

index <- which(rownames(M_batch1_batch2_stage)=="cg02011723")
cg02011723 <- as.data.frame(t(M_batch1_batch2_stage[index,]))
plot(as.factor(cov_stage$Endometriosis_Yes_No_), as.numeric(cg02011723$cg02011723))
index <- which(rownames(M_batch1_batch2_stage)=="cg02623400")
cg02623400 <- as.data.frame(t(M_batch1_batch2_stage[index,]))
plot(as.factor(cov_stage$Endometriosis_Yes_No_), as.numeric(cg02623400$cg02623400))

# DMR
myAnnotation_CC <- DMRcate::cpg.annotate(object = as.matrix(M_batch1_batch2_stage), datatype = "array", what = "M",
                                         analysis.type = "differential", design = as.matrix(mymod_stage_cc),
                                         contrasts = TRUE, cont.matrix = contrasts_case_control_stage,
                                         coef = "Endometriosis_Yes_No_Yes", arraytype = "EPIC", fdr = 0.1)
DMRs_CC <- dmrcate(myAnnotation_CC, lambda=1000, C=2,pcutoff = 0.05)
results.ranges_CC <- extractRanges(DMRs_CC)
results.ranges_CC_select<-results.ranges_CC[results.ranges_CC$Fisher<0.1,]
write.csv(results.ranges_CC_select,"DMRcate_results_select_C_III_IV.csv")

## I-II vs III_IV
index <- which(cov$Stage_grouped=="I_II"|cov$Stage_grouped=="III_IV")
cov_C_I_II_III_IV <- cov[index,]
M_batch1_batch2_C_I_II_III_IV <- M_batch1_batch2[,index]

mod_stage <- model.matrix(~0+Stage_grouped + Cycle_phase_for_Analysis, data=cov_C_I_II_III_IV)
mod_0_stage <- model.matrix(~1, data=cov_C_I_II_III_IV)

Y.r <- t(resid(lm(t(M_batch1_batch2_C_I_II_III_IV) ~ Stage_grouped + Cycle_phase_for_Analysis, data=cov_C_I_II_III_IV)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1

svaout_stage_C_I_II_III_IV <- smartsva.cpp(as.matrix(M_batch1_batch2_C_I_II_III_IV),n.sv=n.sv, mod_stage, mod_0_stage)
svaout_sv_stage_C_I_II_III_IV <- svaout_stage_C_I_II_III_IV$sv
rownames(svaout_sv_stage_C_I_II_III_IV) <- rownames(cov_C_I_II_III_IV)
saveRDS(svaout_sv_stage_C_I_II_III_IV, "svaout_sv_stage_I_II_III_IV.rds")
write.table(svaout_sv_stage_C_I_II_III_IV,"svaout_sv_stage_I_II_III_IV.txt", quote=F, row.names=T, sep="\t")

svaout_sv_stage <- readRDS("svaout_sv_stage_I_II_III_IV.rds")
index <- which(cov$Stage_grouped=="I_II"|cov$Stage_grouped=="III_IV")
cov_stage <- cov[index,]
M_batch1_batch2_stage <- M_batch1_batch2[,index]

regressme_stage <- data.frame(cov_stage$Cycle_phase_for_Analysis,cov_stage$Stage_grouped,svaout_sv_stage)
colnames(regressme_stage)[1] <- "Cycle_phase_for_Analysis"
colnames(regressme_stage)[2] <- "Endometriosis_Stage"

mymod_stage_cc <- data.frame(model.matrix(~Cycle_phase_for_Analysis+Endometriosis_Stage+., regressme_stage))
fit_stage_cc <- lmFit(M_batch1_batch2_stage, mymod_stage_cc)

contrasts_case_control_stage <- makeContrasts(Endometriosis_StageIII_IV,levels = mymod_stage_cc)

fit_case_control_stage <- contrasts.fit(fit_stage_cc, contrasts_case_control_stage)
fit_case_control_stage <- eBayes(fit_case_control_stage)

summary(decideTests(fit_case_control_stage),adjust.method="BH",p.value=0.05)
saveRDS(fit_case_control_stage, "fit_I_II_III_IV.rds")

plot_qq_manhattan(paste("CaseControl_M_values_I_II_III_IV",colnames(fit_case_control_stage)[1], sep=""), fit_case_control_stage[,colnames(fit_case_control_stage)[1]])
fit_case_control_stage <- readRDS("fit_I_II_III_IV.rds")
res <- topTable(fit_case_control_stage[,colnames(fit_case_control_stage)[1]],number = nrow(M_batch1_batch2))
write.csv(res,"I_II_III_IV.csv", row.names=T, quote = F)

# DMR
myAnnotation_CC <- DMRcate::cpg.annotate(object = as.matrix(M_batch1_batch2_stage), datatype = "array", what = "M",
                                         analysis.type = "differential", design = as.matrix(mymod_stage_cc),
                                         contrasts = TRUE, cont.matrix = contrasts_case_control_stage,
                                         coef = "Endometriosis_Yes_No_Yes", arraytype = "EPIC", fdr = 0.1)
DMRs_CC <- dmrcate(myAnnotation_CC, lambda=1000, C=2,pcutoff = 0.05)
results.ranges_CC <- extractRanges(DMRs_CC)
results.ranges_CC_select<-results.ranges_CC[results.ranges_CC$Fisher<0.1,]
write.csv(results.ranges_CC_select,"DMRcate_results_select_I_II_III_IV.csv")


##### PRISTINE CONTROL ANALYSIS #####
head(nupp)
index <- which(annotation_file$Sample.ID %in% nupp$sampleid)
nupp2 <- as.data.frame(annotation_file$Epic_Complete.Bar.code[index])
index <- which(rownames(cov) %in% nupp2$`annotation_file$Epic_Complete.Bar.code[index]` | cov$Endometriosis_Yes_No_=="Yes")
cov_pristine <- cov[index,]
M_batch1_batch2_pristine <- M_batch1_batch2[,index]

PCs_pristine <- prcomp(t(M_batch1_batch2_pristine),center=F, scale.=F)
plot_pca("PCA_presva_pristine.pdf", PCs_pristine, cov_pristine)

pdf("corrplot_presva_pristine.pdf", w = 10, h = 10)
plot_corr(PCs_pristine, cov_pristine[,-6])
dev.off()

mod_pristine <- model.matrix(~0+Endometriosis_Yes_No_ + Cycle_phase_for_Analysis, data=cov_pristine)
mod_0_pristine <- model.matrix(~1, data=cov_pristine)

Y.r_pristine <- t(resid(lm(t(M_batch1_batch2_pristine) ~ Endometriosis_Yes_No_ + Cycle_phase_for_Analysis, data=cov_pristine)))
n.sv_pristine <- EstDimRMT(Y.r_pristine, FALSE)$dim + 1 

svaout_pristine <- smartsva.cpp(as.matrix(M_batch1_batch2_pristine),n.sv=n.sv_pristine, mod_pristine, mod_0_pristine)
svaout_sv_pristine <- svaout_pristine$sv
svaout_sv_pristine <- readRDS("svaout_sv_pristine.rds")
rownames(svaout_sv_pristine)<-colnames(M_batch1_batch2_pristine)

# Corrected Data PCA and Correlation Plot
M_batch1_batch2_pristine_corrected <- get_corrected_data(as.matrix(M_batch1_batch2_pristine),mod_pristine,svaout_pristine)

PCs_corrected_pristine <- prcomp(t(M_batch1_batch2_pristine_corrected),center=F, scale.=F)
plot_pca("PCA_corrected_pristine_postsva.pdf", PCs_corrected_pristine, cov_pristine)

pdf("corrplot_corrected_pristine_postsva.pdf", w = 10, h = 10)
plot_corr(PCs_corrected_pristine, cov_pristine[,-6])
dev.off()

regressme_pristine <- data.frame(cov_pristine$Cycle_phase_for_Analysis,cov_pristine$Endometriosis_Yes_No_,svaout_sv_pristine) 
colnames(regressme_pristine)[1] <- "Cycle_phase_for_Analysis"
colnames(regressme_pristine)[2] <- "Endometriosis_Yes_No_"

caco<-factor(cov_pristine$Endometriosis_Yes_No_)
table(caco)
cyclephase<-factor(cov_pristine$Cycle_phase_for_Analysis)
table(cyclephase)

mymod_pristine <- model.matrix(~0+Cycle_phase_for_Analysis+Endometriosis_Yes_No_+., regressme_pristine)
fit_pristine <- lmFit(M_batch1_batch2_pristine,mymod_pristine)

contrasts_cyclephase_pristine<- makeContrasts(Cycle_phase_for_AnalysisESE-Cycle_phase_for_AnalysisLSE, Cycle_phase_for_AnalysisESE-Cycle_phase_for_AnalysisMSE, Cycle_phase_for_AnalysisESE-Cycle_phase_for_AnalysisPE, Cycle_phase_for_AnalysisLSE-Cycle_phase_for_AnalysisMSE,Cycle_phase_for_AnalysisLSE-Cycle_phase_for_AnalysisPE, Cycle_phase_for_AnalysisMSE-Cycle_phase_for_AnalysisPE, Cycle_phase_for_AnalysisMSE-Cycle_phase_for_AnalysisMenstrual, Cycle_phase_for_AnalysisESE-Cycle_phase_for_AnalysisMenstrual, Cycle_phase_for_AnalysisLSE-Cycle_phase_for_AnalysisMenstrual, Cycle_phase_for_AnalysisPE-Cycle_phase_for_AnalysisMenstrual, levels = mymod_pristine)

fit_cycle_phase_pristine<-contrasts.fit(fit_pristine, contrasts_cyclephase_pristine)
fit_cycle_phase_pristine<- eBayes(fit_cycle_phase_pristine)

summary(decideTests(fit_cycle_phase_pristine),adjust.method="BH",p.value=0.05)

mymod_pristine <- model.matrix(~Cycle_phase_for_Analysis+Endometriosis_Yes_No_+., regressme_pristine)
fit_pristine <- lmFit(M_batch1_batch2_pristine,mymod_pristine)

contrasts_case_control_pristine<- makeContrasts(Endometriosis_Yes_No_Yes,levels = mymod_pristine)

fit_case_control_pristine<-contrasts.fit(fit_pristine, contrasts_case_control_pristine)
head(fit_pristine$coefficients)
fit_case_control_pristine<- eBayes(fit_case_control_pristine)

summary(decideTests(fit_case_control_pristine),adjust.method="BH",p.value=0.05)
saveRDS(fit_case_control_pristine, "fit_case_control_pristine.rds")
fit_case_control_pristine <- readRDS("fit_case_control_pristine.rds")

## CASE CONTROL ##
plot_qq_manhattan(paste("CaseControl_M_values_pristine",colnames(fit_case_control_pristine)[1], sep=""), fit_case_control_pristine[,colnames(fit_case_control_pristine)[1]])

# Single Site
CC_pristine <- topTable(fit_case_control_pristine[,colnames(fit_case_control_pristine)[1]],number = nrow(M_batch1_batch2_pristine))
write.csv(CC_pristine,"CC_pristine.csv", row.names=T, quote = F)

CC_pristine <- rownames(CC_pristine)[CC_pristine$adj.P.Val < 0.05]
plot_heatmap(name = "Case_Control_pristine", data = M_batch1_batch2_pristine,cov = cov_pristine,probes = CC_pristine,samples = rownames(cov_pristine))
CC_pristine <- rownames(CC_all_pristine)[CC_all_pristine$P.Value < 0.05]
write.csv(CC_pristine[CC_pristine$adj.P.Val < 0.05,], "CC_pristine_FDR.csv")

# DMR
myAnnotation_CC <- DMRcate::cpg.annotate(object = as.matrix(M_batch1_batch2_pristine), datatype = "array", what = "M",
                                         analysis.type = "differential", design = as.matrix(mymod_pristine),
                                         contrasts = TRUE, cont.matrix = contrasts_case_control_pristine,
                                         coef = "Endometriosis_Yes_No_Yes", arraytype = "EPIC", fdr = 0.1)
DMRs_CC <- dmrcate(myAnnotation_CC, lambda=1000, C=2,pcutoff = 0.05)
results.ranges_CC <- extractRanges(DMRs_CC)
results.ranges_CC_select<-results.ranges_CC[results.ranges_CC$Fisher<0.1,]
write.csv(results.ranges_CC_select,"DMRcate_results_select_CC_pristine.csv")

## By disease stage ##

# Pristine control vs I_II

index <- which(annotation_file$Sample.ID %in% nupp$sampleid)
nupp2 <- as.data.frame(annotation_file$Epic_Complete.Bar.code[index])
index <- which(rownames(cov) %in% nupp2$`annotation_file$Epic_Complete.Bar.code[index]` | cov$Endometriosis_Yes_No_=="Yes")
cov_pristine <- cov[index,]
M_batch1_batch2_pristine <- M_batch1_batch2[,index]

index <- which(cov_pristine$Stage_grouped=="CONTROL"|cov_pristine$Stage_grouped=="I_II")
cov_pristine_C_I_II <- cov_pristine[index,]
M_batch1_batch2_C_I_II <- M_batch1_batch2_pristine[,index]

mod_stage <- model.matrix(~0+Stage_grouped + Cycle_phase_for_Analysis, data=cov_pristine_C_I_II)
mod_0_stage <- model.matrix(~1, data=cov_pristine_C_I_II)

Y.r <- t(resid(lm(t(M_batch1_batch2_C_I_II) ~ Stage_grouped + Cycle_phase_for_Analysis, data=cov_pristine_C_I_II)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1

svaout_stage_C_I_II <- smartsva.cpp(as.matrix(M_batch1_batch2_C_I_II),n.sv=n.sv, mod_stage, mod_0_stage)
svaout_sv_stage_C_I_II <- svaout_stage_C_I_II$sv
rownames(svaout_sv_stage_C_I_II) <- rownames(cov_pristine_C_I_II)
saveRDS(svaout_sv_stage_C_I_II, "svaout_sv_stage_pristineC_I_II.rds")
write.table(svaout_sv_stage_C_I_II,"svaout_sv_stage_pristineC_I_II.txt", quote=F, row.names=T, sep="\t")

index <- which(cov_pristine$Stage_grouped=="CONTROL"|cov_pristine$Stage_grouped=="I_II")
cov_pristine <- cov_pristine[index,]
M_batch1_batch2_pristine <- M_batch1_batch2_pristine[,index]

svaout_sv_pristine <- readRDS("svaout_sv_stage_pristineC_I_II.rds")
regressme_stage_pristine <- data.frame(cov_pristine$Cycle_phase_for_Analysis,cov_pristine$Stage_grouped,svaout_sv_pristine)
colnames(regressme_stage_pristine)[2] <- "Endometriosis_Yes_No_"
colnames(regressme_stage_pristine)[1] <- "Cycle_phase_for_Analysis"

mymod_pristine <- model.matrix(~Cycle_phase_for_Analysis+Endometriosis_Yes_No_+., regressme_stage_pristine)
fit_pristine <- lmFit(M_batch1_batch2_pristine,mymod_pristine)

contrasts_case_control_pristine<- makeContrasts(Endometriosis_Yes_No_I_II,levels = mymod_pristine)

fit_case_control_pristine<-contrasts.fit(fit_pristine, contrasts_case_control_pristine)
fit_case_control_pristine<- eBayes(fit_case_control_pristine)
summary(decideTests(fit_case_control_pristine),adjust.method="BH",p.value=0.05)
saveRDS(fit_case_control_pristine, "fit_I_II_PristineControl.rds")

plot_qq_manhattan(paste("I_II_Control_M_values_pristine",colnames(fit_case_control_pristine)[1], sep=""), fit_case_control_pristine[,colnames(fit_case_control_pristine)[1]])

# Single Site
CC_pristine <- topTable(fit_case_control_pristine[,colnames(fit_case_control_pristine)[1]],number = nrow(M_batch1_batch2_pristine))
write.csv(CC_pristine,"C_I_II_pristine.csv", row.names=T, quote = F)

# DMR
myAnnotation_CC <- DMRcate::cpg.annotate(object = as.matrix(M_batch1_batch2_pristine), datatype = "array", what = "M",
                                         analysis.type = "differential", design = as.matrix(mymod_pristine),
                                         contrasts = TRUE, cont.matrix = contrasts_case_control_pristine,
                                         coef = "Endometriosis_Yes_No_I_II", arraytype = "EPIC", fdr = 0.1)
DMRs_CC <- dmrcate(myAnnotation_CC, lambda=1000, C=2,pcutoff = 0.05)
results.ranges_CC <- extractRanges(DMRs_CC)
results.ranges_CC_select<-results.ranges_CC[results.ranges_CC$Fisher<0.1,]
write.csv(results.ranges_CC_select,"DMRcate_results_select_I_II_pristine.csv")

# Pristine control vs III_IV

index <- which(annotation_file$Sample.ID %in% nupp$sampleid)
nupp2 <- as.data.frame(annotation_file$Epic_Complete.Bar.code[index])
index <- which(rownames(cov) %in% nupp2$`annotation_file$Epic_Complete.Bar.code[index]` | cov$Endometriosis_Yes_No_=="Yes")
cov_pristine <- cov[index,]
M_batch1_batch2_pristine <- M_batch1_batch2[,index]

index <- which(cov_pristine$Stage_grouped=="CONTROL"|cov_pristine$Stage_grouped=="III_IV")
cov_pristine_C_III_IV <- cov_pristine[index,]
M_batch1_batch2_C_III_IV <- M_batch1_batch2_pristine[,index]

mod_stage <- model.matrix(~0+Stage_grouped + Cycle_phase_for_Analysis, data=cov_pristine_C_III_IV)
mod_0_stage <- model.matrix(~1, data=cov_pristine_C_III_IV)

Y.r <- t(resid(lm(t(M_batch1_batch2_C_III_IV) ~ Stage_grouped + Cycle_phase_for_Analysis, data=cov_pristine_C_III_IV)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1

svaout_stage_C_III_IV <- smartsva.cpp(as.matrix(M_batch1_batch2_C_III_IV),n.sv=n.sv, mod_stage, mod_0_stage)
svaout_sv_stage_C_III_IV <- svaout_stage_C_III_IV$sv
rownames(svaout_sv_stage_C_III_IV) <- rownames(cov_pristine_C_III_IV)
saveRDS(svaout_sv_stage_C_III_IV, "svaout_sv_stage_pristineC_III_IV.rds")
write.table(svaout_sv_stage_C_III_IV,"svaout_sv_stage_pristineC_III_IV.txt", quote=F, row.names=T, sep="\t")

index <- which(cov_pristine$Stage_grouped=="CONTROL"|cov_pristine$Stage_grouped=="III_IV")
cov_pristine <- cov_pristine[index,]
M_batch1_batch2_pristine <- M_batch1_batch2_pristine[,index]

svaout_sv_pristine <- readRDS("svaout_sv_stage_pristineC_III_IV.rds")
regressme_stage_pristine <- data.frame(cov_pristine$Cycle_phase_for_Analysis,cov_pristine$Stage_grouped,svaout_sv_pristine)
colnames(regressme_stage_pristine)[2] <- "Endometriosis_Yes_No_"
colnames(regressme_stage_pristine)[1] <- "Cycle_phase_for_Analysis"

mymod_pristine <- model.matrix(~Cycle_phase_for_Analysis+Endometriosis_Yes_No_+., regressme_stage_pristine)
fit_pristine <- lmFit(M_batch1_batch2_pristine,mymod_pristine)

contrasts_case_control_pristine<- makeContrasts(Endometriosis_Yes_No_III_IV,levels = mymod_pristine)

fit_case_control_pristine<-contrasts.fit(fit_pristine, contrasts_case_control_pristine)
fit_case_control_pristine<- eBayes(fit_case_control_pristine)
summary(decideTests(fit_case_control_pristine),adjust.method="BH",p.value=0.05)
saveRDS(fit_case_control_pristine, "fit_III_IV_PristineControl.rds")

plot_qq_manhattan(paste("III_IV_Control_M_values_pristine",colnames(fit_case_control_pristine)[1], sep=""), fit_case_control_pristine[,colnames(fit_case_control_pristine)[1]])

# Single Site
CC_pristine <- topTable(fit_case_control_pristine[,colnames(fit_case_control_pristine)[1]],number = nrow(M_batch1_batch2_pristine))
write.csv(CC_pristine,"C_III_IV_pristine.csv", row.names=T, quote = F)

CC_pristine$CpG <- rownames(CC_pristine)
res2 <- as.data.frame(merge(CC_pristine,annEpic,by.x="CpG",by.y="Name"))
res2 <- res2[order(res2$P.Value),]
write.csv(res2[res2$adj.P.Val<0.05,],"C_III_IV_pristine_FDR.csv", row.names=T, quote = F)
write.csv(res2,"C_III_IV_pristine_anno.csv", row.names=T, quote = F)

index <- which(rownames(M_batch1_batch2_pristine)=="cg18305031")
cg18305031 <- as.data.frame(t(M_batch1_batch2_pristine[index,]))
plot(as.factor(cov_pristine$Endometriosis_Yes_No_), as.numeric(cg18305031$cg18305031))

# DMR
myAnnotation_CC <- DMRcate::cpg.annotate(object = as.matrix(M_batch1_batch2_pristine), datatype = "array", what = "M",
                                         analysis.type = "differential", design = as.matrix(mymod_pristine),
                                         contrasts = TRUE, cont.matrix = contrasts_case_control_pristine,
                                         coef = "Endometriosis_Yes_No_III_IV", arraytype = "EPIC", fdr = 0.1)
DMRs_CC <- dmrcate(myAnnotation_CC, lambda=1000, C=2,pcutoff = 0.05)
results.ranges_CC <- extractRanges(DMRs_CC)
results.ranges_CC_select<-results.ranges_CC[results.ranges_CC$Fisher<0.1,]
write.csv(results.ranges_CC_select,"DMRcate_results_select_III_IV_pristine.csv")

## WGCNA ###

# 50% most variable probes
# see SVA_for_WGCNA_SM scripts

load("for_sv_calc_WGCNA.RData")

# Uncorrected data
bwnet_raw_50_var = blockwiseModules(t(M_batch1_batch2_50_var), maxBlockSize = 40000,
                             power = 12, TOMType = "signed", minModuleSize = 30,
                             reassignThreshold = 0, mergeCutHeight = 0.25,
                             numericLabels = TRUE,
                             saveTOMs = FALSE,
                             verbose = 10)


bwnet <- bwnet_raw_50_var 
saveRDS(bwnet_raw_50_var, "WGCNA_bwnet.rds")
bwnet_colors<-data.frame(bwnet$colors)
colnames(bwnet_colors)<-c("WGCNA_Block")
write.table(bwnet_colors,paste("WGCNA_Blocks.txt", sep = ""))
bwnet_colors[,"Probes"]<-rownames(bwnet_colors)

bwLabels = bwnet$colors

# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)

pdf(paste("blocks_wgcna.pdf", sep = ""))
for(i  in unique(bwnet$blocks[!is.na(bwnet$blocks)])) {
  # Plot the dendrogram and the module colors underneath for blocks 1 and 2
  plotDendroAndColors(bwnet$dendrograms[[i]], bwModuleColors[bwnet$blockGenes[[i]]],
                      "Module colors", main = paste("Dendrogram and module colors in block", i),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}
dev.off()
bwnet <- readRDS("WGCNA_bwnet.rds")

# get locations
locs <- rep(0, 6*length(names(bwnet$MEs)))
dim(locs) <- c(length(names(bwnet$MEs)), 6)
colnames(locs) = c("Island", "N_Shelf", "N_Shore", "OpenSea", "S_Shelf", "S_Shore")
rownames(locs) <- names(bwnet$MEs)

for (i in 1:length(names(bwnet$MEs))) {
  block <-bwnet_colors[bwnet_colors$WGCNA_Block == i-1,2]
  print(paste("Block:",i, "Length", length(block)))
  block_location <- annEpic[block,]
  locs[i,names(table(block_location$Relation_to_Island))] <- table(block_location$Relation_to_Island) 
}
write.csv(locs, paste("module_locations.csv", sep = ""))


analyze_wgcna(bwnet_raw_50_var, M_batch1_batch2_50_var, "raw_50_var", svaout_sv_50_var)

compute_associations_pval_or <- function(MEs, vars, covariates){
  
  pvals <- (matrix(nrow=dim(MEs)[2], ncol = length(vars)))
  rownames(pvals)<-colnames(MEs)
  colnames(pvals)<-vars
  
  odds_ratio <- (matrix(nrow=dim(MEs)[2], ncol = length(vars)))
  rownames(odds_ratio)<-colnames(MEs)
  colnames(odds_ratio)<-vars
  
    for (j in 1:dim(MEs)[2]){
      
      model <- lm(MEs[,j]~Cycle_phase_for_Analysis+Endometriosis_Yes_No_+., covariates)
      pvals[j,1] <- summary(model)$coefficients[2,4]
      odds_ratio[j,1] <- exp(summary(model)$coefficients[2,1])
      
      model <- lm(MEs[,j]~Endometriosis_Yes_No_+Cycle_phase_for_Analysis+., covariates)
      pvals[j,2] <- summary(model)$coefficients[2,4]
      odds_ratio[j,2] <- exp(summary(model)$coefficients[2,1])
    }

  result <- c()
  result$pvals <- pvals
  result$odds_ratio <- odds_ratio
  return(result)
}

analyze_wgcna <- function(bwnet, data, name, svout) { 

datTraits<-cov[,3:4]
datTraits$Endometriosis_Yes_No_<-as.numeric(as.factor(datTraits$Endometriosis_Yes_No_))
datTraits$Cycle_phase_for_Analysis<-as.numeric(as.factor(datTraits$Cycle_phase_for_Analysis))
colnames(datTraits)<-c("CaseControl","CyclePhase")
moduleTraitCor = cor(bwnet$MEs, datTraits, use= "p")

MEs0 = moduleEigengenes(t(data), bwLabels,nPC = 2)
MEs = MEs0$eigengenes
rownames(MEs)<-rownames(t(data))

regressme <- data.frame(cov[,4], cov[,3],svout)
colnames(regressme)[1] <- "Cycle_phase_for_Analysis"
colnames(regressme)[2] <- "Endometriosis_Yes_No_"

compute_associations_pval_or(MEs, colnames(datTraits), regressme)

#Correlation without covariates
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datTraits))

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)

#display the corelation values with a heatmap plot
pdf(paste("Module_trait_relationships_", name, ".pdf", sep = ""))
par(mar= c(6, 8.5, 3, 3))
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= names(bwnet$MEs),
               ySymbols= names(bwnet$MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.5,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))
dev.off() 

# get locations
locs <- rep(0, 6*length(names(bwnet$MEs)))
dim(locs) <- c(length(names(bwnet$MEs)), 6)
colnames(locs) = c("Island", "N_Shelf", "N_Shore", "OpenSea", "S_Shelf", "S_Shore")
rownames(locs) <- names(bwnet$MEs)

#get pathways and locations
for (i in 1:length(names(bwnet$MEs))) {
  block <-bwnet_colors[bwnet_colors$WGCNA_Block == i-1,2]
  print(paste("Block:",i, "Length", length(block)))
  gst_block <- gometh(sig.cpg=block, all.cpg=rownames(data), collection="KEGG")
  #gst_block <- gst_block[gst_block$FDR<0.1,]
  gst_block <- gst_block[order(gst_block$P.DE,decreasing = F),]
  write.csv(gst_block, paste("me_", i-1, "_pathways_", name, ".csv", sep = ""))
  block_location <- annEpic[block,]
  locs[i,names(table(block_location$Relation_to_Island))] <- table(block_location$Relation_to_Island) 
  #pdf(paste("block_location_", i-1, "_", name, ".pdf", sep = ""))
  #ggplot(data=data.frame(block_location), aes(x=Relation_to_Island)) + geom_bar(width = 0.5,aes(y = (..count..)/sum(..count..))) + theme_bw() 
  #dev.off()
}
write.csv(locs, paste("module_locations_", name, ".csv", sep = ""))
}
