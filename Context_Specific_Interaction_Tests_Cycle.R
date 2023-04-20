################################################################
# This script performs context specific for endometrial mQTLs.

# Author: Sally Mortlock
# July 2022
################################################################

##### Cycle Stage Context Specific Analysis #####

# library
library(snpStats)
library(lme4)

# Read in significant mQTL data
bonf <- readRDS("R01_Euro_cis_mQTL_Endometrium_July2022_Bonf_Anno_1mb_Stage_Endo_SVA.rds")
head(bonf)
bonf2 <- read.table("mQTL_COJO_July2022_V2.txt", header=T)

# Read in methylation data
probe <- read.csv("M_batch1_batch2_July2022_mQTL_bonf.csv", header=T, check.names = F)
probe[1:10,1:10]
sample_info <- read.csv("cov_July2022.csv", header=T)
sva <- readRDS("mQTL_svaout_sv.rds")
sva[1:10,1:10]
sva <- as.data.frame(sva)
probe_info <- read.csv("M_batch1_batch2_July2022_mQTL_pos_bonf.csv", header=T)
probe[1:10,1:10]
rownames(probe) <- probe$id
probe$id <- NULL

probe <- as.data.frame(t(probe))

index <- which(probe_info$geneid %in% colnames(probe))
probe_info <- probe_info[index,]
probe_info <- probe_info[order(probe_info$geneid),]
probe <- probe[,order(colnames(probe))]
index <- which(probe_info$geneid == colnames(probe))
length(index) 
index <- which(sample_info$X %in% rownames(probe))
sample_info <- sample_info[index,]
head(sample_info)
sample_info$Sample.ID <- sample_info$X
sample_info <- sample_info[order(sample_info$Sample.ID),]
index <- which(sample_info$Sample.ID %in% rownames(probe))
length(index)
sample_info <- sample_info[index,]
index <- which(sample_info$Sample.ID == rownames(probe))
length(index) 
index <- which(rownames(sva) %in% rownames(probe))
sva <- sva[order(rownames(sva)),]
index <- which(rownames(sva) == rownames(probe))
length(index) 

variants <- unique(bonf2$SNP)
snplist <- variants
plinkfile <- "R01_May2020_Imputed_R2Filtered_updatedIDs_Euro"

# Read in raw geno data
rawdata <- read.plink(bed=plinkfile, select.snps=snplist)
geno <- apply(rawdata$genotypes@.Data, 2, as.numeric)
geno[1:10,1:10]
geno[geno==0] <- NA
geno <- geno - 1
a1 <- rawdata$map$allele.1
a2 <- rawdata$map$allele.2
rownames(geno) <- rawdata$fam$member
colnames(geno) <- rawdata$map$snp.name
geno <- as.data.frame(geno[order(rownames(geno)),])
index <- which(rownames(geno)==sample_info$Sample.ID)
length(index)
geno[1:10,1:10]

########################
#### Stage of cycle ####
########################

### Apply to all mQTLs using loop ###
setwd("/Context_Specific/")

# filter according to proliferative and secretory
sample_info$PS <- 0
i1 <- which(sample_info$Cycle_phase_for_Analysis == "PE")
sample_info$PS[i1] <- 1 
i2 <- which(sample_info$Cycle_phase_for_Analysis == "SE"|sample_info$Cycle_phase_for_Analysis == "ESE"|sample_info$Cycle_phase_for_Analysis == "MSE"|sample_info$Cycle_phase_for_Analysis == "LSE")
sample_info$PS[i2] <- 2

sample_info2 <- sample_info[which(sample_info$PS!=0),]
index <- which(rownames(geno) %in% sample_info2$Sample.ID)
geno2 <- geno[index,]
probe[1:10,1:10]
index <- which(rownames(probe) %in% sample_info2$Sample.ID)
probe2 <- probe[index,]
index <- which(rownames(sva) %in% sample_info2$Sample.ID)
sva <- sva[index,]

out <- array(0, c(10000,3))
for(i in 1:10000) {
  snplist <- bonf2$SNP[i]
  
  # Read in raw geno data
  y <- as.numeric(probe2[,which(probe_info$geneid==bonf2$CpG[i])])
  Y <- as.data.frame(cbind(as.character(sample_info2$Sample.ID), sample_info2$Endometriosis_Yes_No_,as.character(sample_info2$PS), y,sva[,1:39]))
  names(Y)[1] <- "ID"
  names(Y)[2] <- "Endo"
  names(Y)[3] <- "PS"
  
  Y <- Y[order(Y$ID),]
  
  #combine the coloumn 
  data <- cbind(Y, geno2[,which(colnames(geno2)==snplist)])
  
  # change the number accordingly
  data$geno <- as.factor(data$`geno2[, which(colnames(geno2) == snplist)]`)
  levels(data$geno)[levels(data$geno)=="0"]=paste(a1[which(colnames(geno)==snplist)],a1[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="1"]=paste(a1[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="2"]=paste(a2[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  
  # Linear models testing interaction
  
  mod1 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  mod2 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  
  out[i,1] <- snplist
  out[i,2] <- bonf2$CpG[i]
  
  mod1 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  mod2 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  
  foo <- anova(mod1, mod2)
  out[i,3] <- foo$`Pr(>F)`[2] }

out <- as.data.frame(out)
names(out) <- c("SNP", "Probe_ID", "ANOVA_mod1_mod2_pval")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P1.csv", row.names = F, quote = F)

out <- array(0, c(10000,3))
for(i in 10001:20000) {
  snplist <- bonf2$SNP[i]
  
  # Read in raw geno data
  y <- as.numeric(probe2[,which(probe_info$geneid==bonf2$CpG[i])])
  Y <- as.data.frame(cbind(as.character(sample_info2$Sample.ID), sample_info2$Endometriosis_Yes_No_,as.character(sample_info2$PS), y,sva[,1:39]))
  names(Y)[1] <- "ID"
  names(Y)[2] <- "Endo"
  names(Y)[3] <- "PS"
  
  Y <- Y[order(Y$ID),]
  
  #combine the coloumn 
  data <- cbind(Y, geno2[,which(colnames(geno2)==snplist)])
  

  # change the number accordingly
  data$geno <- as.factor(data$`geno2[, which(colnames(geno2) == snplist)]`)
  levels(data$geno)[levels(data$geno)=="0"]=paste(a1[which(colnames(geno)==snplist)],a1[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="1"]=paste(a1[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="2"]=paste(a2[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  
  # Linear models testing interaction
  
  mod1 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  mod2 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  
  out[i-10000,1] <- snplist
  out[i-10000,2] <- bonf2$CpG[i]
  
  mod1 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  mod2 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  
  foo <- anova(mod1, mod2)
  out[i-10000,3] <- foo$`Pr(>F)`[2] }

out <- as.data.frame(out)
names(out) <- c("SNP", "Probe_ID", "ANOVA_mod1_mod2_pval")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P2.csv", row.names = F, quote = F)

out <- array(0, c(10000,3))
for(i in 20001:30000) {
  snplist <- bonf2$SNP[i]
  
  # Read in raw geno data
  y <- as.numeric(probe2[,which(probe_info$geneid==bonf2$CpG[i])])
  Y <- as.data.frame(cbind(as.character(sample_info2$Sample.ID), sample_info2$Endometriosis_Yes_No_,as.character(sample_info2$PS), y,sva[,1:39]))
  names(Y)[1] <- "ID"
  names(Y)[2] <- "Endo"
  names(Y)[3] <- "PS"
  
  Y <- Y[order(Y$ID),]
  
  #combine the coloumn 
  data <- cbind(Y, geno2[,which(colnames(geno2)==snplist)])
  

  # change the number accordingly
  data$geno <- as.factor(data$`geno2[, which(colnames(geno2) == snplist)]`)
  levels(data$geno)[levels(data$geno)=="0"]=paste(a1[which(colnames(geno)==snplist)],a1[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="1"]=paste(a1[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="2"]=paste(a2[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  
  # Linear models testing interaction
  
  mod1 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  mod2 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  
  out[i-20000,1] <- snplist
  out[i-20000,2] <- bonf2$CpG[i]
  
  mod1 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  mod2 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  
  foo <- anova(mod1, mod2)
  out[i-20000,3] <- foo$`Pr(>F)`[2] }

out <- as.data.frame(out)
names(out) <- c("SNP", "Probe_ID", "ANOVA_mod1_mod2_pval")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P3.csv", row.names = F, quote = F)

out <- array(0, c(10000,3))
for(i in 30001:40000) {
  snplist <- bonf2$SNP[i]
  
  # Read in raw geno data
  y <- as.numeric(probe2[,which(probe_info$geneid==bonf2$CpG[i])])
  Y <- as.data.frame(cbind(as.character(sample_info2$Sample.ID), sample_info2$Endometriosis_Yes_No_,as.character(sample_info2$PS), y,sva[,1:39]))
  names(Y)[1] <- "ID"
  names(Y)[2] <- "Endo"
  names(Y)[3] <- "PS"
  
  Y <- Y[order(Y$ID),]
  
  #combine the coloumn 
  data <- cbind(Y, geno2[,which(colnames(geno2)==snplist)])
  

  # change the number accordingly
  data$geno <- as.factor(data$`geno2[, which(colnames(geno2) == snplist)]`)
  levels(data$geno)[levels(data$geno)=="0"]=paste(a1[which(colnames(geno)==snplist)],a1[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="1"]=paste(a1[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="2"]=paste(a2[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  
  # Linear models testing interaction
  
  mod1 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  mod2 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  
  out[i-30000,1] <- snplist
  out[i-30000,2] <- bonf2$CpG[i]
  
  mod1 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  mod2 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  
  foo <- anova(mod1, mod2)
  out[i-30000,3] <- foo$`Pr(>F)`[2] }

out <- as.data.frame(out)
names(out) <- c("SNP", "Probe_ID", "ANOVA_mod1_mod2_pval")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P4.csv", row.names = F, quote = F)

out <- array(0, c(10000,3))
for(i in 40001:50000) {
  snplist <- bonf2$SNP[i]
  
  # Read in raw geno data
  y <- as.numeric(probe2[,which(probe_info$geneid==bonf2$CpG[i])])
  Y <- as.data.frame(cbind(as.character(sample_info2$Sample.ID), sample_info2$Endometriosis_Yes_No_,as.character(sample_info2$PS), y,sva[,1:39]))
  names(Y)[1] <- "ID"
  names(Y)[2] <- "Endo"
  names(Y)[3] <- "PS"
  
  Y <- Y[order(Y$ID),]
  
  #combine the coloumn 
  data <- cbind(Y, geno2[,which(colnames(geno2)==snplist)])
  

  # change the number accordingly
  data$geno <- as.factor(data$`geno2[, which(colnames(geno2) == snplist)]`)
  levels(data$geno)[levels(data$geno)=="0"]=paste(a1[which(colnames(geno)==snplist)],a1[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="1"]=paste(a1[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="2"]=paste(a2[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  
  # Linear models testing interaction
  
  mod1 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  mod2 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  
  out[i-40000,1] <- snplist
  out[i-40000,2] <- bonf2$CpG[i]
  
  mod1 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  mod2 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  
  foo <- anova(mod1, mod2)
  out[i-40000,3] <- foo$`Pr(>F)`[2] }

out <- as.data.frame(out)
names(out) <- c("SNP", "Probe_ID", "ANOVA_mod1_mod2_pval")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P5.csv", row.names = F, quote = F)

out <- array(0, c(10000,3))
for(i in 50001:60000) {
  snplist <- bonf2$SNP[i]
  
  # Read in raw geno data
  y <- as.numeric(probe2[,which(probe_info$geneid==bonf2$CpG[i])])
  Y <- as.data.frame(cbind(as.character(sample_info2$Sample.ID), sample_info2$Endometriosis_Yes_No_,as.character(sample_info2$PS), y,sva[,1:39]))
  names(Y)[1] <- "ID"
  names(Y)[2] <- "Endo"
  names(Y)[3] <- "PS"
  
  Y <- Y[order(Y$ID),]
  
  #combine the coloumn 
  data <- cbind(Y, geno2[,which(colnames(geno2)==snplist)])
  

  # change the number accordingly
  data$geno <- as.factor(data$`geno2[, which(colnames(geno2) == snplist)]`)
  levels(data$geno)[levels(data$geno)=="0"]=paste(a1[which(colnames(geno)==snplist)],a1[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="1"]=paste(a1[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="2"]=paste(a2[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  
  # Linear models testing interaction
  
  mod1 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  mod2 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  
  out[i-50000,1] <- snplist
  out[i-50000,2] <- bonf2$CpG[i]
  
  mod1 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  mod2 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  
  foo <- anova(mod1, mod2)
  out[i-50000,3] <- foo$`Pr(>F)`[2] }

out <- as.data.frame(out)
names(out) <- c("SNP", "Probe_ID", "ANOVA_mod1_mod2_pval")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P6.csv", row.names = F, quote = F)

out <- array(0, c(10000,3))
for(i in 60001:70000) {
  snplist <- bonf2$SNP[i]
  
  # Read in raw geno data
  y <- as.numeric(probe2[,which(probe_info$geneid==bonf2$CpG[i])])
  Y <- as.data.frame(cbind(as.character(sample_info2$Sample.ID), sample_info2$Endometriosis_Yes_No_,as.character(sample_info2$PS), y,sva[,1:39]))
  names(Y)[1] <- "ID"
  names(Y)[2] <- "Endo"
  names(Y)[3] <- "PS"
  
  Y <- Y[order(Y$ID),]
  
  #combine the coloumn 
  data <- cbind(Y, geno2[,which(colnames(geno2)==snplist)])
  

  # change the number accordingly
  data$geno <- as.factor(data$`geno2[, which(colnames(geno2) == snplist)]`)
  levels(data$geno)[levels(data$geno)=="0"]=paste(a1[which(colnames(geno)==snplist)],a1[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="1"]=paste(a1[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="2"]=paste(a2[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  
  # Linear models testing interaction
  
  mod1 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  mod2 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  
  out[i-60000,1] <- snplist
  out[i-60000,2] <- bonf2$CpG[i]
  
  mod1 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  mod2 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  
  foo <- anova(mod1, mod2)
  out[i-60000,3] <- foo$`Pr(>F)`[2] }

out <- as.data.frame(out)
names(out) <- c("SNP", "Probe_ID", "ANOVA_mod1_mod2_pval")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P7.csv", row.names = F, quote = F)

out <- array(0, c(10000,3))
for(i in 70001:80000) {
  snplist <- bonf2$SNP[i]
  
  # Read in raw geno data
  y <- as.numeric(probe2[,which(probe_info$geneid==bonf2$CpG[i])])
  Y <- as.data.frame(cbind(as.character(sample_info2$Sample.ID), sample_info2$Endometriosis_Yes_No_,as.character(sample_info2$PS), y,sva[,1:39]))
  names(Y)[1] <- "ID"
  names(Y)[2] <- "Endo"
  names(Y)[3] <- "PS"
  
  Y <- Y[order(Y$ID),]
  
  #combine the coloumn 
  data <- cbind(Y, geno2[,which(colnames(geno2)==snplist)])
  

  # change the number accordingly
  data$geno <- as.factor(data$`geno2[, which(colnames(geno2) == snplist)]`)
  levels(data$geno)[levels(data$geno)=="0"]=paste(a1[which(colnames(geno)==snplist)],a1[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="1"]=paste(a1[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="2"]=paste(a2[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  
  # Linear models testing interaction
  
  mod1 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  mod2 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  
  out[i-70000,1] <- snplist
  out[i-70000,2] <- bonf2$CpG[i]
  
  mod1 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  mod2 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  
  foo <- anova(mod1, mod2)
  out[i-70000,3] <- foo$`Pr(>F)`[2] }

out <- as.data.frame(out)
names(out) <- c("SNP", "Probe_ID", "ANOVA_mod1_mod2_pval")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P8.csv", row.names = F, quote = F)

out <- array(0, c(10000,3))
for(i in 80001:90000) {
  snplist <- bonf2$SNP[i]
  
  # Read in raw geno data
  y <- as.numeric(probe2[,which(probe_info$geneid==bonf2$CpG[i])])
  Y <- as.data.frame(cbind(as.character(sample_info2$Sample.ID), sample_info2$Endometriosis_Yes_No_,as.character(sample_info2$PS), y,sva[,1:39]))
  names(Y)[1] <- "ID"
  names(Y)[2] <- "Endo"
  names(Y)[3] <- "PS"
  
  Y <- Y[order(Y$ID),]
  
  #combine the coloumn 
  data <- cbind(Y, geno2[,which(colnames(geno2)==snplist)])
  

  # change the number accordingly
  data$geno <- as.factor(data$`geno2[, which(colnames(geno2) == snplist)]`)
  levels(data$geno)[levels(data$geno)=="0"]=paste(a1[which(colnames(geno)==snplist)],a1[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="1"]=paste(a1[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="2"]=paste(a2[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  
  # Linear models testing interaction
  
  mod1 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  mod2 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  
  out[i-80000,1] <- snplist
  out[i-80000,2] <- bonf2$CpG[i]
  
  mod1 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  mod2 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  
  foo <- anova(mod1, mod2)
  out[i-80000,3] <- foo$`Pr(>F)`[2] }

out <- as.data.frame(out)
names(out) <- c("SNP", "Probe_ID", "ANOVA_mod1_mod2_pval")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P9.csv", row.names = F, quote = F)

out <- array(0, c(10000,3))
for(i in 90001:10000) {
  snplist <- bonf2$SNP[i]
  
  # Read in raw geno data
  y <- as.numeric(probe2[,which(probe_info$geneid==bonf2$CpG[i])])
  Y <- as.data.frame(cbind(as.character(sample_info2$Sample.ID), sample_info2$Endometriosis_Yes_No_,as.character(sample_info2$PS), y,sva[,1:39]))
  names(Y)[1] <- "ID"
  names(Y)[2] <- "Endo"
  names(Y)[3] <- "PS"
  
  Y <- Y[order(Y$ID),]
  
  #combine the coloumn 
  data <- cbind(Y, geno2[,which(colnames(geno2)==snplist)])
  

  # change the number accordingly
  data$geno <- as.factor(data$`geno2[, which(colnames(geno2) == snplist)]`)
  levels(data$geno)[levels(data$geno)=="0"]=paste(a1[which(colnames(geno)==snplist)],a1[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="1"]=paste(a1[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="2"]=paste(a2[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  
  # Linear models testing interaction
  
  mod1 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  mod2 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  
  out[i-90000,1] <- snplist
  out[i-90000,2] <- bonf2$CpG[i]
  
  mod1 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  mod2 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  
  foo <- anova(mod1, mod2)
  out[i-90000,3] <- foo$`Pr(>F)`[2] }

out <- as.data.frame(out)
names(out) <- c("SNP", "Probe_ID", "ANOVA_mod1_mod2_pval")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P10.csv", row.names = F, quote = F)

out <- array(0, c(18185,3))
for(i in 100001:118185) {
  snplist <- bonf2$SNP[i]
  
  # Read in raw geno data
  y <- as.numeric(probe2[,which(probe_info$geneid==bonf2$CpG[i])])
  Y <- as.data.frame(cbind(as.character(sample_info2$Sample.ID), sample_info2$Endometriosis_Yes_No_,as.character(sample_info2$PS), y,sva[,1:39]))
  names(Y)[1] <- "ID"
  names(Y)[2] <- "Endo"
  names(Y)[3] <- "PS"
  
  Y <- Y[order(Y$ID),]
  
  #combine the coloumn 
  data <- cbind(Y, geno2[,which(colnames(geno2)==snplist)])
  

  # change the number accordingly
  data$geno <- as.factor(data$`geno2[, which(colnames(geno2) == snplist)]`)
  levels(data$geno)[levels(data$geno)=="0"]=paste(a1[which(colnames(geno)==snplist)],a1[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="1"]=paste(a1[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  levels(data$geno)[levels(data$geno)=="2"]=paste(a2[which(colnames(geno)==snplist)],a2[which(colnames(geno)==snplist)],sep="")
  
  # Linear models testing interaction
  
  mod1 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  mod2 <- summary(lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45]))
  
  out[i-100000,1] <- snplist
  out[i-100000,2] <- bonf2$CpG[i]
  
  mod1 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  mod2 <- lm(data$y ~ as.factor(data$geno) + as.factor(data$PS) + as.factor(data$geno):as.factor(data$PS) + as.factor(data$Endo) + .,data[,7:45])
  
  foo <- anova(mod1, mod2)
  out[i-100000,3] <- foo$`Pr(>F)`[2] }

out <- as.data.frame(out)
names(out) <- c("SNP", "Probe_ID", "ANOVA_mod1_mod2_pval")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P11.csv", row.names = F, quote = F)

out1 <- read.csv("Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P1.csv", header=T)
out2 <- read.csv("Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P2.csv", header=T)
out3 <- read.csv("Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P3.csv", header=T)
out4 <- read.csv("Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P4.csv", header=T)
out5 <- read.csv("Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P5.csv", header=T)
out6 <- read.csv("Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P6.csv", header=T)
out7 <- read.csv("Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P7.csv", header=T)
out8 <- read.csv("Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P8.csv", header=T)
out9 <- read.csv("Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P9.csv", header=T)
out10 <- read.csv("Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P10.csv", header=T)
out11 <- read.csv("Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_P11.csv", header=T)

out <- rbind(out1,out2,out3,out4,out5,out6,out7,out8,out9,out10,out11)
out$FDR <- p.adjust(out$ANOVA_mod1_mod2_pval, method = "BH")
out$Bonf <- p.adjust(out$ANOVA_mod1_mod2_pval, method = "bonferroni")
write.csv(out,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage.csv", row.names = F, quote = F)

write.table(sample_info2[,c(1,1)],"Endo_SampleInfo2.txt", row.names = F, quote = F)
index <- which(sample_info2$PS == 1)
write.table(sample_info2[index,c(1,1)],"Endo_SampleInfo2_P.txt", row.names = F, quote = F)
index <- which(sample_info2$PS==2)
write.table(sample_info2[index,c(1,1)],"Endo_SampleInfo2_S.txt", row.names = F, quote = F)

# In Plink
#./plink_1.9 --bfile R01_May2020_Imputed_R2Filtered_updatedIDs_Euro --freqx --keep Endo_SampleInfo2.txt --out Endo_Cycle_Sample2_Geno_Aug2022_Freq
#./plink_1.9 --bfile R01_May2020_Imputed_R2Filtered_updatedIDs_Euro --freqx --keep Endo_SampleInfo2_P.txt --out Endo_Sample2_P_Geno_Aug2022_Freq
#./plink_1.9 --bfile R01_May2020_Imputed_R2Filtered_updatedIDs_Euro --freqx --keep Endo_SampleInfo2_S.txt --out Endo_Sample2_S_Geno_Aug2022_Freq

genotypes <- read.delim("Endo_Cycle_Sample2_Geno_Aug2022_Freq.frqx", header=T)
genotypes_case <- read.delim("Endo_Sample2_P_Geno_Aug2022_Freq.frqx", header=T)
genotypes_control <- read.delim("Endo_Sample2_S_Geno_Aug2022_Freq.frqx", header=T)

cs_endo2 <- merge(out, genotypes[,1:7], by.x="SNP",by.y="SNP")
cs_endo2 <- merge(cs_endo2, genotypes_case[,c(2,5:7)], by.x="SNP",by.y="SNP")
cs_endo2 <- merge(cs_endo2, genotypes_control[,c(2,5:7)], by.x="SNP",by.y="SNP")
write.csv(cs_endo2,"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_plusFreq.csv", row.names=F, quote = F)
index <- which(cs_endo2$Bonf<0.05 & cs_endo2$C.HOM.A1.>=10 & cs_endo2$C.HOM.A1..y>=10)
index <- which(cs_endo2$Bonf<0.05)
index <- which(cs_endo2$Bonf<0.05 & cs_endo2$C.HOM.A1.>=10 & cs_endo2$C.HOM.A1..y>=10)
write.csv(cs_endo2[index,],"Bonf_1mb_Stage_Endo_SVA_August2022_Context_Specific_CycleStage_plusFreq_Bonf_Filtered.csv", row.names=F, quote = F)
