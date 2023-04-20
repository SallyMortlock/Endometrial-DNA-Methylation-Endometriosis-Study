
library(WGCNA)

enableWGCNAThreads()
allowWGCNAThreads()

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
