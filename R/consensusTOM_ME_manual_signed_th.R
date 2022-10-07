library(dynamicTreeCut, lib.loc = "~/R")
library(fastcluster, lib.loc = "~/R")
library(WGCNA, lib.loc = "~/R")

enableWGCNAThreads()

treehouse_gene_mat_names <- load(file = 
                                 "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/Consensus-dataInput_Treehouse_tpm.RData")

nSets <- checkSets(treehouse_gene_mat)$nSets

load(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/TOM_th_manual_signed.RData")


# Scale TOMs
# Define the reference percentile
scaleP <- 0.95
set.seed(777)
# Sample sufficiently large number of TOM entries
nSamples <- as.integer(1/(1-scaleP) * 1000)
scaleSample <- sample(nGenes*(nGenes-1)/2, size = nSamples)

TOMScalingSamples <- list()
# These are TOM values at reference percentile
scaleQuant <- rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers <- rep(1, nSets)

# Loop over sets scale all sets after the first one
for (set in 1:nSets){
# Select the sampled TOM entries
TOMScalingSamples[[set]] <- as.dist(TOM[set, , ])[scaleSample]
# Calculate the 90th percentile
scaleQuant[set] <- quantile(TOMScalingSamples[[set]], probs = scaleP, type = 8, na.rm = TRUE)
# Scale other TOMs
    if (set>1){
        scalePowers[set] <- log(scaleQuant[1])/log(scaleQuant[set])
        TOM[set, ,] <- TOM[set, ,]^scalePowers[set]
    }
}
message("All TOMs scaled")

# Plot scaled and unscaled TOMs
# For plotting, also scale the sampled TOM entries
scaledTOMSamples <- list()

for (set in 1:nSets){
    scaledTOMSamples[[set]] <- TOMScalingSamples[[set]]^scalePowers[set]}

pdf(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/TOMScaling-QQPlot_th_manual_signed1_2.pdf", wi = 6, he = 6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()

pdf(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/TOMScaling-QQPlot_th_manual_signed1_3.pdf", wi = 6, he = 6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[3]], plot.it = TRUE, cex = 0.6,
xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[3]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[3]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()

pdf(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/TOMScaling-QQPlot_th_manual_signed1_4.pdf", wi = 6, he = 6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[4]], plot.it = TRUE, cex = 0.6,
xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[4]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[4]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()

pdf(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/TOMScaling-QQPlot_th_manual_signed1_5.pdf", wi = 6, he = 6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[5]], plot.it = TRUE, cex = 0.6,
xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[5]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[5]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()

pdf(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/TOMScaling-QQPlot_th_manual_signed1_6.pdf", wi = 6, he = 6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[6]], plot.it = TRUE, cex = 0.6,
xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[6]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[6]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()

pdf(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/TOMScaling-QQPlot_th_manual_signed1_7.pdf", wi = 6, he = 6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[7]], plot.it = TRUE, cex = 0.6,
xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[7]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[7]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()











# Calculate consensus TOM - cTOM is large only if genes are large in all datasets
## Done on scaled TOMs
consensusTOM <- pmin(TOM[1, , ], TOM[2, , ], TOM[3, , ],TOM[4, , ],TOM[5, , ],TOM[6, , ],TOM[7, , ])

message("consensus TOM done with pmin")

save(consensusTOM, file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/consensusTOM_manual_th.RData")

# Cluster consensus TOMs
consTree <- hclust(as.dist(1-consensusTOM), method = "average")
minModuleSize <- 20
unmergedLabels <- cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                                deepSplit = 2, cutHeight = 0.99,
                                minClusterSize = minModuleSize,
                                pamRespectsDendro = FALSE)
save(unmergedLabels, file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/unmergedLabels_consensus_th_manual_signed_20.RData")

unmergedColors <- labels2colors(unmergedLabels)
message("table of unmerged modules")
print(table(unmergedLabels))

pdf(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/consTree_unmergedColors_th_manual_signed_20.pdf", wi = 20, he = 20)
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Calculate MEs to find similar modules and merge them together
unmergedMEs <- multiSetMEs(treehouse_gene_mat, colors = NULL, universalColors = unmergedColors)
# Consensus dissimilarity
consMEDiss <- consensusMEDissimilarity(unmergedMEs)
# Cluster consensus modules
consMETree <- hclust(as.dist(consMEDiss), method = "average")
pdf(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/clustering_MEs_th_manual_signed_20.pdf", wi = 20, he = 20)
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col = "red")
dev.off()

# Merge those MEs with more than 75% similarity
merge <- mergeCloseModules(treehouse_gene_mat, unmergedLabels, cutHeight = 0.25, verbose = 3)

# Use these merged MEs 
moduleLabels <- merge$colors
moduleColors <- labels2colors(moduleLabels)

message("table of merged modules")
print(table(moduleLabels))

# Eigengenes of the new merged modules
consMEs <- merge$newMEs

pdf(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/consTree_unmerged_merged_th_manual_signed_20.pdf", wi = 20, he = 20)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

save(consMEs, moduleColors, moduleLabels, consTree, 
     file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/Consensus_Network_TH_manual_signed_20.RData")
