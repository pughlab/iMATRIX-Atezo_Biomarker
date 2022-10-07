
library(dynamicTreeCut, lib.loc = "~/R")
library(fastcluster, lib.loc = "~/R")
library(WGCNA, lib.loc = "~/R")

treehouse_gene_mat <- vector(mode = "list", length = 7)
mycancers <- c("Lymphoma","NBL","OS","RMS","SARC","WILMS", "EWS")
setLabels <- mycancers

for(i in 1:7){  
    gene_mat <- read.csv(paste0("/cluster/projects/pughlab/projects/INDICATE/wgcna/th/tpm_th_", mycancers[i], ".csv"), 
                         row.names = 1, header = TRUE, stringsAsFactors = F, check.names = F)

    treehouse_gene_mat[[i]] <- list(data = as.data.frame(t(gene_mat)))
    names(treehouse_gene_mat[[i]]$data) <- rownames(gene_mat) 
}


analysisSize <- checkSets(treehouse_gene_mat)
print(analysisSize)

message("keep good genes")
gsg <- goodSamplesGenesMS(treehouse_gene_mat, verbose = 3)
gsg$allOK

if (!gsg$allOK){
    if (sum(!gsg$goodGenes) > 0)
        printFlush(paste("Removing genes:", paste(names(treehouse_gene_mat[[1]]$data)[!gsg$goodGenes],
                                                  collapse = ", ")))
    for (set in 1:analysisSize$nSets){
        if (sum(!gsg$goodSamples[[set]]))
            printFlush(paste("In set", setLabels[set], "removing samples",
                             paste(rownames(treehouse_gene_mat[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        # Remove the offending genes and samples
        treehouse_gene_mat[[set]]$data = treehouse_gene_mat[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes]
    }
# Update analysisSize
analysisSize <- checkSets(treehouse_gene_mat)
}

print(analysisSize)

message("cluster samples and find 75th percentile of the distance to identify outliers:")
sampleTrees <- list()
quantile75th <- list()
for (i in 1:7){
sampleTrees[[i]] <- hclust(dist(treehouse_gene_mat[[i]]$data), method = "average")
quantile75th[[i]] <- quantile(unlist(dist(treehouse_gene_mat[[i]]$data)), 0.75)
}

# Plot the dendrograms including the cut lines

pdf(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/SampleClustering.pdf", width = 20, height = 20)
for (i in 1:7){
plot(sampleTrees[[i]], main = paste("Sample clustering on all genes in", setLabels[i]),
     xlab="", sub="", cex = 0.7)
    abline(h=quantile75th[i], col = "red")
}
dev.off()


for (i in 1:7){
# Find clusters cut by the line
    clust <- cutreeStatic(sampleTrees[[i]], cutHeight = quantile75th[i], minSize = 10)
# Keep the largest one (labeled by the number 1
    print(clust)
    treehouse_gene_mat[[i]]$data <- treehouse_gene_mat[[i]]$data[(clust==1), ]
}

collectGarbage()
# Check the size of the leftover data
analysisSize <- checkSets(treehouse_gene_mat)
print(analysisSize)

nGenes <- analysisSize$nGenes
nSamples <- analysisSize$nSamples

save(treehouse_gene_mat, nGenes, nSamples, setLabels, analysisSize,
     file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/Consensus-dataInput_Treehouse_tpm.RData")
