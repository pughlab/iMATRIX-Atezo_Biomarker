library(dynamicTreeCut, lib.loc = "~/R")
library(fastcluster, lib.loc = "~/R")
library(WGCNA, lib.loc = "~/R")

enableWGCNAThreads()

treehouse_gene_mat_names <- load(file = 
                                 "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/Consensus-dataInput_Treehouse_tpm.RData")

nSets <- checkSets(treehouse_gene_mat)$nSets

# Determine list of powers for each dataset. I decided each p based on model fit (ideally > 0.8)
## This is to ensure we get the optimal set of gene associations within each dataset
listofpowers <- c(22,20,14,16,12,16,18)

# Calculate TOMs
# Initialize an appropriate array to hold the adjacencies and TOMs

adjacencies <- array(0, dim = c(nSets, nGenes, nGenes))
TOM <- array(0, dim = c(nSets, nGenes, nGenes))

# TOM type is unsigned because network type is signed and negative correlations are considered unconnected.
## Use pearson because input is TPM
## In the tutorial, they used adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower, which I think is unsigned network (abs removes the negative sign). So I used adjacency function to get signed network
## use = p in adjacency function means pairwise.complete.obs, see documentation

for (set in 1:nSets){
adjacencies[set, , ] <- adjacency(treehouse_gene_mat[[set]]$data, 
				corFnc = "cor", corOptions = list(use = "p", method = "pearson"),
                                type = "signed", power = listofpowers[set])
message("adjacency done for:", print(set))
}


for (set in 1:nSets){
    TOM[set, , ] <- TOMsimilarity(adjacencies[set,,], TOMType = "unsigned")
message("TOM done for:", print(set))
}

message("all TOMs calculated")

save(TOM, file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/th/TOM_th_manual_signed.RData")


