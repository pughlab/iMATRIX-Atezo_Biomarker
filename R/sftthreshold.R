library(dynamicTreeCut, lib.loc = "~/R")
library(fastcluster, lib.loc = "~/R")
library(WGCNA, lib.loc = "~/R")

enableWGCNAThreads()

treehouse_gene_mat_names <- load(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/Consensus-dataInput_Treehouse_tpm.RData")
nSets <- checkSets(treehouse_gene_mat)$nSets

powers = seq(from = 6, to=30, by=2)
powerTables <- vector(mode = "list", length = nSets)
# Call the network topology analysis function for each set
for (i in 1:nSets){
    powerTables[[i]] <- list(data = pickSoftThreshold(treehouse_gene_mat[[i]]$data, powerVector=powers,
                                                     networkType = "signed", verbose = 3)[[2]])
}

collectGarbage()

colors <- c("black", "red", "blue", "green", "darkgreen", "lightblue", "yellow")
# Will plot these columns of the returned scale free analysis tables
plotCols <- c(2,5,6,7)
colNames <- c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")
# Get the minima and maxima of the plotted points
ylim <- matrix(NA, nrow = 2, ncol = 4)
for (i in 1:nSets){
    for (col in 1:length(plotCols)){
        ylim[1, col] <- min(ylim[1, col], powerTables[[i]]$data[, plotCols[col]], na.rm = TRUE)
        ylim[2, col] <- max(ylim[2, col], powerTables[[i]]$data[, plotCols[col]], na.rm = TRUE)
    }
}


pdf(file = "/cluster/projects/pughlab/projects/INDICATE/wgcna/softthreshold_plots_tpm.pdf", 
    width = 20, height = 20)

for (col in 1:length(plotCols)) for (set in 1:nSets){
    if (set==1){
        plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
             main = colNames[col])
        addGrid()
    }
    if (col==1){
        text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
             labels=powers,col=colors[set])
    } else
        text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
             labels=powers,col=colors[set])
    if (col==1){
        legend("bottomright", legend = setLabels, col = colors, pch = 20) }
    else
        legend("topright", legend = setLabels, col = colors, pch = 20)
}
    
dev.off()

