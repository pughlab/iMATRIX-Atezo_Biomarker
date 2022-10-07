###############
# Figure 5
###############

library(survminer)
library(survival)
library(Hmisc)
library(fgsea)
library(ggrepel)

mainpath <- "~/git/iMATRIX-Atezo_Biomarker/"
datapath <- paste0(mainpath,"data/")

dir.create(file.path(paste0(mainpath, "out/")))

plotpath <- paste0(mainpath,"out/")

setwd(mainpath)

source(paste0(mainpath, "R/ggplot2_theme.R"))
source(paste0(mainpath, "R/plotting_functions.R"))
source(paste0(mainpath, "R/Heatmap_functions.R"))

#Fig5A

load(file = paste0(datapath, "DESeq2_results_allsamples.RData"))

res_all[which(res_all$padj < 0.1),]

Hs.H <- read.table(paste0(datapath, "h.all.v7.1.symbols.gmt"), 
                   header = F, check.names = F, sep = "\t", fill = T, stringsAsFactors = F)

#clean up V1 for aesthetics
Hs.H$V1 <- gsub("HALLMARK_", "", Hs.H$V1)
Hs.H$V1 <- gsub("_", " ", Hs.H$V1)

set.seed(111)
pathwayplot <- gsea.fx(res_all)

#note that order of a few nonsignificant pathways are random due to ties in stats
pdf(paste0(plotpath, "Fig5A.pdf"), width = 20, height = 20)
pathwayplot
dev.off()

#Fig5B

res_all$threshold <- NA
res_all$threshold[ res_all$log2FoldChange > 1 &res_all$padj < 0.05] <- "Up-regulated"
res_all$threshold[ res_all$log2FoldChange < -1 & res_all$padj < 0.05] <- "Down-regulated"
res_all$threshold[ is.na(res_all$threshold)] <- "not significant"

vp <- volcano.fx(res_all, 1, 0.05, "Differential gene expression\n(PR + SD vs. PD)") +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color = "black") 

pdf(paste0(plotpath, "Fig5B.pdf"),
    width = 10, height = 10)
vp
dev.off()






wgcnagenes_mat <- read.csv(file = paste0(datapath, "th_tpm_wgcna_genes.csv"),
                          header = T, stringsAsFactors = F, row.names = 1, check.names = F) 

gene_module <- read.csv(file = paste0(datapath, "th_gene_module_wgcna.csv"),
                        header = T, stringsAsFactors = F,row.names = 1, check.names = F) 

GO_module <- read.csv(file = paste0(datapath, "th_GO_module_wgcna.csv"),
                      header = T, stringsAsFactors = F, check.names = F) 

#Modules ordered based on hierarchical clustering of consensus module eigengenes
orderedlabels <- c('4','18','17','28','6','7','14','24','23','27','3','13','22','26')

# create a vector of ordered modules with genes as names
modulelabels <- as.character(gene_module$moduleLabel)
names(modulelabels) <- rownames(gene_module)
modulelabels <- modulelabels[order(match(modulelabels, orderedlabels))]

# order tpm matrix based on modules
wgcnagenes_mat <- wgcnagenes_mat[names(modulelabels),]

# correlation matrix for genes from wgcna
cor_mat <- rcorr(t(wgcnagenes_mat), type = "pearson")

# Heatmap for modules and gene correlations
set.seed(8)
myColors <- randomcoloR::distinctColorPalette(14)

hm_modules <- Heatmap(modulelabels,
                      #titles and names
                      name = "Modules",
                      show_row_names = FALSE,
                      show_column_names = FALSE,    
                      #clusters
                      cluster_columns = FALSE,
                      cluster_rows = FALSE,
                      #aesthestics
                      col = myColors,
                      row_names_gp = gpar(fontsize = 20),
                      width = unit(1, "cm"),
                      column_title_gp = gpar(fontsize = 42),
                      row_title = NULL)

col_fun <- colorRamp2(c(-1, 0, 1), c("blue","white", "red"))
hm_cor <- Heatmap(cor_mat$r,
                  #titles and names   
                  name = "Gene correlations",   
                  show_row_names = FALSE,
                  show_column_names = FALSE,     
                  #clusters and orders  
                  cluster_columns = FALSE,
                  cluster_rows =  FALSE,
                  show_column_dend = FALSE,
                  #aesthestics
                  col = col_fun, 
                  column_names_gp = gpar(fontsize = 20),
                  row_names_gp = gpar(fontsize = 20),
                  column_title_gp = gpar(fontsize = 20),
                  column_title = NULL,
                  row_title = NULL, use_raster = TRUE,
                  raster_quality = 5)

pdf(file = paste0(plotpath,"Fig5A.pdf"),
    width = 20, 
    height = 20,
    useDingbats = FALSE)
draw(hm_cor + hm_modules)
dev.off()

#Fig5B

metadata <- read.csv(file.path(datapath,"IND_estimate_metadata_IHC_trb_tmb.csv"),header = T, stringsAsFactors = F, check.names = F)


tpm_mat <- read.csv(file = paste0(datapath, "exp_mat/INDICATE.tpm_hg38_ENSG_HUGO.csv"), 
                    header = T, stringsAsFactors = F,check.names = F) 



tpm_mat_ensg <- tpm_mat[!is.na(tpm_mat$ensembl_id),]

tpm_mat_ensg$HGNC_symbol <- NULL
rownames(tpm_mat_ensg) <- gsub("[.].*", "",tpm_mat_ensg$ensembl_id)
tpm_mat_ensg$ensembl_id <- NULL

tpm_mat_ensg_t <- t(tpm_mat_ensg)

rownames(tpm_mat_ensg_t) <- gsub(".*rnaaccess_", "", rownames(tpm_mat_ensg_t))
rownames(tpm_mat_ensg_t) <- gsub("_.*", "", rownames(tpm_mat_ensg_t))
rownames(tpm_mat_ensg_t) <- toupper(rownames(tpm_mat_ensg_t))

tpm_mat_ensg_t <- log2(tpm_mat_ensg_t + 1)

tpm_mat_ensg_t_matched <- tpm_mat_ensg_t[rownames(tpm_mat_ensg_t) %in% metadata$sample_id,]






mygene_modules <- gene_module_th[ gene_module_th$ensembl_id %in% colnames(tpm_mat_ensg_t_matched),]


#Order same as heatmap
labelorders <- c('TH_4','TH_18','TH_17','TH_28','TH_6','TH_7','TH_14',
                 'TH_24','TH_23','TH_27','TH_3','TH_13','TH_22','TH_26')

module_sample <- matrix(ncol = nrow(tpm_mat_ensg_t_matched), nrow = length(unique(mygene_modules$moduleLabel)))
rownames(module_sample) <- labelorders
colnames(module_sample) <- rownames(tpm_mat_ensg_t_matched)

for(mod in 1:nrow(module_sample)){
  mymod <- rownames(module_sample)[mod]    
  modGenes <- mygene_modules$ensembl_id[which(mygene_modules$moduleLabel == mymod)]
  genes <- tpm_mat_ensg_t_matched[,modGenes]
  if(length(modGenes) > 1){
    averagegenes <- apply(genes,1,mean)
    module_sample[mod,] <- averagegenes}
}


module_sample_t <- as.data.frame(t(module_sample))
module_sample_t$sample_id <- rownames(module_sample_t)

metadata_modules <- merge(metadata, module_sample_t, by = "sample_id")

module_cox <- matrix(nrow = nrow(module_sample), ncol = 7)
rownames(module_cox) <- rownames(module_sample)
colnames(module_cox) <- c('coef', 'exp(coef)', 'se(coef)', 'z', 'Pr(>|z|)', 'lower.95' ,'upper.95')

for(i in rownames(module_cox)){
  # message(i)
  f <- as.formula(paste0("Surv(TRTDUR, progressed)~", i))  
  coxmodel <- coxph(f, data=metadata_modules) 
  sumcox <- summary(coxmodel)
  # print(sumcox)
  module_cox[i,1:5] <- sumcox$coefficients[1,1:5]
  module_cox[i,6:7] <- sumcox$conf.int[1,3:4]
}

fdr_df <- as.matrix(p.adjust(module_cox[,5], method = "fdr"))
colnames(fdr_df)[1] <- "fdr"
module_cox <- cbind(module_cox, fdr_df)

module_cox <- as.data.frame(module_cox)


module_cox$term <- NA

module_cox$term <- GO_modules$myterm[match(rownames(module_cox),GO_modules$moduleLabel)]

module_cox$term[module_cox$term == "antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-dependent" ] <- 
  "antigen processing and presentation"
module_cox$term[module_cox$term == "SRP-dependent cotranslational protein targeting to membrane" ] <- "protein targeting to ER"


hrs <- module_cox[,c(2,6,7)]
hrs <- rbind(NA,hrs)

tabtext <- as.data.frame(module_cox[,c(2,5,8)])
tabtext <- format(round(tabtext, 2))
tabtext <- cbind(module_cox[,9],tabtext)
colnames(tabtext) <- c("Signature", "Hazard Ratio", "p-value", "FDR")
rownames(tabtext) <- NULL
tabtext <- rbind(colnames(tabtext),tabtext)

tabtext$`Hazard Ratio` <- NULL


fplot <- forestplot(tabtext, hrs, new_page = TRUE, xlog = TRUE,
                    title = "Univariable model (iMATRIX)", boxsize = 0.25)


