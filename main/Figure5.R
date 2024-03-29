###############
# Figure 5
###############

library(survminer)
library(survival)
library(Hmisc)
library(fgsea)
library(ggrepel)
library(forestplot)

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

Hs.H <- read.table(paste0(datapath, "h.all.v7.1.symbols.gmt"), header = F, check.names = F, sep = "\t", fill = T, stringsAsFactors = F)

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

pdf(paste0(plotpath, "Fig5B.pdf"),width = 10, height = 10)
vp
dev.off()

#Fig5C
load(paste0(datapath, "Consensus_Network_TH_manual_signed_20.RData"))

mycons <- rbind(consMEs[[1]]$data, consMEs[[2]]$data, consMEs[[3]]$data, consMEs[[4]]$data,
                consMEs[[5]]$data, consMEs[[6]]$data, consMEs[[7]]$data)

# remove unclustered module==ME0
mycons$ME0 <- NULL

# make a correlation matrix 
cor_mat <- rcorr(as.matrix(mycons), type = "pearson")

col_fun <- colorRamp2(c(-1, 0, 1), c("blue","white", "red"))
hm_cor <- Heatmap(cor_mat$r,
                  #titles and names   
                  name = "Gene correlations",   
                  show_row_names = TRUE,
                  show_column_names = TRUE,     
                  #clusters and orders  
                  cluster_columns = TRUE,
                  cluster_rows =  TRUE,
                  show_column_dend = TRUE,
                  row_dend_width = unit(5, "cm"),
                  column_dend_height = unit(5, "cm"),
                  #aesthestics
                  col = col_fun, 
                  column_names_gp = gpar(fontsize = 30),
                  row_names_gp = gpar(fontsize = 30),
                  column_title_gp = gpar(fontsize = 30),
                  column_title = NULL,
                  row_title = NULL,
                  show_heatmap_legend = TRUE)

pdf(file = paste0(plotpath,"Fig5C.pdf"),
    width = 20, 
    height = 20,
    useDingbats = FALSE)
draw(hm_cor)
dev.off()

#Fig5D
metadata <- read.csv(file.path(datapath,"anonymized_iMATRIX_Atezo_metadata_IHC_TRB_TMB_v3.csv"), header = T, stringsAsFactors = F, check.names = F)

## gene_module
gene_module_th <- read.csv(file = paste0(datapath, "gene_module_treehouse_manual_signed_20.csv"),header = T, stringsAsFactors = F,row.names = 1,check.names = F) 
gene_module_th <- gene_module_th[ gene_module_th$moduleColor != "grey",]
gene_module_th$moduleLabel <- paste0("TH_", gene_module_th$moduleLabel)

## GO_module
GO_modules <- read.csv(file.path(datapath,"GO_TH_cons_manual_signed_20.csv"), header = T, stringsAsFactors = F, check.names = F, row.names = 1)  
GO_modules <- GO_modules[ GO_modules$module != "grey",]
GO_modules$moduleLabel <- gene_module_th$moduleLabel[match(GO_modules$module, gene_module_th$moduleColor)]

## Add one MF for TH_28 module
TH28_MF <- GO_modules[ GO_modules$termName == "transcription regulator activity",]
GO_modules <- GO_modules[GO_modules$termOntology == "BP",]
GO_modules <- rbind(GO_modules, TH28_MF)

## just to remove small GO terms
GO_modules <- GO_modules[ GO_modules$nModGenesInTerm >= 9,]

## get one rather simplified GO term for each module by choosing the GO terrm with minimum number of genes (> 9)
GO_modules$myterm <- NA
for(i in unique(GO_modules$moduleLabel)){
  tmp <- GO_modules[ GO_modules$moduleLabel == i,]
  moduleterm <- tmp$termName[tmp$nModGenesInTerm == min(tmp$nModGenesInTerm)][1] #use one term if min returns two terms
  GO_modules$myterm[GO_modules$moduleLabel == i] <- moduleterm
  
}

dim(gene_module_th)
dim(GO_modules)

# expression matrix
tpm_mat <- read.csv(file = paste0(datapath, "anonymized_iMATRIX_Atezo_tpm_matrix_hg38.csv"), header = T, stringsAsFactors = F,check.names = F, row.names = 1) 

hgnc_ensembl_ids <- read.csv(file = paste0(datapath, "hgnc_ensembl_ids.csv"), header = T, stringsAsFactors = F,check.names = F, row.names = 1) 

#Use Ensembl id to compare iMatrix data and TH data. They are more reliable
rownames(tpm_mat) <- hgnc_ensembl_ids$ensembl_id[ match(rownames(tpm_mat), hgnc_ensembl_ids$hgnc_symbol)]

# log2
tpm_mat_log2 <- log2(tpm_mat + 1)

# genes in TH and iMATRIX
mygene_modules <- gene_module_th[ gene_module_th$ensembl_id %in% rownames(tpm_mat_log2),]

#Order same as heatmap
labelorders <- c('TH_7','TH_14','TH_17','TH_24','TH_28','TH_6','TH_18',
                 'TH_4','TH_22','TH_27','TH_3','TH_13','TH_26','TH_23')

module_sample <- matrix(ncol = ncol(tpm_mat_log2), nrow = length(unique(mygene_modules$moduleLabel)))
rownames(module_sample) <- labelorders
colnames(module_sample) <- colnames(tpm_mat_log2)

for(mod in 1:nrow(module_sample)){
  mymod <- rownames(module_sample)[mod]    
  modGenes <- mygene_modules$ensembl_id[which(mygene_modules$moduleLabel == mymod)]
  genes <- tpm_mat_log2[modGenes,]
  if(length(modGenes) > 1){
    averagegenes <- apply(genes,2,mean)
    module_sample[mod,] <- averagegenes}
}

module_sample_t <- as.data.frame(t(module_sample))
module_sample_t$sample_id <- rownames(module_sample_t)

metadata_modules <- merge(metadata, module_sample_t, by.x =  "trunc_anonymized_rnaseq_sample_id", by.y =  "sample_id")

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

# some cleanup for terms
module_cox$term[module_cox$term == "antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-dependent" ] <- 
  "antigen processing and presentation"
module_cox$term[module_cox$term == "SRP-dependent cotranslational protein targeting to membrane" ] <- "protein targeting to ER"

# HR table
hrs <- module_cox[,c(2,6,7)]
hrs <- rbind(NA,hrs)

tabtext <- as.data.frame(module_cox[,c(2,5,8)])
tabtext <- format(round(tabtext, 2))
tabtext <- cbind(module_cox[,9],tabtext)
colnames(tabtext) <- c("Signature", "Hazard Ratio", "p-value", "FDR")
rownames(tabtext) <- NULL
tabtext <- rbind(colnames(tabtext),tabtext)
tabtext$`Hazard Ratio` <- NULL

# from SO
fn <- local({
  i = 0
  b_clrs = c(rep("black", 9), "red", "black", rep("red", 2), "black" )
  l_clrs = c(rep("black", 9), "red", "black", rep("red", 2), "black" )
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
  }
})

pdf(file = paste0(plotpath,"Fig5D.pdf"),
    width = 7, height = 10,
    useDingbats = FALSE, onefile = F)
forestplot(tabtext, fn.ci_norm = fn, hrs, new_page = TRUE, xlog = TRUE,
           title = "Survival analysis (iMATRIX-atezo)(p-value, FDR)", boxsize = 0.25)
dev.off()

# Fig5E
mygen <- c('ENSG00000116824', 'ENSG00000167286', 'ENSG00000198851', 'ENSG00000160654', 
           'ENSG00000186810', 'ENSG00000234745', 'ENSG00000204642', 'ENSG00000240065', 
           'ENSG00000172673', 'ENSG00000160185')

pedges_ssgsea <- GSVA::gsva(as.matrix(tpm_mat_log2), list(mygen), method = "ssgsea", ssgsea.norm= TRUE)

metadata_ges <- cbind(metadata, pedges_ssgsea[,metadata$trunc_anonymized_rnaseq_sample_id])
colnames(metadata_ges)[ colnames(metadata_ges) == "pedges_ssgsea[, metadata$trunc_anonymized_rnaseq_sample_id]"] <- "PedGES_GSEA"

summary(metadata_ges$PedGES_GSEA)

metadata_ges$GES_Group <- NA
metadata_ges$GES_Group[ metadata_ges$PedGES_GSEA >= 1.4272] <- "High"
metadata_ges$GES_Group[ metadata_ges$PedGES_GSEA < 1.4272 &
                          metadata_ges$PedGES_GSEA > 1.0413] <- "Intermediate"
metadata_ges$GES_Group[ metadata_ges$PedGES_GSEA <= 1.0413] <- "Low"

sfit <- survfit(Surv(TRTDUR, progressed)~ GES_Group, data= metadata_ges)

kmplot <- ggsurvplot(sfit, 
                     conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 10, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 20,
                     xlim = c(0,800),break.x.by = 200,
                     legend = c(0.7, 0.95), font.legend = 35, legend.title = "",
                     legend.labs = c("High","Intermediate", "Low"),
                     font.main = 35, font.x = 35, font.y = 35, font.tickslab = 35,
                     ylab = "Progression-free survival") 

kmplot$table <- ggrisktable(sfit, data = metadata_ges, 
                            color = "strata", 
                            palette = c("#ED2024", "#adadad", "#3953A4"),
                            fontsize = 12, 
                            xlim = c(0,800),break.x.by = 200, 
                            tables.theme = theme_cleantable(), font.tickslab = 30, 
                            y.text = TRUE, ylab = "",  xlab = "Time (days)",
                            legend.labs = c("High","Intermediate", "Low"))

kmplot$table <- kmplot$table + theme(plot.title = element_blank(),
                                     axis.text.x = element_text(size = 35), axis.title.x = element_text(size = 35))

kmplot$plot <- kmplot$plot + labs(title = "iMATRIX-atezo\n") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig5E.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()

rm(kmplot)

# Fig5F
#modified from tableS2 riaz et al
load(file = paste0(datapath,"riaz_tabS2.RData")) 


summary(riaz_tabS2$PedGES_GSEA)

riaz_tabS2$GES_group <- NA
riaz_tabS2$GES_group[ riaz_tabS2$PedGES_GSEA >= 0.60221] <- "High"
riaz_tabS2$GES_group[ riaz_tabS2$PedGES_GSEA < 0.60221 &
                      riaz_tabS2$PedGES_GSEA > 0.26879] <- "Intermediate"
riaz_tabS2$GES_group[ riaz_tabS2$PedGES_GSEA <= 0.26879] <- "Low"

sfit <- survfit(Surv(PFS, progressed)~ GES_group, 
                data= riaz_tabS2[riaz_tabS2$Cohort == "NIV3-PROG",])

kmplot <- ggsurvplot(sfit, 
                     conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 10, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 20,
                     legend = c(0.7, 0.95), font.legend = 35, legend.title = "",
                     legend.labs = c("High","Intermediate", "Low"),
                     font.main = 35, font.x = 35, font.y = 35, font.tickslab = 35,
                     xlim = c(0,1000),break.x.by = 250,
                     ylab = "Progression-free survival") 

kmplot$table <- ggrisktable(sfit, data = riaz_tabS2[riaz_tabS2$Cohort == "NIV3-PROG",], 
                            color = "strata", 
                            palette = c("#ED2024", "#adadad", "#3953A4"),
                            fontsize = 12,
                            tables.theme = theme_cleantable(), font.tickslab = 30,
                            xlim = c(0,1000),break.x.by = 250, 
                            y.text = TRUE, ylab = "",  xlab = "Time (days)",
                            legend.labs = c("High","Intermediate", "Low"))

kmplot$table <- kmplot$table + theme(plot.title = element_blank(),
                                     axis.text.x = element_text(size = 35),
                                     axis.title.x = element_text(size = 35),
                                     plot.margin = unit(c(0,1,0,0), "cm"),)

kmplot$plot <- kmplot$plot + labs(title = "Riaz et al 2017\nNivo treatment after Ipi Prog") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0,1,0,0), "cm"),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig5F.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()

rm(kmplot)

# Fig5G
# modified from Cindy Nat Comm paper
load(file = paste0(datapath,"ins.RData")) 

summary(ins$PedGES_GSEA)

ins$GES_group <- NA
ins$GES_group[ ins$PedGES_GSEA >= 1.3693] <- "High"
ins$GES_group[ ins$PedGES_GSEA < 1.3693 &
                 ins$PedGES_GSEA > 1.0031] <- "Intermediate"
ins$GES_group[ ins$PedGES_GSEA <= 1.0031] <- "Low"

sfit <- survfit(Surv(OSTIME_Months, OSevent)~ GES_group, data= ins)

kmplot <- ggsurvplot(sfit, 
                     conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 10, pval.coord = c(30, 0.75),
                     risk.table=TRUE, fontsize = 20,
                     xlim = c(0,40),break.x.by = 10,
                     legend = c(0.7, 0.95), font.legend = 35, legend.title = "",
                     legend.labs = c("High","Intermediate", "Low"),
                     font.main = 35, font.x = 35, font.y = 35, font.tickslab = 35,
                     ylab = "Overall survival") 

kmplot$table <- ggrisktable(sfit, data = ins, 
                            color = "strata", 
                            palette = c("#ED2024", "#adadad", "#3953A4"),
                            fontsize = 12,
                            xlim = c(0,40),break.x.by = 10,
                            tables.theme = theme_cleantable(), font.tickslab = 30, 
                            y.text = TRUE, ylab = "",  xlab = "Time (months)",
                            legend.labs = c("High","Intermediate", "Low"))

kmplot$table <- kmplot$table + theme(plot.title = element_blank(),
                                     axis.text.x = element_text(size = 35), axis.title.x = element_text(size = 35))

kmplot$plot <- kmplot$plot + labs(title = "INSPIRE\n") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig5G.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()

rm(kmplot)
