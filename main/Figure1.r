###############
# Figure 1
###############

mainpath <- "~/git/iMATRIX-Atezo_Biomarker/"
datapath <- paste0(mainpath,"data/")

dir.create(file.path(paste0(mainpath, "out/")))

plotpath <- paste0(mainpath,"out/")

setwd(mainpath)

source(paste0(mainpath, "R/ggplot2_theme.R"))
source(paste0(mainpath, "R/Heatmap_functions.R"))

metadata <- read.csv(file.path(datapath,"anonymized_iMATRIX_Atezo_metadata_IHC_TRB_TMB_v2.csv"), header = T, stringsAsFactors = F, check.names = F)

# Order by cancer group and treatment duration
metadata <- metadata[order(metadata$cancer, metadata$TRTDUR), ]

PFS_ha <- HeatmapAnnotation(`Treatment duration` = anno_barplot(metadata$TRTDUR), height = unit(3, "cm"),
                          annotation_name_gp = gpar(fontsize = 10))


cols <- c("IHC_CD8", "IHC_CD3", "IHC_PDL1TC", "IHC_PDL1IC")
ihc_mat <- metadata[, cols]

rownames(ihc_mat) <- metadata$sample_id

ihc_mat[is.na(ihc_mat)] <- NA
ihc_mat[ihc_mat == ""] <- NA
ihc_mat[ihc_mat == "<1"] <- 0.5


ihc_mat[,cols] <- sapply(ihc_mat[,cols],as.numeric)

ihc_mat <- as.matrix(ihc_mat)

mode(ihc_mat) <- "numeric"

CD8_ha <- HeatmapAnnotation(`%CD8+` = anno_barplot(ihc_mat[,1]), height = unit(3, "cm"),
                           annotation_name_gp = gpar(fontsize = 10))

CD3_ha <- HeatmapAnnotation(`%CD3+` = anno_barplot(ihc_mat[,2]), height = unit(3, "cm"),
                           annotation_name_gp = gpar(fontsize = 10))

PDL1_TC_ha <- HeatmapAnnotation(`%PDL1+ in\n tumour cells` = anno_barplot(ihc_mat[,3]), height = unit(3, "cm"),
                               annotation_name_gp = gpar(fontsize = 10))

PDL1_IC_ha <- HeatmapAnnotation(`%PDL1+ in\n immune cells` = anno_barplot(ihc_mat[,4]), height = unit(3, "cm"),
                               annotation_name_gp = gpar(fontsize = 10))

Div_ha <- HeatmapAnnotation(`TCRb diversity` = anno_barplot(log10(metadata$observed_Shannon), height = unit(3, "cm"),
                                                          axis_param = list(at = c(0,1,2,3),
                                                                            labels = c("1", "10", "100", "1000"))),
                           annotation_name_gp = gpar(fontsize = 10))

#recode 0 to 0.01 for log10 transformation
metadata$TMB <- NA
metadata$TMB <- as.numeric(metadata$TMB_Score)
metadata$TMB[metadata$TMB == 0] <- 0.09

TMB_ha <- HeatmapAnnotation(`Tumour mutation\nburden` = anno_points(log10(metadata$TMB), height = unit(3,"cm"),
                                                                   axis_param = list(at = c(-1,0, 1, 2),
                                                                                     labels = c("<0.1", "1", "10", "100"))),
                           annotation_name_gp = gpar(fontsize = 10))

age_ha <- HeatmapAnnotation(Age = anno_barplot(metadata$AGE, height = unit(2, "cm")),annotation_name_gp = gpar(fontsize = 10))

mycohort <- metadata$disease_group
names(mycohort) <- metadata$sample_id
mycohorts <- t(as.matrix(mycohort))
rownames(mycohorts) <- "Cohort"

cohorts_hm <- cohorts_hm.fx(mycohorts)

myresponse <- metadata$`BCOR-INV`
names(myresponse) <- metadata$sample_id
myresponse <- t(as.matrix(myresponse))
rownames(myresponse) <- "Objective response"

response_hm <- response_hm.fx(myresponse)

mycancer <- metadata$cancer
names(mycancer) <- metadata$sample_id
mycancer <- t(as.matrix(mycancer))
rownames(mycancer) <- "Cancer"

cancer_hm <- cancer_hm.fx(mycancer)

metadata$origin <- NA
metadata$origin[ metadata$sample_origin == "Lymph node"] <- "Lymph node"
metadata$origin[ metadata$sample_origin != "Lymph node"] <- "Other tissue"

myorigin <- metadata$origin
names(myorigin) <- metadata$sample_id
myorigin <- t(as.matrix(myorigin))
rownames(myorigin) <- "Lymph node/other tissue"

origin_hm <- origin_hm.fx(myorigin)

mytype <- metadata$tumor_type
names(mytype) <- metadata$sample_id
mytype <- t(as.matrix(mytype))
rownames(mytype) <- "Tumour type"

type_hm <- type_hm.fx(mytype)

heatmaps <- PFS_ha %v% response_hm %v% origin_hm %v%
  CD8_ha %v% CD3_ha %v% PDL1_IC_ha %v% PDL1_TC_ha %v% TMB_ha %v% Div_ha 

pdf(paste0(plotpath, "Fig1.pdf"),
   width = 10, height = 20, useDingbats = FALSE)

draw(heatmaps,heatmap_legend_side = "bottom")

dev.off()

