###############
# Figure 2
###############

library(survminer)
library(survival)

mainpath <- "~/git/iMATRIX-Atezo_Biomarker/"
datapath <- paste0(mainpath,"data/")

dir.create(file.path(paste0(mainpath, "out/")))

plotpath <- paste0(mainpath,"out/")

setwd(mainpath)

source(paste0(mainpath, "R/ggplot2_theme.R"))
source(paste0(mainpath, "R/Heatmap_functions.R"))
source(paste0(mainpath, "R/Splot_function.R"))

metadata <- read.csv(file.path(datapath,"anonymized_iMATRIX_Atezo_metadata_IHC_TRB_TMB_v3.csv"), header = T, stringsAsFactors = F, check.names = F)

message("summary for CD8 staining (IHC) in the iMATRIX-Atezo:")
summary(metadata$IHC_CD8)
message("summary for CD3 staining (IHC) in the iMATRIX-Atezo:")
summary(metadata$IHC_CD3)

#Fig2A
metadata$CD8group <- NA
metadata$CD8group[metadata$IHC_CD8 >= 1.242] <- "High"
metadata$CD8group[metadata$IHC_CD8 > 0.070 &
                  metadata$IHC_CD8 < 1.242] <- "Intermediate"
metadata$CD8group[metadata$IHC_CD8 <= 0.070] <- "Low"

sfit <- survfit(Surv(TRTDUR, progressed)~ CD8group, data = metadata)

kmplot <- ggsurvplot(sfit, 
                     conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 10, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 20,
                     xlim = c(0,800),break.x.by = 200,
                     legend = c(0.7, 0.95), font.legend = 35, legend.title = "",
                     legend.labs = c("High","Intermediate", "Low"),
                     font.main = 35, font.x = 35, font.y = 35, font.tickslab = 35,
                     ylab = "Progression-free survival") 

kmplot$table <- ggrisktable(sfit, data = metadata, 
                            color = "strata", 
                            palette = c("#ED2024", "#adadad", "#3953A4"),
                            fontsize = 12, 
                            xlim = c(0,800),break.x.by = 200, 
                            tables.theme = theme_cleantable(), font.tickslab = 30, 
                            y.text = TRUE, ylab = "",  xlab = "Time (days)",
                            legend.labs = c("High","Intermediate", "Low"))

kmplot$table <- kmplot$table + theme(plot.title = element_blank(),
                                     axis.text.x = element_text(size = 35), axis.title.x = element_text(size = 35))

kmplot$plot <- kmplot$plot + labs(title = "CD8 protein expression (IHC)\n") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit( c(0.1,1,0,0), "cm"),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig2A.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()

rm(kmplot)

# Fig2B
#remove lymph nodes
sfit <- survfit(Surv(TRTDUR, progressed)~ CD8group, 
                data= metadata[ metadata$sample_origin != "Lymph node",])

kmplot <- ggsurvplot(sfit, 
                     conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 10, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 20,
                     xlim = c(0,800),break.x.by = 200,
                     legend = c(0.7, 0.95), font.legend = 35, legend.title = "",
                     legend.labs = c("High","Intermediate", "Low"),
                     font.main = 35, font.x = 35, font.y = 35, font.tickslab = 35,
                     ylab = "Progression-free survival") 

kmplot$table <- ggrisktable(sfit, data = metadata[ metadata$sample_origin != "Lymph node",],
                            color = "strata", 
                            palette = c("#ED2024", "#adadad", "#3953A4"),
                            fontsize = 12,
                            xlim = c(0,800),break.x.by = 200,
                            tables.theme = theme_cleantable(), font.tickslab = 30, 
                            y.text = TRUE, ylab = "",  xlab = "Time (days)",
                            legend.labs = c("High","Intermediate", "Low"))

kmplot$table <- kmplot$table + theme(plot.title = element_blank(),
                                     axis.text.x = element_text(size = 35), axis.title.x = element_text(size = 35))

kmplot$plot <- kmplot$plot + labs(title = "CD8 protein expression\nexcluding lymph nodes (IHC)") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit( c(0.1,1,0,0), "cm"),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig2B.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()


rm(kmplot)

# Fig2C
metadata$CD3group <- NA
metadata$CD3group[metadata$IHC_CD3 >= 3.05] <- "High"
metadata$CD3group[metadata$IHC_CD3 > 0.16 &
                  metadata$IHC_CD3 < 3.05] <- "Intermediate"
metadata$CD3group[metadata$IHC_CD3 <= 0.16] <- "Low"

sfit <- survfit(Surv(TRTDUR, progressed) ~ CD3group, data = metadata)

kmplot <- ggsurvplot(sfit, 
                     conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 10, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 20,
                     xlim = c(0,800),break.x.by = 200,
                     legend = c(0.7, 0.95), font.legend = 35, legend.title = "",
                     legend.labs = c("High","Intermediate", "Low"),
                     font.main = 35, font.x = 35, font.y = 35, font.tickslab = 35,
                     ylab = "Progression-free survival") 

kmplot$table <- ggrisktable(sfit, data = metadata, 
                            color = "strata", 
                            palette = c("#ED2024", "#adadad", "#3953A4"),
                            fontsize = 12,
                            xlim = c(0,800),break.x.by = 200,
                            tables.theme = theme_cleantable(), font.tickslab = 30, 
                            y.text = TRUE, ylab = "",  xlab = "Time (days)",
                            legend.labs = c("High","Intermediate", "Low"))

kmplot$table <- kmplot$table + theme(plot.title = element_blank(),
                                     axis.text.x = element_text(size = 35), axis.title.x = element_text(size = 35))

kmplot$plot <- kmplot$plot + labs(title = "CD3 protein expression (IHC)\n") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit( c(0.1,1,0,0), "cm"),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig2C.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()

# Fig2D

metadata <- read.csv(file.path(datapath,"anonymized_iMATRIX_Atezo_metadata_IHC_TRB_TMB_v3.csv"), header = T, stringsAsFactors = F, check.names = F)
deconv <- read.csv(file.path(datapath,"anonymized_all_immunedeconv.csv"), header = T, stringsAsFactors = F, check.names = F)

metadata_deconv <- merge(metadata, deconv, by = "trunc_anonymized_rnaseq_sample_id")


message("summary for CD8 estimate from cibersort in the iMATRIX-Atezo:")
summary(metadata_deconv$T.cell.CD8._cibersort_abs)

metadata_deconv$CD8group <- NA
metadata_deconv$CD8group[metadata_deconv$T.cell.CD8._cibersort_abs >= 0.0938328] <- "High"

metadata_deconv$CD8group[metadata_deconv$T.cell.CD8._cibersort_abs < 0.0938328 &
                           metadata_deconv$T.cell.CD8._cibersort_abs > 0.0007922] <- "Intermediate"

metadata_deconv$CD8group[metadata_deconv$T.cell.CD8._cibersort_abs <= 0.0007922 ] <- "Low"

sfit <- survfit(Surv(TRTDUR, progressed)~ CD8group, data= metadata_deconv)

kmplot <- ggsurvplot(sfit, 
                     conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 10, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 20,
                     xlim = c(0,800),break.x.by = 200,
                     legend = c(0.7, 0.95), font.legend = 35, legend.title = "",
                     legend.labs = c("High","Intermediate", "Low"),
                     font.main = 35, font.x = 35, font.y = 35, font.tickslab = 35,
                     ylab = "Progression-free survival") 

kmplot$table <- ggrisktable(sfit, data = metadata_deconv, 
                            color = "strata", 
                            palette = c("#ED2024", "#adadad", "#3953A4"),
                            fontsize = 12,
                            xlim = c(0,800),break.x.by = 200,
                            tables.theme = theme_cleantable(), font.tickslab = 30, 
                            y.text = TRUE, ylab = "",  xlab = "Time (days)",
                            legend.labs = c("High","Intermediate", "Low"))

kmplot$table <- kmplot$table + theme(plot.title = element_blank(),
                                     axis.text.x = element_text(size = 35), axis.title.x = element_text(size = 35))

kmplot$plot <- kmplot$plot + labs(title = "CD8 gene signature (CIBERSORT)\n") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit( c(0.1,1,0,0), "cm"),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig2D.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()

rm(kmplot)

# Fig2E
#remove lymph node
sfit <- survfit(Surv(TRTDUR, progressed)~ CD8group, 
                data= metadata_deconv[metadata_deconv$sample_origin != "Lymph node",])

kmplot <- ggsurvplot(sfit, 
                     conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 10, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 20,
                     xlim = c(0,800),break.x.by = 200,
                     legend = c(0.7, 0.95), font.legend = 35, legend.title = "",
                     legend.labs = c("High","Intermediate", "Low"),
                     font.main = 35, font.x = 35, font.y = 35, font.tickslab = 35,
                     ylab = "Progression-free survival") 

kmplot$table <- ggrisktable(sfit, 
                            data = metadata_deconv[metadata_deconv$sample_origin != "Lymph node",], 
                            color = "strata", 
                            palette = c("#ED2024", "#adadad", "#3953A4"),
                            fontsize = 12,
                            xlim = c(0,800),break.x.by = 200,
                            tables.theme = theme_cleantable(), font.tickslab = 30, 
                            y.text = TRUE, ylab = "",  xlab = "Time (days)",
                            legend.labs = c("High","Intermediate", "Low"))

kmplot$table <- kmplot$table + theme(plot.title = element_blank(),
                                     axis.text.x = element_text(size = 35), axis.title.x = element_text(size = 35))

kmplot$plot <- kmplot$plot + labs(title = "CD8 gene signature\nexcluding lymph nodes (CIBERSORT)") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit( c(0.1,2.5,0,0), "cm"),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig2E.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()

rm(kmplot)

#Fig2F
summary(metadata_deconv$T.cell_mcpcounter)

metadata_deconv$CD8group <- NA
metadata_deconv$CD8group[metadata_deconv$T.cell_mcpcounter >= 4.6172] <- "High"

metadata_deconv$CD8group[metadata_deconv$T.cell_mcpcounter < 4.6172 &
                           metadata_deconv$T.cell_mcpcounter > 0.9173] <- "Intermediate"

metadata_deconv$CD8group[metadata_deconv$T.cell_mcpcounter <= 0.9173 ] <- "Low"

sfit <- survfit(Surv(TRTDUR, progressed)~ CD8group, data= metadata_deconv)

kmplot <- ggsurvplot(sfit, 
                     conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 10, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 20,
                     xlim = c(0,800),break.x.by = 200,
                     legend = c(0.7, 0.95), font.legend = 35, legend.title = "",
                     legend.labs = c("High","Intermediate", "Low"),
                     font.main = 35, font.x = 35, font.y = 35, font.tickslab = 35,
                     ylab = "Progression-free survival") 

kmplot$table <- ggrisktable(sfit, 
                            data = metadata_deconv, 
                            color = "strata", 
                            xlim = c(0,800),break.x.by = 200,
                            palette = c("#ED2024", "#adadad", "#3953A4"),
                            fontsize = 12,
                            tables.theme = theme_cleantable(), font.tickslab = 30, 
                            y.text = TRUE, ylab = "",  xlab = "Time (days)",
                            legend.labs = c("High","Intermediate", "Low"))

kmplot$table <- kmplot$table + theme(plot.title = element_blank(),
                                     axis.text.x = element_text(size = 35), axis.title.x = element_text(size = 35))

kmplot$plot <- kmplot$plot + labs(title = "CD8 gene signature (MCPcounter)\n") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit( c(0.1,1,0,0), "cm"),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig2F.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()

rm(kmplot)




