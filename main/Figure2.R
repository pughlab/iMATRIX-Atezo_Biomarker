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

metadata <- read.csv(file.path(datapath,"metadata.csv"), header = T, stringsAsFactors = F, check.names = F)

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

sfit <- survfit(Surv(TRTDUR, progressed)~ CD8group, data= metadata)

kmplot_CD8 <- ggsurvplot(sfit, conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 8, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 8,
                     legend = c(0.7, 0.95), font.legend = 16, legend.title = "",
                     font.main = 23, font.x = 25,font.y = 25, font.tickslab = 20) 

kmplot_CD8$table <- kmplot_CD8$table + theme(axis.text.x = element_text(size = 20))
kmplot_CD8$plot <- kmplot_CD8$plot + labs(title = "PFS analysis of CD8 protein expression")

pdf(file = paste0(plotpath,"Fig2A.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot_CD8
dev.off()

# Fig2B
metadata$CD3group <- NA
metadata$CD3group[metadata$IHC_CD3 >= 3.05] <- "High"
metadata$CD3group[metadata$IHC_CD3 > 0.16 &
                  metadata$IHC_CD3 < 3.05] <- "Intermediate"
metadata$CD3group[metadata$IHC_CD3 <= 0.16] <- "Low"

sfit <- survfit(Surv(TRTDUR, progressed)~ CD3group, data= metadata)
kmplot_CD3 <- ggsurvplot(sfit, conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 8, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 8,
                     legend = c(0.7, 0.95), font.legend = 16, legend.title = "",
                     font.main = 23, font.x = 25,font.y = 25, font.tickslab = 20) 

kmplot_CD3$table <- kmplot_CD3$table + theme(axis.text.x = element_text(size = 20))
kmplot_CD3$plot <- kmplot_CD3$plot + labs(title = "PFS analysis of CD3 protein expression")

pdf(file = paste0(plotpath,"Fig2B.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot_CD3
dev.off()

# Fig2C
metadata_CD8 <- read.csv(file.path(datapath,"metadata_CD8CIBERSORT.csv"), header = T, stringsAsFactors = F, check.names = F)

message("summary for CD8 estimate from cibersort in the iMATRIX-Atezo:")
summary(metadata_CD8$T_cell_CD8_CIBERSORT_ABS)

metadata_CD8$CD8group <- NA
metadata_CD8$CD8group[metadata_CD8$T_cell_CD8_CIBERSORT_ABS >= 0.10383] <- "High"
metadata_CD8$CD8group[metadata_CD8$T_cell_CD8_CIBERSORT_ABS < 0.10383 &
                      metadata_CD8$T_cell_CD8_CIBERSORT_ABS > 0] <- "Intermediate"
metadata_CD8$CD8group[metadata_CD8$T_cell_CD8_CIBERSORT_ABS == 0 ] <- "Low"

sfit <- survfit(Surv(TRTDUR, progressed)~ CD8group, data= metadata_CD8)
kmplot_ciber <- ggsurvplot(sfit, conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 8, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 8,
                     legend = c(0.7, 0.95), font.legend = 16, legend.title = "",
                     font.main = 23, font.x = 25,font.y = 25, font.tickslab = 20) 

kmplot_ciber$table <- kmplot_ciber$table + theme(axis.text.x = element_text(size = 20))
kmplot_ciber$plot <- kmplot_ciber$plot + labs(title = "PFS analysis of CD8 gene signature (CIBERSORT)")

pdf(file = paste0(plotpath,"Fig2C.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot_ciber
dev.off()

# Fig2D
th_data <- read.csv(file.path(datapath,"treehouse_CD8.csv"), header = T, stringsAsFactors = F, 
                    check.names = F, row.names = NULL)


th_data$CD8Level <- NA
th_data$CD8Level[th_data$CD8_CIBERSORT >= 0.10383] <- "High"
th_data$CD8Level[th_data$CD8_CIBERSORT < 0.10383 & 
                 th_data$CD8_CIBERSORT > 0] <- "Intermediate"
th_data$CD8Level[th_data$CD8_CIBERSORT == 0] <- "Low"

#recode 0 for visualization
th_data$CD8_CIBERSORT[ th_data$CD8_CIBERSORT == 0 ] <- 0.00001

th_data$cohort <- th_data$diseasetype
th_data$group <- "TH"

#Splot
th_mediandf <- median.cohorts.fx(th_data)
list_sorted_df <- sort.cohorts.fx(th_data, th_mediandf)
cd8_Splot <- Splot.fx(list_sorted_df, "Data from Treehouse (n = 1851)") + 
  annotation_logticks(sides = "l") + scale_y_continuous(trans = "log10",
                                                        breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1),
                                                        label = c(0, 0.0001, 0.001, 0.01, 0.1)) 

pdf(file = paste0(plotpath,"Fig2D.pdf"),
    width = 25, 
    height = 12,
    useDingbats = FALSE,
    onefile = FALSE)
cd8_Splot + theme(legend.position = "right")
dev.off()
