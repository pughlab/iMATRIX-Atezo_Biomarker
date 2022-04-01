###############
# Figure 3
###############

library(survminer)
library(survival)
library(forestmodel)

mainpath <- "~/git/iMATRIX-Atezo_Biomarker/"
datapath <- paste0(mainpath,"data/")

dir.create(file.path(paste0(mainpath, "out/")))

plotpath <- paste0(mainpath,"out/")

setwd(mainpath)

source(paste0(mainpath, "R/ggplot2_theme.R"))
source(paste0(mainpath, "R/Heatmap_functions.R"))
source(paste0(mainpath, "R/Splot_function.R"))

metadata <- read.csv(file.path(datapath,"metadata.csv"), header = T, stringsAsFactors = F, check.names = F)

message("summary for PDL1 immune cells (IHC) in the iMATRIX-Atezo:")
summary(metadata$IHC_PDL1IC)
message("summary for PDL1 tumor cells (IHC) in the iMATRIX-Atezo:")
summary(metadata$IHC_PDL1TC)


#Fig3A

metadata$PDL1group <- NA
metadata$PDL1group[metadata$PDL1 == "High"] <- "High"
metadata$PDL1group[metadata$PDL1 == "Low"] <- "Low/No"
metadata$PDL1group[metadata$PDL1 == "No expression"] <- "Low/No"


sfit <- survfit(Surv(TRTDUR, progressed)~ PDL1group, data= metadata)

kmplot <- ggsurvplot(sfit, conf.int=FALSE, palette = c("#ED2024", "#3953A4"),
                     pval = TRUE, pval.size = 8, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 8,
                     legend = c(0.7, 0.95),font.legend = 16, legend.title = "",
                     font.main = 23, font.x = 25,font.y = 25, font.tickslab = 20) 

kmplot$table <- kmplot$table + theme(axis.text.x = element_text(size = 20))
kmplot$plot <- kmplot$plot + labs(title = "PD-L1 protein expression")

pdf(file = paste0(plotpath,"Fig3A.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()

# Fig3B

metadata$PDL1group <- factor(metadata$PDL1group, levels = c("Low/No", "High"))
coxmodel <- coxph(Surv(TRTDUR, progressed)~ PDL1group + cancer + AGE, data= metadata) 

fig3b <- forest_model(coxmodel,exponentiate = TRUE) + 
  labs(title = "PDL1 protein expression")

pdf(file = paste0(plotpath,"Fig3B.pdf"),
    width = 15, 
    height = 10,
    useDingbats = FALSE)
fig3b
dev.off()


