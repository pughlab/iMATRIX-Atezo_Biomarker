###############
# Figure 3
###############
options(scipen = 999)


library(survminer)
library(survival)
library(forestmodel)

library(ggsignif)
library(ggbeeswarm)

mainpath <- "~/git/iMATRIX-Atezo_Biomarker/"
datapath <- paste0(mainpath,"data/")

dir.create(file.path(paste0(mainpath, "out/")))

plotpath <- paste0(mainpath,"out/")

setwd(mainpath)

source(paste0(mainpath, "R/ggplot2_theme.R"))

metadata <- read.csv(file.path(datapath,"IND_metadata_IHC_trb_tmb.csv"), header = T, stringsAsFactors = F, check.names = F)

message("summary for PDL1 immune cells (IHC) in the iMATRIX-Atezo:")
summary(metadata$IHC_PDL1IC)
message("summary for PDL1 tumor cells (IHC) in the iMATRIX-Atezo:")
summary(metadata$IHC_PDL1TC)

#Fig3A

metadata$PDL1 <- factor(metadata$PDL1, levels = c("No expression", "Low", "High"))

metadata$PDL1group <- NA
metadata$PDL1group[metadata$PDL1 == "High"] <- "High"
metadata$PDL1group[metadata$PDL1 == "Low"] <- "Low/No"
metadata$PDL1group[metadata$PDL1 == "No expression"] <- "Low/No"

sfit <- survfit(Surv(TRTDUR, progressed)~ PDL1group, data= metadata)

kmplot <- ggsurvplot(sfit, 
                     conf.int=FALSE, palette = c("#ED2024", "#3953A4"),
                     pval = TRUE, pval.size = 10, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 20,
                     xlim = c(0,800),break.x.by = 200,
                     legend = c(0.7, 0.95), font.legend = 35, legend.title = "",
                     legend.labs = c("High", "Low/No"),
                     font.main = 35, font.x = 35, font.y = 35, font.tickslab = 35,
                     ylab = "Progression-free survival") 

kmplot$table <- ggrisktable(sfit, data = metadata, 
                            color = "strata", 
                            palette = c("#ED2024", "#3953A4"),
                            fontsize = 12, 
                            xlim = c(0,800),break.x.by = 200, 
                            tables.theme = theme_cleantable(), font.tickslab = 30, 
                            y.text = TRUE, ylab = "",  xlab = "Time (days)",
                            legend.labs = c("High", "Low/No"))

kmplot$table <- kmplot$table + theme(plot.title = element_blank(),
                                     axis.text.x = element_text(size = 35), axis.title.x = element_text(size = 35))

kmplot$plot <- kmplot$plot + labs(title = "PD-L1 protein expression (IHC)\n") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.1,1,0,0), "cm"),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig3A.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()

# Fig3B
metadata$PDL1group <- factor(metadata$PDL1group, levels = c("Low/No", "High"))

#cox model
coxmodel <- coxph(Surv(TRTDUR, progressed)~ PDL1group + cancer + AGE, data= metadata) 

summary(coxmodel)

pcox <- forest_model(coxmodel,exponentiate = TRUE) + 
  labs(title = "Cox for PDL1 protein expression adjusting for cancer group")

pdf(file = paste0(plotpath,"Fig3B.pdf"),
    width = 15, 
    height = 10,
    useDingbats = FALSE)
pcox
dev.off()

# Fig 3C

genemat <- read.csv(file.path(datapath,"IND_tpm_hg38_final.csv"),header = T, stringsAsFactors = F, check.names = F, row.names = 1)

pdl1_exp <- genemat[which(rownames(genemat) == "CD274"),, drop = T]

metadata$PDL1_gene <- NA
metadata$PDL1_gene <- pdl1_exp[match(metadata$sample_id, names(pdl1_exp))]

metadata$PDL1_gene <- as.numeric(metadata$PDL1_gene)

summary(metadata$PDL1_gene)

tapply(metadata$PDL1_gene, metadata$PDL1group, summary)

metadata$PDL1g <- NA
metadata$PDL1g[ metadata$PDL1_gene >= 4.817] <- "High"
metadata$PDL1g[ metadata$PDL1_gene < 4.817 &
                  metadata$PDL1_gene > 0.475] <- "Intermediate"
metadata$PDL1g[ metadata$PDL1_gene <= 0.475] <- "Low"


sfit <- survfit(Surv(TRTDUR, progressed) ~ PDL1g, data= metadata)

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

kmplot$plot <- kmplot$plot + labs(title = "PD-L1 gene expression\n") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig3C.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()


# Fig 3D

tmp <- metadata[!is.na(metadata$IHC_PDL1TC),]
tmp<- tmp[ !is.na(tmp$IHC_CD8),]

tmp$PDL1TC <- NA
tmp$PDL1TC[ tmp$IHC_PDL1TC == 0] <- "PD-L1 (TC) = 0"
tmp$PDL1TC[ tmp$IHC_PDL1TC > 0] <- "PD-L1 (TC) > 0"

p <- ggplot(data = tmp, aes( x = PDL1TC, y = IHC_CD8)) + 
  geom_beeswarm(alpha = 1, cex = 4, size = 5, color = "grey") + 
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black") + 
  myplot + myaxis +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 30, angle = 0, hjust = 0.5),
        plot.title = element_text(size = 30, hjust = 0.5)) +
  annotation_logticks(sides = "l") + 
  scale_y_continuous(trans = "log10") +
  labs( y = "CD8 staining (IHC)", title = "PD-L1 protein expression\non tumor cells") +
  geom_signif(comparisons = list(c("PD-L1 (TC) = 0", "PD-L1 (TC) > 0")), y_position = 2,
              map_signif_level=TRUE, textsize = 20, test = "t.test", vjust = 0.5)

pdf(file = paste0(plotpath,"Fig3D.pdf"),
    width = 8, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
p
dev.off()


