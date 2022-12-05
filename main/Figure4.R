###############
# Figure 4
###############

library(survminer)
library(survival)

mainpath <- "~/git/iMATRIX-Atezo_Biomarker/"
datapath <- paste0(mainpath,"data/")

dir.create(file.path(paste0(mainpath, "out/")))

plotpath <- paste0(mainpath,"out/")

setwd(mainpath)

source(paste0(mainpath, "R/ggplot2_theme.R"))
source(paste0(mainpath, "R/plotting_functions.R"))
source(paste0(mainpath, "R/Heatmap_functions.R"))


metadata <- read.csv(file.path(datapath,"anonymized_iMATRIX_Atezo_metadata_IHC_TRB_TMB_v3.csv"),header = T, stringsAsFactors = F, check.names = F)

metadata <- metadata[ !is.na(metadata$observed_Shannon),]
metadata <- metadata[order(metadata$cancer, metadata$observed_Shannon), ]
metadata$trunc_anonymized_tcrseq_sample_id <- factor(metadata$trunc_anonymized_tcrseq_sample_id, levels = metadata$trunc_anonymized_tcrseq_sample_id)

#Fig4A
# diversity plot
summary(metadata$observed_Shannon)

metadata$Div_group <- NA
metadata$Div_group[metadata$observed_Shannon >= 80.269] <- "High"
metadata$Div_group[metadata$observed_Shannon < 80.269 &
                     metadata$observed_Shannon > 9.527] <- "Intermediate"
metadata$Div_group[metadata$observed_Shannon <= 9.527] <- "Low"

colpal <- c("High" = "#ED2024",
            "Intermediate" = "#adadad",
            "Low" = "#3953A4")

divplot <- ggplot(data = metadata, 
                  aes(y = observed_Shannon, 
                      x = trunc_anonymized_tcrseq_sample_id)) + 
  geom_bar(aes(fill = Div_group),colour = "#000000", stat="identity", width = 0.8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20), 
        legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border=element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")) 

divplot <- divplot + annotation_logticks(sides = "l") +
  scale_y_continuous(trans = "log10") + scale_fill_manual(values = colpal) + labs(y = "TCRb diversity")

# clone plot
compldfle <- read.csv(file.path(datapath,"anonymized_TableS2.csv"), header = T, stringsAsFactors = F, check.names = F)

compldfle$trunc_anonymized_sample_id <- as.factor(compldfle$trunc_anonymized_sample_id)
compldfle$cloneno <- as.factor(compldfle$cloneno)
compldfle$trunc_anonymized_sample_id <- factor(compldfle$trunc_anonymized_sample_id, levels = metadata$trunc_anonymized_tcrseq_sample_id)

clonenocol <- rep("#ffffff", nlevels(compldfle$cloneno))

clonpt <- ggplot(data = compldfle, aes(y = cloneFraction, x = trunc_anonymized_sample_id, fill = cloneno)) + 
  geom_bar(colour = "#000000", stat="identity", width = 0.8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20), 
        legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border=element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")) +
  scale_fill_manual(values = clonenocol,guide = FALSE) +
  labs(y = "TCRb clonal\nfraction")

#Gliph2 plot
gliph2 <- read.csv(file.path(datapath,"anonymized_TableS3.csv"), header = T, stringsAsFactors = F, check.names = F)

# Add number of sample specific gliph types
metadata$gliph2_type <- NA
for(i in 1:nrow(metadata)){
  mysample <- metadata$trunc_anonymized_tcrseq_sample_id[i]
  mygliph <- gliph2[ gliph2$trunc_anonymized_sample_id == mysample,]
  gliphtab <- as.data.frame(table(mygliph$type), stringsAsFactors = F)
  #count number of gliph types that have more than one cdr3 in a given sample    
  if(nrow(gliphtab) !=0 ){
    metadata$gliph2_type[i] <- nrow(gliphtab[gliphtab$Freq > 1,])}
}

#convert NAs to zero
metadata$gliph2_type[ is.na(metadata$gliph2_type)] <- 0

gliphplot <- ggplot(data = metadata, aes(y = gliph2_type, x = trunc_anonymized_tcrseq_sample_id)) + 
  geom_bar(colour = "#000000", fill = "#adadad", stat="identity", width = 0.8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20), 
        legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border=element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")) + labs(y = "TCRb specificity\ngroups")

#align and export
pdf(file = paste0(plotpath,"Fig4A_1.pdf"),
    width = 20, 
    height = 12,
    useDingbats = FALSE, 
    onefile = F)
align_plots1(divplot, clonpt, gliphplot)
dev.off()

# Heatmap for sample origin
metadata$origin <- NA
metadata$origin[ metadata$sample_origin == "Lymph node"] <- "Lymph node"
metadata$origin[ metadata$sample_origin != "Lymph node"] <- "Other tissue"
myorigin <- metadata$origin
names(myorigin) <- metadata$trunc_anonymized_tcrseq_sample_id
myorigin <- t(as.matrix(myorigin))
rownames(myorigin) <- "Lymph node\n/Other tissue"

origin_hm <- origin_hm.fx(myorigin)

# Heatmap for tumour type
mytype <- metadata$tumor_type
names(mytype) <- metadata$trunc_anonymized_tcrseq_sample_id
mytype <- t(as.matrix(mytype))
rownames(mytype) <- "Primary/Metastatic"
type_hm <- type_hm.fx(mytype)

#export
pdf(paste0(plotpath, "Fig4A_2.pdf"),
    width = 10, height = 20,
    useDingbats = FALSE)
draw(type_hm %v% origin_hm,heatmap_legend_side = "bottom")
dev.off()

# Fig4B
metadata$Div_group <- NA
metadata$Div_group[metadata$observed_Shannon >= 80.269] <- "High"
metadata$Div_group[metadata$observed_Shannon < 80.269 &
                     metadata$observed_Shannon > 9.527] <- "Intermediate"
metadata$Div_group[metadata$observed_Shannon <= 9.527] <- "Low"

sfit <- survfit(Surv(TRTDUR, progressed)~ Div_group, data= metadata)

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

kmplot$plot <- kmplot$plot + labs(title = "TCR diversity\n") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig4B.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE, onefile = FALSE)
kmplot
dev.off()

rm(kmplot)

#Fig4C
# remove lymph nodes
tmp <- metadata[ metadata$sample_origin != "Lymph node",]

#remove intermediate group to compare High to Low
tmp$Div_group[ tmp$Div_group == "Intermediate"] <- NA

sfit <- survfit(Surv(TRTDUR, progressed)~ Div_group, data= tmp)

kmplot <- ggsurvplot(sfit, 
                     conf.int=FALSE, palette = c("#ED2024", "#3953A4"),
                     pval = TRUE, pval.size = 10, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 20,
                     xlim = c(0,800),break.x.by = 200,
                     legend = c(0.7, 0.95), font.legend = 35, legend.title = "",
                     legend.labs = c("High", "Low"),
                     font.main = 35, font.x = 35, font.y = 35, font.tickslab = 35,
                     ylab = "Progression-free survival") 

kmplot$table <- ggrisktable(sfit, data = tmp, 
                            color = "strata", 
                            palette = c("#ED2024", "#3953A4"),
                            fontsize = 12, 
                            xlim = c(0,800),break.x.by = 200, 
                            tables.theme = theme_cleantable(), font.tickslab = 30, 
                            y.text = TRUE, ylab = "",  xlab = "Time (days)",
                            legend.labs = c("High", "Low"))

kmplot$table <- kmplot$table + theme(plot.title = element_blank(),
                                     axis.text.x = element_text(size = 35), axis.title.x = element_text(size = 35))

kmplot$plot <- kmplot$plot + labs(title = "TCR diversity\nexcluding lymph nodes") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), legend.key.size = unit(1, 'cm'))

pdf(file = paste0(plotpath,"Fig4C.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE,
    onefile = FALSE)
kmplot
dev.off()

rm(kmplot)


