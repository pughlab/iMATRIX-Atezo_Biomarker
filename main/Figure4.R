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
source(paste0(mainpath, "R/Heatmap_functions.R"))
source(paste0(mainpath, "R/Splot_function.R"))

#metadata <- read.csv(file.path(datapath,"metadata.csv"), header = T, stringsAsFactors = F, check.names = F)
#metadata <- metadata[ !is.na(metadata$observed_Shannon),]

#Fig4A

compldfle <- read.csv(file.path(datapath,"anonymized_TableS2.csv"), header = T, stringsAsFactors = F, check.names = F)

compldfle$samplename <- as.factor(compldfle$trunc_anonymized_sample_id)

compldfle$cloneno <- as.factor(compldfle$cloneno)

clonenocol = rep("#ffffff", nlevels(compldfle$cloneno))

clonpt <- ggplot(data = compldfle, aes(y = cloneFraction, x = samplename, fill = cloneno)) + 
  geom_bar(colour = "#000000", stat="identity", width = 0.8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20), 
        legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA),
        panel.border=element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")) +
  scale_fill_manual(values = clonenocol,guide = FALSE) 




metadata_trb$Div_group <- NA
metadata_trb$Div_group[metadata_trb$observed_Shannon >= 80.269] <- "High"
metadata_trb$Div_group[metadata_trb$observed_Shannon < 80.269 &
                         metadata_trb$observed_Shannon > 9.527] <- "Intermediate"
metadata_trb$Div_group[metadata_trb$observed_Shannon <= 9.527] <- "Low"


colpal <- c("High" = "#ED2024",
            "Intermediate" = "#adadad",
            "Low" = "#3953A4")

divplot <- ggplot(data = metadata_trb, aes(y = observed_Shannon, x = sample_id_DNA)) + 
  geom_bar(aes(fill = Div_group),colour = "#000000", stat="identity", width = 0.8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20), 
        legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border=element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")) 

divplot <- divplot + annotation_logticks(sides = "l") +
  scale_y_continuous(trans = "log10") + scale_fill_manual(values = colpal)


# Add number of sample specific gliph types
metadata_trb$gliph2_type <- NA

for(i in 1:nrow(metadata_trb)){
  mysample <- metadata_trb$sample_id_DNA[i]
  mygliph <- gliph2[gliph2$Sample == mysample,]
  gliphtab <- as.data.frame(table(mygliph$type), stringsAsFactors = F)
  #count number of gliph types that have more than one cdr3 in a given sample    
  if(nrow(gliphtab) !=0 ){
    metadata_trb$gliph2_type[i] <- nrow(gliphtab[gliphtab$Freq > 1,])}
}

#convert NAs to zero
metadata_trb$gliph2_type[ is.na(metadata_trb$gliph2_type)] <- 0

gliphplot <- ggplot(data = metadata_trb, aes(y = gliph2_type, x = sample_id_DNA)) + 
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
        plot.margin = unit(c(0,0,0,0),"cm")) 

pdf(file = paste0(plotpath,"TRB_clonplot.pdf"),
    width = 20, 
    height = 12,
    useDingbats = FALSE, 
    onefile = F)
align_plots1(divplot, clonpt, gliphplot)
dev.off()

# Fig4B

#log10
metadata$log10shann <- log10(metadata$observed_Shannon)

message("summary of log10 shannon diversity")
summary(metadata$log10shann)

metadata$Div_group <- NA
metadata$Div_group[metadata$log10shann >= 1.905] <- "High"
metadata$Div_group[metadata$log10shann < 1.905 &
                     metadata$log10shann > 0.979] <- "Intermediate"
metadata$Div_group[metadata$log10shann <= 0.979] <- "Low"

sfit <- survfit(Surv(TRTDUR, progressed)~ Div_group, data= metadata)

kmplot <- ggsurvplot(sfit, conf.int=FALSE, palette = c("#ED2024", "#adadad", "#3953A4"),
                     pval = TRUE, pval.size = 8, pval.coord = c(300, 0.75),
                     risk.table=TRUE, fontsize = 8,
                     legend = c(0.7, 0.95),font.legend = 16, legend.title = "",
                     font.main = 23, font.x = 25,font.y = 25, font.tickslab = 20) 

kmplot$table <- kmplot$table + theme(axis.text.x = element_text(size = 20))
kmplot$plot <- kmplot$plot + labs(title = "PFS analysis of TRB diversity")

pdf(file = paste0(plotpath,"Fig4B.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE, onefile = FALSE)
kmplot
dev.off()


