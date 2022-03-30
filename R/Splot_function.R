library(ggplot2)
library(ComplexHeatmap)
library(circlize)

median.cohorts.fx <- function(df){
  
  Div_df.cohorts <- cbind.data.frame(unique(df$group), unique(df$cohort), NA)
  
  colnames(Div_df.cohorts) <- c("group","cohort", "median") 
  
  for(i in 1:nrow(Div_df.cohorts)){
    Div_df.cohorts$median[i] <- median(na.omit(df$CD8_CIBERSORT[df$cohort == Div_df.cohorts$cohort[i]]))
  } 
  
  Div_df.cohorts <- Div_df.cohorts[order(Div_df.cohorts$median),]
  return(Div_df.cohorts)
}


sort.cohorts.fx <- function(df, median_df){
  disease.width <- (nrow(df)/nrow(median_df)) 
  sorted.df <- df[0,]
  start = 0
  for(i in 1:(nrow(median_df))){
    tmp <- df[df$cohort==median_df$cohort[i],]
    tmp <- tmp[order(tmp$CD8_CIBERSORT),]
    tmp <- tmp[!is.na(tmp$CD8_CIBERSORT),] 
    #create range of x values to squeeze dots into equal widths of the plot for each Disease regardless of the number of samples
    div <- disease.width/nrow(tmp)
    #If there is only one sample, put the dot in the middle of the alloted space
    if(dim(tmp)[1]==1){
      tmp$Xpos<-start+(disease.width/2)
    } 
    else tmp$Xpos <- seq(from = start, 
                         to = start+disease.width, 
                         by = div)[-1]
    sorted.df <- rbind(sorted.df, tmp)  
    median_df$Median.start[i] <- tmp$Xpos[1]
    median_df$Median.stop[i] <- tmp$Xpos[nrow(tmp)]
    median_df$N[i]<-nrow(tmp)
    start <- start+disease.width+30
    
  }
  median_df$medianloc <- median_df$Median.start+
    ((median_df$Median.stop-median_df$Median.start)/2)
  sorted.df$cohort <- factor(sorted.df$cohort,
                             levels = median_df$cohort)
  return(list(sorted.df,median_df))
}

Splot.fx <- function(list.sorted_df.median, plottitle){       
  sorted_df <- as.data.frame(list.sorted_df.median[1])
  median_df <- as.data.frame(list.sorted_df.median[2])
  disease.width <- (nrow(sorted_df)/nrow(median_df)) 
  
  colpal <- c("High" = "Red",
              "Low" = "Blue",
              "Intermediate" = "Light grey")
  
  Splot <- ggplot() +
    geom_point(data = sorted_df, aes(x = Xpos ,y = CD8_CIBERSORT, color = CD8Level), size = 5) +
    geom_crossbar(data = median_df, 
                  aes(x =medianloc, y = median,
                      ymin = median, ymax = median),
                  width = disease.width) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 22),
          axis.title = element_text(size = 22), 
          plot.title = element_text(size=22),
          legend.position = "none") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",
                                          colour = NA),
          panel.border=element_blank(),
          plot.margin = unit(c(1.2,1,0,1),"cm")) +
    scale_color_manual(values = colpal) +
    scale_x_continuous(breaks = seq((disease.width)/2,max(sorted_df$Xpos),
                                    disease.width+30),
                       labels = median_df$cohort,
                       expand = c(0,20)) + 
    labs(y = "CD8 gene signature (CIBERSORT)", title = plottitle) 
  return(Splot)  
}




