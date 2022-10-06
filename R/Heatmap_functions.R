library(ggplot2)
library(ComplexHeatmap)
library(circlize)

age_hm.fx <- function(age_mat){
  
  col_fun= colorRamp2(c(0, 18), c("#ffffff", "#000000"))
  
  age_hm = Heatmap(age_mat,
                   #titles and names
                   name = "Age",
                   show_row_names = TRUE,
                   show_column_names = FALSE,    
                   #clusters
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   #aesthestics
                   col = col_fun,
                   row_names_gp = gpar(fontsize = 10),
                   height = unit(1, "cm"),
                   column_title_gp = gpar(fontsize = 10),
                   row_title = NULL)
  return(age_hm)
}

type_hm.fx <- function(tumortype_mat){
  
  colpal <- c("Primary" = "#d8daeb",
              "Metastatic" = "#542788")
  
  type_hm = Heatmap(tumortype_mat,
                    #titles and names
                    name = "Tumour type",
                    show_row_names = TRUE,
                    show_column_names = FALSE,    
                    #clusters
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    #aesthestics
                    col = colpal,
                    row_names_gp = gpar(fontsize = 10),
                    height = unit(1, "cm"),
                    column_title_gp = gpar(fontsize = 10),
                    row_title = NULL)
  return(type_hm)
}


cancer_hm.fx <- function(cancer_mat){
  
  colpal <- c("Lymphoma" = "#1f78b4",
              "Solid tumour" = "#a6cee3")
  
  cancer_hm = Heatmap(cancer_mat,
                      #titles and names
                      name = "Cancer group",
                      show_row_names = TRUE,
                      show_column_names = FALSE,    
                      #clusters
                      cluster_columns = FALSE,
                      cluster_rows = FALSE,
                      #aesthestics
                      col = colpal,
                      row_names_gp = gpar(fontsize = 10),
                      height = unit(1, "cm"),
                      column_title_gp = gpar(fontsize = 10),
                      row_title = NULL)
  return(cancer_hm)
}



origin_hm.fx <- function(origin_mat){
  
  colpal <- c("Lymph node" = "#b15928",
              "Other tissue" = "#ffff99")
  
  origin_hm = Heatmap(origin_mat,
                      #titles and names
                      name = "Sample origin",
                      show_row_names = TRUE,
                      show_column_names = FALSE,    
                      #clusters
                      cluster_columns = FALSE,
                      cluster_rows = FALSE,
                      #aesthestics
                      col = colpal,
                      row_names_gp = gpar(fontsize = 10),
                      height = unit(1, "cm"),
                      column_title_gp = gpar(fontsize = 10),
                      row_title = NULL)
  return(origin_hm)
}





response_hm.fx <- function(response_mat){
  
  colpal <- c("PD" = "#e41a1c",
              "PR" = "#4daf4a",
              "SD" = "#ffff33",
              "NE" = "grey")
  
  response_hm = Heatmap(response_mat,
                        #titles and names
                        name = "Response",
                        show_row_names = TRUE,
                        show_column_names = FALSE,    
                        #clusters
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        #aesthestics
                        col = colpal,
                        row_names_gp = gpar(fontsize = 10),
                        height = unit(1, "cm"),
                        column_title_gp = gpar(fontsize = 10),
                        row_title = NULL)
  return(response_hm)
  
}



cohorts_hm.fx <- function(cohorts_mat){
  
  colpal <- c("Lymphoma" = "#006d2c", #green
              "NBL" = "#cccc00", #yellow
              "OS" = "#f1eef6",   # blue      
              "RMS" = "#bdc9e1",  # blue
              "SARC" = "#74a9cf", # blue
              "EWS" = "#045a8d",  # blue
              "WILMS" = "#9F1214", #red
              "RT" = "#884692", #purple
              "Other_PDL1pos" = "#757172", #grey
              "Other_PDL1neg" = "#b0aeae") #grey
  
  cohort_hm = Heatmap(mycohorts,
                      #titles and names
                      name = "Cohort",
                      show_row_names = TRUE,
                      show_column_names = FALSE,    
                      #clusters
                      cluster_columns = FALSE,
                      cluster_rows = FALSE,
                      #aesthestics
                      col = colpal,
                      column_names_gp = gpar(fontsize = 10),
                      height = unit(1, "cm"),
                      row_names_gp = gpar(fontsize = 10))
  return(cohort_hm)   
}




ihc_hm.fx <- function(ihc_mat){
  
  col_fun= colorRamp2(c(0, 50), c("#fed976", "#f03b20"))
  
  ihc_hm = Heatmap(ihc_mat,
                   #titles and names   
                   name = "IHC % cells",   
                   show_row_names = TRUE,
                   show_column_names = FALSE,     
                   #clusters and orders  
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   show_column_dend = TRUE,
                   #row_dend_width = unit(8, "cm"),
                   #aesthestics
                   col = col_fun,
                   column_names_gp = gpar(fontsize = 10),
                   row_names_gp = gpar(fontsize = 10),
                   height = unit(nrow(ihc_mat), "cm"),
                   column_title_gp = gpar(fontsize = 10),
                   column_title = NULL,
                   row_title = NULL)
  return(ihc_hm)   
}


























