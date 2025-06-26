rm(list = ls())
library(BayesSpace)
library(Matrix)
library(SingleCellExperiment)
library(dplyr)
`%notin%` <- Negate(`%in%`)

sample_id <- list(c("151507", "151509"),
                  c("151669", "151671"),
                  c("151673", "151675"))
data_path <- "./"
load(paste0(data_path, 'Combined_simulated.rda'))

for(i in c(1,3,4,5,6)){
  sce <- sce_list[[i]]
  set.seed(123)
  dec <- scran::modelGeneVar(sce)
  n=2000
  top <- scran::getTopHVGs(dec, n=n)
  
  sce <- scater::runPCA(sce, subset_row=top)
  platform <- "Visium"
  sce <- spatialPreprocess(sce, platform=platform, skip.PCA=TRUE)
  q=5
  colData(sce)$col <- colData(sce)$imagecol
  colData(sce)$row <- colData(sce)$imagerow
  
  sce <- spatialCluster(sce, q=q,platform=platform,init.method = c("mclust"),
                        model = c("t"), save.chain=TRUE)
  save(sce, file = paste0("./simulated_sample_", i, "_BayesSpace.rda"))
  
}

rm(sce_list)

file_paths <- sprintf("./simulated_sample_%d_BayesSpace.rda", 1:6)
# Create a list to store the loaded data
sce_list <- vector("list", length(file_paths))

# Use lapply to load the files and store them in sce_list
sce_list <- lapply(seq_along(file_paths), function(i) {
  load(file_paths[i]) 
  sce 
})

### manually match BayesSpace cluters across samples
###################################### plot ################################
library(ggplot2)
library(SingleCellExperiment) 

create_ggplot <- function(sce) {
  cell_meta <- colData(sce)
  ggplot(cell_meta, aes(x = row, y = col, color = factor(spatial.cluster))) +
    geom_point() +
    labs(x = "X", y = "Y", color = "Pattern") +
    theme_minimal() + scale_x_reverse()
}

plots <- lapply(sce_list, create_ggplot)

for (i in seq_along(plots)) {
  plot_name <- paste0("BayesSpace_", i, ".png")
  ggsave(plot = plots[[i]], filename = plot_name, width = 5, height = 5)
}

############################## match cluster ###########################
# Define different mappings for each element in sce_list
cluster_mappings <- list(
  c("1" = 4, "2" = 2, "3" = 1, "4" = 5, "5" = 3),
  c("1" = 1, "2" = 5, "3" = 2, "4" = 3, "5" = 4),
  c("1" = 1, "2" = 3, "3" = 5, "4" = 2, "5" = 4),
  c("1" = 5, "2" = 3, "3" = 1, "4" = 2, "5" = 4),
  c("1" = 2, "2" = 5, "3" = 1, "4" = 4, "5" = 3),
  c("1" = 1, "2" = 5, "3" = 3, "4" = 2, "5" = 4)
)

apply_cluster_mapping <- function(sce, mapping) {
  df <- colData(sce)  
    df$matched_cluster <- as.numeric(recode(as.character(df$spatial.cluster), !!!mapping))
    colData(sce) <- df
  return(sce)
}
sce_list_matched <- Map(apply_cluster_mapping, sce_list, cluster_mappings)

save(sce_list_matched, file = paste0("./simulated_BayesSpace_matched.rda"))