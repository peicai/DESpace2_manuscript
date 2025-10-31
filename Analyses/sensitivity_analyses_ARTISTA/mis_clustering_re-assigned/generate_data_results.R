library(limma)
library(BiocParallel)
library(edgeR)
library(scuttle)
library(DESpace)
path <- "./"
spe <- readRDS("./Simulated_data/ARTISTA/highly_abundant/combined_simulated.rda")

`%notin%` <- Negate(`%in%`)

set.seed(123)
# Function to randomly relabel n% of cells
relabel_clusters <- function(labels, frac = 0.1) {
  n <- length(labels)
  k <- unique(labels)
  # number of cells to relabel
  n_change <- ceiling(frac * n)    
  # pick random cells
  idx <- sample(seq_len(n), n_change)   
  new_labels <- labels
  
  for (i in idx) {
    # exclude original cluster
    possible <- setdiff(k, labels[i])  
    new_labels[i] <- sample(possible, 1)
  }
  
  return(new_labels)
}


spe$shuffle_001 <- relabel_clusters(spe$Banksy_recomputed, frac = 0.01)
spe$shuffle_005 <- relabel_clusters(spe$Banksy_recomputed, frac = 0.05)
spe$shuffle_01 <- relabel_clusters(spe$Banksy_recomputed, frac = 0.1)
spe$shuffle_02 <- relabel_clusters(spe$Banksy_recomputed, frac = 0.2)
save(spe, file = paste0(path, "combined_simulated_shuffled.rda"))

## DESpace
for(i in paste0("shuffle_0", c("01", "05", "1", "2"))){
  res_edgeR <- dsp_test(spe,
                        cluster_col = i,
                        sample_col = "Sample",
                        condition_col = "Condition",
                        filter_gene = FALSE,
                        filter_cluster = FALSE,
                        verbose = FALSE)
  save(res_edgeR, file = paste0(path, "DESpace_global_", i, ".rda"))
  
  results_res1 <- individual_dsp(spe,
                                 cluster_col = i,
                                 sample_col = "Sample",
                                 condition_col = "Condition",
                                 filter_gene = FALSE,
                                 filter_cluster = TRUE)
  
  save(results_res1, file = paste0(path, "DESpace_individual_", i, ".rda"))
}

## competitors
source("./sensitive_analyses/01_competitors_wrapper.R", echo=TRUE)
library(SingleCellExperiment)
library(dplyr)
library(edgeR)
library(muscat)
`%notin%` <- Negate(`%in%`)
spe$Condition <- ifelse(spe$Condition == "Condition_1", "Condition1", "Condition2")

for(i in paste0("shuffle_0", c("01", "05", "1", "2"))){
  run_seurat_FindMarkers(spe,
                         output_path = path,
                         cluster_method = i,
                         cluster_col = i)
  
  run_seurat_PseudoBulk_FindMarkers(spe,
                                    output_path = path,
                                    cluster_method = i,
                                    cluster_col = i)
  run_spatialLIBD_pairwise(spe,
                           output_path = path,
                           cluster_method = i,
                           cluster_col = i)
  
  run_scran_findMarkers(spe,
                        output_path = path,
                        cluster_method = i,
                        cluster_col = i)
  
}
