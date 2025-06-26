rm(list = ls())
path <- "./Simulated_data/ARTISTA/lowly_abundant/"
sce_combined <- readRDS(paste0(path, "combined_simulated.rda"))
source("01_competitors_wrapper.R", echo=TRUE)
##################################### GT ######################################
sce_combined$Condition <- ifelse(sce_combined$Condition == "Condition_1", "Condition1", "Condition2")
run_seurat_FindMarkers(sce_combined,
                          output_path = path,
                          cluster_method = "GT",
                          cluster_col = "Banksy_smooth")

run_seurat_PseudoBulk_FindMarkers(sce_combined,
                                  output_path = path,
                                  cluster_method = "GT",
                                  cluster_col = "Banksy_smooth")

run_spatialLIBD_pairwise(sce_combined,
                       output_path = path,
                       cluster_method = "GT",
                       cluster_col = "Banksy_smooth")

run_scran_findMarkers(sce_combined,
                      output_path = path,
                      cluster_method = "GT",
                      cluster_col = "Banksy_smooth")

#########################################################################
rm(list = ls())
path <- "./Simulated_data/ARTISTA/lowly_abundant/"
spe <- readRDS(file = paste0(path, "combined_simulated_updated.rda")) # see `05_run_DESpace.R`
source("01_competitors_wrapper.R", echo=TRUE)
library(SingleCellExperiment)
library(dplyr)
library(edgeR)
library(muscat)
`%notin%` <- Negate(`%in%`)
spe$Condition <- ifelse(spe$Condition == "Condition_1", "Condition1", "Condition2")

run_seurat_FindMarkers(spe,
                          output_path = path,
                          cluster_method = "Banksy",
                          cluster_col = "Banksy_recomputed")
  
run_seurat_PseudoBulk_FindMarkers(spe,
                                  output_path = path,
                                  cluster_method = "Banksy",
                                  cluster_col = "Banksy_recomputed")
run_spatialLIBD_pairwise(spe,
                       output_path = path,
                       cluster_method = "Banksy",
                       cluster_col = "Banksy_recomputed")

run_scran_findMarkers(spe,
                      output_path = path,
                      cluster_method = "Banksy",
                      cluster_col = "Banksy_recomputed")