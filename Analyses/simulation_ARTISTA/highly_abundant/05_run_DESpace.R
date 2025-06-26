rm(list = ls())
library(SingleCellExperiment)
library(dplyr)
library(SeuratObject)
library(DESpace) # Need to use R 4.5 or manually load source code
library(Seurat)
library(limma)
library(edgeR)
`%notin%` <- Negate(`%in%`)
path <- "./Simulated_data/ARTISTA/highly_abundant/"
############################################ GT #######################################
sce_combined <- readRDS(paste0(path, "combined_simulated.rda"))
table(sce_combined$Banksy_smooth)
res_edgeR <- dsp_test(spe = sce_combined,
                          cluster_col = "Banksy_smooth",
                          sample_col = "Sample",
                          condition_col = "Condition",
                          filter_gene = FALSE,
                          filter_cluster = FALSE,
                          verbose = FALSE)
save(res_edgeR, file = paste0(path, "DESpace_global_GT.rda"))

######################### individual test #############################
results_res1 <- individual_dsp(sce_combined,
                                        cluster_col = "Banksy_smooth",
                                        sample_col = "Sample",
                                        condition_col = "Condition",
                                        filter_gene = FALSE,
                                        filter_cluster = TRUE)

save(results_res1, file = paste0(path, "DESpace_individual_GT.rda"))

############################################ Banksy cluster #######################################
# after recomuting spatial clusters (in `03_Banksy_cluster.R`)
rm(list = ls())
library(SingleCellExperiment)
library(dplyr)
library(SeuratObject)
library(DESpace)
library(Seurat)
library(limma)
library(edgeR)
library(ggplot2)
`%notin%` <- Negate(`%in%`)
path <- "./Simulated_data/ARTISTA/highly_abundant/"
load(paste0(path, "Banksy_harmony_smooth.rda"))
spe$Banksy_recomputed <- gsub("Layer", "", spe$Banksy_recomputed)
sce_combined <- readRDS(paste0(path, "combined_simulated.rda"))
sce_combined$Truth <- sce_combined$Banksy_smooth
sce_combined$Banksy_recomputed <- spe$Banksy_recomputed
sce_combined$Sample <- spe$sample_id
sce_combined$Condition <- spe$subject_id
saveRDS(sce_combined, file = paste0(path, "combined_simulated_updated.rda"))

spe <- readRDS(file = paste0(path, "combined_simulated_updated.rda"))
library(SpatialExperiment)
# -------------------------------------------------------------------------
res_edgeR <- dsp_test(spe = spe,
                          cluster_col = "Banksy_recomputed",
                          sample_col = "Sample",
                          condition_col = "Condition",
                          filter_gene = FALSE,
                          filter_cluster = FALSE,
                          verbose = FALSE)
save(res_edgeR, file = paste0(path, "DESpace_global_Banksy.rda"))
######################### individual test #############################
results_res1 <- individual_dsp(sce_combined,
                               cluster_col = "Banksy_recomputed",
                               sample_col = "Sample",
                               condition_col = "Condition",
                               filter_gene = FALSE,
                               filter_cluster = TRUE)

save(results_res1, file = paste0(path, "DESpace_individual_Banksy.rda"))

