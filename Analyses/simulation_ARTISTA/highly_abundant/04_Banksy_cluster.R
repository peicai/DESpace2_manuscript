rm(list = ls())
library(Matrix)
library(SingleCellExperiment)
library(dplyr)
`%notin%` <- Negate(`%in%`)
library(Banksy)
library(SummarizedExperiment)
library(SpatialExperiment)
library(harmony)
library(data.table)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)

SEED <- 202406
path <- "./Simulated_data/ARTISTA/highly_abundant/"
sce_combined <- readRDS(paste0(path, "combined_simulated.rda"))
######################### pre process data ####################################
spe <- SpatialExperiment::SpatialExperiment(
  rowData = rowData(sce_combined),
  colData = colData(sce_combined),
  assays = assays(sce_combined),
  reducedDims = NULL,
  sample_id = sce_combined$Sample,
  spatialCoordsNames = c("sdimx", "sdimy"),
  imgData = NULL
)

assay(spe, "logcounts") <- NULL
reducedDims(spe) <- NULL
rowData(spe) <- NULL
colData(spe) <- DataFrame(
  sample_id = spe$Sample,
  subject_id = factor(spe$Condition),
  row.names = colnames(spe)
)
invisible(gc())


#' Stagger spatial coordinates
#' stagger the spatial coordinates across the samples so that spots from different samples do not overlap.
locs <- data.frame(colData(sce_combined)[, c("sdimx", "sdimy")])
locs <- cbind(locs, sample_id = as.numeric(factor(spe$sample_id)))
locs_dt <- data.table(locs)
colnames(locs_dt) <- c("sdimx", "sdimy", "group")
locs_dt[, sdimx := sdimx - min(sdimx), by = group]
global_max <- max(locs_dt$sdimx) * 1.5
locs_dt[, sdimx := sdimx + group * global_max]
locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(spe)
spatialCoords(spe) <- locs


#' Get HVGs
library(Seurat)
sample_ids <- unique(spe$sample_id)
spe_list <- lapply(sample_ids, function(x) spe[, spe$sample_id == x])

#' Normalize data
seu_list <- lapply(spe_list, function(x) {
  x <- as.Seurat(x, data = NULL)
  NormalizeData(x, scale.factor = 5000, normalization.method = 'RC')
})

#' Compute HVGs
hvgs <- lapply(seu_list, function(x) {
  VariableFeatures(FindVariableFeatures(x, nfeatures = 2000))
})
hvgs <- Reduce(union, hvgs)

#' Add data to SpatialExperiment and subset to HVGs
aname <- "normcounts"
spe_list <- Map(function(spe, seu) {
  assay(spe, aname, withDimnames=FALSE) <- GetAssayData(seu)
  spe[hvgs,]
}, spe_list, seu_list)

invisible(gc())

compute_agf <- FALSE
k_geom <- 18
spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, 
                   compute_agf = compute_agf, k_geom = k_geom)

spe_joint <- do.call(cbind, spe_list)
seu_joint <- as.Seurat(spe_joint, data = NULL)
invisible(gc())
group_name = c("sample_id")
dimx = "sdimx"
dimy = "sdimy"
col_data <- cbind(colData(spe_joint), spatialCoords(spe_joint))
seu_joint <- AddMetaData(seu_joint, metadata = as.data.frame(col_data))
head(seu_joint@meta.data)
seu_joint = RunBanksy(
  seu_joint, lambda = 0.8, assay = 'originalexp', slot = 'data',
  group = group_name, dimx = dimx, dimy = dimy,features = "all",
  k_geom = 18, split.scale = TRUE)
seu_joint = RunPCA(seu_joint, features = rownames(spe_joint), 
                   npcs = 20)
seu_joint = RunHarmony(seu_joint, group.by.vars=group_name,
                       reduction.name='pca', reduction.save='harmony')

A1 <- plot_grid(
  DimPlot(seu_joint, reduction = 'pca'),
  DimPlot(seu_joint, reduction = 'harmony')
)

# Cluster BANKSY-Harmony embedding
seu_joint = FindNeighbors(seu_joint, dims = 1:20, reduction = 'harmony')
seu_joint = FindClusters(seu_joint, resolution = 0.01, graph.name = 'BANKSY_snn')

spe_joint = spe
spe_joint$BANKSY_snn_res0.3 <- seu_joint@meta.data$BANKSY_snn_res.0.3
spe_list <- lapply(sample_ids, function(x) spe_joint[, spe_joint$sample_id == x])
spe_list <- lapply(spe_list, smoothLabels, cluster_names = "BANKSY_snn_res0.3", k = 6L, verbose = FALSE)
names(spe_list) <- paste0("sample_", sample_ids)
table(colData(spe_list[[1]])$BANKSY_snn_res0.3_smooth)
new_labels <- c("Layer3", "Layer2", "Layer1", "Layer1", "Layer0", "Layer0", "Layer4", "Layer4")
spe <- lapply(spe_list, function(df) {
  # Check if 'BANKSY_snn_res0.3_smooth' exists in colData
  if("BANKSY_snn_res0.3_smooth" %in% colnames(colData(df))) {
    colData(df)$Banksy_recomputed <- factor(
      colData(df)$BANKSY_snn_res0.3_smooth,
      levels = 0:7,
      labels = new_labels
    )
  }
  return(df)
})

spe_joint <- do.call(cbind, spe)

seu_joint@meta.data$Banksy_recomputed <- spe_joint$Banksy_recomputed

spe <- spe_joint
save(spe, file = paste0(output_dir, "Banksy_harmony_smooth.rda"))