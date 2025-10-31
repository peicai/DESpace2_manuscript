rm(list = ls())
library(Matrix)
library(SingleCellExperiment)
library(dplyr)
library(Banksy)
library(SummarizedExperiment)
library(SpatialExperiment)
library(harmony)
library(data.table)
library(ggplot2)
library(cowplot)
library(Seurat)
library(SeuratWrappers)
`%notin%` <- Negate(`%in%`)
SEED <- 202406
cluster_cols <- paste0("nclust_", c(6:10,3,4)) 
output_path <- "./Simulated_data/mis_cluster_v2/"
source("./sensitive_analyses/search_res.R", echo=TRUE)

extract_nclust <- function(seu_object, 
                           cluster_col, 
                           sample_col = "sample_id") {
  print(colnames(seu_object@meta.data))
  print(cluster_col)
  df <- seu_object@meta.data[, c(cluster_col, sample_col), drop = FALSE]
  nclust_per_sample <- tapply(df[[cluster_col]], df[[sample_col]], function(x) length(unique(x)))
  min(nclust_per_sample)
}
custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                   "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", 
                   "#FC8D62", "#8DA0CB", "#E78AC3")
sample_order <- c("2DPI_2", "5DPI_1", "10DPI_2", "20DPI_1",
                  "2DPI_1", "5DPI_2", "10DPI_1", "20DPI_2")

sce_combined <- readRDS("./Simulated_data/ARTISTA/highly_abundant/combined_simulated.rda")

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
dimx = "sdimx"; dimy = "sdimy"
col_data <- cbind(colData(spe_joint), spatialCoords(spe_joint))
seu_joint <- AddMetaData(seu_joint, metadata = as.data.frame(col_data))
seu_joint = RunBanksy(
  seu_joint, lambda = 0.8, assay = 'originalexp', slot = 'data',
  group = group_name, dimx = dimx, dimy = dimy,features = "all",
  k_geom = 18, split.scale = TRUE)
seu_joint = RunPCA(seu_joint, features = rownames(spe_joint), 
                   npcs = 50, assay = 'BANKSY')
seu_joint = RunHarmony(seu_joint, group.by.vars=group_name,
                       reduction.name='BANKSY', reduction.save='harmony')
seu_joint = RunUMAP(seu_joint, dims = 1:30, reduction = "harmony")
# Cluster BANKSY-Harmony embedding
seu_joint = FindNeighbors(seu_joint, dims = 1:30, reduction = 'harmony')
for(cluster_col in cluster_cols){
  cluster_path <- file.path(output_path, cluster_col)
  file_path <- file.path(output_path, cluster_col, "combined_simulated_Banksy.rds")
  if (!dir.exists(cluster_path)) {
    dir.create(cluster_path, recursive = TRUE)
  }
  
  ######################### pre process data ####################################
  # Loop over target cluster numbers
  n_clust_target <- as.numeric(sub("nclust_", "", as.character(cluster_col)))
  # search across resolution
  result <- tryCatch({
    binary_search(
      seu_joint, 
      n_clust_target = n_clust_target, 
      extract_nclust = extract_nclust,
      do_clustering = FindClusters, 
      graph.name = "BANKSY_snn"
    )
  }, warning = function(w) {
    message("Warning for n_clust_target = ", n_clust_target, ": ", conditionMessage(w))
    invokeRestart("muffleWarning")
  }, error = function(e) {
    message("Error for n_clust_target = ", n_clust_target, ": ", conditionMessage(e))
    return(NULL)  # skip to next target cluster numbers
  })
  
  if (is.null(result)) next  # skip if failed
  
  # Extract last column from metadata
  last_col <- tail(colnames(result@meta.data), 1)
  sce_combined$Banksy_cluster <- result@meta.data[[last_col]]
  saveRDS(sce_combined, file = file.path(output_path, cluster_col, "combined_simulated_Banksy.rds"))
  #}
  ## visualization
  df <- data.frame(colData(sce_combined)) %>%
    dplyr::mutate(Sample = factor(Sample, levels = sample_order))
  
  A1 <- ggplot(df, aes(x = sdimx, y = sdimy, col = Banksy_cluster)) + 
    geom_point(size = 0.1) +
    scale_color_manual(values = custom_colors) +
    theme_classic() +
    facet_wrap(~ Sample, scales = "free", ncol = 4, nrow = 2) +
    guides(colour = guide_legend(override.aes = list(size = 3))) 
  
  ggsave(file.path(output_path, cluster_col,  "Banksy.png"), plot = A1, 
         width = 15, height = 8)
}
saveRDS(seu_joint, file = file.path(output_path, "seu_joint.rds"))

### For cluster_col %in% c('nclust_3', 'nclust_4')
seu_joint = RunBanksy(
  seu_joint, lambda = 0.8, assay = 'originalexp', slot = 'data',
  group = group_name, dimx = dimx, dimy = dimy,features = "all",
  k_geom = 18, split.scale = TRUE)
seu_joint = RunPCA(seu_joint, features = rownames(spe_joint), 
                   npcs = 50, assay = 'BANKSY')
seu_joint = RunHarmony(seu_joint, group.by.vars=group_name,
                       reduction.name='BANKSY', reduction.save='harmony')
seu_joint = RunUMAP(seu_joint, dims = 1:20, reduction = "harmony")
seu_joint = FindNeighbors(seu_joint, dims = 1:20, reduction = 'harmony',
                          k.param = 90)

library(dplyr)

## cluster_col == 'nclust_3'
result <- binary_search(
  seu_joint, 
  n_clust_target = 5, 
  extract_nclust = extract_nclust,
  do_clustering = FindClusters, 
  graph.name = "BANKSY_snn"
)

# last_col <- tail(colnames(result@meta.data), 1)
# sce_combined$Banksy_cluster <- result@meta.data[[last_col]]

clu <- result@meta.data[["BANKSY_snn_res.0.0625"]]
centroids <- as.data.frame(emb) %>%
  mutate(cluster = as.character(clu)) %>%
  group_by(cluster) %>%
  summarize(across(everything(), mean), .groups = "drop")

rownames(centroids) <- centroids$cluster
mat <- as.matrix(centroids[ , setdiff(colnames(centroids), "cluster")])

hc <- hclust(dist(mat), method = "average")
merged_k <- 4
merge_map <- cutree(hc, k = merged_k)

sce_combined$Banksy_cluster <- factor(merge_map[clu])
table(original = clu, merged = sce_combined$Banksy_cluster)

saveRDS(sce_combined, file = file.path(output_path, cluster_col, "combined_simulated_Banksy.rds"))


## cluster_col == 'nclust_4'
result <- binary_search(
  seu_joint, 
  n_clust_target = 6, 
  extract_nclust = extract_nclust,
  do_clustering = FindClusters, 
  graph.name = "BANKSY_snn"
)

# last_col <- tail(colnames(result@meta.data), 1)
# sce_combined$Banksy_cluster <- result@meta.data[[last_col]]

clu <- result@meta.data[["BANKSY_snn_res.0.5"]]
centroids <- as.data.frame(emb) %>%
  mutate(cluster = as.character(clu)) %>%
  group_by(cluster) %>%
  summarize(across(everything(), mean), .groups = "drop")

rownames(centroids) <- centroids$cluster
mat <- as.matrix(centroids[ , setdiff(colnames(centroids), "cluster")])

hc <- hclust(dist(mat), method = "average")
merged_k <- 5
merge_map <- cutree(hc, k = merged_k)

sce_combined$Banksy_cluster <- factor(merge_map[clu])
table(original = clu, merged = sce_combined$Banksy_cluster)

df <- colData(sce_combined)[, c("Banksy_cluster", "Sample"), drop = FALSE]
(nclust_per_sample <- tapply(df[["Banksy_cluster"]], df[["Sample"]], function(x) length(unique(x))))

saveRDS(sce_combined, file = file.path(output_path, cluster_col, "combined_simulated_Banksy.rds"))

## manually match clusters across samples if needed