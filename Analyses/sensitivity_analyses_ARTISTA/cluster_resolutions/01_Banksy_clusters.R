rm(list = ls())
load("./muSpaData_Data/Wei22_full.rda")
output_path <- "./Simulated_data/cluster_res/"
library(data.table)
library(SpatialExperiment)
library(Seurat)
library(Banksy)
library(SeuratWrappers)
library(harmony)

spe <- SpatialExperiment(
  rowData = NULL,
  colData = colData(spe),
  assays = assays(spe),
  reducedDims = NULL,
  sample_id = as.character(spe$sample),
  spatialCoordsNames = c("imagecol", "imagerow"),
  imgData = NULL
)
reducedDims(spe) <- NULL
colData(spe) <- DataFrame(
  sample_id = spe$sample,
  group = factor(spe$time),
  row.names = colnames(spe)
)

group1 <- c("10DPI_2", "20DPI_1",  "2DPI_2",  "5DPI_1")
group2 <- c("10DPI_1", "20DPI_2",  "2DPI_1",  "5DPI_2" )

spe <- spe[, spe$sample_id %in% c(group1, group2)]
spe$sample_id <- droplevels(spe$sample_id)
spe$subject_id <- ifelse(spe$sample_id %in% group1, "Condition_1", "Condition_2")
# Stagger spatial coordinates
# stagger the spatial coordinates across the samples so that spots from different samples do not overlap.
locs <- spatialCoords(spe)
locs <- cbind(locs, sample_id = as.numeric(factor(spe$sample_id)))
locs_dt <- data.table(locs)
colnames(locs_dt) <- c("sdimx", "sdimy", "group")
locs_dt[, sdimx := sdimx - min(sdimx), by = group]
global_max <- max(locs_dt$sdimx) * 1.5
locs_dt[, sdimx := sdimx + group * global_max]
locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(spe)
spatialCoords(spe) <- locs

# Get HVGs
sample_ids <- unique(spe$sample_id)
spe_list <- lapply(sample_ids, function(x) spe[, spe$sample_id == x])

# Normalize data
seu_list <- lapply(spe_list, function(x) {
  x <- as.Seurat(x, data = NULL)
  NormalizeData(x, scale.factor = 5000, normalization.method = 'RC')
})

# Compute HVGs
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

# Banksy cluster
compute_agf <- FALSE
k_geom <- 18
spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, 
                   compute_agf = compute_agf, k_geom = k_geom)

spe_joint <- do.call(cbind, spe_list)
# colnames(spe_joint) <- paste0(colnames(spe_joint), ".", 
#                               spe_joint$sample_id)
seu_joint <- as.Seurat(spe_joint, data = NULL)

group_name = c("sample_id")
dimx = "sdimx"; dimy = "sdimy"
col_data <- cbind(colData(spe_joint), spatialCoords(spe_joint))
seu_joint <- AddMetaData(seu_joint, metadata = as.data.frame(col_data))

# Fit Banksy
seu_joint$sample_id <- droplevels(factor(seu_joint$sample_id))
seu_joint = RunBanksy(
  seu_joint, lambda = 0.8, assay = 'originalexp', slot = 'data',
  group = group_name, dimx = dimx, dimy = dimy,features = "all",
  k_geom = 18, split.scale = TRUE)
seu_joint = RunPCA(seu_joint, features = rownames(spe_joint), 
                   npcs = 20)
seu_joint = RunHarmony(seu_joint, group.by.vars=group_name,
                       reduction.name='pca', reduction.save='harmony')
seu_joint = FindNeighbors(seu_joint, dims = 1:20, reduction = 'harmony')
seu_joint = FindClusters(seu_joint, resolution = 0.08, graph.name = 'BANKSY_snn')


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

cluster_list <- list()
resolutions <- numeric()

# Loop over target cluster numbers
for(n_clust_target in c(2, 4, 6, 8, 10, 12)) {
  
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
    return(NULL)  
  })
  
  if (is.null(result)) next 
  
  # Extract last column from metadata
  last_col <- tail(colnames(result@meta.data), 1)
  cluster_vec <- result@meta.data[[last_col]]
  
  # Extract the corresponding resolution 
  resolution <- sub(".*res\\.", "", last_col)
  resolutions <- c(resolutions, resolution)
  cluster_list[[paste0("nclust_", n_clust_target)]] <- cluster_vec
}

cluster_df <- as.data.frame(cluster_list)
cluster_df$cell <- rownames(result@meta.data)
cluster_df <- cluster_df[, c("cell", paste0("nclust_", c(2,4,6,8,10,12))[paste0("nclust_", c(2,4,6,8,10,12)) %in% names(cluster_list)])]

seu_joint@meta.data$barcode <- rownames(seu_joint@meta.data)

merged_meta <- merge(
  seu_joint@meta.data,
  cluster_df,
  by.x = "barcode",  
  by.y = "cell",     
  all.x = TRUE       
)

rownames(merged_meta) <- merged_meta$barcode
merged_meta$barcode <- NULL
seu_joint@meta.data <- merged_meta

cluster_cols <- grep("^nclust_", colnames(seu_joint@meta.data), value = TRUE)

cluster_df <- seu_joint@meta.data
cluster_df$cell <- rownames(cluster_df)
cluster_df_ordered <- cluster_df[colnames(spe_joint), ]
cluster_df_ordered$cell <- NULL
colData(spe_joint) <- DataFrame(cluster_df_ordered)
rotate_coord <- function(spe) {
  for (sample in unique(spe$sample_id)) {
    idx <- spe$sample_id == sample
    sdimx <- spe[["sdimx"]][idx]
    sdimy <- spe[["sdimy"]][idx]
    
    if (sample %in% c('2DPI_1', '5DPI_3')) {
      # Rotate 180 degrees: Reverse both x and y
      sdimx <- -sdimx
      sdimy <- -sdimy
    } else if (sample %in% c('5DPI_1', '5DPI_2', '10DPI_2', 
                             '15DPI_1', '15DPI_2', '15DPI_3',
                             '15DPI_4', '20DPI_1', '20DPI_2')) {
      # Rotate 90 degrees counterclockwise: Swap x and y
      sdimx <- -sdimy
      sdimy <- spe[["sdimx"]][idx]
    } else if (sample %in% c('10DPI_3')) {
      # Rotate 90 degrees clockwise: Swap x and y, reverse new y
      sdimx <- sdimy
      sdimy <- -spe[["sdimx"]][idx]
    }
    
    # Scale to ensure all coordinates are positive **within this sample**
    spe$sdimx[idx] <- sdimx - min(sdimx) + 1
    spe$sdimy[idx] <- sdimy - min(sdimy) + 1
  }
  return(spe)
}

spe <- rotate_coord(spe_joint)
cluster_cols <- paste0("nclust_", c(2,4,6,8,10,12))
colData(spe)[, cluster_cols] <- lapply(colData(spe)[, cluster_cols], as.factor)
saveRDS(spe, file = paste0(output_path, "spe_nclust_rotated.rds"))