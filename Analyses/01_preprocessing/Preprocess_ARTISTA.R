################################################################################
##                     1.  Load and build datasets
################################################################################
# Manually download the data:
# Download all RDS files containing "DPI" from the following link: https://db.cngb.org/stomics/artista/download/
# Alternatively (this method is much slower):
# Define the DPI groups and the number of files in each group
# dpi_groups <- list(
#   "2DPI" = 3, "5DPI" = 3, "10DPI" = 3, "15DPI" = 4, "20DPI" = 3
# )
# # Generate the file names
# file_names <- unlist(lapply(names(dpi_groups), function(dpi) {
#   paste0(dpi, "_", seq_len(dpi_groups[[dpi]]))
# }))
# urls <- paste0("https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000056/RDS/", file_names, ".rds")
# download_dir <- "./Data"
# dir.create(download_dir, showWarnings = FALSE) 
# # Download and load each .rds file
# for (url in urls) {
#   # Extract the file name from the URL
#   file_name <- basename(url)
#   # Define the local file path
#   local_file <- file.path(download_dir, file_name)
#   download.file(url, local_file, mode = "wb", method = "curl") 
#   output_file <- file.path(download_dir, basename(url))
# }
library(Seurat)
library(SingleCellExperiment)
library(scuttle)
library(scater)
# Put all downloaded RDS data into `folder_path`
folder_path <- "./Data"
files <- list.files(path = folder_path, pattern = "^[0-9].*[0-9]\\.rds$", full.names = TRUE)

sce_list <- list()
i <- 1
for(file in files){
  file_name <- tools::file_path_sans_ext(basename(file))
  seurat <- readRDS(file)
  sce <- as.SingleCellExperiment(seurat)
  matches <- regmatches(file_name, gregexpr("[0-9]+", file_name))
  
  if(identical(rownames(colData(sce)), rownames(seurat@images[[1]]@coordinates))){
    colData(sce) <- cbind(colData(sce), seurat@images[[1]]@coordinates)
  }else{
    print("match cell names!")
  }
  
  sce$time <- paste0(matches[[1]][1], "DPI")
  sce$rep <- paste0("rep", matches[[1]][2])
  sce$sample <- file_name
  sce_new <- SingleCellExperiment(
    assays = list(counts = counts(sce),
                  logcounts = logcounts(sce)),
    colData = colData(sce)
  )
  sce_list[[i]] <- sce_new
  i = i+1
}

common_gene <- Reduce(intersect, lapply(sce_list, rownames))
sce_all <- do.call(cbind, lapply(sce_list, function(sce) sce[common_gene,]))
save(sce_all, file = paste0("./Data/injury_DPI_all.rda"))

################################################################################
##                      2. Quality control
################################################################################
# Load the generated singleCellExperiment object
load("./Data/injury_DPI_all.rda")

library(SpatialExperiment)
sample_ids <- sce_all$sample |> unique()

# Discard low-quality cells
spe <- do.call(cbind, lapply(sample_ids, function(sample){
  # Subset the spe object for the given sample
  spe_subset <- sce_all[, sce_all$sample == sample]
  thres_sum <- quantile(colData(spe_subset)$sum, probs = 0.01)
  thres_detected <- quantile(colData(spe_subset)$detected, probs = 0.01)
  # Filter cells based on library size and number of expressed
  filter_condition <- with(colData(spe_subset), 
                           sum > thres_sum & 
                             detected > thres_detected)
  spe_subset <- spe_subset[, filter_condition]  
  return(spe_subset)
}))

# Discard lowly abundant genes, which were detected in less than 20 spots.
qc_low_genes <- lapply(sample_ids, function(sample){
  # Subset the spe object for the given sample
  spe_subset <- spe[, spe$sample == sample]
  # Select QC threshold for lowly expressed genes: at least 20 non-zero cells:
  qc_low_gene <- rowSums(assays(spe_subset)$counts > 0) >= 20
  return(rownames(spe_subset)[qc_low_gene])
})

# Find common genes across all samples
common_genes <- Reduce(intersect, qc_low_genes)
spe <- spe[common_genes, ]
save(spe, file = paste0("./Data/clean_data.rda"))

################################################################################
##                      3. Banksy spatial cluster
################################################################################

load("./Data/clean_data.rda")
library(data.table)
library(SpatialExperiment)
library(Seurat)
library(Banksy)
library(SeuratWrappers)
library(harmony)
# Pre-process data
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
seu_joint <- as.Seurat(spe_joint, data = NULL)

group_name = c("sample_id")
dimx = "sdimx"; dimy = "sdimy"
col_data <- cbind(colData(spe_joint), spatialCoords(spe_joint))
seu_joint <- AddMetaData(seu_joint, metadata = as.data.frame(col_data))

# Fit Banksy
seu_joint = RunBanksy(
  seu_joint, lambda = 0.8, assay = 'originalexp', slot = 'data',
  group = group_name, dimx = dimx, dimy = dimy,features = "all",
  k_geom = 18, split.scale = TRUE)
seu_joint = RunPCA(seu_joint, features = rownames(spe_joint), 
                   npcs = 20)
seu_joint = RunHarmony(seu_joint, group.by.vars=group_name,
                       reduction.name='pca', reduction.save='harmony')
seu_joint = FindNeighbors(seu_joint, dims = 1:20, reduction = 'harmony')
seu_joint = FindClusters(seu_joint, resolution = 0.1, graph.name = 'BANKSY_snn')

# Banksy smooth
spe_joint = spe
spe_joint$BANKSY_snn_res0.1 <- seu_joint@meta.data$BANKSY_snn_res.0.1
spe_list <- lapply(sample_ids, function(x) spe_joint[, spe_joint$sample_id == x])
spe_list <- lapply(spe_list, smoothLabels, cluster_names = "BANKSY_snn_res0.1", k = 6L, verbose = FALSE)
names(spe_list) <- paste0("sample_", sample_ids)
spe_Banksy <- do.call(cbind, spe_list)

save(spe_Banksy, file = paste0("./Data/data_Banksy.rda"))

################################################################################
##                      4. Gene name annotation
################################################################################
gene_ids <- read.csv("./Data/gene_name_list.csv")
# R scripts to generate 'gene_name_list.csv' file
# Obtain gene names and matched gene IDs from .h5ad files

# Download a `.h5ad` file, e.g., Regeneration.h5ad`, from the following link: 
# https://db.cngb.org/stomics/artista/download/
# library(zellkonverter)
# library(tidyr)
# setwd("./Data")
# ## Convert data from anndatato to sce
# sce <- readH5AD("Regeneration.h5ad")
# df <- colData(sce)
# df <- data.frame(row_name = rownames(sce))
# # Function to split row names into gene and identifier columns
# split_row_names <- function(row_name) {
#   if (grepl(" \\| ", row_name)) {
#     parts <- strsplit(row_name, " \\| ")[[1]]
#     gene <- parts[1]
#     identifier <- parts[2]
#   } else {
#     gene <- NA
#     identifier <- row_name
#   }
#   return(data.frame(gene_id = gene, gene_name = identifier, stringsAsFactors = FALSE))
# }
# 
# # Apply function to each row name
# split_df <- df %>%
#   rowwise() %>%
#   mutate(split = list(split_row_names(row_name))) %>%
#   unnest(cols = c(split))
# write.csv(split_df, "gene_name_list.csv", row.names = FALSE)

################################################################################
##                      5. Build final SPE objects
################################################################################
load("./Data/data_Banksy.rda")
library(dplyr)
# Add matched gene ids to rowData
gene_name <- rownames(spe_Banksy)
gene_id <- lapply(gene_name, function(gene) {
  gene_id <- gene_ids[which(gene_ids$gene_name == gene), "gene_id"]
  if(is.na(gene_id)){gene_id <- gene}
  gene_id
}) %>% unlist()
rowData(spe_Banksy) <- data.frame(gene_name, gene_id)
# Avoid duplicate cell names
colnames(spe_Banksy) <- paste0(colnames(spe_Banksy), ".", 
                               spe_Banksy$sample_id)
################################################################################
##                      5.1 Wei2022_full
################################################################################
colData(spe_Banksy) <- cbind(colData(spe_Banksy), 
                             spatialCoords(spe_Banksy))
colnames(colData(spe_Banksy)) <- c("sample_id", "condition",
                                   "Banksy", "Banksy_smooth", 
                                   "sdimx", "sdimy")
colData(spe_Banksy) <- subset(colData(spe_Banksy), select = -Banksy)
spe <- spe_Banksy; rm(spe_Banksy)

sample_ids <- unique(spe$sample_id)
# rotate tissue coordinates to ensure all sections look consistent.
rotate_coord <- function(spe, sample) {
  spe_object <- spe[, spe$sample_id == sample]
  sdimx <- spe_object[["sdimx"]]
  sdimy <- spe_object[["sdimy"]]
  if (sample %in% c('2DPI_1', '5DPI_3')) {
    # Rotate 180 degrees: Reverse both x and y
    sdimx <- -sdimx
    sdimy <- -sdimy
  } else if (sample %in% c('5DPI_1', '5DPI_2', '10DPI_2', 
                              '15DPI_1', '15DPI_2', '15DPI_3',
                              '15DPI_4', '20DPI_1', '20DPI_2')) {
    # Rotate 90 degrees counterclockwise: Swap x and y
    sdimx <- -sdimy
    sdimy <- spe_object[["sdimx"]]
  } else if (sample %in% c('10DPI_3')) {
    # Rotate 90 degrees clockwise: Swap x and y, reverse new y
    sdimx <- sdimy
    sdimy <- -spe_object[["sdimx"]]
  }
  
  # Scale to ensure all coordinates are positive
  spe_object$sdimx <- sdimx - min(sdimx) + 1
  spe_object$sdimy <- sdimy - min(sdimy) + 1
  return(spe_object)
}
all_spe <- lapply(sample_ids, function(x) rotate_coord(spe, x))
spe <- do.call(cbind, all_spe)
save(spe, file = paste0("./Data/Wei22_full.rda"))
