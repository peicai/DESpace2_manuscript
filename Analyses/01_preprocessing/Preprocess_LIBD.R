rm(list=ls())
library("spatialLIBD")
library(SpatialExperiment)
library(scater)
`%notin%` <- Negate(`%in%`)
setwd("./Data/LIBD/")
#load("human_sce.rda")
ehub <- ExperimentHub::ExperimentHub()
# Download the full real data (about 2.1 GB in RAM) use:
sce <- spatialLIBD::fetch_data(type = "sce", eh = ehub)
rm(ehub)

(spe <- SpatialExperiment(
  assays = list(counts = assay(sce)),
  colData = colData(sce)
))
# subset to keep only spots over tissue
spatialCoords(spe) <- as.matrix(colData(sce)[,6:7])
spatialCoordsNames(spe) <- c("x_coord","y_coord")

spatialData(spe) <- cbind(spatialData(spe), 
                          in_tissue = colData(sce)$sample_name)

# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe,)
sample_names <- unique(colData(sce)$sample_name)

qc_lib_size_all <- c(250, 250, 250, 250, 300, 300, 
                     300, 250, 350, 700, 500, 500)
qc_detected_all <- c(250, 250, 250, 250, 300, 300, 
                     300, 250, 350, 700, 350, 500)
qc_mito_all <- c(0.3, 0.3, 0.3, 0.3, 0.35, 0.35, 
                 0.35, 0.3, 0.3, 0.27, 0.34, 0.31)
qc_cells_per_spot_all <- c(10, 12, 12, 10, 8, 16, 
                           8, 12, 15, 14, 13, 14)

for(i in seq_along(sample_names)){
    #### sample size:
    sample_id <- sample_names[i]
    spe1 <- spe[, spatialData(spe)$in_tissue == sample_id]
    dim(spe1)
    
    ## Library size
    # select QC threshold for library size
    qc_lib_size <- colData(spe1)$sum < qc_lib_size_all[i]
    colData(spe1)$qc_lib_size <- qc_lib_size
    
    ## Number of expressed features
    # select QC threshold for expressed genes
    qc_detected <- colData(spe1)$detected < qc_detected_all[i]
    colData(spe1)$qc_detected <- qc_detected
    
    ## Proportion of mitochondrial read
    # select QC threshold for mitochondrial read proportion
    qc_mito <- colData(spe1)$expr_chrM_ratio > qc_mito_all[i]
    colData(spe1)$qc_mito <- qc_mito
    
    ## Number of cells per spot
    # select QC threshold for number of cells per spot
    qc_cell_count <- colData(spe1)$cell_count > qc_cells_per_spot_all[i]
    colData(spe1)$qc_cell_count <- qc_cell_count
    
    ## Remove low-quality spots
    # combined set of discarded spots
    discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count

    # store in object
    colData(spe1)$discard <- discard
    
    # remove combined set of low-quality spots
    spe1 <- spe1[, !colData(spe1)$discard]
    
    ## Gene level
    ### basic filter: remove undetected genes
    # calculate logcounts (log-transformed normalized counts) and store in object
    spe1 <- logNormCounts(spe1)
    qc_undetected_gene <- rowSums(assays(spe1)$logcounts > 0) > 0
    
    ### basic filter: remove lowly expressed genes: at least xx non-zero cells:
    qc_low_gene <- rowSums(assays(spe1)$counts > 0) >=20
    ## remove combined set of low-quality spots
    
    spe1 <- spe1[qc_low_gene & qc_undetected_gene,]
    dim(spe1)
    
    ## store the filtered dataset
    sce_one = spe1
    path = './Data/LIBD/'
    save(sce_one, file = paste0(path,sample_id,"_sce_LIBD.rda"))
}