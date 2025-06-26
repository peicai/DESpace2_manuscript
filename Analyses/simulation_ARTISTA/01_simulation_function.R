suppressMessages({library(SingleCellExperiment)
  library(BiocParallel)
  library(scater)
  library(limma)
  library(scuttle)
  library(edgeR)
  #library(DEXSeq)
  library(DESeq2)
  library(devtools)  # if not installed: install.packages('devtools')
  library(remotes)
  library(data.table)
  library(rlang)
  #library(qvalue)
  library(Seurat)  
  library(dplyr)
  library(MAST)
  library(scran)
})
`%notin%` <- Negate(`%in%`)

addPerCellQCMetrics <- function(x, ...) {
  colData(x) <- cbind(colData(x), scuttle::perCellQCMetrics(x, ...))
  x
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# Define a function to extract raw counts and cell metadata from an SCE object
process_subset <- function(subset_sce, cluster_col = "clust_Harmony_BANKSY_k50_res0.15",
                           coord_col = c("imagerow", "imagecol")) {
  # Get raw counts from subset
  raw_counts <- counts(subset_sce)
  # Get log norm counts
  if ("logcounts" %notin% names(assays(subset_sce))){
    # Normalize counts
    subset_sce <- scuttle::logNormCounts(subset_sce)
  }
  lognorm_counts <- logcounts(subset_sce)
  # Extract metadata and convert to data.table
  cell_meta <- as.data.table(
    colData(subset_sce)[, c(coord_col, "barcode", cluster_col,
                            "Scenario", "Condition", "Sample")]
  )
  
  # Return a list with counts and metadata
  list(raw_counts = raw_counts, lognorm_counts = lognorm_counts, cell_meta = cell_meta)
}

# Determine regions based on the scenario
determine_regions <- function(scenario, region_all) {
  if (scenario == "DSP1") {
    # add mixture: Unif vs. SV with randomly chosen group
    return(c(region_all[1], sample(region_all, 1)))
  } else if (scenario == "DSP2") {
    # add mixture: SV vs. SV, clusters in both groups randomly chosen, but different
    return(sample(region_all, 2))
  } else if (scenario == "NULL1") {
    # no mixture
    return(c(region_all[1], region_all[1]))
  } else if (scenario == "NULL2") {
    # add mixture: both groups SV in the same cluster, chosen at random among the n clusters
    return(rep(sample(region_all, 1), 2))
  }
}

# Function to process each gene
process_gene <- function(gene, gene_ind, 
                         gene_list, scenario, 
                         spatial_probs, para, pattern_name,
                         processed_sce_list, region_gene,
                         pattern_shape = "regular_mix",
                         cluster_col = "clust_Harmony_BANKSY_k50_res0.15"){
  print(paste0(round(gene_ind/length(gene_list)*100),"%"))
  # for every SV gene, randomly draw p from a beta
  para2 <- estBetaParams(mu = 0.5, var = 0.025^2)
  # for every SV gene, randomly draw p from a beta
  spatial_probs_updated <- sapply(spatial_probs, function(x) {
    ifelse(x == 0.5, rbeta(1, shape1 = para2$alpha, shape2 = para2$beta), {
      value <- rbeta(1, shape1 = para$alpha, shape2 = para$beta)
      while (value <= 0.5) {
        value <- rbeta(1, shape1 = para$alpha, shape2 = para$beta)
      }
      value
    })
  })

  # Prepare result for each group
  results <- list()
  
  for(grp in seq_along(region_gene)){
    # get the cell metadata and define the pattern
    cell_meta <- processed_sce_list[[grp]]$cell_meta
    pattern <- cell_meta[cell_meta[[cluster_col]] == region_gene[grp],]
    pattern_ids <- pattern$barcode
    if(pattern_shape == "regular_mix"){
      cell_meta[, (pattern_name) := ifelse(barcode %in% pattern_ids, 'in', 'out')]
    }else if(pattern_shape == "inverted_mix"){
      cell_meta[, (pattern_name) := ifelse(barcode %in% pattern_ids, 'out', 'in')]
    }
    
    colnames(cell_meta) <- c("x", "y", "cell_ID", "label", "Scenario", 
                             "Condition", "Sample", "pattern")
    
    # get gene counts
    gene_raw_counts <- processed_sce_list[[grp]]$raw_counts[gene, ]
    gene_lognorm_counts <- processed_sce_list[[grp]]$lognorm_counts[gene, ]
    
    # simulate the pattern
    gg <- simulateOneGenePatternObject(cell_meta = cell_meta,
                                       gene_raw_counts = gene_raw_counts,
                                       gene_lognorm_counts = gene_lognorm_counts,
                                       spatial_prob = spatial_probs_updated[grp])
    names(gg) <- cell_meta$cell_ID
    
    # define metadata
    df <- DataFrame(gene_name = gene,
                    in_pattern = region_gene[grp],
                    spatial_prob = spatial_probs_updated[grp],
                    group = grp,
                    scenario = scenario
    )
    results[[grp]] <- list(gg, df)
  }
  return(results)
}

# Function to generate SingleCellExperiment objects given the index of 'results_all'
create_sce_object <- function(index, gene_list, results_all, processed_sce_list) {
  # Extract the required data from results_all
  extracted_data <- lapply(results_all, function(x) x[[index]])
  
  # Prepare raw counts and row data
  new_raw_counts <- do.call(rbind, lapply(extracted_data, `[[`, 1))
  new_rowData <- do.call(rbind, lapply(extracted_data, `[[`, 2))
  
  # Set row names
  rownames(new_raw_counts) <- gene_list
  rownames(new_rowData) <- gene_list
  
  # Create SingleCellExperiment object with normalized counts
  sce_object <- SingleCellExperiment(
    assays = list(counts = new_raw_counts),
    colData = processed_sce_list[[index]]$cell_meta,
    rowData = new_rowData
  )
  print(unique(sce_object$Sample))
  print(rowData(sce_object)[1,])
  # Normalize counts
  sce_object <- scuttle::logNormCounts(sce_object)
  
  return(sce_object)
}

## Simulate single-gene spatial patterns ####
simulateOneGenePatternObject = function(cell_meta,
                                        #coordinate_name = c("row", "col"),
                                        pattern_name = 'pattern',
                                        gene_raw_counts,
                                        gene_lognorm_counts,
                                        spatial_prob = 0.95,
                                        # show_pattern = TRUE,
                                        # pattern_colors = c('in' = 'green', 'out' = 'red'),
                                        # save_raw = TRUE,
                                        # save_dir = '~',
                                        # save_name = 'counts',
                                        ...) {
  # data.table variables
  cell_ID = sdimx_y = sdimx = sdimy = NULL
  ## get number of cells within pattern
  cell_number = nrow(cell_meta[get(pattern_name) == 'in'])
  
  ## create the spatial expression pattern for the specified gene
  # 1. rank all gene values from the cells from high to low
  # 2. move the highest expressing values to the spatial pattern using a probability
  #     - 0.5 is the control = random
  #     - 1 is perfection: all the highest values go to the pattern
  #     - 0.5 to 1 is decreasing noise levels
  
  # rank genes
  gene_vector = gene_lognorm_counts
  sort_expr_gene = sort(gene_vector, decreasing = T)
  
  # number of cells in and out the pattern
  total_cell_number = length(sort_expr_gene)
  remaining_cell_number = total_cell_number - cell_number
  # print(paste0("total_cell_number: ", total_cell_number))
  # print(paste0("cell_number: ", cell_number))
  # print(paste0("remaining cell number: ", remaining_cell_number))
  # calculate outside probability
  outside_prob = 1 - spatial_prob
  prob_vector = c(rep(spatial_prob, cell_number), rep(outside_prob, remaining_cell_number))
  
  # first get the 'in' pattern sample values randomly
  sample_values = sample(sort_expr_gene, replace = F, size = cell_number, prob = prob_vector)
  
  # then take the remaining 'out' pattern values randomly
  remain_values = sort_expr_gene[!names(sort_expr_gene) %in% names(sample_values)]
  remain_values = sample(remain_values, size = length(remain_values))
  
  
  ## A. within pattern ##
  # ------------------- #
  in_cell_meta = cell_meta[get(pattern_name) == 'in']
  in_ids = in_cell_meta$cell_ID
  
  # preparation for raw matrix
  sample_values_id_vector = names(sample_values)
  names(sample_values_id_vector) = in_ids
  
  ## B. outside pattern ##
  # -------------------- #
  out_ids = cell_meta[get(pattern_name) == 'out']$cell_ID
  
  # preparation for raw matrix
  remain_values_id_vector = names(remain_values)
  names(remain_values_id_vector) = out_ids
  
  
  ## raw matrix
  # swap the cell ids #
  raw_gene_vector = gene_raw_counts
  
  raw_new_sample_vector = raw_gene_vector[sample_values_id_vector]
  names(raw_new_sample_vector) = names(sample_values_id_vector)
  
  raw_new_remain_vector = raw_gene_vector[remain_values_id_vector]
  names(raw_new_remain_vector) = names(remain_values_id_vector)
  
  new_sim_raw_values = c(raw_new_sample_vector, raw_new_remain_vector)
  new_sim_raw_values = new_sim_raw_values[names(raw_gene_vector)] # new row counts for the given gene

  return(new_sim_raw_values)
  
}


CombineSimulatedObject = function(processed_sce_list,
                                  #coordinate_name = c("row", "col"),
                                  pattern_name = 'pattern',
                                  spatial_prob = 0.95,
                                  scenario = "DSP1",
                                  regions = "WM",
                                  mu = 0.90,
                                  var = 0.025^2,  
                                  pattern_shape = "regular_mix",
                                  cluster_col = "clust_Harmony_BANKSY_k50_res0.15",
                                  ...) {
  set.seed(123)
  # we use all genes as 'SVGs' now; could be adapted later
  gene_list = rownames(processed_sce_list[[1]]$raw_counts)
  para <- estBetaParams(mu = mu, var = var)
  # Main loop over genes
  results_all <- lapply(seq_along(gene_list), function(gene_ind) {
    gene <- gene_list[gene_ind]
    spatial_probs <- spatial_prob[[gene_ind]]
    region_gene <- regions[[gene_ind]]
    process_gene(gene, gene_ind, gene_list, scenario, spatial_probs, para, 
                 pattern_name, processed_sce_list, region_gene, pattern_shape, cluster_col)
  })
  print(seq_along(spatial_prob[[1]]))
  sce_object_list <- lapply(seq_along(spatial_prob[[1]]), function(i) {
    create_sce_object(i, gene_list, results_all, processed_sce_list)
  }) 
  
  return( sce_object_list)
}


