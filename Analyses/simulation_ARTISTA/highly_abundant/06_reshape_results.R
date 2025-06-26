###########################################################################################
### transform the data: aggregate resutls across domains and calculate q value ############
###########################################################################################
source("01_competitors_wrapper.R", echo=TRUE)
`%notin%` <- Negate(`%in%`)
generalize_transformation <- function(dir, file_prefix, gene_order=gene_order,
                                      pval_col = "p_val", 
                                      padj_col = "p_val_adj", 
                                      cluster_col = "cluster", 
                                      condition_filter = "Condition1") {
  # Create a temporary environment
  temp_env <- new.env(parent = globalenv())
  # Copy function arguments to temp_env
  assign("dir", dir, envir = temp_env)
  assign("file_prefix", file_prefix, envir = temp_env)
  assign("pval_col", pval_col, envir = temp_env)
  assign("padj_col", padj_col, envir = temp_env)
  assign("cluster_col", cluster_col, envir = temp_env)
  assign("condition_filter", condition_filter, envir = temp_env)
  
  # Run the function logic
  evalq({
    # Load the data
    load(paste0(dir, file_prefix, ".rda"))
    results_var <- ls()[sapply(ls(), function(x) is.list(get(x)) && length(get(x)) > 0)][1]
    results_res <- get(results_var)
    # Combine data across all clusters
    filtered_list <- lapply(names(results_res), function(name) {
      df <- results_res[[name]]
      gene_names <- if ("gene" %in% colnames(df)) df$gene else rownames(df)
      if (any(duplicated(gene_names)) & cluster_col %in% colnames(df)) {
        df <- df[df[[cluster_col]] == condition_filter, ]
        gene_names <- if ("gene" %in% colnames(df)) df$gene else rownames(df)
      }
      # Extract gene names from either row names or the 'gene' column
      df <- df[match(gene_order, gene_names), ]
      colnames(df) <- paste0(colnames(df), "_", name)
      return(df)
    })
    combined_df <- do.call(cbind, filtered_list)
    
    # Extract p-values and adjusted p-values
    p_val_columns <- grep(paste0("^", pval_col, "_\\d+"), colnames(combined_df), value = TRUE)
    p_val_adj_columns <- grep(paste0("^", padj_col, "_\\w+"), colnames(combined_df), value = TRUE)
    
    # Compute minimum p-values for each gene
    min_p_val <- apply(combined_df[, p_val_columns], 1, min)
    min_p_val_adj <- apply(combined_df[, p_val_adj_columns], 1, min)
    res <- cbind(combined_df, min_p_val, min_p_val_adj)
    if(file_prefix == "spatialLIBD_anova_Banksy_all5") 
      rownames(res) <- res$gene_0
    # Prepare long-format dataframe for perGeneQValue
    filtered_list2 <- lapply(names(results_res), function(name) {
      df <- results_res[[name]]
      gene_names <- if ("gene" %in% colnames(df)) df$gene else rownames(df)
      if (any(duplicated(gene_names)) & cluster_col %in% colnames(df)) {
        df <- df[df[[cluster_col]] == condition_filter, ]
        gene_names <- if ("gene" %in% colnames(df)) df$gene else rownames(df)
      }else{
        if(file_prefix == "spatialLIBD_anova_Banksy_all5") rownames(df) <- df$gene
        df$gene <- rownames(df)
      }
      df$featureID <- name
      df <- df[match(gene_order, gene_names), ]
      return(df)
    })
    if(file_prefix == "DESpace_individual_test_all5"){
      filtered_list2 <- lapply(filtered_list2, function(df) {
        df[, !grepl("logFC\\.", colnames(df))]
      })
    }
    res_long <- do.call(rbind, filtered_list2)
    
    # Apply perGeneQValue_modified
    res_long <- perGeneQValue_modified(res_long, pval = pval_col, padj = padj_col)
    
    # Merge long format results
    res <- merge(res, res_long, by = "row.names")
    rownames(res) <- res$Row.names
    res <- res[, -1]
    print(head(res))
    # Save the final transformed data
    save(res, file = paste0(dir, file_prefix, "_transformed.rda"))
  }, envir = temp_env)

  #print(ls(temp_env))
}

path <- "./Simulated_data/ARTISTA/highly_abundant/"

load("./Simulated_data/ARTISTA/highly_abundant/DESpace_global_GT.rda")
gene_order <- res_edgeR$gene_id; rm(res_edgeR)
#cluster <- "GT"
cluster <- "Banksy"
generalize_transformation(
  dir = path, 
  file_prefix = paste0("seurat_FindMarkers_", cluster)
)

generalize_transformation(
  dir = path, 
  file_prefix = paste0("scran_findMarkers_", cluster),
  pval_col = "p.value", 
  padj_col = "FDR"
)

generalize_transformation(
  dir = path, 
  file_prefix = paste0("seurat_PseudoBulk_FindMarkers_", cluster)
)

generalize_transformation(
  dir = path, 
  file_prefix = paste0("spatialLIBD_pairwise_", cluster),
  pval_col = "p_value_Condition1-Condition2",
  padj_col = "fdr_Condition1-Condition2"
)
generalize_transformation(
  dir = path, 
  file_prefix = paste0("DESpace_individual_", cluster),
  pval_col = "PValue", 
  padj_col = "FDR"
)

######################################### Load and transform results ##########################
###########################################################################################
rm(list = ls())
source("01_competitors_wrapper.R", echo=TRUE)
dir <- "./Simulated_data/ARTISTA/highly_abundant/"

process_method_results <- function(method_name_file, method_name_script, dir, 
                                   dir_output, data_name, 
                                   pval_col_name, padj_col_name,
                                   logFC_col_name) {
  cat("Processing method:", method_name_script, "\n")
  load(paste0(dir, method_name_file, ".rda"))
  results_res <- get(data_name)
  gene_order <- rownames(results_res[[1]])
  
  # Function to filter and format data
  filter_and_format <- function(df, gene_order, method_name_script, method_name_file, domain_name, pval_col_name, padj_col_name) {
    filtered_df <- data.frame(
      pval = df[, pval_col_name],  
      padj = df[, padj_col_name], 
      logFC = df[, logFC_col_name],
      gene = rownames(df),
      domain = domain_name,
      method = method_name_script
    )
    return(filtered_df)
  }
  
  if(is.null(names(results_res))) names(results_res) <- seq_len(length(results_res))
  # Prepare filtered data frame list for each domain
  filtered_list <- lapply(names(results_res), function(domain_name) {
    filtered_df <- filter_and_format(results_res[[domain_name]], gene_order, method_name_script, method_name_file, domain_name, pval_col_name, padj_col_name)
    return(filtered_df)
  })
  combined_df <- do.call(rbind, filtered_list)
  save(combined_df, file = paste0(dir_output, method_name_file, "_domains.rda"))
  cat("Processed method:", method_name_script, "\n")
}

# Mapping between method names in file paths and method names in R script
method_mapping <- list(
  "DESpace_individual_Banksy" = "DESpace",
  "seurat_FindMarkers_Banksy" = "FindMarkers",
  "scran_findMarkers_Banksy" = "findMarkers",
  "seurat_PseudoBulk_FindMarkers_Banksy" = "pseudo.bulk_FindMarkers",
  "spatialLIBD_pairwise_Banksy" = "spatialLIBD_pairwise"
  
)

# Define column names for pval and padj
pval_col_names <- c("PValue",
                    "p_val", 
                    "p.value", 
                    "p_val", 
                    "p_value_Condition1-Condition2")
padj_col_names <- c("FDR","p_val_adj",
                    "FDR", "p_val_adj", 
                    "fdr_Condition1-Condition2")
logFC_col_names <- c("logFC","avg_log2FC",
                     "summary.logFC","avg_log2FC", 
                     "logFC_Condition1-Condition2")

data_names <- c("results_res1","results_res1", 
                rep("results_res", 3))
# Process each method using the defined function and method_mapping
for (i in seq_along(method_mapping)) {
  method_file <- names(method_mapping)[i]
  method_script <- method_mapping[[method_file]]
  print(method_file)
  dir_sub = dir
  process_method_results(method_name_file = method_file, 
                         method_name_script = method_script, 
                         dir = dir_sub, dir_output = dir_sub, data_name = data_names[i], 
                         pval_col_name = pval_col_names[i], 
                         padj_col_name = padj_col_names[i],
                         logFC_col_name = logFC_col_names[i])
}

