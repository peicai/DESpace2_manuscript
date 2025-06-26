rm(list=ls())
input_path <- "./"
method_mapping <- list(
  "DESpace_global_test_2vs20" = "DESpace",
  "DESpace_global_test_all5" = "DESpace (multi-group)",
  "DESpace_individual_test_2vs20_transformed" = "DESpace (individual)",
  "DESpace_individual_test_all5_transformed" = "DESpace (multi-group with individual)",
  "DESpace_spline_all5" = "DESpace (multi-group with spline)",
  "seurat_FindMarkers_Banksy_2vs20_transformed" = "FindMarkers",
  "scran_findMarkers_Banksy_2vs20_transformed" = "findMarkers",
  "scran_findMarkers_Banksy_all5_transformed" = "findMarkers (multi-group)",
  "seurat_PseudoBulk_FindMarkers_Banksy_2vs20_transformed" = "pseudo.bulk_FindMarkers",
  "spatialLIBD_pairwise_Banksy_2vs20_transformed" = "spatialLIBD_pairwise",
  "spatialLIBD_anova_Banksy_all5_transformed" = "spatialLIBD_anova (multi-group)"
)

combine_results <- function(input_path, method_mapping, i) {
  print(names(method_mapping)[i])
  load(paste0(input_path, names(method_mapping)[i], ".rda"))
  if(i %in% c(1,2,5)){
    result <- res_edgeR$gene_results[, c("PValue", "FDR", "gene_id")]
  }else{
      if(class(res) != "data.frame") res <- res |> as.data.frame()
      result <- res[, c("min_p_val", "qval")]
      result$gene <- rownames(result)
  }
  colnames(result) <- c("p_value", "FDR", "gene")
  result$method <- method_mapping[[i]][1]
  return(result)
}

# Apply function over each comparison, then bind and reshape results
results_list <- lapply(seq_along(method_mapping), function(i) {
  combine_results(input_path, method_mapping, i)
})

# Combine all results into one data frame
library(dplyr)
combined_results <- bind_rows(results_list)

############################ gene lists from MSigDB ################################
gene_folder = "./data/ARTISTA/gene_list_MSigDB/"
BUZZ_WORDS = list.files(gene_folder)
# remove .json
BUZZ_WORDS = substring(BUZZ_WORDS, 1,  nchar(BUZZ_WORDS) - 5)
library(rjson)
GENES <- lapply(BUZZ_WORDS, function(WORDS){
  print(WORDS)
  a = unlist(fromJSON(file = paste0(gene_folder, WORDS,".json")))
  sel = grepl("geneSymbols|Gene", names(a))
  return(a[sel] |> as.vector())
})
str(GENES)
names(GENES) <- BUZZ_WORDS # healing and regeneration

gene_list <- read.csv("./data/ARTISTA/gene_name_list.csv", stringsAsFactors = FALSE)
library(stringr)
#install.packages("fuzzyjoin")
library(fuzzyjoin)
# Ensure both columns are lowercase for case-insensitive matching
gene_list <- gene_list %>%
  mutate(extracted_name = tolower(gene_id))
library(dtplyr)
gene_list$extracted_name <- as.character(gene_list$extracted_name)

gene_list_cleaned <- gene_list %>%
  mutate(
    # Extract the [nr] part (removing [nr] tag)
    nr_name = str_extract(extracted_name, "(^|\\|)[^\\|\\[]+\\[nr\\]") %>%
      str_remove_all("\\[nr\\]|\\|"),
    
    # Extract the [hs] part (removing [hs] tag)
    hs_name = str_extract(extracted_name, "(^|\\|)[^\\|\\[]+\\[hs\\]") %>%
      str_remove_all("\\[hs\\]|\\|"),
    
    # Clean up empty strings to NA
    nr_name = ifelse(nr_name == "", NA, nr_name),
    hs_name = ifelse(hs_name == "", NA, hs_name),
    raw_name = ifelse(
      !str_detect(extracted_name, "\\[nr\\]") & !str_detect(extracted_name, "\\[hs\\]"),
      extracted_name,
      NA
    )
  )

# Define row names for metrics
row_names <- c("n_genes", # number of genes associated in this GENE[[i]]
               "n_gene_ids_annotated", # number of gene ids matched
               "n_gene_names_annotated", # number of gene names matched
               "n_gene_ids_analyzed",    # number of gene ids analyzed 
               "n_gene_ids_data",       # number of gene ids included in n_total_genes
               "n_top_genes_in_list", # number of significant genes in this GENE[[i]]
               "n_sig_genes", # number of significant gene
               "n_top_genes", # top n genes
               "n_total_genes",       # number of genes for this sample comparison
               "n_total_analyzed_genes", # number of genes analyzed by DESpace
               "pct_top_genes_in_list", # number of sig genes in this GENE[[i]]/number of genes in GENE are analyzed by DESpace
               "pct_genes_analyzed", # number of sig genes/number of genes analyzed by DESpace
               "fraction_top_genes_in_list", # pct_sig_genes_in_list/pct_genes_analyzed
               "avg_rank_in_list1" ,           # average ranking in list (of all genes in the list that have been analyzed
               "avg_rank_in_list2",
               "avg_rank_in_list3" 
)
get_rank <- function(matched_name, gene_list) {
  matched_ranks <- gene_list %>%              
    mutate(FDR_rank = dense_rank(FDR)) %>%   
    filter(gene %in% matched_name) %>%
    dplyr::select(gene, method, FDR, FDR_rank)  # Select relevant columns
  round(mean(matched_ranks[["FDR_rank"]]))
}

################ to calculate average ranking -> issue: n_total_analyzed_genes are not comparable across methods #####
################# idea 1: only count overlapped genes ###################
# Count occurrences of each gene across methods
gene_counts <- combined_results %>%
  group_by(gene) %>%
  summarize(method_count = n_distinct(method))

# Find genes that are present in all methods
total_methods <- n_distinct(combined_results$method)
overlapped_genes <- gene_counts %>%
  filter(method_count == total_methods) %>%
  pull(gene)

# Result: vector of genes that appear in all methods
overlapped_genes
############ idea 2: FDR = 1 for those genes that are filtered out ########
# Create a data frame of all combinations of genes and methods
all_combinations <- expand.grid(
  gene = unique(combined_results$gene),
  method = unique(combined_results$method)
)

# Left join with combined_results to identify missing entries
filled_results <- all_combinations %>%
  left_join(combined_results, by = c("gene", "method")) %>%
  mutate(
    p_value = ifelse(is.na(p_value), 1, p_value),
    FDR = ifelse(is.na(FDR), 1, FDR)
  )

# filled_results now contains the original data with dummy rows for missing gene-method pairs
head(filled_results)
table(filled_results$method)
####################### create the table #############
# Function to create a table for each comparison and gene list
load("./Simulated_data/ARTISTA/Banksy_harmony.rda")
n=200
create_comparison_table <- function(com) {
  print(com)
  # Precompute values for this comparison
  n_total_analyzed_genes <- sum(combined_results$method == com)
  # numebr of genes for 5-group and 2-group comparison
  n_total_genes <- 13890  
  genes_analyzed <- combined_results[combined_results$method == com, ]
  sig_genes <- genes_analyzed[genes_analyzed$FDR <= 0.05, ]
  top_genes <- genes_analyzed %>%
    arrange(FDR) %>%     # Extract the gene IDs after sorting
    head(n = n)
  # Create the table for each gene list
  table_for_genes <- lapply(GENES, function(gene_list) {
    gene_names <- tolower(gene_list)
    n_genes <- length(gene_names)
    
    # Match genes in gene_list_cleaned with current gene list
    matched_genes <- gene_list_cleaned %>%
      mutate(across(c(hs_name, raw_name), tolower)) %>%
      filter(hs_name %in% gene_names | raw_name %in% gene_names) %>%
      distinct(gene_id, gene_name) 
    
    avg_rank_in_list1 <- get_rank(matched_genes$gene_name, genes_analyzed)
    avg_rank_in_list2 <- get_rank(matched_genes$gene_name, 
                                  genes_analyzed[genes_analyzed$gene %in% overlapped_genes,])
    avg_rank_in_list3 <- get_rank(matched_genes$gene_name, 
                                  filled_results[filled_results$method == com,])
    
    # Calculate metrics
    n_gene_ids_annotated <- length(unique(matched_genes$gene_name))
    n_gene_names_annotated <- length(unique(matched_genes$gene_id))
    n_gene_ids_data <- length(intersect(matched_genes$gene_name, rownames(spe)))
    n_gene_ids_analyzed <- length(intersect(matched_genes$gene_name, genes_analyzed$gene))
    n_sig_genes <- nrow(sig_genes)
    n_top_genes <- nrow(top_genes)
    n_top_genes_in_list <- length(intersect(matched_genes$gene_name, top_genes$gene))
    pct_top_genes_in_list <- n_top_genes_in_list/n_gene_ids_analyzed
    pct_genes_analyzed <- n_top_genes/n_total_analyzed_genes
    fraction_top_genes_in_list <- round(pct_top_genes_in_list/pct_genes_analyzed,0)
    print(intersect(matched_genes$gene_name, top_genes$gene))
    # Compile the metrics as a named vector
    c(n_genes, n_gene_ids_annotated, n_gene_names_annotated, 
      n_gene_ids_analyzed, n_gene_ids_data, n_top_genes_in_list, 
      n_sig_genes,n_top_genes, n_total_genes, n_total_analyzed_genes,
      round(pct_top_genes_in_list,3)*100,round(pct_genes_analyzed,3)*100,fraction_top_genes_in_list,
      avg_rank_in_list1, avg_rank_in_list2, avg_rank_in_list3)

  })
  
  # Convert to dataframe, set row names, and return
  comparison_table <- as.data.frame(do.call(cbind, table_for_genes))
  rownames(comparison_table) <- row_names
  colnames(comparison_table) <- names(GENES)
  
  comparison_table
}

# Apply the function to each unique comparison and store in a list
comparison_tables <- lapply(unique(combined_results$method), create_comparison_table)

# Set names of the list by comparison
names(comparison_tables) <- unique(combined_results$method)

# Display the list of tables for each comparison
comparison_tables
#install.packages("writexl")
library(writexl)
comparison_tables_with_row_names <- lapply(comparison_tables, function(df) {
  df <- cbind(Row_Names = rownames(df), df)  # Add row names as a new column
  rownames(df) <- NULL  # Remove row names to avoid duplication
  return(df)
})
names(comparison_tables_with_row_names) <- c(
  "DESpace", 
  "DESpace_multi_group",
  "DESpace_individual",
  "DESpace_individual_multi_group",
  "DESpace_spline_multi_group", 
  "seurat_FindMarkers", 
  "scran_findMarkers", 
  "findMarkers_multi_group", 
  "pseudo_bulk_FindMarkers", 
  "spatialLIBD_pairwise", 
  "spatialLIBD_anova_multi_group"
)
# Save the list of tables to an Excel file with each comparison as a separate sheet
write_xlsx(comparison_tables_with_row_names, path = paste0("./",
                                                           "comparison_tables_top", n, "_MSigDB_v2.xlsx"))

summary_list <- list()
# 
# # Loop through each method in the list
for (method_name in names(comparison_tables_with_row_names)) {

  # Get the table
  tbl <- comparison_tables_with_row_names[[method_name]]

  # Ensure it's a data frame and extract the rows of interest
  top_genes_row <- tbl[tbl$Row_Names == "n_top_genes_in_list", ]
  avg_rank_row  <- tbl[tbl$Row_Names == "avg_rank_in_list2", ]
  n_sig_genes  <- tbl[tbl$Row_Names == "n_sig_genes", ]
  # Extract values from both columns (assuming column order is known or fixed)
  regen_top_genes <- top_genes_row[[2]]
  wound_top_genes <- top_genes_row[[3]]

  regen_avg_rank  <- avg_rank_row[[2]]
  wound_avg_rank  <- avg_rank_row[[3]]

  # Combine into a data frame row
  summary_list[[method_name]] <- data.frame(
    Method = method_name,
    Regeneration_TopGenes = regen_top_genes,
    Wounding_TopGenes = wound_top_genes,
    Total = regen_top_genes + wound_top_genes,
    Regeneration_AvgRank = regen_avg_rank,
    Wounding_AvgRank = wound_avg_rank,
    num_sig_genes =  as.numeric(n_sig_genes[2])
  )
}
summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL 

summary_df <- summary_df[order(-summary_df$Total), ]
print(summary_df)

write.csv(summary_df, file = paste0("./",
                                    "table_top", n, "_MSigDB.csv"),
          row.names = FALSE)
