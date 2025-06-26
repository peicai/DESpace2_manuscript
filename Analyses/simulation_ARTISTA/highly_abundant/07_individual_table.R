rm(list=ls())
library(SpatialExperiment)
library(dplyr)
library(purrr)
library(tidyr)
dataset <- "ARTISTA"
pattern <- "mixture"
cluster <- "Banksy"
path <- "./Simulated_data/ARTISTA/highly_abundant/"
sce.combined <- readRDS(paste0(path, "/combined_simulated.rda"))
gene_list_all <- rowData(sce.combined)
gene_order <- rownames(sce.combined)
null_gene_df <- gene_list_all %>%
  as.data.frame() %>%
  filter(scenario %in% c("NULL1", "NULL2"))
null_gene_df$status <- 0

# Step 1: Filter the data frame for scenarios "DSP1" or "DSP2"
filtered_df <- gene_list_all %>%
  as.data.frame() %>%
  filter(scenario %in% c("DSP1", "DSP2"))

df1 <- filtered_df %>%
  dplyr::select(gene_name, in_pattern_1, spatial_prob_1, scenario) 
colnames(df1) <- c("gene_name", "in_pattern", "spatial_prob", "scenario")
df1$group <- "1"
df2 <- filtered_df %>%
  dplyr::select(gene_name, in_pattern_2, spatial_prob_2, scenario)
colnames(df2) <- c("gene_name", "in_pattern", "spatial_prob", "scenario")
df2$group <- "2"
# Combine the two dataframes, handling row names by resetting them
final_df <- rbind(df1, df2)

# Step 2: For each level in "in_pattern", filter out rows where spatial_prob is not equal to 0.5
filtered_df <- final_df %>%
  group_by(in_pattern) %>%
  filter(round(as.numeric(spatial_prob),1) != 0.5) 

filtered_df$in_pattern <- paste0("Layer", filtered_df$in_pattern)

# Relabel the Truth column
# filtered_df$in_pattern <- recode(filtered_df$in_pattern, !!!relabel_map)
# Step 3: Split the data frame into a list of data frames based on in_pattern
result_list <- split(filtered_df, filtered_df$in_pattern)
result_list <- do.call(rbind, result_list)

files <- c(
  paste0("spatialLIBD_pairwise_", cluster, "_domains.rda"),
  paste0("seurat_PseudoBulk_FindMarkers_", cluster, "_domains.rda"),
  paste0("seurat_FindMarkers_", cluster, "_domains.rda"),
  paste0("scran_findMarkers_", cluster, "_domains.rda"),
  paste0("DESpace_individual_", cluster, "_domains.rda")
)

# Define the variable names for storing results
var_names <- c(
  "spatialLIBD",
  "PseudoBulk_FindMarkers",
  "seurat_FindMarkers",
  "scran_findMarkers",
  "DESpace"
)
`%notin%` <- Negate(`%in%`)
# Load files, reorder results, and assign to variables
null_gene_order <- gene_order[gene_order %in% null_gene_df$gene_name]
DSP_gene_order <- gene_order[gene_order %in% result_list$gene_name]
reorder_results <- function(gene_order = gene_order, results_df, 
                            only_pval = TRUE, only_padj = FALSE){
  if("pval" %notin% colnames(results_df)) colnames(results_df) <- c("pval", "padj", 
                                                                    "logFC", "gene", 
                                                                    "domain", "method")
  reorder_df <- lapply(levels(factor(results_df$domain)), function(dom) {
    results_df_sub <- results_df[results_df$domain == dom, ]
    results_df_sub$gene <- sub("\\..*", "", results_df_sub$gene)
    reordered_df_sub <- data.frame(
      do.call(rbind, lapply(gene_order, function(gene) {
        if (gene %in% results_df_sub$gene) {
          results_df_sub[results_df_sub$gene == gene, , drop = FALSE]
        } else {
          data.frame(pval = 1, padj = 1, logFC = 0, gene = gene,
                     domain = dom, method = unique(results_df$method)) 
        }
      }))
    )
  })
  
  # Combine into a single data frame and sort by domain index
  combined_df_sorted <- do.call(rbind, reorder_df)
  if(only_pval){
    top2 <- combined_df_sorted %>%
      group_by(gene) %>%
      # First, arrange by pval (ascending) and logFC (descending)
      arrange(pval, desc(abs(logFC))) %>%  
      # Select the two rows with the lowest pval
      slice_min(order_by = pval, n = 2, with_ties = FALSE)  

    top1 <- top2 %>%
      group_by(gene) %>%
      arrange(pval, desc(abs(logFC))) %>%  
      slice_min(order_by = pval, with_ties = FALSE)
  }
  if(only_padj){
    top2 <- combined_df_sorted %>%
      group_by(gene) %>%
      arrange(pval, desc(abs(logFC))) %>%  
      slice_min(order_by = padj, n=2, with_ties = FALSE)
    
    top1 <- top2 %>%
      group_by(gene) %>%
      arrange(pval, desc(abs(logFC))) %>%  
      slice_min(order_by = padj, with_ties = FALSE)
  }
  return(list(top1, top2))
}

for (i in seq_along(files)) {
  print(i)
  load(paste0(path, files[i]))
  assign(var_names[i], reorder_results(gene_order = DSP_gene_order, combined_df, only_pval = FALSE, only_padj = TRUE))
  rm(combined_df)
}

list_all <- list(
  spatialLIBD,
  PseudoBulk_FindMarkers,
  seurat_FindMarkers,
  scran_findMarkers,
  DESpace
)
extract_and_recode <- function(data_list, index) {
  map_dfr(data_list, ~ .x[[index]]) 
}

# Apply for top1 and top2
res_all_top2 <- extract_and_recode(list_all, 2)
res_all_top1 <- extract_and_recode(list_all, 1)

DSP1_gene <- result_list %>%
  filter(scenario == "DSP1") %>%
  filter(!(spatial_prob - 0.5 < 0.1))
res_all4_DSP1 <- res_all_top1 %>% filter(gene %in% DSP1_gene$gene_name)
# Merge the dataframes on the gene column
merged_df <- merge(res_all4_DSP1, DSP1_gene, by.x = "gene", by.y = "gene_name")

(match_percentages <- merged_df %>%
    mutate(match = in_pattern == paste0("Layer", domain)) %>%
    group_by(method) %>%
    summarise(
      total = n(),
      matches = sum(match),
      percentage = (matches / total) * 100
    ))

# Filtering result_list by scenario
filtered_df <- result_list %>% filter(scenario == "DSP2")
res_all4_DSP2 <- res_all_top2 %>% filter(gene %in% filtered_df$gene_name)
df_wide <- res_all4_DSP2 %>%
  group_by(gene, method) %>%
  mutate(row = row_number()) %>%  
  pivot_wider(names_from = row, values_from = domain, names_prefix = "domain.") %>%
  summarise(
    domain.1 = coalesce(first(domain.1), last(domain.1)),  # Combine non-NA values
    domain.2 = coalesce(first(domain.2), last(domain.2))   # Combine non-NA values
  ) %>%
  ungroup() 

# Split by group and rename columns
grouped_dfs <- filtered_df %>%
  group_by(group) %>%
  group_split()

grouped_dfs <- merge(grouped_dfs[[1]], grouped_dfs[[2]], by = "gene_name")
merged_df <- merge(df_wide, grouped_dfs, by.x = "gene", by.y = "gene_name")

(match_percentages <- merged_df %>%
    rowwise() %>%  
    mutate(
      match = sum(c(paste0("Layer", domain.1), paste0("Layer", domain.2))
                  %in% c(in_pattern.x, in_pattern.y)) == 2  # Check if both domain.1 and domain.2 are in in_pattern.x or in_pattern.y
    ) %>%
    group_by(method) %>%
    summarise(total = n(),
              matches = sum(match),
              percentage = (matches / total) * 100
    ))  # Count how many rows have both matches


