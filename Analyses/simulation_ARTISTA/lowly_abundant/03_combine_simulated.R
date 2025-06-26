rm(list = ls())
library(SingleCellExperiment)
path <- "./Simulated_data/ARTISTA/lowly_abundant/"
# Load the simulated datasets and store them in a list
file_paths <- c(
  paste0(path,"DSP1_inverted_simulated.rda"),
  paste0(path,"DSP2_inverted_simulated.rda"),
  paste0(path,"NULL1_inverted_simulated.rda"),
  paste0(path,"NULL2_inverted_simulated.rda")
)

# Create an empty list to store the loaded data
simulated_sce_list <- vector("list", length(file_paths))
simulated_list <- list()
# Load and assign the data
for (i in seq_along(file_paths)) {
  load(file_paths[i])
  simulated_list[[i]] <- simulated_sce_list
  rm(simulated_sce_list)  # Clear the temporary variable
}

##### merge rowData #####
# Function to process rowData for a single simulated dataset
process_rowData <- function(sce_list) {
  do.call(cbind, lapply(sce_list[1:2], function(sce_one) {
    sample <- unique(sce_one$Sample)
    group <- ifelse(sample %in% c( "10DPI_2", "20DPI_1",  "2DPI_2",  "5DPI_1"), "1", "2")

    df <- rowData(sce_one)
    
    # Rename columns
    names(df)[names(df) == 'in_pattern'] <- paste0('in_pattern_', group)
    names(df)[names(df) == 'spatial_prob'] <- paste0('spatial_prob_', group)
    
    # Remove the 'group' column if it exists
    df <- df[, !(colnames(df) %in% "group")]
    
    return(df)
  }))
}

# Apply the rowData processing function to all simulated lists
simulated_rowData <- lapply(simulated_list, process_rowData)

# Combine all results into a single data frame
rowData_df <- do.call(rbind, simulated_rowData)
rowData_df <- rowData_df[,c(1:3,6:8)]
dim(rowData_df)
##### merge colData #####
head(colData(simulated_list[[1]][[1]]))
meta_all <- lapply(simulated_list, function(scenario_list){
  meta_df <- do.call(rbind,
                     lapply(scenario_list, 
                            function(sce) colData(sce))
  )
  names(meta_df)[names(meta_df) == 'Pattern'] <- paste0('Pattern_', 
                                                        unique(meta_df$Scenario))
  meta_df <- meta_df[, !(colnames(meta_df) %in% c("Scenario",
                                                  "sizeFactor"))]
  meta_df
})
meta_all_combined <- do.call(cbind,meta_all)
meta_all_combined <- meta_all_combined[, !duplicated(colnames(meta_all_combined))]
dim(meta_all_combined)
##### merge counts #####
counts_list <- lapply(simulated_list, function(sce){
  do.call(cbind, lapply(sce, counts))
})
gc()
raw_counts <- do.call(rbind,
                      counts_list)
dim(raw_counts)
sce_object <- SingleCellExperiment(assays = list(counts = raw_counts),
                                   rowData =  rowData_df,
                                   colData = meta_all_combined)
# Log normalize counts
sce_object <- scuttle::logNormCounts(sce_object)
saveRDS(sce_object, file = paste0(path, "combined_simulated.rda"))