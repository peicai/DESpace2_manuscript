rm(list = ls())
source("./simulation_ARTISTA/01_simulation_function.R", echo=TRUE)
output_path <- "./Simulated_data/cluster_res/"
samples <- c("2DPI_1", "2DPI_2", 
             "5DPI_1", "5DPI_2", 
             "10DPI_1", "10DPI_2",
             "20DPI_1", "20DPI_2" )
dpi <- sub("_.*", "", samples)

set.seed(17) 
dpi_groups <- split(samples, dpi)
selected_samples <- unlist(lapply(dpi_groups, function(group) {
  sample(group, 2)
}))
(group1_samples <- selected_samples[seq(1,8,2)])
# "10DPI_2" "20DPI_1"  "2DPI_2"  "5DPI_1" 
(group2_samples <- selected_samples[seq(2,8,2)])
# "10DPI_1" "20DPI_2"  "2DPI_1"  "5DPI_2" 

spe <- readRDS("./Simulated_data/cluster_res/spe_nclust_rotated.rds")
spe_filtered = spe[, spe$sample_id %in% c(group1_samples,
                                          group2_samples)];
dim(spe_filtered); rm(spe)
table(spe_filtered$nclust_2, spe_filtered$sample_id)
common_gene <- rownames(spe_filtered)
split_length <- ceiling(length(common_gene) / 4)

set.seed(123) 
shuffled_gene <- sample(common_gene)  
split_gene <- split(shuffled_gene, rep(1:4, length.out = length(shuffled_gene)))
scenarios <-c("DSP1", "DSP2", "NULL1", "NULL2")
## DSP1: unif vs. SV
## DSP2: SV1 vs. SV2
## NULL1: unif vs. unif
## NULL2: SV vs. SV

probs <- list(c(0.5, 1),
              c(1, 1),
              c(0.5, 0.5),
              c(1, 1))
cluster_cols <- paste0("nclust_", c(2,4,6,8,10,12))
colData(spe_filtered)$barcode <- rownames(colData(spe_filtered))

create_subset_sce <- function(sce_objects, SVG, sample, condition, scenario) {
  subset_sce <- sce_objects[SVG, sce_objects$sample_id == sample]
  subset_sce$Scenario <- scenario
  subset_sce$Condition <- condition
  subset_sce$Sample <- sample
  subset_sce[, colSums(counts(subset_sce)) > 0]
  print(dim(subset_sce))
  return(subset_sce)
}
for(cluster_col in cluster_cols){
  message("Processing cluster: ", cluster_col)
  
  cluster_path <- file.path(output_path, cluster_col)
  if (!dir.exists(cluster_path)) dir.create(cluster_path, recursive = TRUE)
  
  region_all <- unique(colData(spe_filtered)[[cluster_col]])
  
  for (s in 1){
    SVG <- split_gene[[s]]
    spatial_probs <- probs[[s]]
    scenario_each <- scenarios[s]
    
    inf <- lapply(SVG, function(g){
      regions <- determine_regions(scenario_each, region_all)
      if (scenario_each == "DSP1") {
        ind <- sample.int(2,2)
        spatial_probs <- spatial_probs[ind]
        regions <- regions[ind]
      }
      return(list(rep(regions, length(group1_samples)), rep(spatial_probs,length(group1_samples))))
    })
    names(inf) <- SVG
    
    subset_list <- mapply(function(g1, g2) {
      list(
        create_subset_sce(spe_filtered, SVG, g1, "Condition_1", scenario_each),
        create_subset_sce(spe_filtered, SVG, g2, "Condition_2", scenario_each)
      )
    }, group1_samples, group2_samples, SIMPLIFY = FALSE)
    
    # Flatten the list
    subset_list <- do.call(c, subset_list)
    processed_sce <- lapply(subset_list, function(x) process_subset(subset_sce = x,
                                                                    cluster_col = cluster_col,
                                                                    coord_col = c("sdimx", "sdimy")))
    simulated_sce_list <- CombineSimulatedObject(processed_sce_list = processed_sce,
                                                 pattern_name = 'Pattern',
                                                 spatial_prob = lapply(inf, `[[`, 2),
                                                 scenario = scenario_each,
                                                 regions = lapply(inf, `[[`, 1),
                                                 cluster_col = cluster_col)
    save(simulated_sce_list, 
         file = file.path(cluster_path, paste0(scenario_each, "_simulated.rda")))
    save(inf, 
         file = file.path(cluster_path, paste0("RegionsGT_", scenario_each, "_simulated.rda")))
  }
}
