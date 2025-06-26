rm(list = ls())
source("01_simulation_function.R", echo=TRUE)

data_path <- "./Real_data/LIBD_filtered/"
output_path <- "./"

sample_id <- list(c("151507", "151509"),
                  c("151669", "151671"),
                    c("151673", "151675"))

sce_objects <- lapply(unlist(sample_id), function(id) {
  print(id)
  load(paste0(data_path, id, "_sce_LIBD.rda"))
  sce_one$layer_guess_reordered_droplevel <- sce_one$layer_guess_reordered
  ind <- sce_one$layer_guess_reordered_droplevel %in% c("Layer1", "Layer2")
  if (sum(ind) != 0) sce_one[,ind]$layer_guess_reordered_droplevel <- "Layer6"
  sce_one$layer_guess_reordered_droplevel <- droplevels(sce_one$layer_guess_reordered_droplevel)
  sce_one <- sce_one[,!is.na(sce_one$layer_guess_reordered_droplevel)]
  print(levels(sce_one$layer_guess_reordered_droplevel))
  return(sce_one)
})

common_gene <- Reduce(intersect, lapply(sce_objects, rownames))

# Calculate the length of each split
split_length <- ceiling(length(common_gene) / 4)

# Split the vector into approximately equal-length parts
split_gene <- split(common_gene, rep(1:4, each = split_length, 
                                    length.out = length(common_gene)))
scenarios <-c("DSP1", "DSP2", "NULL1", "NULL2")
## DSP1: unif vs. SV
## DSP2: SV1 vs. SV2
## NULL1: unif vs. unif
## NULL2: SV vs. SV

# First pair of samples: 151507 & 151509
## for scenario s
probs <- list(c(0.5, 1),
              c(1, 1),
              c(0.5, 0.5),
              c(1, 1))
region_all <- c("WM", "Layer6", "Layer5", "Layer4", "Layer3")

for (s in 1:4){
  SVG <- split_gene[[s]]
  # in DSP_1 (Unif vs. SV), for each gene, randomly choose the group which is Unif/SV.
  spatial_probs <- probs[[s]]
  scenario_each <- scenarios[s]
  
  inf <- lapply(SVG, function(g){
    # Every time we simulate a SV pattern, 
    # the “in” pattern randomly chosen among the possible clusters: 'region_all'
    regions <- determine_regions(scenario_each, region_all)
    # In DSP1 (Unif vs. SV), for each gene, randomly choose the group which is Unif/SV.
    if (scenario_each == "DSP1") {
      ind <- sample.int(2,2)
      spatial_probs <- spatial_probs[ind]
      regions <- regions[ind]
    }
    return(list(regions, spatial_probs))
  })
  names(inf) <- SVG
  
  for (i in 1:3){
    #paired_sample <- sample_id[[i]]
    ### Condition 1: sample 1, 3, 5 = 2*i-1
    subset_sce1 <- sce_objects[[(2*i-1)]][SVG, ]
    subset_sce1$Scenario <- scenario_each
    subset_sce1$Condition <- "Condition_1"
    subset_sce1$Sample <- unlist(sample_id)[(2*i-1)]
    ### Condition 2: sample 2, 4, 6 = 2*i
    subset_sce2 <- sce_objects[[(2*i)]][SVG, ]
    subset_sce2$Scenario <- scenario_each
    subset_sce2$Condition <- "Condition_2"
    subset_sce2$Sample <- unlist(sample_id)[(2*i)]

    processed_sce <- lapply(list(subset_sce1,subset_sce2), process_subset)
    simulated_sce_list <- CombineSimulatedObject(processed_sce_list = processed_sce,
                                  pattern_name = 'Pattern',
                                  spatial_prob = lapply(inf, `[[`, 2),
                                  scenario = scenario_each,
                                  regions = lapply(inf, `[[`, 1))
    save(simulated_sce_list, file = paste0(output_path, 
                                           paste(sample_id[[i]], collapse = "_"),
                                           "_",
                                           scenario_each, 
                                           '_simulated.rda'))
    
  }
  save(inf, file = paste0(output_path, "RegionsGT_",
                                         scenario_each, 
                                         '_simulated.rda'))
}