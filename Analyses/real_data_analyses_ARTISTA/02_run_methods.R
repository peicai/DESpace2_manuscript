library(devtools)
# library(DESpace) # manually source code if R version < 4.5
path <- "./"
source("01_competitors_wrapper.R")
########################################################################################
################################ 2 vs. 20 DPI ##########################################
########################################################################################
load("./muSpaData_Data/Wei22_full.rda")
spe <- spe[, spe$condition %in% c("2DPI", "20DPI")]
############################ DESpace 2 conditions ###################
A = system.time({res_edgeR <- dsp_test(
  spe,
  cluster_col = "Banksy_smooth",
  sample_col = "sample_id",
  condition_col = "condition",
  filter_gene = FALSE,
  filter_cluster = FALSE,
  verbose = TRUE
)})

save(res_edgeR, file = paste0(path, "DESpace_global_test_2vs20.rda")); rm(res_edgeR)
B = system.time({res_edgeR <- individual_dsp(
  spe,
  cluster_col = "Banksy_smooth",
  sample_col = "sample_id",
  condition_col = "condition",
  filter_gene = FALSE,
  filter_cluster = FALSE
)
})
save(res_edgeR, file = paste0(path, "DESpace_individual_test_2vs20.rda")); rm(res_edgeR)

############################ competitors 2vs20 ########################################
spe$Sample <- spe$sample_id
spe$Condition <- spe$condition
spe$Condition <- ifelse(spe$Condition == "2DPI", "Condition1", "Condition2")
spe$barcode <- colnames(spe)
D = system.time({run_seurat_FindAllMarkers(spe,
                                           output_path = path,
                                           cluster_method = "Banksy_2vs20",
                                           cluster_col = "Banksy_smooth")})

D2 = system.time({run_seurat_FindMarkers(spe,
                                         output_path = path,
                                         cluster_method = "Banksy_2vs20",
                                         cluster_col = "Banksy_smooth")})
E = system.time({run_seurat_PseudoBulk_FindMarkers(spe,
                                                   output_path = path,
                                                   cluster_method = "Banksy_2vs20",
                                                   cluster_col = "Banksy_smooth")})
G = system.time({ run_spatialLIBD_pairwise(spe,
                                         output_path = path,
                                         cluster_method = "Banksy_2vs20_v2",
                                         cluster_col = "Banksy_smooth")})
H = system.time({ run_scran_findMarkers(spe,
                                        output_path = path,
                                        cluster_method = "Banksy_2vs20",
                                        cluster_col = "Banksy_smooth")})

system_times <- mget(ls(), envir = .GlobalEnv, inherits = FALSE)
system_times <- system_times[sapply(system_times, inherits, "proc_time")]

computational_cost <- do.call(rbind, system_times) %>% as.data.frame()
computational_cost$method <- c("DESpace_global","DESpace_individual",
                               "DESpace_global","DESpace_smooth_spline",
                               "DESpace_individual","Seurat_FindAllMarkers",
                               "Seurat_FindMarkers","Seurat_PseudoBulk_FindMarkers", 
                               "spatialLIBD_pairwise","scran_findMarkers",
                               "spatialLIBD_pairwise","spatialLIBD_anova", "spatialLIBD_enrichment",
                               "seurat_FindAllMarkers","scran_findMarkers")
computational_cost$group <- c(2,2,5,5,5,2,2,2,2,2,5,5,5,5,5)
write.csv(computational_cost, file = paste0(path, "system_times.csv"), row.names = TRUE)