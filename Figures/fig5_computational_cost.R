rm(list = ls())
library(ggplot2)
library(dplyr)
library(forcats)
library(readr)
library(lubridate)

path = "./data/"
path_save = "./figure/"

color_bar = c(   "DESpace_global" = "#A6CEE3",
                 "DESpace_individual" = "#08306B",
                 "DESpace_smooth_spline" = "blue",
                 "Seurat_FindAllMarkers" = "#B2DF8A",#"#33A02C",
                 "seurat_FindAllMarkers" = "#B2DF8A",
                 "Seurat_FindMarkers" = "#B2DF8A",
                 "FindAllMarkers" = "#B2DF8A",
                 "Seurat_PseudoBulk_FindMarkers" = "#F0E68C",
                "spatialLIBD_pairwise" = "mediumpurple",
                "scran_findMarkers" = "#B22222", 
                "findMarkers" = "#B22222",
                "spatialLIBD_enrichment" = "mediumpurple",
                "spatialLIBD_anova" = "mediumpurple") 
            

df <- read.csv(paste0(path, "timing_results.csv"))
df <- df %>%
  mutate(
    Time_Minutes = round(elapsed / 60, 1)# Convert elapsed to minutes
    #color_all = color_bar  
    ) %>%
    mutate(Method_clean = str_remove(method, "^(Seurat_|seurat_|scran_)")) %>%  
    mutate(Method = if_else(group == 5,
                            paste0(Method_clean, " (5 groups)"),
                            Method_clean)) %>%
  mutate(Method = recode(Method, 
                         "PseudoBulk_FindMarkers" = "pseudo-bulk_FindMarkers")) %>%  # Add more replacements as needed
  mutate(Method = fct_reorder(Method, elapsed, .fun = sum)) %>%  # Reorder x-axis by elapsed
  group_by(Method) %>%
  mutate(label_ypos = cumsum(elapsed) - (elapsed / 2)  ) # Compute label positions

legend_labels <- c(
  "DESpace_global" = "DESpace2 global",
  "DESpace_individual" = "DESpace2 individual",
  "FindMarkers" = "FindMarkers (Seurat)",
  "pseudo-bulk_FindMarkers" = "pseudo-bulk FindMarkers",
  "findMarkers" = "findMarkers (scran)",
  "spatialLIBD_pairwise" = "spatialLIBD"
)

gg_Time_2groups <- df %>%
  filter(Method %in% c("DESpace_global", "DESpace_individual", "FindMarkers", 
                       "pseudo-bulk_FindMarkers", "findMarkers", "spatialLIBD_pairwise")) %>%
  mutate(Method = recode(Method, !!!legend_labels)) %>%
  ggplot(aes(x = Method, y = elapsed, fill = method)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # No fill here
  geom_text(aes(label = round(elapsed), y = elapsed, hjust = -0.1), 
  size = 3.5, color = "black") +
  theme_bw() +
  scale_fill_manual(values = color_bar)+
  xlab("") +
  ylab("Time (Seconds)") +
  #coord_trans(y = "sqrt") +
  #scale_y_sqrt(breaks = c(10, 50, 5000, 9500), limits = c(0, 9500)) +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000),
    limits = c(1, 10000)
  ) +
  theme(
    axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1),
    axis.text.y = element_text(size = rel(1.2), angle = 0, hjust = 1),
    axis.title = element_text(size = rel(1.2)),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(size = 0.2, color = "grey", linetype = 2),
    #aspect.ratio = 0.5,
    legend.position = "none"
  ) +
  coord_flip()

ggsave(filename = paste0(path_save, "computation_cost.png"),
       plot = gg_Time_2groups,
       width = 8.7, height = 4.1, 
       units = "in",
       dpi = 600, bg = "white")

#### 5-groups comparison
legend_labels <- c(
  "DESpace_global (5 groups)" = "DESpace2 global",
  "DESpace_individual (5 groups)" = "DESpace2 individual",
  "DESpace_smooth_spline (5 groups)" = "spline-based DESpace2 ",
  "FindAllMarkers (5 groups)" = "FindMarkers (Seurat)",
  "findMarkers (5 groups)" = "findMarkers (scran)",
  "spatialLIBD_anova (5 groups)" = "spatialLIBD"
)
(gg_Time_5groups <- df %>%
  filter(Method %in% c("DESpace_global (5 groups)", "DESpace_individual (5 groups)", 
                       "DESpace_smooth_spline (5 groups)", 
                       "FindAllMarkers (5 groups)",
                       "findMarkers (5 groups)", 
                       "spatialLIBD_anova (5 groups)")) %>%
    mutate(Method_rename = recode(Method, !!!legend_labels)) %>%
  ggplot(aes(x = Method_rename, y = elapsed, fill = sub(" \\(5 groups\\)", "", Method))) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = round(elapsed), y = elapsed, hjust = -0.1
                ), 
              size = 2.5, color = "black") +
  theme_bw() +
  scale_fill_manual(values = color_bar) +
  xlab("") +
  ylab("Time (Seconds)") +
  coord_trans(y = "sqrt") +
  #scale_y_sqrt(breaks = c(10, 100, 1000, 45000), limits = c(0, 47000)) +
    scale_y_log10(
      breaks = c(1, 10, 100, 1000, 10000, 45000),
      limits = c(1, 50000)) +
  theme(
    axis.text.x = element_text(size = rel(0.8), angle = 90, hjust = 1),
    axis.text.y = element_text(size = rel(0.8), angle = 0, hjust = 1),
    axis.title = element_text(size = rel(0.8)),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(size = 0.2, color = "grey", linetype = 2),
    #aspect.ratio = 0.5,
    legend.position = "none"
  ) +
  coord_flip())

ggsave(filename = paste0(path_save, "supp_computation_cost.png"),
       plot = gg_Time_5groups,
       width = 8, height = 3, 
       units = "in",
       dpi = 600, bg = "white")
