library("spatialLIBD")
library("SingleCellExperiment")
library(Seurat)
library(data.table)
registration_pseudobulk_modified <-
  function(
    sce,
    var_registration,
    var_sample_id,
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL) {
    ## Check that inputs are correct
    stopifnot(is(sce, "SingleCellExperiment"))
    stopifnot(var_registration %in% colnames(colData(sce)))
    stopifnot(var_sample_id %in% colnames(colData(sce)))
    stopifnot(all(
      !c("registration_sample_id", "registration_variable") %in% colnames(colData(sce))
    ))
    
    ## Avoid any incorrect inputs that are otherwise hard to detect
    stopifnot(!var_registration %in% covars)
    stopifnot(!var_sample_id %in% covars)
    stopifnot(var_registration != var_sample_id)
    
    ## Check that the values in the registration variable are ok
    uniq_var_regis <- unique(sce[[var_registration]])
    if (any(grepl("\\+|\\-", uniq_var_regis))) {
      stop(
        "Remove the + and - signs in colData(sce)[, '",
        var_registration,
        "'] to avoid downstream issues.",
        call. = FALSE
      )
    }
    
    ## Pseudo-bulk for our current BayesSpace cluster results
    message(Sys.time(), " make pseudobulk object")
    ## I think this needs counts assay
    sce_pseudo <- scuttle::aggregateAcrossCells(
      sce,
      DataFrame(
        registration_variable = sce[[var_registration]],
        registration_sample_id = sce[[var_sample_id]]
      )
    )
    colnames(sce_pseudo) <-
      paste0(
        sce_pseudo$registration_sample_id,
        "_",
        sce_pseudo$registration_variable
      )
    
    ## Check that the covariates are present
    if (!is.null(covars)) {
      for (covariate_i in covars) {
        if (sum(is.na(sce_pseudo[[covariate_i]])) == ncol(sce_pseudo)) {
          stop(
            "Covariate '",
            covariate_i,
            "' has all NAs after pseudo-bulking. Might be due to not being a sample-level covariate.",
            call. = FALSE
          )
        }
      }
    }
    
    ## Drop pseudo-bulked samples that had low initial contribution
    ## of raw-samples. That is, pseudo-bulked samples that are not
    ## benefiting from the pseudo-bulking process to obtain higher counts.
    if (!is.null(min_ncells)) {
      message(
        Sys.time(),
        " dropping ",
        sum(sce_pseudo$ncells < min_ncells),
        " pseudo-bulked samples that are below 'min_ncells'."
      )
      sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= min_ncells]
    }
    
    if (is.factor(sce_pseudo$registration_variable)) {
      ## Drop unused var_registration levels if we had to drop some due
      ## to min_ncells
      sce_pseudo$registration_variable <- droplevels(sce_pseudo$registration_variable)
    }
    
    ## Drop lowly-expressed genes
    # message(Sys.time(), " drop lowly expressed genes")
    # keep_expr <-
    #   edgeR::filterByExpr(sce_pseudo, group = sce_pseudo$registration_variable)
    # sce_pseudo <- sce_pseudo[which(keep_expr), ]
    # 
    ## Compute the logcounts
    message(Sys.time(), " normalize expression")
    logcounts(sce_pseudo) <-
      edgeR::cpm(edgeR::calcNormFactors(sce_pseudo),
                 log = TRUE,
                 prior.count = 1
      )
    
    if (is(sce_pseudo, "SpatialExperiment")) {
      ## Drop things we don't need
      spatialCoords(sce_pseudo) <- NULL
      imgData(sce_pseudo) <- NULL
    }
    if (!is.null(pseudobulk_rds_file)) {
      message(Sys.time(), " saving sce_pseudo to ", pseudobulk_rds_file)
      saveRDS(sce_pseudo, file = pseudobulk_rds_file)
    }
    return(sce_pseudo)
  }
run_LIBD_registration1 <- function(sce.combined,
                                   output_path = "",
                                   cluster_method = "BayesSpace",
                                   cluster_col = "matched_cluster"){
  ## Add gene-level information
  rowData(sce.combined)$gene_name <- rownames(sce.combined)
  sce.combined$Cluster <- as.factor(as.character(sce.combined[[cluster_col]]))

  results_res <- list()
  for (each_layer in levels(sce.combined$Cluster)){
    sce_subset = sce.combined[, colData(sce.combined)[, "Cluster"] == each_layer]
    sce_pseudo <- registration_pseudobulk_modified(
      sce = sce_subset,
      var_registration = "Condition",
      var_sample_id = "Sample",
      min_ncells = 0)
    
    registration_mod <- registration_model(sce_pseudo)
    head(registration_mod)
    block_cor <- registration_block_cor(sce_pseudo, registration_mod)
  
  results_pairwise <- registration_stats_pairwise(sce_pseudo,
                                                  registration_mod,
                                                  block_cor,
                                                  gene_name = "gene_name",
                                                  gene_ensembl = "gene_name"
  )
  results_res[[each_layer]] <- results_pairwise
   }
  save(results_res, file = paste0(output_path, "registration_pairwiseNaN_", 
                                   cluster_method, ".rda"))
}

# This will instruct findMarkers() to perform pairwise t
# -tests between clusters using only cells on the same level of the blocking factor. It will then combine p
# -values from different plates using Stouffer’s Z method to obtain a single p
# -value per gene. http://129.217.206.11/packages/3.9/workflows/vignettes/simpleSingleCell/inst/doc/de.html
run_scran_findMarkers <- function(sce.combined,
                                  output_path = "",
                                  cluster_method = "BayesSpace",
                                  cluster_col = "matched_cluster"){
  expr <- assays(sce.combined)$counts
  libsizes <- colSums(expr)
  size.factors <- libsizes/mean(libsizes)
  logcounts(sce.combined) <- log2(t(t(expr)/size.factors) + 1)
  sce.combined$Cluster <- as.factor(as.character(sce.combined[[cluster_col]]))
  sce.combined$Donor <- ifelse(sce.combined$Sample %in% c("151507", "151509"), "Donor1", 
                               ifelse(sce.combined$Sample %in% c("151669", "151671"), "Donor2", "Donor3"))
  results_res <- list()
  for (each_layer in levels(sce.combined$Cluster)){
    sce_subset = sce.combined[, colData(sce.combined)[, "Cluster"] == each_layer]
    layer = as.factor(colData(sce_subset)[["Condition"]])
    markers_block <- scran::findMarkers(
      sce_subset, groups = layer,
    block = sce_subset$Donor, 
    # effectively performs the cluster comparisons in each batch, 
    # and subsequently combines the results into a single p-value
    pval.type = "all"
  )
    results_res[[each_layer]] <- markers_block[[1]] # Results for condition 1 and condition 2 are same, only have opposite logFC
  }
  save(results_res, file = paste0(output_path, "scran_findMarkers_", 
                                  cluster_method, ".rda"))
}

# https://satijalab.org/seurat/articles/de_vignette.html
run_seurat_PseudoBulk_FindMarkers <- function(sce.combined,
                                              output_path = "",
                                              cluster_method = "BayesSpace",
                                              cluster_col = "matched_cluster"){
  sce.combined$Cluster <- as.character(sce.combined[[cluster_col]])
  sce.combined$Cluster <- droplevels(factor( sce.combined$Cluster))
  counts <- assays(sce.combined)[[1]]
  colnames(counts) <- paste0(colData(sce.combined)$barcode, "_", colData(sce.combined)$Sample)
  seurat <- CreateSeuratObject(counts = counts)
  t2 <- seurat@meta.data
  t <- data.frame(colData(sce.combined))
  rownames(t) <- paste0(t$barcode, "_", t$Sample)
  t3 <- data.frame(t,t2)
  seurat@meta.data <- t3
  layer = as.factor(colData(sce.combined)[["Cluster"]])
  Idents(seurat) <- layer
  # pseudobulk the counts based on donor-condition-celltype
  pseudo_seurat <- AggregateExpression(seurat, assays = "RNA", return.seurat = T, 
                                     group.by = c("Condition", "Sample","Cluster"))
  
  # each 'cell' is a donor-condition-celltype pseudobulk profile
  tail(Cells(pseudo_seurat))
  
  pseudo_seurat$celltype.stim <- paste(pseudo_seurat$Condition, pseudo_seurat$Cluster, sep = "_")
  Idents(pseudo_seurat) <- "celltype.stim"
  n_cluster <-  nlevels(sce.combined$Cluster)
  grp1 <- paste0("Condition-1_",levels(sce.combined$Cluster))
  grp2 <- paste0("Condition-2_",levels(sce.combined$Cluster))
  results_res <- list()
  for(i in seq_len(n_cluster)){
    domain_test <- levels(sce.combined$Cluster)[i]
    res <- Seurat::FindMarkers(pseudo_seurat,
                        ident.1 = grp1[i],  
                        ident.2 = grp2[i], 
                        #grouping.var = "Sample",
                        #logfc.threshold = -Inf,
                        #min.pct = 0,
                        test.use = "DESeq2"# used in the vignette; use DESeq2 for pesudobulk/multiple
                        #min.diff.pt = -Inf,
                        #min.cells.feature = 0,
                        #min.cells.group = 0
                        # return.thresh = 2#, min.pct = 0.25, logfc.threshold = 0.25
    )
    results_res[[domain_test]] <- res
  }
  save(results_res, file = paste0(output_path, "seurat_PseudoBulk_FindMarkers_", 
                                   cluster_method, ".rda"))
}

run_seurat_FindMarkers <- function(sce.combined,
                                      output_path = "",
                                      cluster_method = "BayesSpace",
                                      cluster_col = "matched_cluster"){
  sce.combined$Cluster <- as.character(sce.combined[[cluster_col]])
  counts = assays(sce.combined)[[1]]
  colnames(counts) <- paste0(colData(sce.combined)$barcode, "_", colData(sce.combined)$Sample)
  seurat_object <- CreateSeuratObject(counts = counts                                       )
  t2 <- seurat_object@meta.data
  t <- data.frame(colData(sce.combined))
  rownames(t) <- paste0(t$barcode, "_", t$Sample)
  t3 <- data.frame(t,t2)
  seurat_object@meta.data <- t3
  results_res1 <- list()
  results_res2 <- list()
  # For each cluster:
  for (each_layer in levels(factor(sce.combined$Cluster))){
    seurat_subset = seurat_object[, seurat_object@meta.data[, "Cluster"] == each_layer]
    layer = as.factor(seurat_subset@meta.data[["Condition"]])
    Idents(seurat_subset) <- layer
    seurat_norm = NormalizeData(seurat_subset, 
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000)
    res <- Seurat::FindMarkers(seurat_norm,
                                ident.1 = "Condition_1",
                                grouping.var = "Sample",
                                verbose = T, only.pos = F,
                                logfc.threshold = -Inf,
                                min.pct = 0,
                                test.use = "wilcox",
                                min.diff.pt = -Inf,
                                min.cells.feature = 0,
                                min.cells.group = 0,
                                return.thresh = 2)
    results_res1[[each_layer]] <- res
  }

  save(results_res1, file = paste0(output_path, "seurat_FindMarkers_", 
                                   cluster_method, ".rda"))

}
# We compute a per-gene adjusted p-value, 
# using the perGeneQValue function, 
# which aggregates evidence from multiple tests within a gene to a single p-value for the gene and 
# then corrects for multiple testing across genes
# the use case: given per-cluster p-values for null hypothesis H0 (not significant difference across conditions), 
# we can determine the number of genes in which at least for one cluster H0 is rejected
# what is the associated false discovery rate?

## MODIFY: only to skip "stopifnot( is(object, "DEXSeqResults"))"
## in the original function
perGeneQValue_modified <- function(res_long, 
                                   pval = "p_val",
                                   padj = "p_val_adj",
                                   gene_id = "gene",
                                   ...
                                   ){
  wTest <- which( !is.na( res_long[[padj]] ) )
  ## use only those exons that were testable
  pvals  = res_long[[pval]][wTest]
  ## 'factor' removes ununsed levels
  geneID    = factor(res_long[[gene_id]][wTest])
  geneSplit = split(seq(along=geneID), geneID)
  ## summarise p-values of exons/clusters for one gene: take the minimum
  pGene = sapply(geneSplit, function(i) min(pvals[i]))
  stopifnot(all(is.finite(pGene)))
  ## Determine the thetas to be used
  theta = unique(sort(pGene))
  
  ## compute q-values associated with each theta
  q = perGeneQValueExact(pGene, theta, geneSplit)
  ## return a named vector of q-values per gene
  res_long        = rep(NA_real_, length(pGene))
  res_long        = q[match(pGene, theta)]
  res_long = pmin(1, res_long)
  names(res_long) = names(geneSplit)
  # ## compute q-values associated with each theta
  # q = perGeneQValueExact(pGene, theta, geneSplit)
  res_long = as.data.frame(res_long)
  colnames(res_long) <- "qval"
  return(res_long)
}

perGeneQValueExact <- function(pGene, theta, geneSplit) {
  stopifnot(length(pGene)==length(geneSplit))
  
  ## Compute the numerator \sum_{i=1}^M 1-(1-theta)^{n_i}
  ## Below we first identify the summands which are the same
  ## (because they have the same n_i), then do the sum via the
  ## mapply
  numExons     = listLen(geneSplit)
  tab          = tabulate(numExons)
  notZero      = (tab>0)
  numerator    = mapply(function(m, n) m * (1 - (1-theta)^n),
                        m = tab[notZero],
                        n = which(notZero))
  numerator    = rowSums(numerator)
  
  ## Compute the denominator: for each value of theta, the number
  ## of genes with pGene <= theta[i].
  ## Note that in cut(..., right=TRUE), the intervals are
  ## right-closed (left open) intervals.
  bins   = cut(pGene, breaks=c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)
  counts = tabulate(bins, nbins = nlevels(bins))
  denom  = cumsum(counts)
  stopifnot(denom[length(denom)]==length(pGene))
  
  return(numerator/denom)
}
