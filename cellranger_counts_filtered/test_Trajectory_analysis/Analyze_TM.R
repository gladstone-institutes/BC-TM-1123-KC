setwd("~/Dropbox (Gladstone)/Gladstone chapters/project_files/Tamer_TM-1123/cellranger_counts_filtered/test_Trajectory_analysis")
library(monocle3)
library(magrittr)

load("../dat_with_TM_defined_clustering_02_26_20.RData")

seurat <- dat

gene_annotation <- seurat@reductions[["pca"]]@feature.loadings %>%
  rownames() %>% 
  as.data.frame(., row.names = .) %>%
  set_colnames(., "gene_short_name")

# gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
# colnames(gene_annotation) <- "gene_short_name"

# part two, cell information

cell_metadata <- seurat@assays[["RNA"]]@counts@Dimnames[[2]] %>%
  as.data.frame(., row.names = .) %>%
  set_colnames(., "barcode")
cell_metadata <- cbind(cell_metadata, 
                       seurat@meta.data[rownames(cell_metadata) %>% as.character(), ])

# part three, counts sparse matrix

expression_matrix <- seurat@assays[["RNA"]]@counts
expression_matrix <- expression_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]

### Construct the basic cds object

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

reducedDims(cds_from_seurat)$UMAP <- seurat@reductions$umap@cell.embeddings
reducedDims(cds_from_seurat)$PCA <- seurat@reductions$pca@cell.embeddings

plot_cells(cds_from_seurat, reduction_method = "UMAP",
           color_cells_by = "manual_clustering_by_Tamer")


list_cluster <- seurat@meta.data$manual_clustering_by_Tamer
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

#-------------------------
### Could be a space-holder, but essentially fills out louvain parameters
#Need to check
# cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

### Assign feature loading for downstream module analysis

# cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings

#-----------------------------

cds_from_seurat <- learn_graph(cds_from_seurat)

pdf("Trajectory_from_monocle_for_all_TM_clusters.pdf",
    width = 10, height =10)
plot_cells(cds_from_seurat,
           reduction_method = "UMAP",
           color_cells_by = "manual_clustering_by_Tamer",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
dev.off()

#-----------------------
#-----------------------
#Try other clustering approaches
# clust <- "Cell_sample_ident"
# list_cluster <- seurat@meta.data[[clust]]
# names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]
# 
# cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
# 
# #-------------------------
# ### Could be a space-holder, but essentially fills out louvain parameters
# #Need to check
# # cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
# 
# recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
# names(recreate.partition) <- cds_from_seurat@colData@rownames
# recreate.partition <- as.factor(recreate.partition)
# 
# cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

### Assign feature loading for downstream module analysis

# cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings

#-----------------------------

# cds_from_seurat <- learn_graph(cds_from_seurat)
# 
# plot_cells(cds_from_seurat,
#            reduction_method = "UMAP",
#            color_cells_by = "manual_clustering_by_Tamer",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE)

#-------------------------
#-------------------------
#Remove TM10
to_remove <- list(c("TM10"), 
                  c("TM7", "TM8", "TM9", "TM10"),
                  c("TM1", "TM2", "TM7", "TM8", "TM9", "TM10"))
for (x in to_remove) {
  filtered_cell_metadata <- subset(cell_metadata, 
                                   !(manual_clustering_by_Tamer %in% x))
  keep <- colnames(expression_matrix) %in% (filtered_cell_metadata$barcode %>% 
                                              as.character())
  filtered_expression_matrix <- expression_matrix[, keep]
  filtered_cds_from_seurat <- new_cell_data_set(filtered_expression_matrix,
                                                cell_metadata = filtered_cell_metadata,
                                                gene_metadata = gene_annotation)
  
  keep <- rownames(seurat@reductions$umap@cell.embeddings) %in% (filtered_cell_metadata$barcode %>% 
                                                                   as.character())
  reducedDims(filtered_cds_from_seurat)$UMAP <- seurat@reductions$umap@cell.embeddings[keep, ]
  keep <- rownames(seurat@reductions$pca@cell.embeddings) %in% (filtered_cell_metadata$barcode %>% 
                                                                  as.character())
  reducedDims(filtered_cds_from_seurat)$PCA <- seurat@reductions$pca@cell.embeddings[keep, ]
  
  # plot_cells(filtered_cds_from_seurat, reduction_method = "UMAP",
  #            color_cells_by = "manual_clustering_by_Tamer")
  
  
  list_cluster <- filtered_cell_metadata$manual_clustering_by_Tamer
  names(list_cluster) <- rownames(filtered_cell_metadata)
  
  filtered_cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
  
  #-------------------------
  ### Could be a space-holder, but essentially fills out louvain parameters
  #Need to check
  # cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
  
  recreate.partition <- c(rep(1, length(filtered_cds_from_seurat@colData@rownames)))
  names(recreate.partition) <- filtered_cds_from_seurat@colData@rownames
  recreate.partition <- as.factor(recreate.partition)
  
  filtered_cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
  
  ### Assign feature loading for downstream module analysis
  
  # cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings
  
  #-----------------------------
  
  filtered_cds_from_seurat <- learn_graph(filtered_cds_from_seurat)
  
  pdf(paste0("Trajectory_from_monocle_for_TM_clusters_excluding_",
      paste(x, collapse = "_"), ".pdf"),
      width = 5, height =5)
  print(plot_cells(filtered_cds_from_seurat,
                   reduction_method = "UMAP",
                   color_cells_by = "manual_clustering_by_Tamer",
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE))
  dev.off()
}

# plot_cells(filtered_cds_from_seurat,
#            color_cells_by = "manual_clustering_by_Tamer",
#            label_cell_groups=FALSE,
#            label_leaves=TRUE,
#            label_branch_points=TRUE,
#            graph_label_size=1.5)
# 
# filtered_cds_from_seurat <- order_cells(filtered_cds_from_seurat)

#-----------------------
#-----------------------
library(scmap)
tm1_or_tm2 <- cell_metadata %>% 
  subset(., manual_clustering_by_Tamer %in% c("TM1", "TM2"), 
         select = "barcode") %>% 
  rownames()
atlas_data <- expression_matrix[, tm1_or_tm2] %>% 
  apply(., 2, function(x) x*10000/sum(x)) #Normalize using the Seurat 3 method.
atlas_ann <- cell_metadata[tm1_or_tm2, "manual_clustering_by_Tamer", drop = F]

tm3 <- cell_metadata %>% 
  subset(., manual_clustering_by_Tamer %in% c("TM3"), select = "barcode") %>% 
  rownames()
query_data <- expression_matrix[, tm3] %>%
  apply(., 2, function(x) x*10000/sum(x)) #Normalize using the Seurat 3 method.
query_ann <- cell_metadata[tm3, "manual_clustering_by_Tamer", drop = F]

objLabels <- c(query = "Exp", atlas = "Atlas")
for (i in c("query", "atlas")) {
  
  #Create the SingleCellExperiment object.
  
  sce  <-  SingleCellExperiment(assays = i %>% 
                                  paste0(., "_data") %>% 
                                  get() %>% 
                                  as.matrix() %>%
                                  list(normcounts = .),
                                colData = i %>%
                                  paste0(., "_ann") %>%
                                  get()
  )
  
  #Log transform the normalized counts.
  logcounts(sce) <- sce %>%
    normcounts() %>%
    `+`(., 1) %>%
    log2()
  
  rowData(sce)$feature_symbol <- rownames(sce)
  
  #Assign this object to scExp if i is equal to "query" or, to scAtlas if i is equal to "ref". 
  this_name <- objLabels[i] %>% paste0("sc", .)
  assign(this_name, sce)
}

scAtlas <- scAtlas %>% 
  selectFeatures(.,
                 suppress_plot = FALSE)

scAtlas <- scAtlas %>% indexCluster(., cluster_col = "manual_clustering_by_Tamer")
scmapCluster_results <- scmapCluster(
  projection = scExp, 
  index_list = scAtlas %>%
    metadata() %>%
    `[[`(j = "scmap_cluster_index") %>%
    list(Atlas = .)
)

plot(
  getSankey(
    scExp %>%
      colData() %>%
      `[[`(i = "manual_clustering_by_Tamer"), 
    scmapCluster_results %>%
      `[[`(j = "scmap_cluster_labs") %>%
      as.character(),
    plot_height = 400
  )
)


