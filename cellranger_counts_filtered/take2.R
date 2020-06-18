setwd("~/Dropbox (Gladstone)/Gladstone chapters/project_files/Tamer_TM-1123/cellranger_counts_filtered")
library(Seurat)
library(magrittr)
library(tidyverse)
library(plyr)

dat <- Read10X("S1234")  %>% 
  CreateSeuratObject(., min.cells = 3, min.features = 200)

dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")

pdf("Distributions_of_nGenes_nUMIs_nMitoCounts.1234.pdf")
VlnPlot(dat, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"), 
        ncol = 3)
dev.off()

plot1 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

ggplot(data.frame(x=dat$nFeature_RNA, y = dat$nCount_RNA), aes(x)) +
  geom_histogram()

# lower_thres_nFeature <- quantile(dat$nFeature_RNA, 0.05)
# upper_thres_nFeature <- quantile(dat$nFeature_RNA, 0.95)
dat <- subset(dat, subset = percent.mt < 30 #& 
              # nFeature_RNA > lower_thres_nFeature &
              # nFeature_RNA < upper_thres_nFeature
)

dat <- NormalizeData(dat)

dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dat), 10)

# plot variable features with and without labels
# plot1 <- VariableFeaturePlot(dat)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(dat)
dat <- ScaleData(dat, features = all.genes)
dat <- RunPCA(dat, features = VariableFeatures(object = dat))
DimPlot(dat, reduction = "pca")
ElbowPlot(dat)

dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.5)
dat <- RunUMAP(dat, dims = 1:10)

pdf("UMAP_showing_clusters.1234.pdf", width = 10, height =10)
DimPlot(dat, reduction = "umap")
dev.off()

cell_sample_matching <- dat@meta.data %>%
  rownames() %>% strsplit(., split = "-") %>% 
  laply(., extract2, i = 2) %>% 
  set_names(., colnames(dat))

dat <- AddMetaData(object = dat, 
                   metadata = cell_sample_matching, 
                   col.name = "Cell_sample_ident")

pdf("UMAP_showing_sample_number.1234.pdf", width = 10, height =10)
DimPlot(dat, reduction = "umap", group.by = "Cell_sample_ident")
dev.off()

dat <- CellCycleScoring(dat, s.features = cc.genes$s.genes, 
                        g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

pdf("UMAP_showing_cell_cycle_phase_from_Seurat.1234.pdf", width = 10, height =10)
DimPlot(dat, reduction = "umap", group.by = "Phase")
dev.off()

interesting_cc_genes <- c("MKI67", 
                          "AURKB", "AURKA", "PCNA", "CDC20", "TK1", 
                          "CCNA2", "BRCA", "PLK1", "ANLN", "E2F1", "CDC25")

interesting_sarcomeric_genes <- c("TNNT2", "TNNI1", "MYH6", "MYH7", "ACTN2", 
                                  "TNNC1", "ACTA2", "MYL2", "MYL4", "PLN", "ATP1B1", 
                                  "SCN5A", "TNNI3")

overexpressed_genes <- c("CDK1", "CDK4", "CCND1", "CCNB1")

pdf("UMAP_showing_expression_of_Cell_cycle_genes.1234.pdf", width = 10, height =10)
FeaturePlot(dat, 
            features = interesting_cc_genes)
dev.off()

pdf("UMAP_showing_expression_of_sarcomeric_genes.1234.pdf", width = 10, height =10)
FeaturePlot(dat, 
            features = interesting_sarcomeric_genes)
dev.off()

pdf("UMAP_showing_expression_of_overexpressed_genes.1234.pdf", width = 10, height =10)
FeaturePlot(dat, 
            features = overexpressed_genes)
dev.off()
# Visualize the distribution of cell cycle markers across
# RidgePlot(dat, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

#---------------------
#Extract cell ids associated with clusters as manually annotated by Tamer.

cluster1 <-  dat@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  subset(., UMAP_2 < -4.5 & UMAP_1 < 2.5) %>%
  rownames() %>%
  intersect(., cell_sample_matching %>%
              subset(., equals(., 1)) %>%
              names()) %>% 
  set_names(rep("TM1", length(.)), .) 

cluster2 <- dat@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  subset(., UMAP_2 < 6 & UMAP_1 < 5.7 & UMAP_2 > -3.8 &
           UMAP_1 > -2.4 &
           (4.13*UMAP_1 + UMAP_2 > -5.01) &
           (-1.09*UMAP_1 + UMAP_2 > -5.91)) %>%
  rownames() %>%
  intersect(., cell_sample_matching %>%
              subset(., equals(., 1)) %>%
              names()) %>%
  set_names(rep("TM2", length(.)), .)

cluster3 <- dat@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  subset(., UMAP_2 < 0 & UMAP_1 < -5.2 & UMAP_2 > -4 &
           (2*UMAP_1 + UMAP_2 < -14.2)) %>%
  rownames() %>%
  intersect(., cell_sample_matching %>%
              subset(., equals(., 2)) %>%
              names()) %>%
  set_names(rep("TM3", length(.)), .)

cluster4 <- dat@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  subset(., UMAP_2 < 5.9 & UMAP_1 > -6.65 & UMAP_2 > 1.71 &
           UMAP_1 < -3.9 &
           (-2.11*UMAP_1 + UMAP_2 > 13.949)) %>%
  rownames() %>%
  intersect(., cell_sample_matching %>%
              subset(., equals(., 3)) %>%
              names()) %>%
  set_names(rep("TM4", length(.)), .)

cluster5 <- dat@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% 
  subset(., UMAP_2 < 1.71 & UMAP_2 > -3.2 & UMAP_1 < -3.7 &
           (2*UMAP_1 + UMAP_2 > -14.2) &
           (2*UMAP_1 + UMAP_2 < -8.76))  %>%
  rownames() %>%
  intersect(., cell_sample_matching %>%
              subset(., equals(., 3)) %>%
              names()) %>%
  set_names(rep("TM5", length(.)), .)

cluster6 <- dat@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% 
  subset(., UMAP_2 < 4.33 & UMAP_2 > -1.36 & UMAP_1 < -3.1 &
           (-2.11*UMAP_1 + UMAP_2 < 13.949) &
           (2*UMAP_1 + UMAP_2 > -8.76))  %>%
  rownames() %>%
  intersect(., cell_sample_matching %>%
              subset(., equals(., 4)) %>%
              names()) %>%
  set_names(rep("TM6", length(.)), .)

cluster7 <- dat@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  subset(., UMAP_2 < 6 & UMAP_1 < 5.7 & UMAP_2 > -3.8 &
           UMAP_1 > -2.4 &
           (4.13*UMAP_1 + UMAP_2 > -5.01) &
           (-1.09*UMAP_1 + UMAP_2 > -5.91)) %>%
  rownames() %>%
  intersect(., cell_sample_matching %>%
              subset(., equals(., 2)) %>%
              names()) %>%
  set_names(rep("TM7", length(.)), .)

cluster8 <- dat@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  subset(., UMAP_2 < 6 & UMAP_1 < 5.7 & UMAP_2 > -3.8 &
           UMAP_1 > -2.4 &
           (4.13*UMAP_1 + UMAP_2 > -5.01) &
           (-1.09*UMAP_1 + UMAP_2 > -5.91)) %>%
  rownames() %>%
  intersect(., cell_sample_matching %>%
              subset(., equals(., 3)) %>%
              names()) %>%
  set_names(rep("TM8", length(.)), .)

cluster9 <- dat@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  subset(., UMAP_2 < 6 & UMAP_1 < 5.7 & UMAP_2 > -3.8 &
           UMAP_1 > -2.4 &
           (4.13*UMAP_1 + UMAP_2 > -5.01) &
           (-1.09*UMAP_1 + UMAP_2 > -5.91)) %>%
  rownames() %>%
  intersect(., cell_sample_matching %>%
              subset(., equals(., 4)) %>%
              names()) %>%
  set_names(rep("TM9", length(.)), .)

cluster10 <- setdiff(names(cell_sample_matching), unique(c(cluster1, cluster2, 
                                                           cluster3, cluster4,
                                                           cluster5, cluster6,
                                                           cluster7, cluster8, 
                                                           cluster9) %>% names())) %>%
  set_names(rep("TM10", length(.)), .)


#Cells that are ignored
cluster_all <- c(cluster1,
                 cluster2,
                 cluster3,
                 cluster4,
                 cluster5,
                 cluster6,
                 cluster7, cluster8, cluster9, cluster10)


dat <- AddMetaData(object = dat, 
                   metadata = cluster_all, 
                   col.name = "manual_clustering_by_Tamer")

pdf("UMAP_showing_clusters_manually_annotated_by_TM.1234.pdf", width = 10, height =10)
DimPlot(dat, reduction = "umap", 
        group.by = "manual_clustering_by_Tamer")
dev.off()

#------------------------------------------
#Pairwise differential analysis.

diff_comparisons <- list(
  comp1 = c("TM1", "TM2"),
  comp2 = c("TM1", "TM3"),
  comp3 = c("TM3", "TM7"),
  comp4 = c("TM4", "TM5"),
  comp5 = c("TM4", "TM8"),
  comp6 = c("TM5", "TM8"),
  comp7 = c("TM4", "TM6"),
  comp8 = c("TM5", "TM6"),
  comp9 = c("TM6", "TM9"),
  comp10 = c("TM2", "TM7"),
  comp11 = c("TM2", "TM8"),
  comp12 = c("TM2", "TM9")
)

diff_results <- lapply(diff_comparisons, 
                       function(x) FindMarkers(dat, ident.1 = x[1], 
                                               ident.2 = x[2], 
                                               group.by = "manual_clustering_by_Tamer"))

diff_res_log2 <- lapply(diff_results, function(x) {
  x$avg_logFC <- log2(exp(1))*x$avg_logFC
  x
})


mapply(function(x, y) {
  write.table(x %>% 
                cbind(., gene = rownames(.)), file = paste0("Diff_results_", y[1], "_", y[2], ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}, diff_res_log2, diff_comparisons)

table(cluster_all) %>% 
  data.frame() %>% 
  set_colnames(., c("Cluster_label", "N_cells")) %>%
  write.table(., file = "Number_of_cells_in_clusters_labelled_by_Tamer.txt",
              quote = FALSE, sep = "\t", 
              row.names = FALSE)

#---------------------------------
#Compare cluster TM1 with all other LacZ cells.
red1 <- cluster1
red2 <- cell_sample_matching %>%
  subset(., equals(., "1"))
red2 <- red2[setdiff(names(red2), names(red1))]
red2 <- rep("LacZ_not_in_TM1", length(red2)) %>% 
  set_names(., names(red2))
others <- cell_sample_matching %>% 
  subset(., equals(., "1") %>% not())

cluster_red <- c(red1, red2, others)

dat <- AddMetaData(object = dat, 
                   metadata = cluster_red, 
                   col.name = "LacZ_cells_grouping_by_Tamer")

diff_results <- FindMarkers(dat, ident.1 = "TM1", 
                            ident.2 = "LacZ_not_in_TM1", 
                            group.by = "LacZ_cells_grouping_by_Tamer")
diff_results$avg_logFC <- log2(exp(1))*diff_results$avg_logFC
colnames(diff_results) <- c("p_val", 
                            "avg_log2FC_TM1_by_LacZ_not_in_TM1",
                            "pct.TM1",
                            "pct.LacZ_not_in_TM1", 
                            "p_val_adj")
write.table(diff_results %>% 
              cbind(., gene = rownames(.)), 
            file = "Diff_results_TM1_vs_LacZ_not_in_TM1.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)


#----------------------------
#Feature plots for surface marker proteins identified by TM (02-18-20)
surface_marker_genes <- c("CD36", "IGF2R", "JAG1", "FLNA", "ANK3", "EPHA4")

pdf("UMAP_showing_expression_of_surface_marker_genes_split_by_time_point.1234.pdf", width = 10, height =3.5)
map(surface_marker_genes, function(x) 
  FeaturePlot(dat, 
              features = x, 
              split.by = "Cell_sample_ident"))
dev.off()

pdf("UMAP_showing_expression_of_surface_marker_genes.1234.pdf", width = 10, height =10)
FeaturePlot(dat, 
            features = surface_marker_genes)
dev.off()

