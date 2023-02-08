library(SingleCellExperiment)
library(scDblFinder)
library(tidyverse)
library(assertthat)
library(SeuratDisk)
library(Seurat)
library(biomaRt)
library(assertthat)
library(patchwork)
library(scCustomize)
# to do : cell cycle scoring
all <- RenameIdents(all,"P2RY12high" = "P2RY12high MΦ", "P2RY12low" = "P2RY12low MΦ")


human <- Read10X_h5(file.path("data","counts","XENO_human_post_xenocell","count","raw_feature_bc_matrix.h5"))[["Gene Expression"]][,read.csv(file.path("data","counts","XENO_human_post_xenocell","multiplexing_analysis","assignment_confidence_table.csv"))$Barcodes]
mouse <- Read10X_h5(file.path("data","counts","XENO_mouse_post_xenocell","count","raw_feature_bc_matrix.h5"))[["Gene Expression"]][,read.csv(file.path("data","counts","XENO_mouse_post_xenocell","multiplexing_analysis","assignment_confidence_table.csv"))$Barcodes] 
mouse_meta <- readRDS(file.path("data","Niklas_project","mouse.metadata.rds"))

## load orthologue mapping and rename the gene names of the mouse object
## found this from: https://support.bioconductor.org/p/129636/
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

mouse_genes <- mouse_human_genes %>% 
  filter(Common.Organism.Name=="mouse, laboratory") %>% 
  distinct(DB.Class.Key, Symbol)

human_genes <- mouse_human_genes %>% 
  filter(Common.Organism.Name=="human") %>% 
  distinct(DB.Class.Key, Symbol) %>% 
  rename(human_Symbol = Symbol)

orthologs <- mouse_genes %>% 
  inner_join(human_genes) %>% 
  distinct(Symbol, .keep_all = T) %>% 
  distinct(human_Symbol, .keep_all = T) %>% 
  as.data.frame()

rownames(orthologs) <- orthologs$Symbol

mouse <- mouse[rownames(mouse) %in% rownames((orthologs)),]
orthologs <- orthologs[rownames(mouse),]

## sanity check
assert_that(identical(rownames(mouse),rownames(orthologs)))

## rename
rownames(mouse) <- orthologs$human_Symbol

## create seurat object and merge both objects
mouse <- CreateSeuratObject(mouse) %>% 
  AddMetaData("mouse", col.name = "species") %>% 
  AddMetaData(mouse_meta$cell_types,col.name = "celltype")

sce_mouse <- SingleCellExperiment(assays=list(counts=mouse@assays$RNA@counts))
sce_mouse <- scDblFinder(sce_mouse)
mouse <- mouse[,sce_mouse@colData$scDblFinder.class == "singlet"]

## human data
human <- CreateSeuratObject(human) %>% 
  AddMetaData("human", col.name = "species") %>% 
  AddMetaData("tumor",col.name = "celltype")

sce_human <- SingleCellExperiment(assays=list(counts=human@assays$RNA@counts))
sce_human <- scDblFinder(sce_human)
human <- human[,sce_human@colData$scDblFinder.class == "singlet"]

## remove duplicate cells
mouse <- mouse[,!colnames(mouse) %in% colnames(human)]
human <- human[,!colnames(human) %in% colnames(mouse)]

## build seurat object ( res 0.5, trying 02.5)
all <- merge(human, mouse) %>% 
  subset(subset = nFeature_RNA > 400 & nFeature_RNA < 5000) %>% 
  SCTransform(variable.features.n = 10000) %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:30) %>% 
  FindNeighbors(dims=1:30) %>% 
  FindClusters(resolution = .25)

## recluster
all <- all %>% 
  FindSubCluster(cluster=c("7"),
                  resolution = .15,
                  graph.name = "SCT_snn")
#immune cells
all <- all %>% 
  FindSubCluster(cluster=c("3"),
                 resolution = .15,
                 subcluster.name = "sub.cluster",
                 graph.name = "SCT_snn")


## exploratory plots
DimPlot(all, label=T)
p1 <- DimPlot(all, label=T, group.by = "species") + NoLegend() + ggtitle("Species") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(all, label = T, group.by = "cluster_final") + NoLegend() + ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
p3 <- DimPlot(all)
ggsave(file.path("PlotsFinal","Featureplot.jpeg"))
DimPlot(all, label=T, group.by = "sub.cluster")
DimPlot(all, label=T, group.by = "sub.cluster2")
DimPlot(all, label = T, group.by = "celltype") + NoLegend()
FeaturePlot(all, features = features)
P2RY12 <- WhichCells(object = all, expression = P2RY12 > 1)
MRC1 <- WhichCells(object = all, expression = MRC1 > 1)
cells <- list(P2RY12 = P2RY12, MRC1=MRC1)
Cell_Highlight_Plot(seurat_object = all, cells_highlight = cells)

## consolidate cluster naming
all$cluster_final <- case_when(
  as.character(all$seurat_clusters) == "3" ~ as.character(all$sub.cluster),
  T ~ as.character(all$seurat_clusters)
)
Idents(all) <- all$cluster_final
## DimPlot
DimPlot(all, group.by = "cluster_final", label = T)

#SaveObject
saveRDS(all,file.path("data","allcells.rds"))

#Exploring marker genes
markers <- FindAllMarkers(all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- markers %>% 
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
genes <- as.list(rownames(markers))

immunecells <- subset(all, idents = c("P2RY12high MΦ","P2RY12low MΦ","NK"))
features = c("P2RY12","MRC1","COL1A1","ABCC9","MECOM","DCN","COL3A1")
DotPlot(all, features = features) + RotatedAxis()
DoHeatmap(subset(all, downsample = 100), features = features, size = 3)
DoHeatmap(immunecells, features = VariableFeatures(all)[1:30], cells = 1:500, size = 4,
          angle = 90) + NoLegend()

#Renaming Idents
all <- RenameIdents(all, c("Tumor1","Tumor2","Tumor3","Tumor4","Tumor5","Tumor6","Tumor7","Tumor8","Tumor9","Tumor10","Tumor11") == "Tumor")
all <- RenameIdents(all, c("Tumor1","Tumor2","Tumor3","Tumor4","Tumor5","Tumor6","Tumor7","Tumor1","Tumor1","Tumor1","Tumor1"))

#Fibroblasts COL1A1 COL3A1, LUM
all <- RenameIdents(all, "6" = "FB")

#Macrophages P2RY12 MRC1 MSR1
all <- RenameIdents(all, "3_1" = "P2RY12hi", "3_0" = "P2RY12lo")


#Endothelia cells FLT1
all <- RenameIdents(all, "9" = "EC")

#Tumor cells
all <- RenameIdents(all, "1" = "Tumor", "2" = "Tumor")
all <- RenameIdents(all, "4" = "Tumor", "8" = "Tumor", "7" = "Tumor", "5" = "Tumor")


#NK cells NCR1, PRF1
all <- RenameIdents(all, "3_2" = "NK")

features = c("MRC1","P2RY12","MSR1","COL1A1","PRF1","FLT1")
FeaturePlot_scCustom(all,features = features)
DotPlot(all, features = features) + RotatedAxis()

#Composition
#Immune cells pTPRC which encodes for CD45
immune <- FeaturePlot_scCustom(all, "PTPRC")
#mesenchymal cells  (COL1A1, COL1A2, LUM, DCN, ACTA2, RGS5)
mesen <- FeaturePlot_scCustom(all, "RGS5")

#Differential expression testing
p2ry12.de.markers <- FindMarkers(all, ident.1 = "P2RY12high MΦ", ident.2 = "P2RY12low MΦ", test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
DEgenes <- rownames(p2ry12.de.markers)

EnhancedVolcano(p2ry12.de.markers,lab = rownames(p2ry12.de.markers), x = "avg_log2FC", y = "p_val", pCutoff = 10e-6,
                title = "P2RY12 high vs low")

library(enrichR)
DEenrichRPlot(all, ident.1 = "P2RY12high MΦ", ident.2 = "P2RY12low MΦ", max.genes = 200, enrich.database = "GO_Biological_Process_2021")
  



## consolidating celltypes

table(subset.all$species)
## save object
SaveH5Seurat(all, file=file.path("data","all.H5Seurat"))
