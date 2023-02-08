library(Seurat)
library(SeuratDisk)
library(scCustomize)
library(viridis)
library(patchwork)
library(enrichR)

all <- LoadH5Seurat(file=file.path("data","all.H5Seurat"))
pal <- viridis(n = 10, option = "D")
FeaturePlot_scCustom(all, "COL1A1")
all <- RenameIdents(all,"P2RY12hi" = "P2RY12high", "P2RY12lo" = "P2RY12low")

p2 <- DimPlot_scCustom(all, group.by = "species") + Move_Legend("bottom")
p3 <- DimPlot_scCustom(all, group.by = "sub.cluster", label = T) + NoLegend()
p1 <- DimPlot_scCustom(all) + NoLegend()

p1 / (p2 | p3)

#markers used
gene_list_plot <- c("COL1A1", "MRC1", "MSR1", "P2RY12", "FLT1", "PRF1", "PTPRC")
colors_list <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "darkorchid3", "orchid",
                                   "orange", "gold", "gray")

Stacked_VlnPlot(seurat_object = all, features = gene_list_plot, x_lab_rotate = TRUE,
                colors_use = colors_list)

DEenrichRPlot(all, ident.1 = "P2RY12high", ident.2 = "P2RY12low", max.genes = 200, enrich.database = "ChEA_2022")






