#Integrated object
library(SeuratDisk)
library(Seurat)
library(Signac)
meningioma <- LoadH5Seurat(file.path("data","integrated_object.H5Seurat"))

SaveH5Seurat(meningioma,file.path("data","multi_meningioma.H5Seurat"))
DefaultAssay(meningioma) <- "RNA"
meningioma <- SCTransform(meningioma, verbose = FALSE, variable.features.n = 10000) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


DefaultAssay(meningioma) <- "ATAC"
meningioma <- RunTFIDF(meningioma)
meningioma <- FindTopFeatures(meningioma, min.cutoff = 'q0')
meningioma <- RunSVD(meningioma)
meningioma <- RunUMAP(meningioma, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DimPlot(meningioma, reduction = "umap.atac")

meningioma <- FindMultiModalNeighbors(meningioma, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
meningioma <- RunUMAP(meningioma, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
meningioma <- FindClusters(meningioma, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
DimPlot(meningioma,reduction = "wnn.umap")


p1 <- DimPlot(meningioma, reduction = "umap.rna", label = TRUE, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(meningioma, reduction = "umap.atac", label = TRUE, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(meningioma, reduction = "wnn.umap", label = TRUE, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

DefaultAssay(meningioma) <- "SCT"
FeaturePlot(meningioma, reduction = "wnn.umap", label = TRUE, "PTPRC") 
meningioma <- FindSubCluster(meningioma, cluster = "4_0", graph.name = "wsnn", algorithm = 3, resolution = 1)
Idents(meningioma) <- "sub.cluster"
FeaturePlot(meningioma,"ITGAM",label = T, reduction = "wnn.umap")


#ImmuneCells PTPRC SPI1

FeaturePlot(meningioma, reduction = "wnn.umap", label = TRUE, "PTPRC")
FeaturePlot(meningioma, reduction = "wnn.umap", label = TRUE, "P2RY12")

#meningioma <- FindSubCluster(meningioma, cluster = 5, graph.name = "wsnn", algorithm = 3, resolution = 1)
Idents(meningioma) <- "sub.cluster"

#Macrophages MRC1 MSR1 P2RY12
FeaturePlot(meningioma, reduction = "wnn.umap", label = TRUE, c("MRC1","MSR1","P2RY12"))
FeaturePlot(meningioma, reduction = "wnn.umap", label = TRUE, "P2RY12")
FeaturePlot(meningioma, reduction = "wnn.umap", label = TRUE, "PLVAP")


meningioma <- RenameIdents(meningioma, "4_0_0" = "P2RY12high","4_1" = "P2RY12high","4_0_1" = "P2RY12low")




#Endothelia FLT1
FeaturePlot(meningioma, reduction = "wnn.umap", label = TRUE, "FLT1")
FeaturePlot(meningioma, reduction = "wnn.umap", label = TRUE, "PECAM1")



#mesenchymal (Fibroblasts)
FeaturePlot(meningioma, reduction = "wnn.umap", label = TRUE, c("COL1A1","COL1A2","COL3A1","DCN"))
meningioma <- RenameIdents(meningioma, "1" = "FB","8" = "FB","3" = "FB")
levels(meningioma)
meningioma$celltype <- Idents(meningioma)

markers <- FindAllMarkers(meningioma)
markers
p2.markers.multi <- FindMarkers(meningioma, ident.1 = "P2RY12high", ident.2 = "P2RY12low", min.pct = 0.1, logfc.threshold = 0.1, only.pos = F)
write.table(p2ry12.de, file.path("data","de_highvslow.tsv"), sep='\t')


EnhancedVolcano(p2.markers.multi,lab = rownames(p2.markers.multi), x = "avg_log2FC", y = "p_val", pCutoff = 10e-6,
                title = "P2RY12 high vs low Multiome")
ggsave(file.path("PlotsFinal","DEhighlow_multi.jpeg"))

DotPlot(meningioma,features = features)




