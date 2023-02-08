library(NICHES)
library(Seurat)
library(EnhancedVolcano)

xenograft <- LoadH5Seurat(file.path("data","all.H5Seurat"))
xenograft@meta.data$celltype <- Idents(xenograft)
xenograft.niche <- RunNICHES(xenograft,
                        assay = 'SCT',
                        species = 'human',
                        LR.database = 'fantom5',
                        cell_types = "celltype",
                        SystemToCell = T,
                        CellToCell = F) 


# Here we use a previously-created imputed data slot. The same computation may be performed on any data slot including standard log-normalized RNA values.
niche.object <- xenograft.niche$SystemToCell
niche.object <- NormalizeData(niche.object)
niche.object <- ScaleData(niche.object) #Scale
niche.object <- FindVariableFeatures(niche.object,selection.method="disp") #Identify variable features
niche.object <- RunPCA(niche.object) #RunPCA
ElbowPlot(niche.object,ndims=100) #Choose PCs to use for embedding


niche.object <- RunUMAP(niche.object,dims = 1:40) #Embed
umap <- DimPlot(niche.object,reduction = 'umap', label = T)  + ggtitle('SystemToCell') + theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path("PlotsFinal","systemtocell_xeno.jpeg"))

#Markers
Idents(niche.object) <- niche.object[['ReceivingType']]
mark <- FindAllMarkers(niche.object,logfc.threshold = 0.1,min.pct = 0.1,only.pos = F,test.use = 'roc')
writ
# Pull markers of interest to plot
mark$ratio <- mark$pct.1/mark$pct.2
marker.list <- mark %>% group_by(cluster) %>% top_n(20,avg_log2FC)
write.table(mark, file.path("data","de_nichep2.tsv"), sep='\t')

#Plot in Heatmap form
heatmap <- DoHeatmap(niche.object,features = marker.list$gene,cells = WhichCells(niche.object,downsample = 100))
umap + heatmap
## differential markers of the macrophages
diffgenes <- FindMarkers(niche.object, ident.1 = "P2RY12high", ident.2 ="P2RY12low")
diffgenes %>% 
  EnhancedVolcano(lab = rownames(.),
                  x = 'avg_log2FC',
                  y = 'p_val',
                  title = 'System to cell interactions P2RY12 high vs low',
                  pCutoff = 10e-6,
                  FCcutoff = 0.5,
                  pointSize = 3.0,
                  labSize = 6.0)
ggsave(file.path("PlotsFinal","p2DE_xeno.jpeg"))

#cell cell interactions of immune niche
immune.sub <- subset(xenograft, idents = c("GC","P2RY12low","P2RY12high"))

immune.cc <- RunNICHES(immune.sub,
                     assay = 'SCT',
                     species = 'human',
                     LR.database = 'fantom5',
                     cell_types = "celltype",
                     CellToCell = T)
cc.object <- immune.cc$CellToCell
cc.object <- cc.object %>% ScaleData() %>% FindVariableFeatures(method = "disp") %>% RunPCA()
ElbowPlot(cc.object,ndims=100) #Choose PCs to use for embedding
cc.object <- RunUMAP(cc.object,dims = 1:50)
cc.object <- FindNeighbors(cc.object,dims=1:50)
cc.object <- FindClusters(cc.object,resolution = 0.5)
#Embed
DimPlot(cc.object,reduction = "umap", label = F) + NoAxes()+ ggtitle('Xenograft Cell-Cell Signaling')+
  guides(colour=guide_legend(ncol=2,override.aes = list(size=6)))
Idents(cc.object) <- cc.object[['ReceivingType']] 
marker.immun <- FindAllMarkers(cc.object,
                               logfc.threshold = 0.1,
                               min.pct = 0.1,
                               only.pos = T,
                               test.use = 'roc')

return_celltypes(all)


