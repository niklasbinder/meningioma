xenograft.markers <- FindAllMarkers(xenograft, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
xenograft.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

##################################################################
#Subsetting top 100 markers with adjusted p values lower than .05#
##################################################################
top100 <- xenograft.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)

library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)

df <- top100pval[,7:6]
dfsample <- split(df$gene,df$cluster)
length(dfsample)
dfsample$P2RY12high = bitr(dfsample$P2RY12high, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$P2RY12low = bitr(dfsample$P2RY12low, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$GC = bitr(dfsample$GC, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$Tumor = bitr(dfsample$Tumor, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$EC= bitr(dfsample$EC, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$FB= bitr(dfsample$FB, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
genelist <- list("P2RY12high" = dfsample$P2RY12high$ENTREZID, 
                 "P2RY12low" = dfsample$P2RY12low$ENTREZID,
                 "GC" = dfsample$GC$ENTREZID,
                 "Tumor" = dfsample$Tumor$ENTREZID,
                 "EC" = dfsample$EC$ENTREZID,
                 "FB" = dfsample$FB$ENTREZID)
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
dotplot(GOclusterplot) + ggtitle('Gene Ontology Clusterplot') + theme(plot.title = element_text(hjust = 0.5))
KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
dotplot(KEGGclusterplot) + ggtitle('KEGG Clusterplot') + theme(plot.title = element_text(hjust = 0.5))
Diseaseclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDO")
dotplot(Diseaseclusterplot) + ggtitle('Disease Ontology Clusterplot') + theme(plot.title = element_text(hjust = 0.5))
