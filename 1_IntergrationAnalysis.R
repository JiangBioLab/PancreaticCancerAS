rm(list = ls())

rm(list = ls())
set.seed(123)
setwd('/mywork/path')
system("type R") 
file.path(R.home("bin"), "R")

ps <- c('Seurat', 'ggplot2', 'stringr', 'reshape2','ggpubr','dplyr',
        "shiny","harmony",'tidyverse','ggpubr','ggplot2','parallel',
        'hdf5r','RColorBrewer','pillar')
#resolution check
sce = readRDS('sce.10XHarmony.rds')
# Check clustering stability at given resolution  
# Set different resolutions 
res.used <- seq(3,5,by=1)
res.used <- seq(5,15,by=2)

# Loop over and perform clustering of different resolutions 
sce <- FindNeighbors(sce, dims = 1:20,reduction = "harmony")
for(i in res.used){
  print(i)
  sce <- FindClusters(object = sce,
                      verbose = T, resolution = res.used)
  #sce <- FindNeighbors(sce, dims = 1:i)
}
# Make plot 
library(clustree)
clus.tree.out <- clustree(sce) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")

clus.tree.out
sce_clustree_raw= sce

#findCluster
ElbowPlot(sce,'harmony')
sce <- FindNeighbors(sce, dims = 1:20,reduction = 'harmony')
sce <- FindClusters(sce, reduction = "harmony",resolution = 4)
sce <- RunUMAP(object = sce, dims = 1:25, reduction = "harmony",do.fast = T )

DimPlot(sce,reduction = "umap",group.by = "seurat_clusters",label=T)

##Remove clusterers with cell content less than 500####
table(sce@meta.data$seurat_clusters)
#filter cluster(cells <500)
sce_change = subset(sce,seurat_clusters %in% c(0:34))
sce_change

sce_change <- RunUMAP(object = sce_change, dims = 1:25, reduction = "harmony" )
table(sce_change$RNA_snn_res.4)

#filter samples(cells<125)
#save(sce_change,file = 'sce.cluster34.res4Dim25.Rdata')
#load('sce.cluster34.res4Dim25.Rdata')


#Cell manual annotation####
Epithelial <-c('EPCAM','KRT8')
Immune <-c('CD3D',"CD3E",'IL7R') 
Monocyte <-c('CD14','CD163')
Dendritic <-c('GPR183')
Endothelial <-c( 'CDH5')
Fibroblastes <-c('ACTA2','COL1A1')
Mast <-c('TPSAB1', 'KIT')
Acinar <-c('CPA1')
BetaCell <-c('INS')
AlphaCell <-c('GCG')

new_cluster_id = c(
  'Epithelial' ,'Immune' , 'Monocyte' , 'Endothelial' ,
  'Fibroblastes' ,'Mast', 'Acinar','Beta-Cells','Alpha-Cells')
genes_to_check = c(Epithelial, Immune,
                   Monocyte, Endothelial,
                   Fibroblastes,Mast,Acinar  , AlphaCell, BetaCell)
#Rename clusters
celltype=data.frame(ClusterID=0:34,
                    celltype='undefined') 
celltype[celltype$ClusterID %in% c(6,8,10,12:15,18,19,24,26,28,29,33),2]='Epithelial1'
celltype[celltype$ClusterID %in% c(5),2]='Epithelial2'
celltype[celltype$ClusterID %in% c(33),2]='Beta'
celltype[celltype$ClusterID %in% c(28),2]='Alpha'
celltype[celltype$ClusterID %in% c(0,2,7,20,23,25),2]="Immune"
celltype[celltype$ClusterID %in% c(3,4,9,11,16,21,32),2]='Monocyte' 
celltype[celltype$ClusterID %in% c(27,30),2]="Immune"
celltype[celltype$ClusterID %in% c(30),2]='Endothelial'
celltype[celltype$ClusterID %in% c(17,22,31),2]='Fibroblastes'
celltype[celltype$ClusterID %in% c(1),2]='Mast'
celltype[celltype$ClusterID %in% c(34),2]='Acinar'

table(celltype$celltype=='undefined',celltype$ClusterID)
table(celltype$celltype)

cluster = sce_change@meta.data$seurat_clusters
sce_change@meta.data$selfCelltype = celltype[match(cluster,celltype$ClusterID),'celltype']
sce_change@meta.data$selfCelltype = as.factor(sce_change@meta.data$selfCelltype )

saveRDS(sce_change,file = 'scefinal.cluster34.res4.rds')#dim25
save(sce_change,file = 'scefinal.cluster34.res4.Rdata')

sce_change  = readRDS('scefinal.cluster34.res4.rds')
table(sce_change$selfCelltype)

#Remove samples with less than 10% epithelial content####
#Calculate the percentage of cells in a sample
orig.ident.info = as.data.frame(table(sce_change$orig.ident))
colnames(orig.ident.info) = c('Sample ID','Number of cells')
orig.ident.info$DataSource='GSE156405'
orig.ident.info$DataSource[grep('TISSUE',orig.ident.info$`Sample ID`)]='GSE155698'
orig.ident.info$DataSource[grep('SmartSeq2',orig.ident.info$`Sample ID`)]='GSE81547'
orig.ident.info$Technology = '10X'
orig.ident.info$Technology[grep('SmartSeq2',orig.ident.info$`Sample ID`)]='SmartSeq2'
table(orig.ident.info$DataSource,orig.ident.info$Technology)
write.csv(orig.ident.info,'sampleInfoFinal.csv')

df = data.frame(selfCelltype = sce_change$selfCelltype,
                patientID =sce_change$orig.ident)
dftbl  = as.data.frame(table(df))
perDf = dftbl %>%
  group_by(patientID) %>%
  mutate(percent = Freq/sum(Freq))
write.csv(perDf,'Celltype_percent.csv')

dftbl_epi1_Filter = subset(dftbl_epi1,percent>=0.1)
sce_change_patientFilt = subset(sce_change,orig.ident%in%dftbl_epi1_Filter$patientID )
table(sce_change_patientFilt$DataSource,sce_change_patientFilt$orig.ident)
table(sce_change_patientFilt$orig.ident)
saveRDS(sce_change_patientFilt,'scefinal.patientFilterd.rds')

#Extraction of epithelial and immune cells for CNV analysis####
rm(list = ls())
sce = readRDS('scefinal.patientFilterd.rds')
sce.epi = subset(sce,selfCelltype=='Epithelial1')
#saveRDS(sce.epi,'10xsmart_epi_harmony.rds')
sce.epi.norm = subset(sce.epi,cell_type =='Norm')
#saveRDS(sce.epi.norm,'10xsmart_epi_norm.rds')
sce.epi.cancer = subset(sce.epi,cell_type !='Norm')
sce.imm = subset(sce,selfCelltype =='Immune')
sce.imm.epi = merge(sce.epi.cancer,sce.imm)
GetAssay(sce.imm.epi,assay = "RNA")
sce.imm.epi <- FindVariableFeatures(sce.imm.epi, 
                                    selection.method = "vst", 
                                    nfeatures = 2000) 
sce.imm.epi <- ScaleData(sce.imm.epi)
sce.imm.epi <- RunPCA(object = sce.imm.epi, pc.genes = VariableFeatures(sce.imm.epi)) 
system.time({sce.imm.epi <- RunHarmony(sce.imm.epi,lambd=0.6,
                                       group.by.vars = c("orig.ident"), 
                                       max.iter.harmony=50,
                                       max.iter.cluster = 30)})
saveRDS(sce.imm.epi,'10xsmart_epiImm_harmony.rds')

###############################
sce = readRDS('sce.10XHarmony.rds')

pdf('../../Figures/Integration_before_after_umap.pdf',width = 10,height = 3)
after_integration = DimPlot(sce,reduction = "harmony",group.by = "DataSource",label=F,
                            pt.size=0.1,order = c('GSE81547','GSE156405','GSE155698'))+
  scale_color_manual(values = c('GSE156405'='#036EB8','GSE155698' = '#F0CD77','GSE81547'='#E26E42'))
before_integration = DimPlot(sce,reduction = "pca",group.by = "DataSource",label=F,
                             pt.size=0.1,order = c('GSE81547','GSE156405','GSE155698'))+
  scale_color_manual(values = c('GSE156405'='#036EB8','GSE155698' = '#F0CD77','GSE81547'='#E26E42'))
print(before_integration+after_integration)
dev.off()


sce_harmony <- FindNeighbors(sce, dims = 1:10,reduction = 'harmony')
sce_harmony <- FindClusters(sce_harmony, reduction = "harmony",resolution = 0.1)
sce_harmony <- RunUMAP(object = sce_harmony, dims = 1:4, reduction = "harmony",do.fast = T )

table(sce_harmony$seurat_clusters)
DimPlot(sce_harmony,reduction = "harmony", group.by= 'seurat_clusters', split.by= "DataSource",label=F,
        pt.size=0.1,order = c('GSE81547','GSE156405','GSE155698'))

DimPlot(sce_harmony,reduction = "harmony", split.by= 'seurat_clusters',label=F,pt.size=0.1)

table(sce_harmony$seurat_clusters,sce_harmony$DataSource)
sce_harmony$cluster_self =ifelse(sce_harmony$seurat_clusters %in% c(3,6),'cluster1','cluster2')

pdf('../../Figures/harmony_cluster_umap.pdf',width = 5,height = 3)
pclu = DimPlot(sce_harmony,reduction = "harmony", group.by= 'cluster_self',label=F,pt.size=0.1)+
  scale_color_manual(values = c('cluster1'='#EDB2B9','cluster2' = '#8FBD7A'))
print(pclu)
dev.off()

clusterInfo = sce_harmony@meta.data
colnames(clusterInfo)
table(clusterInfo$DataSource,clusterInfo$cluster_self)

pda = ggplot(clusterInfo, aes(cluster_self))+
  geom_bar(aes(fill=DataSource), position = 'fill')+
  scale_fill_manual(values = c('GSE156405'='#036EB8','GSE155698' = '#F0CD77','GSE81547'='#E26E42'))+
  theme_classic()

load('scefinal.cluster34.res4.Rdata')
sce_ant = sce_change@meta.data
sce_ant$cluster_harmony = sce_harmony$cluster_self[match(rownames(sce_ant),
                                                         names(sce_harmony$orig.ident))]

pbar = ggplot(sce_ant, aes(cluster_harmony))+
  geom_bar(aes(fill=selfCelltype), position = 'fill')+
  scale_fill_manual(values = mycol1)+
  theme_classic()

pdf('../../Figures/harmony_clusterInfo_umap.pdf',width = 6,height = 3)
print(pda+pbar)
dev.off()
