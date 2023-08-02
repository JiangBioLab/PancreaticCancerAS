rm(list = ls())
library(SeuratData)
library(SeuratObject)
library(Seurat)
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape)
library(harmony)
library(shiny)

setwd('/home/lxxiao/xiaolixing/pdac/10x/')
#10x datasets
setwd('/home/lxxiao/xiaolixing/pdac/10x/')
system("type R") 
file.path(R.home("bin"), "R")
#multiprocess####
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 1000000 * 1024^2)


if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)

if(length(getOption("BioC_mirror"))==0) options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

options(stringsAsFactors = F)
library(Seurat)

grep("filtered_feature_bc_matrix",list.files('GSE156405_pri_met/'),value = T)
fs_1=grep("filtered_feature_bc_matrix",list.files('GSE156405_pri_met/'),value = T)[c(1:7,19)]
samples_1=str_split(fs_1,'_',simplify = T)[,1]
data_path1 = paste0('GSE156405_pri_met/',fs_1)

grep("_TISSUE_",list.files('No9_hd/'),value = T)
fs_2=grep("_TISSUE_",list.files('No9_hd/'),value = T)
fs_2 = fs_2[-10]#-c("PDAC_TISSUE_14")
samples_2=fs_2
data_path2 = paste0('No9_hd/',fs_2,'/filtered_feature_bc_matrix')

fs = c(fs_1,fs_2)
samples = c(samples_1,samples_2)
data_path = c(data_path1,data_path2)

library(stringr)
folder = './No9_hd/AdjNorm_TISSUE_1/filtered_feature_bc_matrix'
pro = 'PDAC_TISSUE_1'

if(F){ 
  
  library(Seurat)
  sceList_1 = lapply(samples_1,function(pro){ 
    print(pro)
    folder=file.path('./',grep(pro,data_path,value = T)) 
    
    CreateSeuratObject(counts = Read10X(folder), 
                       project = pro )
  })
pro = 'PDAC_TISSUE_16'  
  sceList_2 = lapply(samples_2,function(pro){ 
    print(pro)
    folder=file.path('./',grep(paste0(pro,'/'),data_path,value = T)) 
    
    CreateSeuratObject(counts = Read10X(folder), 
                       project = pro )
  })
  sceList = c(sceList_1,sceList_2)
  sce.big <- merge(sceList[[1]], 
                   y = c(sceList[[2]],
                         sceList[[3]],
                         sceList[[4]],
                         sceList[[5]],
                         sceList[[6]],
                         sceList[[7]],
                         sceList[[8]],
                         sceList[[9]],
                         sceList[[10]],
                         sceList[[11]],
                         sceList[[12]],
                         sceList[[13]],
                         sceList[[14]],
                         sceList[[15]],
                         sceList[[16]],
                         sceList[[17]],
                         sceList[[18]],
                         sceList[[19]],
                         sceList[[20]],
                         sceList[[21]],
                         sceList[[22]],
                         sceList[[23]],
                         sceList[[24]],
                         sceList[[25]],
                         sceList[[26]],
                         sceList[[27]]), 
                   add.cell.ids = samples, 
                   project = "PDAC")
  sce.big
  table(sce.big$orig.ident)
  save(sce.big,file = 'sce.big.merge.dataset.Rdata')
}

load('sce.big.merge.dataset.Rdata')
table(sce.big$orig.ident)
raw_sce = sce.big.Filt
raw_sce[["percent.mt"]] <- PercentageFeatureSet(raw_sce, pattern = "^MT-")
fivenum(raw_sce[["percent.mt"]][,1])
raw_sce[["percent.ribo"]] <- PercentageFeatureSet(raw_sce, pattern = "^RP[SL]")
CombinePlots(plots = list(plot1, plot2))

pro='merge_all'

raw_sce1 <- subset(raw_sce, subset = nFeature_RNA > 1500 & nCount_RNA > 2000 & percent.mt < 20)
raw_sce1
selected_f <- rownames(raw_sce1)[Matrix::rowSums(raw_sce1@assays$RNA@counts > 0 ) > 30]
raw_sce1 <- subset(raw_sce1, features = selected_f)
raw_sce1
sce=raw_sce1
#sce=raw_sce1[, colnames(raw_sce1)%in% phe$cell]
sce
sce <- NormalizeData(sce, normalization.method =  "LogNormalize", 
                     scale.factor = 10000)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000) 
VariableFeaturePlot(sce)
sce <- ScaleData(sce)
scaled_sce = sce

sce_runfinished = sce
sce = sce_runfinished

sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 
PCA_sce = sce
DimPlot(sce,reduction = 'pca',group.by = 'orig.ident',pt.size = 0.1)
###hrmony####
library(harmony)
library(parallel)
system.time({sce <- RunHarmony(sce,lambd=0.6,
                                    group.by.vars = c("orig.ident"), 
                                    max.iter.harmony=50,
                                    max.iter.cluster = 30)})
DimPlot(sce,reduction = 'harmony',group.by = 'orig.ident',pt.size = 0.1)
saveRDS(sce,"all_raw_harmony_PDAC.rds")

sce = readRDS('all_raw_harmony_PDAC.rds')
# Check clustering stability at given resolution  
# Set different resolutions 
res.used <- seq(3,5,by=1)
res.used
res.used <- seq(5,15,by=2)
res.used
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
#RUN TSNE UMAP######################################

#DimHeatmap(sce, dims = 1:20, cells = 100, balanced = TRUE)
ElbowPlot(sce,'harmony')
sce <- RunTSNE(object = sce, dims = 1:10, reduction = "harmony",do.fast = T )
sce <- RunUMAP(object = sce, dims = 1:13, reduction = "harmony",do.fast = T )

sce <- FindNeighbors(sce, dims = 1:20,reduction = 'harmony')
sce <- FindClusters(sce, reduction = "harmony",resolution = 4)
table(sce@meta.data$RNA_snn_res.3)
DimPlot(sce,reduction = 'umap',group.by = 'selfCelltype')
DimPlot(sce,reduction = 'umap',label = T)

saveRDS(sce,"all_harmony_PDAC.rds")

table(sce@meta.data$RNA_snn_res.4) 
# tsne####
set.seed(123)
DimPlot(sce,reduction = "tsne",group.by = "seurat_clusters",label=T)
DimPlot(sce,reduction = "tsne",group.by = "RNA_snn_res.0.1",label=T)
sce@meta.data$orig.ident
DimPlot(sce,reduction = "tsne",group.by = 'orig.ident',label=T)

DimPlot(sce,reduction = "umap",group.by = "seurat_clusters",label=T)
DimPlot(sce,reduction = "umap",group.by = "RNA_snn_res.4",label=T)
DimPlot(sce,reduction = "umap",group.by = "orig.ident",label=T)
saveRDS(sce,"all_RunTSNEUMAP_PDAC.rds")

##Filt cluster##################################
table(sce@meta.data$seurat_clusters)
sce_change = sce[, sce@meta.data$seurat_clusters %in% c(0:65)]
sce_change


table(sce_change@meta.data$seurat_clusters,sce_change@meta.data$orig.ident)

phe=data.frame(cell=rownames(sce_change@meta.data),
               cluster =sce_change@meta.data$seurat_clusters)
dim(phe)
head(phe)
table(phe$cluster)
#save TSNE UMAP position#################################

sce = sce_change
table(sce@meta.data$seurat_clusters,sce@meta.data$orig.ident)

phe=data.frame(cell=rownames(sce@meta.data),
               cluster =sce@meta.data$seurat_clusters)
dim(phe)
head(phe)
table(phe$cluster)

sce_change_befor=sce
sce = sce_change
#sce = sce[, sce@meta.data$seurat_clusters %in% c(0:36)]
#saveRDS(sce,'sce_cluster_filted.rds')
tsne_pos=Embeddings(sce,'tsne')
phe=data.frame(cell=rownames(sce@meta.data),
               cluster =sce@meta.data$seurat_clusters)
DimPlot(sce,reduction = "tsne",label=T,split.by ='orig.ident')

head(tsne_pos) 
dat=cbind(tsne_pos,phe)
save(dat,file=paste0(pro,'_for_tSNE.pos.Rdata')) 

tsne_pos=Embeddings(sce,'umap')
DimPlot(sce,reduction = "umap",label=T,split.by ='orig.ident')

head(tsne_pos) 
dat=cbind(tsne_pos,phe)
save(dat,file=paste0(pro,'_for_UMAP.pos.Rdata')) 

saveRDS(sce,"all_RunTSNEUMAP_PDAC_clusterFilt.rds")

#!/usr/local/bin/R --vanilla
setwd('/home/lxxiao/xiaolixing/pdac/10x/')

.libPaths(c("/home/lxxiao/R/x86_64-pc-linux-gnu-library/4.0","/usr/local/lib64/R/library"))
.libPaths()
#library packages
library(parallel)
library(future)
library(stringr)
library(Seurat)
library(ggplot2)

clu = clusters[10]
clusters_marker = function(clu){
  library(parallel)
  library(stringr)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  sce<-readRDS('all_harmony_PDAC.rds')
  sce_file<-sce
  markers_df <- FindMarkers(object = sce_file, ident.1 = clu, min.pct = 0.25)
  print(x = head(markers_df))
  markers_genes =  rownames(head(x = markers_df, n = 5))
  DotPlot(object = sce_file, features =markers_genes )+coord_flip()
  ggsave(filename=paste0('figures/DotPlot_subcluster_',clu,'_sce.markers_heatmap.pdf'))
  markers_df$cluster = paste0('cluster',clu)
  print(head(markers_df))
  markers_df2 = top_n(markers_df,5, avg_logFC)
  return(markers_df2)
  #FeaturePlot(object = sce, features=markers_genes )
  #ggsave(filename=paste0(pro,'_FeaturePlot_subcluster_',i,'_sce.markers_heatmap.pdf'))
}
#results <- lapply(clusters,clusters_marker)

cl <- makeCluster(5)
clusters = unique(sce@meta.data$seurat_clusters)
clusters2 = clusters[c(5,6)]
results <- parLapply(cl,clusters2,clusters_marker)
saveRDS(results,'result/cluster_top5_marker.rds')

res.df <- do.call('rbind',results) 
stopCluster(cl)
write.csv(res.df,'result/cluster_top5_marker.csv')



sample_data  = str_split(sce$orig.ident,'_',simplify = T)[,1]
sce@meta.data$tissue = ifelse(sample_data=='AdjNorm','AdjNorm',
                              ifelse(str_detect(sample_data,'^P'),'locPDAC',
                                     'metPDAC'))
table(sce@meta.data$tissue)
saveRDS(sce,'all_RunTSNEUMAP_PDAC.rds')
sce = readRDS('all_RunTSNEUMAP_PDAC.rds')
sce.all.list.self <- SplitObject(sce , split.by = "selfCelltype")


#normal cell####
table(sce.epi@meta.data$patient_id)
sce.epi.list = SplitObject(sce.epi.t.cnv,split.by = 'selfCelltype')
sce.epi.cnvonly = subset(sce.epi.t.cnv,selfCelltype=='Epithelial')
table(sce.epi.cnvonly@meta.data$copykat.pred)
saveRDS(sce.epi.cnvonly,'result/Epithelial_only_CNV.rds')

colnames(sce@meta.data)
sce.all = readRDS('all_RunTSNEUMAP_PDAC.rds')
sce.all.list <- SplitObject(sce.all , split.by = "patient_id")
sce.all.list 
names(sce.all.list)
i="LiM"
for (i in names(sce.all.list)) {
  print(i)
  epi_mat = sce.all.list[[i]]@assays$RNA@counts
  #epi_phe = sce.all.list[[i]]@meta.data
  exp.rawdata <- as.matrix(epi_mat)
  norm_cell = colnames(sce.all.list[[i]])[which(sce.all.list[[i]]@meta.data$selfCelltype=='T_lymphocytes')]
  epi_cell  = colnames(sce.all.list[[i]])[which(sce.all.list[[i]]@meta.data$selfCelltype=='Epithelial')]
  exp.rawdata.epi = exp.rawdata[,epi_cell]
  norm_cell2 = intersect(norm_cell,epi_cell)
  sce_r=CreateSeuratObject(counts = epi_mat, 
                         meta.data = epi_phe )
 
  write.table(exp.rawdata,file = paste0('patient/',i,'_counts.txt'),sep = '\t',
              row.names = T,col.names = T,quote = F)
  write.csv(norm_cell,file = paste0('patient/',i,'_norm_name.csv'))
  saveRDS(sce_r,file = paste0('patient/',i,'.rds'))
}
exp.rawdata = read.table(file = paste0('patient/',i,'_counts.txt'),
                         header = T,row.names = 1,sep = '\t')
norm_cell =as.character(read.csv(paste0('patient/',i,'_norm_name.csv'),header = T,row.names = 1)[,1])

sce.big.Filt <- merge(sce.all.list[[1]], 
                 y = c(sce.all.list[[2]],
                       #sce.all.list[[3]],
                       sce.all.list[[4]],
                       sce.all.list[[5]],
                       sce.all.list[[6]],
                       sce.all.list[[7]],
                       sce.all.list[[8]],
                       sce.all.list[[9]],
                       sce.all.list[[10]],
                       sce.all.list[[11]],
                       sce.all.list[[12]],
                       sce.all.list[[13]],
                       #sce.all.list[[14]],
                       sce.all.list[[15]],
                       sce.all.list[[16]],
                       sce.all.list[[17]],
                       sce.all.list[[18]],
                       sce.all.list[[19]],
                       sce.all.list[[20]],
                       sce.all.list[[21]],
                       sce.all.list[[22]],
                       sce.all.list[[23]],
                       sce.all.list[[24]],
                       sce.all.list[[25]],
                       sce.all.list[[26]],
                       sce.all.list[[27]]), 
                 add.cell.ids = patient_id, 
                 project = "PDAC")
sce_all
sce.big.Filt
saveRDS(sce.big.Filt,'patient_filt.rds')

#smart datasets
smart_h = read.csv('No10_health_smart/all_count_matrix.csv',header = T,row.names = 1)
smart_phe = read.csv('No10_health_smart/SraRunTable.txt',header = T,row.names = 1)
colnames(smart_h)[1:10]
gsmid =str_split(colnames(smart_h),'[_]',simplify = T)[,1]
all(gsmid==smart_phe$Sample.Name)
celltype = smart_phe$inferred_cell_type
table(celltype)
smart_pro =paste0('smart_',celltype)
smart_sce = seurat.object=CreateSeuratObject(count = smart_h,project = smart_pro,
                                             min.cell = 10,min.gene = 500)
smart_sce@meta.data[1:3,]
table(smart_sce$orig.ident)
x10_sce = readRDS('patient_filt.rds')
table(x10_sce@meta.data$orig.ident)

merge_sce = merge(smart_sce,x10_sce)

saveRDS(merge_sce,'patient_filt_10xsmart_merge.rds')#34590
raw_sce = merge_sce
raw_sce@meta.data$dataset_source = ifelse(str_split(raw_sce$orig.ident,'[_]',simplify = T)[,1]=='smart',
                                          'smart_Seq',
                                          '10x')

raw_sce[["percent.mt"]] <- PercentageFeatureSet(raw_sce, pattern = "^MT-")
fivenum(raw_sce[["percent.mt"]][,1])
raw_sce[["percent.ribo"]] <- PercentageFeatureSet(raw_sce, pattern = "^RP[SL]")
plot1 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

raw_sce@meta.data[1:4,]
raw_sce1 <- subset(raw_sce, subset = nFeature_RNA > 1000 & nCount_RNA > 2000 & percent.mt < 20)
raw_sce1
selected_f <- rownames(raw_sce1)[Matrix::rowSums(raw_sce1@assays$RNA@counts > 0 ) > 30]
raw_sce1 <- subset(raw_sce1, features = selected_f)
raw_sce1
table(raw_sce1@meta.data$dataset_source)

sce=raw_sce1
sce <- NormalizeData(sce, normalization.method =  "LogNormalize", 
                     scale.factor = 10000)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000) 
sce <- ScaleData(sce)
scaled_sce = sce
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 
PCA_sce = sce
table(sce@meta.data$orig.ident)

table(sce$orig.ident)
sce@meta.data$cell_from = sce@meta.data$orig.ident
sce@meta.data$cell_from[grep('smart',sce@meta.data$cell_from)]='Smart_Norm'
table(sce@meta.data$cell_from)
sce@meta.data$cell_from[grep('PDAC_TISSUE',sce@meta.data$cell_from)]='GSE155698_locPDAC'
sce@meta.data$cell_from[grep('^P',sce@meta.data$cell_from)]='GSE156405_locPDAC'
sce@meta.data$cell_from[grep('AdjNorm',sce@meta.data$cell_from)]='GSE155698_Norm'
sce@meta.data$cell_from[sce@meta.data$cell_from %in% c('LiM','VM','LuM')]='GSE156405_metPDAC'
table(sce@meta.data$cell_from)
sce@meta.data$cell_from2 = str_split(sce@meta.data$cell_from,'[_]',simplify = T)[,2]
table(sce@meta.data$cell_from2)
sce@meta.data$cell_from3 = str_split(sce@meta.data$cell_from,'[_]',simplify = T)[,1]
table(sce@meta.data$cell_from3)
###hrmony####
library(harmony)
library(parallel)
system.time({sce <- RunHarmony(sce,lambd=0.3,
                               group.by.vars = c("orig.ident"), 
                               max.iter.harmony=50,
                               max.iter.cluster = 50)})
table(sce$dataset_source)
DimPlot(sce,reduction = 'harmony',split.by = 'cell_from3',group.by = 'cell_from2',
        label = F,label.size = 5,
        label.box = F,pt.size = 0.1,repel=F)+
  scale_color_simpsons()
#ggsave('figures/figures1_renew2/dataset3_harmony_result.pdf')

saveRDS(sce,'10xsmart_merge_harmony.rds')
#RUN TSNE UMAP######################################
sce = readRDS('10xsmart_merge_harmony.rds')
#DimHeatmap(sce, dims = 1:20, cells = 100, balanced = TRUE)
ElbowPlot(sce,'harmony')
sce <- FindNeighbors(sce, dims = 1:20,reduction = 'harmony')
sce <- FindClusters(sce, reduction = "harmony",resolution = 4)

sce <- RunUMAP(object = sce, dims = 1:29, reduction = "harmony",do.fast = T )
#sce@meta.data$selfCelltype[is.na(sce@meta.data$selfCelltype)]='smart'
DimPlot(sce,reduction = 'umap',group.by = 'selfCelltype')+
  scale_color_manual(values = c('smart'='red'))
DimPlot(sce,reduction = 'umap',label = T)

##Filt cluster##################################
table(sce@meta.data$seurat_clusters)
sce_change = sce[, sce@meta.data$seurat_clusters %in% c(0:59)]
sce_change
table(sce_change@meta.data$seurat_clusters,sce_change@meta.data$orig.ident)

phe=data.frame(cell=rownames(sce_change@meta.data),
               cluster =sce_change@meta.data$seurat_clusters)
dim(phe)
head(phe)
table(phe$cluster)

#save TSNE UMAP position#################################
sce = sce_change
table(sce@meta.data$seurat_clusters,sce@meta.data$orig.ident)

phe=data.frame(cell=rownames(sce@meta.data),
               cluster =sce@meta.data$seurat_clusters)
sce_change_befor=sce
sce = sce_change
#markers####
#cell manual annotation####
Epithelial       <-c('EPCAM','KRT8')
#Epithelial2       <-c('EPCAM','KRT8','UBE2C')
T_lymphocytes    <-c('CD3D',"CD3E",'IL7R') #'PTPRC' also called 'CD45',
Monocyte           <-c('CD14','CD163')#,'LYZ','CD68'
NK                <-c( "NKG7", "GNLY")#'FCGR3A'
B_lymphocytes     <-c("CD79A","MZB1")#"MS4A1"
Dendritic         <-c('GPR183')#"TSPAN13",
Endothelial       <-c( 'CDH5')#"PECAM1",
Fibroblastes      <-c('ACTA2','COL1A1')
#Monocyte          <-c( 'CD163')#'CD68',
Mast              <-c('TPSAB1', 'KIT')
Acinar_cells      <-c('CPA1')


new_cluster_id = c(
  'Epithelial'    ,
  #'Epithelial2'    ,
  'T_lymphocytes'   ,
  'Monocyte'          ,
  'NK'               ,
  'B_lymphocytes'    ,
  #'Dendritic'       ,
  'Endothelial'      ,
  'Fibroblastes'     ,
  'Mast_cells',
  'Acinar_cells'
)


#VlnPlot(sce,features = Fibroblasts,pt.size = 0,ncol = 2)
# Specify genes  
genes_to_check = c(Epithelial    ,
                   #Epithelial2    ,
                   T_lymphocytes  ,
                   NK             ,
                   B_lymphocytes  ,
                   Monocyte       ,
                   #Dendritic      ,
                   Endothelial    ,
                   Fibroblastes   ,
                   #Monocyte      ,
                   Mast           ,
                   Acinar_cells   )
# All on Dotplot 
#ggthemr(palette = "dust", set_theme = FALSE)

library(ggplot2)

p_all_markers <- DotPlot(sce, features = unique(genes_to_check))+ coord_flip()

p_umap=DimPlot(sce, reduction = "umap",
               group.by = "RNA_snn_res.4",label = T) 
#detach("package:ggthemr")
library(patchwork)
p_all_markers+p_umap
##rename clusters####
celltype=data.frame(ClusterID=0:59,
                    celltype=0:59) 
celltype[celltype$ClusterID %in% c(11:15,17,20:22,30,32,33,35,57,
                                   39,40,45,54,55),2]='Epithelial'
celltype[celltype$ClusterID %in% c(0,1,3:5,23,29,37,38,41,44,46:50,52),2]='immune'#'T_lymphocytes' 

celltype[celltype$ClusterID %in% c(1,2,7,8,9,10,18,19,43,56,58),2]='Monocyte' #MYeloid
celltype[celltype$ClusterID %in% c(),2]='immune'#'NK' 
celltype[celltype$ClusterID %in% c(36),2]='immune'#'B_lymphocytes' 
#celltype[celltype$ClusterID %in% c(24,27),2]='Dendritic' 
celltype[celltype$ClusterID %in% c(28),2]='Endothelial'
celltype[celltype$ClusterID %in% c(16,25,26,31),2]='Fibroblastes'
celltype[celltype$ClusterID %in% c(6,34),2]='Mast_cells'
#celltype[celltype$ClusterID %in% c(10),2]='Monocyte'
celltype[celltype$ClusterID %in% c(27,42,51,53,59),2]='Acinar_cells'
celltype[celltype$ClusterID %in% c(24),2]='undefined'

table(celltype$celltype)
head(celltype)
cluster = sce@meta.data$seurat_clusters
sce@meta.data$selfCelltype = celltype[match(cluster,celltype$ClusterID),'celltype']
sce <- RunUMAP(object = sce, dims = 1:15, reduction = "harmony",do.fast = T )
DimPlot(sce, reduction = "umap", group.by = "selfCelltype",label = F)  
DimPlot(sce, reduction = "umap", group.by = "dataset_source",label = T)  
DimPlot(sce, reduction = "umap", label = T)+scale_color_manual(values = c('44'='red'))  

DotPlot(sce, group.by="selfCelltype",
        features = unique(genes_to_check))
saveRDS(sce,'10xsmart_merge_harmony.rds')
sce = readRDS('10xsmart_merge_harmony.rds')
colnames(sce@meta.data)
#get epi####
sce.epi = subset(sce,selfCelltype=='Epithelial')
system.time({sce.epi <- RunHarmony(sce.epi,lambd=0.6,
                               group.by.vars = c("orig.ident"), 
                               max.iter.harmony=50,
                               max.iter.cluster = 30)})
DimPlot(sce.epi,reduction = 'harmony',group.by = 'cell_from')
sce.epi <- FindNeighbors(sce.epi, dims = 1:20,reduction = 'harmony')
sce.epi <- FindClusters(sce.epi, reduction = "harmony",resolution = 1)
sce.epi <- RunUMAP(object = sce.epi, dims = 1:20, reduction = "harmony",do.fast = T )
DimPlot(sce.epi,reduction = 'umap',label =F,group.by = 'cell_from2',
        label.size = 2,
        label.box = F,
        pt.size = 0.1,
        repel=F)+scale_color_manual(values = c('locPDAC' = '#006834',
                                                'metPDAC'='#036EB8',
                                              'Norm' = '#EB6000'))

saveRDS(sce.epi,'10xsmart_epi_harmony.rds')

sce.epi = readRDS('10xsmart_epi_harmony.rds')
#get epi norm/cancer####
table(sce.epi$cell_from2)
sce.epi.norm = subset(sce.epi,cell_from2 =='Norm')
colnames(sce.epi.norm@meta.data)
saveRDS(sce.epi.norm,'10xsmart_epi_norm.rds')

sce.epi.cancer = subset(sce.epi,cell_from2 !='Norm')
sce.imm = subset(sce,selfCelltype =='immune')

#get epi imm####
table(sce.epi.cancer@meta.data$orig.ident)

sce.imm.epi = merge(sce.epi.cancer,
                    sce.imm)
table(sce.imm.epi@meta.data$selfCelltype)
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
table(sce.imm.epi$orig.ident,
      sce.imm.epi$selfCelltype)
saveRDS(sce.imm.epi,'10xsmart_epiImm_harmony.rds')
sce.imm.epi = readRDS('10xsmart_epiImm_harmony.rds')

#got sce4cnv####
#cnv_patient2.R####
{
names(table(sce.imm.epi$orig.ident))
sce.imm.epi.4cnv = subset(sce.imm.epi,orig.ident %in% c("LiM","LuM","VM",
                                                        "P2","P3","P4" ,
                                                        "PDAC_TISSUE_1","PDAC_TISSUE_11B" ,
                                                        # "PDAC_TISSUE_10","PDAC_TISSUE_12","P5",
                                                        "PDAC_TISSUE_13","PDAC_TISSUE_15", 
                                                        "PDAC_TISSUE_16","PDAC_TISSUE_2","PDAC_TISSUE_3"   ,  
                                                        "PDAC_TISSUE_4","PDAC_TISSUE_5","PDAC_TISSUE_6",
                                                        "PDAC_TISSUE_7" ,"PDAC_TISSUE_8" , "PDAC_TISSUE_9" ))
saveRDS(sce.imm.epi.4cnv,'10xsmart_epiImm4cnv_harmony.rds')
}
#imm.epi figure plot####
sce.imm.epi <- FindNeighbors(sce.imm.epi, dims = 1:5,reduction = 'harmony')
sce.imm.epi <- FindClusters(sce.imm.epi, reduction = "harmony",resolution = 0.7)
sce.imm.epi <- RunUMAP(object = sce.imm.epi, dims = 1:3, reduction = "harmony",
                       do.fast = T )
#table(sce.imm.epi@meta.data$RNA_snn_res.0.3)
DimPlot(sce.imm.epi,reduction = 'umap',group.by = 'selfCelltype')

saveRDS(sce.imm.epi,'10xsmart_epiImm_harmony.rds')
table(sce.imm.epi@meta.data$patient_id)

#cnv result####
{#/usr/local/bin/Rscript /home/lxxiao/xiaolixing/pdac/result_dir/code/cnv_result_sum.R
sce.imm.epi.4cnv.list=readRDS('/home/lxxiao/xiaolixing/pdac/10x/10xsmart_epiImm4cnv_harmony.rds')
sce.imm.epi.4cnv.list = SplitObject(sce.imm.epi.4cnv,split.by = 'orig.ident')
pro = names(sce.imm.epi.4cnv.list)
pred.all = c()
CNA.mat = c()
for(i in pro[1:length(pro)]){
  print(i)
  pred.test = read.table(paste0("cnv_res/",i,"_copykat_prediction.txt"),sep = '\t', header = T)
  pred.test = pred.test[-which(pred.test$copykat.pred=="not.defined"),] 
  pred.test$patient_id = i
  pred.all.list = c(pred.all,pred.test)
  CNA.test<-read.table(paste0('cnv_res/',i,'_copykat_CNA_results.txt'),sep = '\t', header = T)
  
  #CNA.test[1:4,1:4]
  cell_id = gsub("-",".",pred.test$cell.names)
  CNA.test = CNA.test[,cell_id]
  #CNA.mat=CNA.test
  CNA.mat.list = c(CNA.mat,CNA.test)
}

CNA.mat = do.call(cbind,CNA.mat.list,)
pred.all = do.call(rbind,pred.all)
colnames(CNA.mat) <-gsub("[.]","-",colnames(CNA.mat))
# table(colnames(CNA.mat))
pred.all[1:4,]
dim(pred.all)
table(pred.all$patient_id)
table(pred.all$copykat.pred)
write.csv(pred.all,'cnv_res/all_prediction.csv')
write.csv(CNA.mat,'cnv_res/all_CNA_mat.csv')
}

#sce get cnv result####
pred.all = read.csv('cnv_res/all_prediction.csv',header = T,row.names = 1)
CNA.mat  = read.csv('cnv_res/all_CNA_mat.csv',header = T,row.names = 1)
sce.all=readRDS('10xsmart_merge_harmony.rds')

##remove undefined cells
copy_data = data.frame(id = names(sce.all$orig.ident),
                       copykat.pred = "not.detected")

copy_data[match(pred.all$cell.names,copy_data$id),2] = pred.all$copykat.pred
table(copy_data$copykat.pred)
copy_data[copy_data$copykat.pred=="not.detected",]
sce.all@meta.data$copykat.pred <- copy_data$copykat.pred
saveRDS(sce.all,'10xsmart_merge_harmony.rds')
DimPlot(sce.all, reduction = "umap",group.by = "copykat.pred")

table(pred.all$patient_id)
table(pred.all$copykat.pred)
##remove undefined cells
copy_data = data.frame(id = names(sce.imm.epi$orig.ident),
                       copykat.pred = "not.detected")

copy_data[match(pred.all$cell.names,copy_data$id),2] = pred.all$copykat.pred
table(copy_data$copykat.pred)
copy_data[copy_data$copykat.pred=="not.detected",]
sce.imm.epi@meta.data$copykat.pred <- copy_data$copykat.pred
saveRDS(sce.imm.epi,'10xsmart_epiImm_CNV.rds')
#cna filt####
#cnv_ms.R####
sce.imm.epi=readRDS('10xsmart_epiImm_CNV.rds')
table(sce.imm.epi$tn_cell)
sce.imm.epi$tn_cell = factor(sce.imm.epi$tn_cell,
                             levels = c('Tumor',
                                        'Immune_cells',
                                        'Nontumor',
                                        'untest'))
#cnv result plot####
{
  genes_to_check = c(Epithelial,T_lymphocytes,
                     NK,B_lymphocytes,Monocyte , #Dendritic,
                     Endothelial,Fibroblastes,
                     Mast,Acinar_cells)
  # All on Dotplot ####
  library(ggplot2)
  
  p_umap1=DimPlot(sce.imm.epi, reduction = "umap",
                  label = F,
                  label.size = 5,
                  label.box = F,
                  pt.size = 0.1,
                  repel=F)+scale_color_simpsons()
  p_umap2=DimPlot(sce.imm.epi, reduction = "umap",group.by = "selfCelltype",
                  label = F,
                  label.size = 5,
                  label.box = F,
                  pt.size = 0.1,
                  repel=F) +scale_color_nejm()
  p_all_markers <- DotPlot(sce.imm.epi, features = unique(genes_to_check))+ 
    coord_flip()+scale_color_gradient(low = "white",high = 'red',
                                      guide = guide_colorbar(ticks.colour = "black",
                                                             frame.colour = "black"),
                                      name = "Average\nexpression")
  p_cnv=DimPlot(sce.imm.epi, reduction = "umap",group.by = "tn_cell",
                label = F,
                label.size = 5,
                label.box = F,
                pt.size = 0.1,
                repel=F)+scale_color_nejm()
  
  p_umap2+p_cnv+p_umap1#+p_all_markers
  #ggsave('figures/figures1_renew2/re_harmony_epi_t.pdf')
}

