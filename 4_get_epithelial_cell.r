setwd('/home/lxxiao/xiaolixing/pdac/10x/')
library(ggplot2)
library(stringr)
library(Seurat)
library(dplyr)
library(reshape2)
library(pheatmap)
##start ##################################
#cnv_ms.R####
#epi_imm####
sce.imm.epi <- readRDS('10xsmart_epiImm_CNV.rds')#epi with non cancdrf
sce <- RunUMAP(object = sce.imm.epi, dims = 1:10, reduction = "harmony",do.fast = T )
table(sce.imm.epi$tn_cell,sce.imm.epi$selfCelltype)

DimPlot(sce,reduction = 'umap',group.by = 'tn_cell',pt.size = 0.01)+
  scale_color_manual(values = c(
    'Immune_cells'='#93c47d',
    'Nontumor' = "gray70",
    'Tumor' = "#C80813FF",
    'untest' = "#47A0DB"
  ))
#ggsave('figures/figures1_renew2/epi_cnvRES_umap.pdf')


#get epi_cancer_cnv####
#epi_cancer####
sce.epi.cancer = subset(sce.imm.epi,selfCelltype=='Epithelial')
sce<-sce.epi.cancer
table(sce$tn_cell,sce$tissue)
cellIdent = ifelse(sce@meta.data$tn_cell=='Nontumor'&sce@meta.data$tissue=="locPDAC","Epi_loc_Norm",
                   ifelse(sce@meta.data$tn_cell=='Nontumor'&sce@meta.data$tissue=="metPDAC","Epi_met_Norm",
                          ifelse(sce@meta.data$tn_cell=='Tumor'&sce@meta.data$tissue=="locPDAC","Epi_loc_Tumor",
                                 ifelse(sce@meta.data$tn_cell=='Tumor'&sce@meta.data$tissue=="metPDAC","Epi_met_Tumor",
                                        'untest'))))
sce@meta.data$cellIdent=cellIdent
table(sce@meta.data$cell_from2)
table(sce@meta.data$cell_from2,sce@meta.data$tn_cell)

#epi.cancer####
sce <- NormalizeData(sce, normalization.method =  "LogNormalize", 
                     scale.factor = 10000)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000) 
sce <- ScaleData(sce)
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 
system.time({sce <- RunHarmony(sce,lambd=0.3,
                               group.by.vars = c("orig.ident"), 
                               max.iter.harmony=50,
                               max.iter.cluster = 50)})

DimPlot(sce, label = T,group.by = 'tissue',reduction = 'harmony')
sce <- FindNeighbors(sce, dims = 1:5,reduction = 'harmony')
sce <- FindClusters(sce, reduction = "harmony",resolution = 0.3)
sce <- RunUMAP(object = sce, dims = 1:4, reduction = "harmony",do.fast = T )

colnames(sce@meta.data)
sce = readRDS('10xsmart_epi_cancer_harmony.rds')
DimPlot(sce,reduction = 'umap',label =F,
        split.by = 'tissue',group.by = 'tn_cell',
        label.size = 5,
        label.box = F,
        pt.size = 0.5,
        repel=F)#+scale_color_nejm(alpha = 0.7)
#saveRDS(sce,'10xsmart_epi_cancer_harmony.rds')
#load epi_norm_cells####
sce.epi.norm <-readRDS('10xsmart_epi_norm.rds')
colnames(sce.epi.norm@meta.data)
table(sce.epi.norm$tissue)
table(is.na(sce.epi.norm$tissue))
sce.epi.norm$tissue[is.na(sce.epi.norm$tissue)] = 'SmartNorm'
sce.epi.norm <- NormalizeData(sce.epi.norm, normalization.method =  "LogNormalize", 
                     scale.factor = 10000)
GetAssay(sce.epi.norm,assay = "RNA")
sce.epi.norm <- FindVariableFeatures(sce.epi.norm, 
                            selection.method = "vst", 
                            nfeatures = 2000) 
sce.epi.norm <- ScaleData(sce.epi.norm)
sce.epi.norm <- RunPCA(object = sce.epi.norm, 
                       pc.genes = VariableFeatures(sce.epi.norm)) 
system.time({sce.epi.norm <- RunHarmony(sce.epi.norm,lambd=0.0001,
                               group.by.vars = c("orig.ident"), 
                               max.iter.harmony=10,
                               max.iter.cluster = 50)})
DimPlot(sce.epi.norm, label = T,group.by = 'tissue',
        reduction = 'harmony')
sce.epi.norm <- FindNeighbors(sce.epi.norm, dims = 1:3,reduction = 'harmony')
sce.epi.norm <- FindClusters(sce.epi.norm, reduction = "harmony",
                             resolution = 0.3)
sce.epi.norm <- RunUMAP(object = sce.epi.norm, dims = 1:4, 
                        reduction = "harmony",do.fast = T )
table(sce.epi.norm$seurat_clusters)
DimPlot(sce.epi.norm,reduction = 'umap')
sce.epi.norm@meta.data$selfcluster = ifelse(sce.epi.norm$seurat_clusters %in% c(2,4),
                                            'Cluster1','Cluster2')
table(sce.epi.norm$tissue,sce.epi.norm@meta.data$selfcluster)

DimPlot(sce.epi.norm,reduction = 'umap',label =F,
        split.by = 'selfcluster',
        group.by = 'tissue',
        label.size = 2,
        label.box = F,
        pt.size = 0.1,
        repel=F)+scale_color_aaas(alpha=0.5)

saveRDS(sce.epi.norm,'10xsmart_epi_norm_harmony.rds')
sce.epi.norm = readRDS('10xsmart_epi_norm_harmony.rds')

#find norm clusters markers####
sce.epi.norm.Markers = read.csv('10xsmart_epi_norm_top10_selfcluster_marker.csv',
                                header = T,row.names = 1)
head(sce.epi.norm.Markers)
rownames(sce.epi.norm.Markers)

#DotPlot norm ####
endocrine_cell = c("GCG","INS","SST")
duct_cells = c('KRT19')
genes_to_check = c(duct_cells,endocrine_cell)
p2 = DotPlot(sce.epi.norm, features = unique(genes_to_check))+ 
  coord_flip()+
  scale_color_gradient(low = "white",high = 'red',
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Average\nexpression")
sce.epi.norm@meta.data$selfcluster2 = ifelse(sce.epi.norm$seurat_clusters %in% c(2,4),
                                            'duct_cells','endocrine_cell')
mycolor = c('duct_cells' ='#EB6000',
            'endocrine_cell' = '#F3E094')
p1= DimPlot(sce.epi.norm,reduction = 'umap',label =F,
        group.by = 'selfcluster2',
        label.size = 4,
        label.box = F,
        pt.size = 0.1,
        repel=F)+scale_color_manual(values = mycolor)
p1
#ggsave('figures/figures1_renew2/epi_norm_clusters.pdf')
p2 = DotPlot(sce.epi.norm, features = unique(genes_to_check),
             group.by = 'selfcluster2')+
  scale_color_gradient(low = "white",high = 'red',
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Average\nexpression")+
  theme_bw()
p1+p2
FeaturePlot(sce.epi.norm,features = genes_to_check,pt.size = 0.1)+
  scale_color_gradientn(low = "grey",high = 'red',
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Average\nexpression")


#ggsave('figures/figures1_renew2/epi_norm_clusters.pdf')
saveRDS(sce.epi.norm,'10xsmart_epi_norm_harmony.rds')

#merge epi_cancerCNV_norm####
sce.epi.norm = readRDS('10xsmart_epi_norm_harmony.rds')
sce.epi.cancer = readRDS('10xsmart_epi_cancer_harmony.rds')
sce.epi.norm.duct = subset(sce.epi.norm,selfcluster2=='duct_cells')
sce.epi.cancer.tumor = subset(sce.epi.cancer,)
sce <- merge(sce.epi.cancer,sce.epi.norm.duct)
table(sce$selfcluster2)
table(sce@meta.data$tissue)
sce@meta.data$tissue2 = ifelse(sce$tissue %in% c('AdjNorm',
                                                  'SmartNorm'),'Norm','PDAC')
table(is.na(sce$tn_cell))
sce$tn_cell[is.na(sce$tn_cell)]='Norm'
table(sce$tn_cell)

table(sce@meta.data$tn_cell,sce@meta.data$tissue2)
table(sce@meta.data$tn_cell)

table(sce@meta.data$tn_cell,sce@meta.data$tissue)
cellIdent = ifelse(sce@meta.data$tn_cell=='Norm'& sce@meta.data$tissue =="SmartNorm" ,'Epi_Norm',
                   ifelse(sce@meta.data$tn_cell=='Norm'& sce@meta.data$tissue =="AdjNorm" ,'Epi_Norm',
                   ifelse(sce@meta.data$tn_cell=='Nontumor'&sce@meta.data$tissue=="locPDAC","Epi_loc_Norm",
                          ifelse(sce@meta.data$tn_cell=='Nontumor'&sce@meta.data$tissue=="metPDAC","Epi_met_Norm",
                                 ifelse(sce@meta.data$tn_cell=='Tumor'&sce@meta.data$tissue=="locPDAC","Epi_loc_Tumor",
                                               ifelse(sce@meta.data$tn_cell=='Tumor'&sce@meta.data$tissue=="metPDAC","Epi_met_Tumor",
                                                      'untest'))))))
table(cellIdent)
cellIdent_raw = ifelse(cellIdent %in% c('Epi_Norm'),'Epi_Norm',
                       ifelse(cellIdent %in% c('Epi_loc_Tumor','Epi_met_Tumor'),'Epi_Tumor',
                              ifelse(cellIdent %in% c('Epi_loc_Norm','Epi_met_Norm'),'Epi_Nontumor','untest')))
table(cellIdent_raw)
sce@meta.data$cellIdent = cellIdent
sce@meta.data$cellIdent_raw = cellIdent_raw
saveRDS(sce,'10xsmart_epi_cancerNorm_CNVRES.rds')
#epi cancer/norm harmony####
sce<-readRDS('10xsmart_epi_cancerNorm_CNVRES.rds')
table(sce$tn_cell,sce$cellIdent_raw)
table(sce$cellIdent_raw)
sce <-subset(sce,cellIdent_raw !='Epi_Nontumor' )
sce <-subset(sce,cellIdent_raw !='untest' )
table(sce$tn_cell,sce$cellIdent_raw)
sce <- NormalizeData(sce, normalization.method =  "LogNormalize", 
                     scale.factor = 10000)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000) 
sce <- ScaleData(sce)
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 
system.time({sce <- RunHarmony(sce,lambd=2,
                               group.by.vars = c("orig.ident"), 
                               max.iter.harmony=100,
                               max.iter.cluster = 50)})
DimPlot(sce, label = T,group.by = 'cellIdent',reduction = 'harmony')
sce <- FindNeighbors(sce, dims = 1:3,reduction = 'harmony')
sce <- FindClusters(sce, reduction = "harmony",resolution = 0.4)
sce <- RunUMAP(object = sce, dims = 1:4, reduction = "harmony",do.fast = T )

#table(is.na(sce$cell_from3))
#sce$cell_from3[is.na(sce$cell_from3)]='Epi_Norm'
#table(sce$tn_cell)
#table(sce$cellIdent_raw)
#table(sce$tissue2)
table(sce$cellIdent)
library(paletteer)
library(ggsci)
library(scales)
pal <- paletteer_d("ggsci::nrc_npg",n=10)

#ggsave('figures/figures1_renew2/cnv_res_cancer_vs_norm.pdf')
saveRDS(sce,'10xsmart_cancer_norm_CNV.rds')

#tumor norm####
sce.epi.tumor.norm = readRDS('10xsmart_cancer_norm_CNV.rds')
sce.epi.tumor.norm <- FindClusters(sce.epi.tumor.norm, reduction = "harmony",resolution = 1)
sce.epi.tumor.norm <- RunUMAP(object = sce.epi.tumor.norm, dims = 1:4, 
                              reduction = "harmony",do.fast = T )
table(sce.epi.tumor.norm$seurat_clusters)
library(RColorBrewer)
library(paletteer)
c_palettes<- palettes_c_names
d_palettes<- palettes_d_names
mycol1<-paletteer_d( "ggthemes::Classic_20",n=20)
my_pal = c('#33658a','#86bbd8','#f6ae2d','#f26419','#d88c9a','#fe6d73','#17c3b2',
           '#55a630','#004b23','#8e443d','#f4442e','#9d6b53','#5a189a','#aeb8fe',
           '#db4c40','#89bd9e','#279af1','#f0c987','#b0f2b4','#a7754d','##d1b3c4')
DimPlot(sce.epi.tumor.norm, label = T,group.by = 'seurat_clusters',split.by = 'cellIdent')+
  scale_color_manual(values = mycol1)
saveRDS(sce.epi.tumor.norm,'10xsmart_cancer_norm_CNV.rds')
sce.epi.tumor.norm = readRDS('10xsmart_cancer_norm_CNV.rds')
#subset Loc/Met/Nom####
table(sce.epi.tumor.norm$cellIdent)

loc = subset(sce.epi.tumor.norm,cellIdent=='Epi_loc_Tumor')
met = subset(sce.epi.tumor.norm,cellIdent=='Epi_met_Tumor')
Norm = subset(sce.epi.tumor.norm,cellIdent=='Epi_Norm')
table(Norm$orig.ident)

Norm = subset(Norm,orig.ident %in% c('AdjNorm_TISSUE_1',
                                     'AdjNorm_TISSUE_2',
                                     'AdjNorm_TISSUE_3',
                                     'smart_ductal'))
saveRDS(loc,'10xsmart_cancer_norm_loc.rds')
saveRDS(met,'10xsmart_cancer_norm_met.rds')
saveRDS(Norm,'10xsmart_cancer_norm_Norm.rds')
orig.ident = names(table(sce.epi.tumor.norm$orig.ident))
orig.ident2 = orig.ident[-c(22,23,24,26)]
sce.epi.tumor.norm2 = subset(sce.epi.tumor.norm,orig.ident %in% orig.ident2)
saveRDS(sce.epi.tumor.norm2,'10xsmart_cancer_norm_CNV_rmSMART.rds')


#cnv cancerNorm RES####
sce.epi.tumor.norm2 = readRDS('10xsmart_cancer_norm_CNV_rmSMART.rds')
table(sce.epi.tumor.norm2$cellIdent_raw)
sce =sce.epi.tumor.norm2
sce = RunUMAP(object = sce, dims = 1:10, 
              reduction = "harmony",do.fast = T )
table(sce$cellIdent)

table(sce$patient_id)

patient_id = sce$orig.ident
table(patient_id)
patient_id[grep('smart_',patient_id)] = 'Norm'
table(patient_id)

smart_sample_num = str_split(names(patient_id)[grep('^Norm',patient_id,value = F)],'[.]',simplify = T)[,2]
smart_sample_num = str_sub(smart_sample_num,1,1)

patient_id[grep('^Norm',patient_id,value = F)] = paste0(patient_id[grep('^Norm',patient_id,value = F)],smart_sample_num)
sce@meta.data$patient_id =  patient_id
table(sce$patient_id)


patient_id2 = patient_id
table(patient_id2)
patient_id2[grep('^Norm',patient_id2,value = F)] = paste0(patient_id2[grep('^Norm',patient_id2,value = F)],
                                                          '_GSE81547')
patient_id2[grep('^AdjNorm_TISSUE_',patient_id2,value = F)] = paste0('Norm',
                                                                     str_split(patient_id2[grep('^AdjNorm_TISSUE_',patient_id2,value = F)],
                                                                              '[_]',simplify = T)[,3],
                                                                     '_GSE155698')

patient_id2[grep('^PDAC_TISSUE_',patient_id2,value = F)] = paste0('P',
                                                                     str_split(patient_id2[grep('^PDAC_TISSUE_',patient_id2,value = F)],
                                                                               '[_]',simplify = T)[,3],
                                                                     '_GSE155698')
patient_id2[-grep('_',patient_id2,value = F)] = paste0(patient_id2[-grep('_',patient_id2,value = F)],
                                                                  '_GSE156405')

sce@meta.data$patient_id2 =  patient_id2
#saveRDS(sce,'10xsmart_cancer_norm_CNV_rmSMART.rds') #rewrite 20230214
#test = readRDS('10xsmart_cancer_norm_CNV_rmSMART.rds')

sample_id = names(table(sce@meta.data$patient_id2))
levels = c(
  sample_id[-c(grep('Norm',sample_id,value = F),1,2,30)],#loc
  'LiM_GSE156405','LuM_GSE156405','VM_GSE156405',#met
  grep('Norm',sample_id,value = T))
write.table(levels,row.names = F,col.names = F,eol = ',')
levels_v2 = c("P1_GSE155698",  "P2_GSE155698", "P3_GSE155698","P4_GSE155698","P5_GSE155698","P6_GSE155698","P7_GSE155698",
             "P8_GSE155698","P9_GSE155698","P11B_GSE155698","P13_GSE155698","P15_GSE155698","P16_GSE155698",
              "P2_GSE156405","P3_GSE156405","P4_GSE156405",
              
              "LiM_GSE156405","LuM_GSE156405","VM_GSE156405",
              
              "Norm1_GSE155698","Norm2_GSE155698","Norm3_GSE155698",
              "NormA_GSE81547","NormB_GSE81547","NormC_GSE81547","NormD_GSE81547","NormE_GSE81547","NormF_GSE81547","NormG_GSE81547","NormH_GSE81547")
sce@meta.data$patient_id2 = factor(sce@meta.data$patient_id2,
                                   levels = levels_v2
                                   )

mycolor=c('#CE2020','#228B22','#1F78B4',#adj
          '#CE2020','#228B22',#met,
          
          '#FDB462','#8B658B', 
          '#4876FF','#00BFFF', '#EE82EE',
          "#E01516" ,"#e98686" ,"#ef6a50",'#b6d7a8',"#7f9775" ,"#ffbb66",
          "#b7911e" ,"#864A68","#09813e" ,"#846dbe", "#B839D8",
          #'#8B8682','#CDC9C9',
          
          '#FDB462',#smart
          '#1F78B4')#met

library(ColorBrewer)
display.brewer.all()
brewer.pal(8,"Set1")
cols<-brewer.pal(8, "Set1")
locpal<-colorRampPalette(cols[1:4])
loc_color<-pal(16)

normpal<-colorRampPalette(cols[5:8])
norm_color<-normpal(11)

brewer.pal(5,"Set3")
metcols<-brewer.pal(5, "Set3")
metpal<-colorRampPalette(metcols)
met_color<-metpal(3)
met_color<-c("#8DD3C7" ,"#BEBADA" ,"#FB8072")

mycolors = c(loc_color,met_color,norm_color)


DimPlot(sce,reduction = 'umap',label =F,
        group.by = 'patient_id2',
        split.by = 'cellIdent',
        label.size = 2,
        label.box = F,
        pt.size = 0.01,
        repel=F)+scale_color_manual(values = mycolors)


sce = RunUMAP(object = sce, dims = 1:7, 
              reduction = "harmony",do.fast = T )
#ggsave('figures/figures1_renew2/replot_epi_patientID_umap.pdf')#4:12

DimPlot(sce,reduction = 'umap',label =F,
        group.by = 'cellIdent',
        #split.by = 'cellIdent',
        label.size = 2,
        label.box = F,
        pt.size = 0.01,
        repel=F)+
  scale_color_manual(values  = c('Epi_met_Tumor'='#5a9740',
                                 'Epi_loc_Tumor'='#549acb',
                                 'Epi_Norm'='#f76205'))
table(sce$cellIdent)
#ggsave('figures/figures1_renew2/epi_finalRES_umap.pdf')
#sce.epi.tumor.norm.filt ####
sce.epi.tumor.norm = readRDS('10xsmart_cancer_norm_CNV_rmSMART.rds')
genes.tumor.norm = rownames(sce.epi.tumor.norm)
genes.tumor.norm.raw = str_split(genes.tumor.norm,'[.]',simplify = T)[,1]
genes.tumor.norm.fdup = genes.tumor.norm.raw[!duplicated(genes.tumor.norm.raw)]
length(genes.tumor.norm.fdup)
genes.tumor.norm.fdup2 = genes.tumor.norm.fdup[-grep('^ERCC',genes.tumor.norm.fdup)]
sce.epi.tumor.norm.fdup = sce.epi.tumor.norm[genes.tumor.norm.fdup2,]
saveRDS(sce.epi.tumor.norm.fdup,'10xsmart_cancer_norm_CNV_rmSMART_fdup.rds')
#cluster_marker_top5.R
sce.epi.tumor.norm.fdup= readRDS('10xsmart_cancer_norm_CNV_rmSMART_fdup.rds')
emt_signature = read.csv('emt_gene_signature.csv',header = T,row.names = 1)
########
#find pdac_vs_hd markers####
sce_file = sce.epi.tumor.norm.fdup[emt_signature$Gene,]
sce_file
#find celltype markers
sce_file %>% dplyr::glimpse()
Idents(sce_file) <-sce_file$cellIdent
#Levels(sinlgecelldata) 

markers_df_loc <- FindMarkers(object = sce_file, 
                              ident.1 = 'Epi_loc_Tumor', 
                              ident.2 = 'Epi_Norm',
                              logfc.threshold = 1,
                              slot = "counts",
                              min.pct = 0.25)
print(x = head(markers_df_loc))
markers_df_loc$emt_gene = rownames(markers_df_loc)
markers_df_loc$cluster = paste0('loc_vs_hd')

markers_df_met <- FindMarkers(object = sce_file, 
                              ident.1 = 'Epi_met_Tumor', 
                              ident.2 = 'Epi_Norm',
                              logfc.threshold = 1,
                              slot = "counts",
                              min.pct = 0.25)
print(x = head(markers_df_met))
markers_df_met$emt_gene = rownames(markers_df_met)
markers_df_met$cluster = paste0('met_vs_hd')

markers_df = rbind(markers_df_loc,markers_df_met)
write.csv(markers_df,'/home/lxxiao/xiaolixing/pdac/10x/pdac_vs_hd_degs.csv')



sce = RunUMAP(object = sce, dims = 1:50, reduction = "harmony",do.fast = T )
DimPlot(sce,reduction = 'umap',group.by = 'cellIdent',pt.size = 0.01)+
  scale_color_manual(values  = c('Epi_loc_Tumor'='#006834',
                                 'Epi_met_Tumor'='#F81919',
                                 'Epi_Norm'='#EB6000'))
#ggsave('figures/figures1_renew2/epi_tissuetype_umap')
DimPlot(sce,reduction = 'umap',group.by = 'orig.ident',pt.size = 0.01)

#run cellmarkers####
#cluster_marker_top5.R####
cellmarkers = read.csv('10xsmart_cancer_norm_CNV_rmSMART_cellIdent_marker.csv',
                       header = T,row.names = 1)
cellmarkers$gene = str_split(rownames(cellmarkers),'[.]',simplify = T)[,1]
cellmarkers_fdup = read.csv('10xsmart_cancer_norm_CNV_rmSMART_fdup_cellIdent_marker.csv',
                       header = T,row.names = 1)
cellmarkers_fdup$gene = rownames(cellmarkers_fdup)
head(cellmarkers)
cellmarkers = cellmarkers[-grep('^ERCC',rownames(cellmarkers)),]
cellmarkers_fdup = cellmarkers
cellmarkers = cellmarkers_fdup

cellmarkers = cellmarkers[-which(rownames(cellmarkers) %in% c('alignment-not-unique',
                             "no-feature",
                             "ambiguous")),]

cellmarkers$gene = rownames(cellmarkers)
cellmarkers_top5 = cellmarkers %>% group_by(cluster) %>% top_n(50,avg_logFC)
head(cellmarkers_top5)

cellmarkers_top5$cluster2 = str_split(cellmarkers_top5$cluster,'cluster',simplify = T)[,2]




#sub_harmony.R
loc_h=readRDS('10xsmart_cancer_norm_loc_harmony.rds')
met_h=readRDS('10xsmart_cancer_norm_met_harmony.rds')
Norm_h=readRDS('10xsmart_cancer_norm_Norm_harmony.rds')
loc_h <- RunUMAP(object = loc_h, dims = 1:10, reduction = "harmony",do.fast = T )

met_h <- RunUMAP(object = met_h, dims = 1:10, reduction = "harmony",do.fast = T )

Norm_h <- FindClusters(Norm_h, reduction = "harmony",resolution = 0.3)
Norm_h <- RunUMAP(object = Norm_h, dims = 1:7, reduction = "harmony",do.fast = T )


DimPlot(met_h,reduction = 'umap',pt.size = 0.1,repel=F)
DimPlot(Norm_h,reduction = 'umap',pt.size = 0.1,repel=F)

saveRDS(loc_h,'10xsmart_cancer_norm_loc_harmony.rds')
saveRDS(met_h,'10xsmart_cancer_norm_met_harmony.rds')
saveRDS(Norm_h,'10xsmart_cancer_norm_Norm_harmony.rds')

#barplot####
library(reshape2)
library(ggplot2)
sce = loc_h
sce = met_h
sce = Norm_h
data_meta = data.frame(tissue = as.character(sce@meta.data$orig.ident),
                       cluster = as.factor(paste0('C',sce@meta.data$seurat_clusters)))
data = data_meta

##long to wide
library(reshape2)
data0 = dcast(data,cluster~tissue)
head(data0)
dim(data0)
sum(data0[,2])+sum(data0[,3])+sum(data0[,4])
#计算百分比
per = function(colindex){
  col_res = (colindex*100)/sum(colindex)
  return(col_res)
}
data_res = data0[,-1]
data1 = lapply(data_res,per)
data1 = as.data.frame(do.call('rbind',data1))
dim(data1)
colnames(data1)=data0$cluster
data1 = data1[order(data1$C0,decreasing = F),]
data1$tissue = as.factor(rownames(data1))

head(data1)
#wide to long dataframe
data2 <- melt(data1, id.vars = c('tissue'))
head(data2)
table(data2$variable)

data2$tissue = factor(data2$tissue,levels = data1$tissue)
library(paletteer)
d_palettes<- palettes_d_names
table(palettes_d_names$package)
palettes_d_names$palette[grep('colorBlindness',palettes_d_names$package)]

mycol1<-paletteer_d('colorBlindness::ModifiedSpectralScheme11Steps',n =(ncol(data1)-1))
p2=ggplot(data2, aes(x = tissue,y=value,fill = variable ))+
  geom_col(position = 'stack', width = 1)+
  geom_bar(position = "stack", stat = "identity", width = 0.6,) +
  scale_fill_manual(values =mycol1)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'bottom')
#ggsave('figures/figures1_renew2/loc_patient_barplot.pdf')
p1=DimPlot(sce,reduction = 'umap',pt.size = 0.1, repel=F,label = F)+
  scale_color_manual(values = mycol1)+
 theme(legend.position = 'none')
p1/p2
#ggsave('figures/figures1_renew2/loc_patient_UMAPBar.pdf')
