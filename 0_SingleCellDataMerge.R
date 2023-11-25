#Seurat3.2.3
rm(list = ls())
set.seed(123)
setwd('/mywork/path')
system("type R") 
file.path(R.home("bin"), "R")

ps <- c('Seurat', 'ggplot2', 'stringr', 'reshape2','ggpubr','dplyr',
        "shiny","harmony",'tidyverse','ggpubr','ggplot2','parallel',
        'hdf5r','RColorBrewer','pillar')

#multiprocess####
#10x datasets
#GSE156405 10X
grep("filtered_feature_bc_matrix",list.files('GSE156405_pri_met/'),value = T)
fs_1=grep("filtered_feature_bc_matrix",list.files('GSE156405_pri_met/'),value = T)[c(1:7,19)]
samples_1=str_split(fs_1,'_',simplify = T)[,1]
data_path1 = paste0('GSE156405_pri_met/',fs_1)

#GSE155698 10X
grep("_TISSUE_",list.files('No9_hd/'),value = T)
fs_2=grep("_TISSUE_",list.files('No9_hd/'),value = T)
fs_2 = fs_2[-10]#-c("PDAC_TISSUE_14")
samples_2=fs_2
data_path2 = paste0('No9_hd/',fs_2,'/filtered_feature_bc_matrix')

fs = c(fs_1,fs_2)
samples = c(samples_1,samples_2)
data_path = c(data_path1,data_path2)
#integration 10x
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
  save(sce.big,file = '../revision_data/data/sce.big.merge.dataset.Rdata')
}
load('../revision_data/data/sce.big.merge.dataset.Rdata')

orig.ident.info = as.data.frame(table(sce.big$orig.ident))

colnames(orig.ident.info) = c('Sample ID','Number of cells')
orig.ident.info$DataSource='GSE156405'
orig.ident.info$DataSource[grep('TISSUE',orig.ident.info$`Sample ID`)]='GSE155698'
orig.ident.info$Technology = '10X'
table(orig.ident.info$DataSource=='GSE156405')#19
table(orig.ident.info$DataSource=='GSE155698')#8


sce.big@meta.data$DataSource='GSE156405'
sce.big@meta.data$DataSource[grep('TISSUE',sce.big$orig.ident)]='GSE155698'
table(table(sce.big$DataSource=='GSE156405',sce.big$orig.ident)[2,]!=0)#8
table(table(sce.big$DataSource=='GSE155698',sce.big$orig.ident)[2,]!=0)#19
table(sce.big$DataSource)

load('sce.big.merge.dataset.Rdata')
raw_sce = sce.big
raw_sce[["percent.mt"]] <- PercentageFeatureSet(raw_sce, pattern = "^MT-")
fivenum(raw_sce[["percent.mt"]][,1])
raw_sce[["percent.ribo"]] <- PercentageFeatureSet(raw_sce, pattern = "^RP[SL]")
raw_sce1 <- subset(raw_sce, subset = nFeature_RNA > 1500 & nCount_RNA > 2000 & percent.mt < 20)
selected_f <- rownames(raw_sce1)[Matrix::rowSums(raw_sce1@assays$RNA@counts > 0 ) > 30]
raw_sce1 <- subset(raw_sce1, features = selected_f)
raw_sce1
table(raw_sce1$orig.ident)#filter samples (cell<125)
samples  = names(table(raw_sce1$orig.ident))
samples_filtered = samples[which(as.numeric(table(raw_sce1$orig.ident))>125)]
raw_sce1 = subset(raw_sce1,orig.ident %in% samples_filtered)
table(raw_sce1$orig.ident)
#save(raw_sce1,file = '../revision_data/data2/sce.10XFiltered.dataset.Rdata')
load('../revision_data/data2/sce.10XFiltered.dataset.Rdata')

#GSE81547 SMART-SEQ2
smart_h = read.csv('No10_health_smart/all_count_matrix.csv',header = T,row.names = 1)
smart_phe = read.csv('No10_health_smart/SraRunTable.txt',header = T,row.names = 1)
colnames(smart_h)[1:10]
dim(smart_h)
head(smart_phe)
gsmid =str_split(colnames(smart_h),'[_]',simplify = T)[,1]
all(gsmid==smart_phe$Sample.Name)
celltype = smart_phe$inferred_cell_type
table(celltype)
smart_pro =paste0('SmartSeq2_',
                  substr(str_split(colnames(smart_h),'[.]',simplify = T)[,2],1,1),
                  ".",celltype)
table(smart_pro)
smart_sce = CreateSeuratObject(count = smart_h,project = smart_pro,
                               min.cell = 10,min.gene = 500)

merge_sce = merge(smart_sce,raw_sce1)#43264 samples
table(merge_sce$orig.ident)
merge_sce@meta.data$Technology = ifelse(str_split(merge_sce$orig.ident,'[_]',simplify = T)[,1]=='SmartSeq2',
                                        'SmartSeq2',
                                        '10X')
merge_sce@meta.data$DataSource='GSE156405'
merge_sce@meta.data$DataSource[grep('TISSUE',merge_sce$orig.ident)]='GSE155698'
merge_sce@meta.data$DataSource[grep('SmartSeq2',merge_sce$orig.ident)]='GSE81547'

merge_sce@meta.data$cell_type = 'TYPE'
merge_sce@meta.data$cell_type[grep('PDAC_TISSUE',merge_sce@meta.data$orig.ident)]='locPDAC'
merge_sce@meta.data$cell_type[grep('^P',merge_sce@meta.data$orig.ident)]='locPDAC'
merge_sce@meta.data$cell_type[grep('AdjNorm',merge_sce@meta.data$orig.ident)]='Norm'
merge_sce@meta.data$cell_type[grep('SmartSeq2',merge_sce@meta.data$orig.ident)]='Norm'
merge_sce@meta.data$cell_type[merge_sce@meta.data$orig.ident %in% c('LiM','VM','LuM')]='metPDAC'
table(merge_sce@meta.data$cell_type)

table(merge_sce@meta.data$Technology)
merge_sce@meta.data$cell_type2 = 'SMARTSEQ2'
samrt_celltype = str_split(merge_sce$orig.ident[grep('SmartSeq2_',merge_sce@meta.data$orig.ident)],
                           '[.]',simplify = T)[,2]
length(samrt_celltype)
merge_sce@meta.data$cell_type2[grep('SmartSeq2',merge_sce@meta.data$orig.ident)]= samrt_celltype

merge_sce@meta.data$cell_type2[grep('PDAC_TISSUE',merge_sce@meta.data$orig.ident)]='locPDAC'
merge_sce@meta.data$cell_type2[grep('^P',merge_sce@meta.data$orig.ident)]='locPDAC'
merge_sce@meta.data$cell_type2[grep('AdjNorm',merge_sce@meta.data$orig.ident)]='Norm'

merge_sce@meta.data$cell_type2[merge_sce@meta.data$orig.ident %in% c('LiM','VM','LuM')]='metPDAC'
table(merge_sce@meta.data$cell_type2)

table(merge_sce$orig.ident)
merge_sce$orig.ident = str_split(merge_sce$orig.ident,'[.]',simplify = T)[,1]

orig.ident.info = as.data.frame(table(merge_sce$orig.ident))
colnames(orig.ident.info) = c('Sample ID','Number of cells')
orig.ident.info$DataSource='GSE156405'
orig.ident.info$DataSource[grep('TISSUE',orig.ident.info$`Sample ID`)]='GSE155698'
orig.ident.info$DataSource[grep('SmartSeq2',orig.ident.info$`Sample ID`)]='GSE81547'
orig.ident.info$Technology = '10X'
orig.ident.info$Technology[grep('SmartSeq2',orig.ident.info$`Sample ID`)]='SmartSeq2'
table(orig.ident.info$DataSource,orig.ident.info$Technology)

table(merge_sce$DataSource)
table(merge_sce$DataSource,merge_sce$orig.ident)#43022
table(table(merge_sce$DataSource=='GSE155698',merge_sce$orig.ident)[2,]!=0)#18
table(table(merge_sce$DataSource=='GSE156405',merge_sce$orig.ident)[2,]!=0)#7
table(table(merge_sce$DataSource=='GSE81547',merge_sce$orig.ident)[2,]!=0)#8
#save(merge_sce,file = '../revision_data/data/sce.big.all.Rdata')#43022 samples
#write.csv(orig.ident.info,'../revision_data/data/orig_ident_info.csv')
#rm(list = ls())
#cell filter####
load('../revision_data/data/sce.big.all.Rdata')

raw_sce = merge_sce
raw_sce@active.ident =as.factor(raw_sce$orig.ident) 
length(table(raw_sce@active.ident))
#rm('merge_sce')
table(raw_sce@active.ident)
raw_sce@active.ident <- as.factor(raw_sce$orig.ident)
raw_sce[["percent.mt"]] <- PercentageFeatureSet(raw_sce, pattern = "^MT-")
fivenum(raw_sce[["percent.mt"]][,1])
raw_sce[["percent.ribo"]] <- PercentageFeatureSet(raw_sce, pattern = "^RP[SL]")
raw_sce1 <- subset(raw_sce, subset = nFeature_RNA > 1000 & nCount_RNA > 2000 & percent.mt < 20)
selected_f <- rownames(raw_sce1)[Matrix::rowSums(raw_sce1@assays$RNA@counts > 0 ) > 30]
raw_sce1 <- subset(raw_sce1, features = selected_f)
raw_sce1
samples = names(table(raw_sce1@meta.data$orig.ident))
samples_filtered = samples[which(as.numeric(table(raw_sce1@meta.data$orig.ident))>125)]
raw_sce1 = subset(raw_sce1,orig.ident%in% samples_filtered)
table(raw_sce1@meta.data$DataSource)
#save(raw_sce1,file = '../revision_data/data/sce.big.Filt.Rdata')#42855 samples
load('../revision_data/data/sce.big.Filt.Rdata')
table(raw_sce1$DataSource)
sce=raw_sce1
sce <- NormalizeData(sce, normalization.method =  "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000)
VariableFeaturePlot(sce)
sce <- ScaleData(sce)
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 

system.time({sce <- RunHarmony(sce,lambd=0.3,
                               group.by.vars = c("orig.ident"), 
                               max.iter.harmony=50,
                               max.iter.cluster = 30)})

pdf('../../Figures/10x_before_integration_umap.pdf',width = 8,height = 3)
after_integration = DimPlot(sce,reduction = "harmony",group.by = "DataSource",label=F,
                            pt.size=0.1,order = c('GSE81547','GSE156405','GSE155698'))+
  scale_color_manual(values = c('GSE156405'='#036EB8','GSE155698' = '#F0CD77','GSE81547'='#E26E42'))
before_integration = DimPlot(sce,reduction = "pca",group.by = "DataSource",label=F,
                             pt.size=0.1,order = c('GSE81547','GSE156405','GSE155698'))+
  scale_color_manual(values = c('GSE156405'='#036EB8','GSE155698' = '#F0CD77','GSE81547'='#E26E42'))
print(before_integration+after_integration)
dev.off()
saveRDS(sce,"../revision_data/data/sce.10XHarmony.rds")
