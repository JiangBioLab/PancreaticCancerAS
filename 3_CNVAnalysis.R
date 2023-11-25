rm(list = ls())
library(Seurat)
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(harmony)
library(shiny)
library(RColorBrewer)
library(pillar)
library(parallel)
library(ggplot2)
library(devtools)
library(infercnv)
display.brewer.all()
brewer.pal.info
display.brewer.pal(5,"Set1")
brewer.pal(5,"Set1")
mycol = colorRampPalette(c("#8fbD7A", "#A99BF9","#E82D43"))(25)
scales::show_col(mycol)
set.seed(123)

#Analyzing CNVs with CopyCat####
#cnv_result.R
#! /usr/bin/env Rscript
library(future)
library(future.apply)
plan(multisession, workers=3)
options(future.globals.maxSize = 1000000 * 1024^5)
library(parallel)
options("sp_evolution_status"=2)
library(sp)
library(sp)
library(Seurat)
library(copykat)
setwd('/data_new/xiaolixing/pdac_revision/data/copykat')
sce =readRDS('/data_new/xiaolixing/pdac_revision/data/10xsmart_epiImm_harmony.rds')
sce.all.list <- SplitObject(sce , split.by = "orig.ident")
print(names(sce.all.list)[1:18])
for (i in names(sce.all.list))  {
  print(i)
  epi_mat = sce.all.list[[i]]@assays$RNA@counts
  exp.rawdata <- as.matrix(epi_mat)
  
  norm_cell = colnames(sce.all.list[[i]])[which(sce.all.list[[i]]@meta.data$selfCelltype=='Immune')]
  print(length(norm_cell))
  
  copykat.test <- copykat(rawmat=exp.rawdata, 
                          id.type="S", 
                          ngene.chr=1, 
                          win.size=25, 
                          KS.cut=0.1, 
                          sam.name=i, 
                          distance="euclidean", 
                          norm.cell.names=norm_cell,
                          genome="hg20",
                          n.cores=5)
  print('test')
}

print('copykat finished')


#! /usr/bin/env Rscript
setwd('/data_new/xiaolixing/pdac_revision/data/')
library(future)
library(future.apply)
plan(multisession, workers=5)
options("sp_evolution_status"=2)
library(sp)
library(Seurat)
sce=readRDS('/data_new/xiaolixing/pdac_revision/data/10xsmart_epiImm_harmony.rds')
sce = SplitObject(sce,split.by = 'orig.ident')
pro = names(sce)
print(pro)
pred.all.list = list()
CNA.mat.list =list()

for(x in 1:18){
  i = names(sce)[x]
  print(i)
  pred.test0 = read.table(paste0("copykat/",i,"_copykat_prediction.txt"),sep = '\t', header = T)
  if(length(as.numeric(table(pred.test0$copykat.pred)))==3){
    pred.test = pred.test0[-which(pred.test0$copykat.pred=="not.defined"),]
  }else{
    pred.test = pred.test0
  }
  pred.test$patient_id = i
  pred.all.list[[x]] = pred.test
  names(pred.all.list)[x] =i
  
  CNA.test<-read.table(paste0('copykat/',i,'_copykat_CNA_results.txt'),sep = '\t', header = T)
  
  #CNA.test[1:4,1:4]
  cell_id = gsub("-",".",pred.test$cell.names)
  CNA.test = CNA.test[,cell_id]
  #CNA.mat=CNA.test
  CNA.mat.list[[x]] = CNA.test
  names(CNA.mat.list)[x] =i
}
saveRDS(CNA.mat.list,'cnvcopykatCNAMat_Result.rds')
saveRDS(pred.all.list,'cnvcopykatpredAll_Result.rds')

CNA.mat.list = readRDS('cnvcopykatCNAMat_Result.rds')
pred.all.list= readRDS('cnvcopykatpredAll_Result.rds')

CNA.mat = as.data.frame(do.call(cbind,CNA.mat.list))
pred.all = as.data.frame(do.call(rbind,pred.all.list))

write.csv(pred.all,'cnvcopykat_all_prediction.csv')
write.csv(CNA.mat,'cnvcopykat_all_CNA_mat.csv')
print('result saved')

sce.all=readRDS('scefinal.patientFilterd.rds')
copy_data = data.frame(id = names(sce.all$orig.ident),
                       copykat.pred = "not.detected")
copy_data[match(pred.all$cell.names,copy_data$id),2] = pred.all$copykat.pred
all(names(sce.all$orig.ident) == copy_data$id)
sce.all@meta.data$copykat.pred <- copy_data$copykat.pred
saveRDS(sce.all,'scefinal.copyKatAllPatient.Rdata')

#epi.imm
sce.all =readRDS('10xsmart_epiImm_harmony.rds')
copy_data = data.frame(id = names(sce.all$orig.ident),copykat.pred = "not.detected")
copy_data[match(pred.all$cell.names,copy_data$id),2] = pred.all$copykat.pred
all(names(sce.all$orig.ident) == copy_data$id)
sce.all@meta.data$copykat.pred <- copy_data$copykat.pred
saveRDS(sce.all,'copyKat_epiImm.rds')

sce.imm.epi = readRDS('copyKat_epiImm.rds')
sce.imm.epi@active.ident = as.factor(sce.imm.epi$orig.ident)
system.time({sce.imm.epi <- RunHarmony(sce.imm.epi,lambd=0.6,
                                       group.by.vars = c("orig.ident"), 
                                       max.iter.harmony=50,
                                       max.iter.cluster = 30)})
sce.imm.epi <- FindNeighbors(sce.imm.epi, dims = 1:5,reduction = 'harmony')
sce.imm.epi <- FindClusters(sce.imm.epi, reduction = "harmony",resolution = 0.7)
sce.imm.epi <- RunUMAP(object = sce.imm.epi, dims = 1:4, reduction = "harmony" )
saveRDS(sce.imm.epi,'/data_new/xiaolixing/pdac_revision/data/copyKat_epiImm_umap.rds')


#Analyzing CNVs with infercnv####
#cd /data/xiaolixing/ref/
#python gtf_to_position.py --attribute_name "gene_name"  gencode.v38.annotation.gtf gene_pos.txt
#114 server /data/xiaolixing/ref/gene_pos.txt
#! /usr/bin/env Rscript
#library(future)
#plan("multiprocess", workers = 2)
#options(future.globals.maxSize = 1000 * 1024^3)
sce =readRDS('/data_new/xiaolixing/pdac_revision/data/10xsmart_epiImm_harmony.rds')
matrix<-as.matrix(sce@assays$RNA@counts)
cellAnn <- data.frame(row.names = names(sce$orig.ident),
                      cellname = names(sce$orig.ident),
                      celltype =sce$selfCelltype )
gene_pos = read.table('/data/xiaolixing/ref/gene_pos.txt',header = T,sep = '\t')
#write.table(matrix,'/data_new/xiaolixing/pdac_revision/data/infercnv/raw_counts_matrixFilt.txt',sep = '\t',row.names = T,col.names = T)
#write.table(cellAnn,'/data_new/xiaolixing/pdac_revision/data/infercnv/cellAnnFilt.txt',sep = '\t',row.names = F,col.names = F)

#!/usr/bin/env Rscript
#infercnv.R
library(future)
library(future.apply)
plan(multisession, workers=8)
options("sp_evolution_status"=2)
library(sp)
library(Seurat)
library(infercnv)

infercnv_obj =CreateInfercnvObject(
  raw_counts_matrix='/data_new/xiaolixing/pdac_revision/data/infercnv_Filt//raw_counts_matrixFilt.txt',
  gene_order_file=  '/data/xiaolixing/ref/gene_pos.txt',
  annotations_file='/data_new/xiaolixing/pdac_revision/data/infercnv_Filt/cellAnnFilt.txt',
  ref_group_names='Immune',
  delim = "\t",
  max_cells_per_group = NULL,
  min_max_counts_per_cell = c(100, +Inf),
  chr_exclude = c("chrX", "chrY", "chrM")
)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='/data_new/xiaolixing/pdac_revision/data/infercnv_Filt/', 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,num_threads=8,
                             no_plot = T,no_prelim_plot = T,
                             HMM=TRUE)
print('finished')
save(infercnv_obj,file = 'infercnv_obj_Filt.Rdata')

#cnvScore result compare(CopyCat vs infercnv)####
setwd('/data_new/xiaolixing/pdac_revision/data/')
cnvScore <- function(data){
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    rescale(to=c(-1, 1)) %>% 
    t()
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

#copykat
sce.all=readRDS('scefinal.patientFilterd.rds')
CNA.mat.list = readRDS('cnvcopykatCNAMat_Result.rds')
pred.all.list= readRDS('cnvcopykatpredAll_Result.rds')

CNA.mat = as.data.frame(do.call(cbind,CNA.mat.list))
pred.all = as.data.frame(do.call(rbind,pred.all.list))
cna_colname = colnames(CNA.mat)
cna_colname = paste0(str_split(cna_colname,'[.]',simplify = T)[,2],'-',
                     str_split(cna_colname,'[.]',simplify = T)[,3])
colnames(CNA.mat) = cna_colname
cnv_score_copykat <- cnvScore(CNA.mat)
colnames(cnv_score_copykat) = 'cnv_score_copykat'
cnv_score_copykat$patientID = rownames(cnv_score_copykat)
write.table(cnv_score_copykat,'cnv_score_copykat.txt',sep = '\t',col.names = T,row.names = T)

#infercnv
infercnv_obj = readRDS("infercnv_Filt/run.final.infercnv_obj")
expr.data = infercnv_obj@expr.data
cnv_score_infercnv <- cnvScore(expr.data)
colnames(cnv_score_infercnv) = 'cnv_score_infercnv'
cnv_score_infercnv$patientID = rownames(cnv_score_infercnv)
write.table(cnv_score_infercnv,'cnv_score_infercnv.txt',sep = '\t',col.names = T,row.names = T)

sce.all@meta.data$cnv_score_copykat = NA
sce.all@meta.data$cnv_score_copykat[match(cnv_score_copykat$patientID,
                                          names(sce.all$orig.ident))] = cnv_score_copykat$cnv_score_copykat
sce.all@meta.data$cnv_score_infercnv= NA
sce.all@meta.data$cnv_score_infercnv[match(cnv_score_infercnv$patientID,
                                           names(sce.all$orig.ident))] = cnv_score_infercnv$cnv_score_infercnv
#saveRDS(sce.all,'scefinal.patientFilterd.CNVscore.rds')
sce.all =readRDS('scefinal.patientFilterd.CNVscore.rds')
sce.all_mata = sce.all@meta.data
sce.all_mata_epiImm = subset(sce.all_mata,selfCelltype %in%c('Epithelial1','Immune'))

copykat = sce.all_mata_epiImm$cnv_score_copykat
infercnv = sce.all_mata_epiImm$cnv_score_infercnv

#CNV Score Range scaling
sce.all_mata_epiImm$copykat_rangeScaling = (copykat-min(copykat,na.rm = T))/max(copykat,na.rm = T)
sce.all_mata_epiImm$infercnv_rangeScaling = (infercnv-min(infercnv,na.rm = T))/max(infercnv,na.rm = T)

pdf('../../Figures/cnv_compare.pdf',width = 6,height = 3)
p1 = ggplot(sce.all_mata_epiImm,aes(x = selfCelltype,
                                    y = copykat_rangeScaling))+
  geom_violin(aes(fill=selfCelltype),color="NA")+
  ylim(0,0.15)+stat_compare_means(aes(group = selfCelltype), label = "p.signif",method = "t.test",)+
  scale_fill_manual(values = c('Immune' = '#8FBD7A','Epithelial1'='#85C7E8'))+
  theme_classic()

p2 = ggplot(sce.all_mata_epiImm,aes(x = selfCelltype,
                                    y = infercnv_rangeScaling))+
  geom_violin(aes(fill=selfCelltype),color="NA")+
  ylim(0,0.15)+stat_compare_means(aes(group = selfCelltype), label = "p.signif",method = "t.test",)+
  scale_fill_manual(values = c('Immune' = '#8FBD7A','Epithelial1'='#85C7E8'))+
  theme_classic()
print(p1+p2)
dev.off()

med_copykat = quantile(sce.all_mata_epiImm$copykat_rangeScaling,na.rm = T)[2]
med_infercnv = quantile(sce.all_mata_epiImm$infercnv_rangeScaling,na.rm = T)[2]

sce.all_mata_epiImm$copykat_zscore_type = ifelse(sce.all_mata_epiImm$copykat_rangeScaling > med_copykat,'up','down')
sce.all_mata_epiImm$infercnv_zscore_type= ifelse(sce.all_mata_epiImm$infercnv_rangeScaling > med_infercnv,'up','down')
Epithelial1_imm_copykat = as.data.frame(table(sce.all_mata_epiImm$copykat_zscore_type,sce.all_mata_epiImm$selfCelltype))
Epithelial1_imm_infercnv= as.data.frame(table(sce.all_mata_epiImm$infercnv_zscore_type,sce.all_mata_epiImm$selfCelltype))


pdf('../../Figures/cnv_barplot.pdf',width = 6,height = 2)
p1 = ggplot(sce.all_mata_epiImm, aes(selfCelltype))+
  geom_bar(aes(fill=copykat_zscore_type), position = 'fill')+
  scale_fill_manual(values = c("up"  =  '#0B68AB',"down"  =  "#8fbd7a"))+
  theme_classic()

p2 = ggplot(sce.all_mata_epiImm, aes(selfCelltype))+
  geom_bar(aes(fill=infercnv_zscore_type), position = 'fill')+
  scale_fill_manual(values = c("up"  =  '#0B68AB',"down"  =  "#8fbd7a"))+
  theme_classic()
print(p1+p2)
dev.off()

#ms####
library(plyr)
library(stringr)
cnv_score_copykat = read.csv('cnv_score_copykat.txt',sep = '\t',header = T,row.names = 1)

CNA.mat.list = readRDS('cnvcopykatCNAMat_Result.rds')
pred.all.list= readRDS('cnvcopykatpredAll_Result.rds')

CNA.mat = as.data.frame(do.call(cbind,CNA.mat.list))
pred.all = as.data.frame(do.call(rbind,pred.all.list))

cna_colname = colnames(CNA.mat)
cna_colname = paste0(str_split(cna_colname,'[.]',simplify = T)[,2],'-',
                     str_split(cna_colname,'[.]',simplify = T)[,3])
cna_colname[1:2]
colnames(CNA.mat) = cna_colname
sh2 = CNA.mat
dim(sh2)
dim(cnv_score_copykat)
CNV_score <- data.frame(CNV_score = cnv_score_copykat$cnv_score_copykat,
                        cellname = cnv_score_copykat$patientID,
                        row.names = rownames(cnv_score_copykat))


## MS top 5% cells
top_MS_cells <- arrange(CNV_score, desc(CNV_score))[1:round(dim(CNV_score)[1]*0.05),]$cellname  # Top 5%

## calculate correlation : corr using 1 cell vs. top_MS_cells
tmp <- data.frame(Ave_tumor = rowMeans(sh2[,top_MS_cells]))
for(i in 1:dim(CNV_score)[1]){
  print(i)
  CNV_score$COR[i] <-  cor(sh2[,i], data.frame(Ave_tumor = rowMeans(sh2[,top_MS_cells])))
}
cutoff.score=0.02
cutoff.corr=0.2
target.celltypes="aneuploid"

cell_info = sce.all_mata_copykat
cell_info$cellname = rownames(sce.all_mata_copykat)
cell_info <- plyr::join(cell_info, CNV_score, by="cellname") # boxplot for celltype
rownames(cell_info)=cell_info$cellname

tumorcells <- filter(cell_info, ((CNV_score > cutoff.score | COR > cutoff.corr) & copykat.pred %in% target.celltypes))$cellname
nontumorcells <- cell_info$cellname[!(cell_info$cellname %in% c(tumorcells))]
immunecells <- cell_info$cellname[cell_info$selfCelltype == "Immune"]

## only classified tumor vs. non-tumor ##
cell_info$cell_index <- rep("X", dim(cell_info)[1])
cell_info[tumorcells,]$cell_index <- "Tumor"
cell_info[nontumorcells,]$cell_index <- "Nontumor"
cell_info[immunecells,]$cell_index <- "Immune_cells"
table(cell_info$cell_index)
#Immune_cells     Nontumor        Tumor
#5083         7808         4072
write.csv(cell_info,'epiImm4cnv_cancer_detected.csv')

sce.all = readRDS('scefinal.patientFilterd.CNVscore.rds')
cell_info = read.csv('epiImm4cnv_cancer_detected.csv',header = T,row.names = 1)
sce.all@meta.data$tn_celltype = 'untest'

sce.all@meta.data$tn_celltype[match(cell_info$cellname,
                                    colnames(sce.all))]=cell_info$cell_index
saveRDS(sce.all,'epiImm_CNVfinal_type.rds')

sce.all = readRDS('epiImm_CNVfinal_type.rds')
table(cell_info$orig.ident)
cell_info$patient_id = factor(cell_info$orig.ident,
                              levels = c( "LiM","LuM","VM","P2","P3","P4" ,"P5",
                                          "PDAC_TISSUE_1","PDAC_TISSUE_2","PDAC_TISSUE_3" , 
                                          "PDAC_TISSUE_5","PDAC_TISSUE_6",
                                          "PDAC_TISSUE_8","PDAC_TISSUE_9",
                                          "PDAC_TISSUE_11B", "PDAC_TISSUE_13",  "PDAC_TISSUE_15" ,
                                          "PDAC_TISSUE_16"))
#2D plot of MS score and correlation 
pdf('../../Figures/cnv_dotFinalplot.pdf',width = 9.6,height = 5)
expos<-ggplot(cell_info, aes(x=CNV_score, y= COR)) + 
  geom_point(aes(fill=cell_index), size=1, alpha=.8, shape=21, colour="black") +
  scale_fill_manual(values = c("Tumor"="red",
                               "Immune_cells" = "gray70",
                               "Nontumor"="dodgerblue1")) +
  facet_wrap(~patient_id,nrow = 3)+xlim(0,0.2)+
  geom_vline(xintercept = cutoff.score, colour="black", size=0.05, linetype = "longdash") + 
  geom_hline(yintercept = cutoff.corr, colour="black", size=0.05, linetype = "longdash") +
  xlab("CNV score") + ylab("CNV correlation") + theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=12), axis.text.x  = element_text(size=12)) +
  theme(axis.title.y = element_text(face="bold", size=12), axis.text.y  = element_text(size=12)) +
  theme(panel.border=element_rect(fill=NA, colour="black", size=1), legend.position = 'right')
print(expos)
dev.off()

sce.imm.epi = subset(sce.all ,selfCelltype%in% c('Epithelial1','Immune'))
sce.imm.epi@active.ident = as.factor(sce.imm.epi$orig.ident)
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
sce.imm.epi <- FindNeighbors(sce.imm.epi, dims = 1:5,reduction = 'harmony')
sce.imm.epi <- FindClusters(sce.imm.epi, reduction = "harmony",resolution = 0.7)
sce.imm.epi <- RunUMAP(object = sce.imm.epi, dims = 1:10, reduction = "harmony" )
saveRDS(sce.imm.epi,'/data_new/xiaolixing/pdac_revision/data/copyKat_epiImm_umap_926.rds')

#epi-Cancer and Norm-Duct####
sce.imm.epi = readRDS('copyKat_epiImm_umap_926.rds')
sce.epi.norm = readRDS('10xsmart_epi_norm.rds')

sce.imm.epi = subset(sce.imm.epi,orig.ident!='P5')
tumor_epi = subset(sce.imm.epi,tn_celltype=='Tumor')
saveRDS(tumor_epi,'copyKat_tumorepi.rds')
sce.tumor_norm = merge(tumor_epi,sce.epi.norm)
sce.tumor_norm$tn_celltype=ifelse(is.na(sce.tumor_norm$tn_celltype),'Norm','Tumor')
table(sce.tumor_norm$orig.ident)
sce.tumor_norm@active.ident = as.factor(sce.tumor_norm$orig.ident)
GetAssay(sce.tumor_norm,assay = "RNA")
sce.tumor_norm <- FindVariableFeatures(sce.tumor_norm, 
                                       selection.method = "vst", 
                                       nfeatures = 2000) 
sce.tumor_norm <- ScaleData(sce.tumor_norm)
sce.tumor_norm <- RunPCA(object = sce.tumor_norm, pc.genes = VariableFeatures(sce.tumor_norm)) 
system.time({sce.tumor_norm <- RunHarmony(sce.tumor_norm,lambd=0.1,
                                          group.by.vars = c("orig.ident"), 
                                          max.iter.harmony=50,
                                          max.iter.cluster = 30)})

pdf('../../Figures/epiTumorNorm_harmony.pdf',height = 4,width = 12)
p1 = DimPlot(sce.tumor_norm,reduction = 'pca',group.by = 'DataSource')+
  scale_color_manual(values = c('GSE156405'='#036EB8','GSE155698' = '#F0CD77','GSE81547'='#E26E42'))+
  theme_classic()
p2 = DimPlot(sce.tumor_norm,reduction = 'harmony',group.by = 'DataSource')+
  scale_color_manual(values = c('GSE156405'='#036EB8','GSE155698' = '#F0CD77','GSE81547'='#E26E42'))+
  theme_classic()
print(p1+p2)
dev.off()

sce.tumor_norm <- FindNeighbors(sce.tumor_norm, dims = 1:25,reduction = 'harmony')
sce.tumor_norm <- FindClusters(sce.tumor_norm, reduction = "harmony",resolution = 0.7)
sce.tumor_norm <- RunUMAP(object = sce.tumor_norm, dims = 1:18, reduction = "harmony" )
saveRDS(sce.tumor_norm,'epi.tumor_norm.rds')