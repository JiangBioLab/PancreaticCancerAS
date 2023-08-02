#!/usr/local/bin/R --vanilla
#setwd('/home/lxxiao/xiaolixing/pdac/10x/patient/epi')
setwd('/home/lxxiao/xiaolixing/pdac/10x/cnv_res/')
.libPaths(c("/home/lxxiao/R/x86_64-pc-linux-gnu-library/4.0","/usr/local/lib64/R/library"))
.libPaths()
library(Seurat)
library(copykat)
					
#sce = readRDS('/home/lxxiao/xiaolixing/pdac/10x/all_RunTSNEUMAP_PDAC.rds')
sce = readRDS('/home/lxxiao/xiaolixing/pdac/10x/10xsmart_epiImm4cnv_harmony.rds')
sce.all.list <- SplitObject(sce , split.by = "orig.ident")
print(names(sce.all.list))
#relist = c('LiM','LuM','P2','P3','P5')
#norm_cell = as.character(read.csv(paste0(i,'_norm_name.csv'),header = T,row.names = 1)[,1])
for (i in names(sce.all.list)[13:14]) {
#for (i in relist){
	#norm_cell = as.character(read.csv(paste0(i,'_norm_name.csv'),header = T,row.names = 1)[,1])  
	print(i)
  	library(Seurat)
  	library(copykat)
  	epi_mat = sce.all.list[[i]]@assays$RNA@counts
        exp.rawdata <- as.matrix(epi_mat)
        #exp.rawdata = read.table(file = paste0(i,'_counts.txt'),
        #                 header = T,row.names = 1,sep = '\t')
  	#norm_cell =as.character(read.csv(paste0(i,'_norm_name.csv'),header = T,row.names = 1)[,1])
        norm_cell = colnames(sce.all.list[[i]])[which(sce.all.list[[i]]@meta.data$selfCelltype=='immune')]
        #epi_cell  = colnames(sce.all.list[[i]])[which(sce.all.list[[i]]@meta.data$selfCelltype=='Epithelial')]  
        
        print(length(norm_cell))
        #exp.rawdata_epi = exp.rawdata[,c(epi_cell,norm_cell)]
        #print(dim(exp.rawdata))
        #print(length(norm_cell))
        copykat.test <- copykat(rawmat=exp.rawdata, 
  						id.type="S", 
  						ngene.chr=1, 
  						win.size=25, 
  						KS.cut=0.1, 
  						sam.name=i, 
  						distance="euclidean", 
  						norm.cell.names=norm_cell,
  						#output.seg="TRUE", 
  						#plot.genes="TRUE", 
  						genome="hg20",
  						n.cores=5)
}

