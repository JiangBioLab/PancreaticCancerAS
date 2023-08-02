pred.all<-read.csv('cnv_res/all_prediction.csv')
CNA.mat.raw<-read.csv('cnv_res/all_CNA_mat.csv')
sdata <- readRDS('10xsmart_epiImm4cnv_harmony.rds')

sdata = sce.imm.epi
sce.imm.epi =sdata 
CNA.mat =CNA.mat.raw[,-1]
CNA.mat[1:30,1:3]
pred.all[1:3,]
names(sdata$orig.ident)[1:3]
rownames(pred.all) = pred.all$cell.names
cell_info =pred.all

cna_colname = colnames(CNA.mat.raw[,-1])
cna_colname = paste0(str_split(cna_colname,'[.]',simplify = T)[,2],'-',
                     str_split(cna_colname,'[.]',simplify = T)[,3])
cna_colname[1:2]
colnames(CNA.mat) = cna_colname

CNA.epi<-CNA.mat[,intersect(colnames(CNA.mat),names(sdata$orig.ident))]

all(colnames(CNA.epi)==rownames(cell_info))#T
sh2 = CNA.epi
sh2[1:3,1:3]
dim(sh2)
length(SD)
CNV_score <- data.frame(MS = colMeans(sh2^2), 
                        SS = apply(sh2^2, 2, sum), 
                        SD = apply(sh2, 2, sd))
##
CNV_score$Row.names <- colnames(sh2)
cell_info$Row.names <- cell_info$cell.names

cell_info2 <- plyr::join(cell_info, CNV_score, by="Row.names") # boxplot for celltype
head(cell_info2)
## MS top 5% cells
top_MS_cells <- arrange(cell_info2, desc(MS))[1:round(dim(cell_info2)[1]*0.05),]$Row.names  # Top 5%

## calculate correlation : corr using 1 cell vs. top_MS_cells
tmp <- data.frame(Ave_tumor = rowMeans(sh2[,top_MS_cells]))
for(i in 1:dim(cell_info2)[1]){
  print(i)
  cell_info2$COR[i] <-  cor(sh2[,i], 
                            data.frame(Ave_tumor = rowMeans(sh2[,top_MS_cells])))
}


##
table(is.na(cell_info2$MS) )
cutoff.score=0.02
cutoff.corr=0.2
target.celltypes="aneuploid"
cell_info3 <- cell_info2
cell_info3$celltype = sdata@meta.data$selfCelltype[match(cell_info3$cell.names,
                                                         names(sdata$orig.ident))]
table(cell_info3$celltype)
table(cell_info3$copykat.pred)
rownames(cell_info3) <- cell_info3$Row.names
head(cell_info3)

tumorcells <- filter(cell_info3, ((MS > cutoff.score | COR > cutoff.corr) & copykat.pred %in% target.celltypes))$Row.names
nontumorcells <- cell_info3$Row.names[!(cell_info3$Row.names %in% c(tumorcells))]
immunecells <- cell_info3$Row.names[cell_info3$celltype == "immune"]

## only classified tumor vs. non-tumor ##
cell_info3$cell_index <- rep("X", dim(cell_info3)[1])
cell_info3[tumorcells,]$cell_index <- "Tumor"
cell_info3[nontumorcells,]$cell_index <- "Nontumor"
cell_info3[immunecells,]$cell_index <- "Immune_cells"
head(cell_info3)
table(cell_info3$cell_index)
table(cell_info3$MS>0.04)
cell_info3$MS2 = cell_info3$MS
cell_info3$MS2[cell_info3$MS2>0.04]=0.04
table(cell_info3$copykat.pred,cell_info3$cell_index)
table(cell_info3$celltype,cell_info3$cell_index)

write.csv(cell_info3,'10xsmart_epiImm4cnv_cancer_detected.csv')
cell_info3 = read.csv('10xsmart_epiImm4cnv_cancer_detected.csv',
                      header = T,row.names = 1)

table(sce.imm.epi$patient_id)
sce.imm.epi@meta.data$tn_cell = 'untest'
table(is.na(cell_info3$cell_index))
length(na.omit(match(cell_info3$cell.names,colnames(sce.imm.epi))))
#cell_info3_t = subset.data.frame(cell_info3,cell_info3$cell_index=="T_lymphocytes")
cell_info4 = cell_info3[na.omit(match(colnames(sce.imm.epi),cell_info3$cell.names)),]
dim(cell_info4)
sce.imm.epi@meta.data$tn_cell[match(cell_info4$cell.names,
                                    colnames(sce.imm.epi))]=cell_info4$cell_index
table(sce.imm.epi@meta.data$selfCelltype, sce.imm.epi@meta.data$tn_cell)
table(sce.all@meta.data$orig.ident, sce.all@meta.data$selfCelltype)
table(sce.imm.epi$tn_cell)


saveRDS(sce.imm.epi,'10xsmart_epiImm_CNV.rds')

names(table(cell_info3$patient_id))

cell_info3$patient_id = factor(cell_info3$patient_id,
                               levels = c( "LiM","LuM","VM","P2","P3","P4" ,
                                           "PDAC_TISSUE_1","PDAC_TISSUE_2","PDAC_TISSUE_3" , 
                                          "PDAC_TISSUE_4" , "PDAC_TISSUE_5","PDAC_TISSUE_6",
                                          "PDAC_TISSUE_7","PDAC_TISSUE_8","PDAC_TISSUE_9",
                                          "PDAC_TISSUE_11B", "PDAC_TISSUE_13",  "PDAC_TISSUE_15" ,
                                          "PDAC_TISSUE_16"))
## 2D plot of MS score and correlation ####
expos<-ggplot(cell_info3, aes(x=MS2, y= COR)) + 
  geom_point(aes(fill=cell_index), size=1, alpha=.8, shape=21, colour="black") +
  scale_fill_manual(values = c("Tumor"="red",
                               "Immune_cells" = "gray70",
                               "Nontumor"="dodgerblue1")) +
  facet_wrap(~patient_id)+
  geom_vline(xintercept = cutoff.score, colour="black", size=0.1, linetype = "longdash") + 
  geom_hline(yintercept = cutoff.corr, colour="black", size=0.1, linetype = "longdash") +
  xlab("MS score") + ylab("CNV correlation") + theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=12), axis.text.x  = element_text(size=12)) +
  theme(axis.title.y = element_text(face="bold", size=12), axis.text.y  = element_text(size=12)) +
  theme(panel.border=element_rect(fill=NA, colour="black", size=1), legend.position = 'right')
expos
#ggsave('figures/figures1_renew2/cnv_ms_epi_Tcells.pdf')