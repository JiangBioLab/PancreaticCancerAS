library(scales)
library(dplyr)
library(BiocManager)
library(stringr)
library(stringi)
library(reshape2)
library(GEOquery)
library(shiny)
library(RColorBrewer)
library(pillar)
library(limma)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(devtools)
library(ggpubr)
library(ggvenn)
library(tidyverse)
library(tibble)
library(stats)
library(SingleCellExperiment)
#install.packages("pacman")
#library(pacman)
#pacman::p_load(apeglm)
library(apeglm)
library(clusterProfiler)
library(org.Hs.eg.db)
library(R.utils)
library(pacman)
library(KEGG.db)
library(estimate)
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text( hjust = 1 ),#angle = 45,
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }

#load data####
pdac1_expr = read.table('/GSE144561_rawCountsAllsamples.txt.gz',header = T,row.names = 1,sep = '\t')
gseid = 'GSE144561'
glpid = 'GSE144561'
gset <- getGEO(gseid, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep(glpfile, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
phe = pData(gset)
colnames(phe)
pdac1_phe = data.frame(title = phe$title,
                       sample_id = phe$geo_accession,
                       description =  gsub('[-]','.',phe$description),
                       patient = phe$`patient diagnosis:ch1`)
match(colnames(pdac1_expr),pdac1_phe$description)
pdac1_phe$description[35]= gsub('[,]','.',pdac1_phe$description[35])

#healthy donor  localized PDAC metastatic PDAC
#21              42              18
#write.csv(pdac1_phe,'CTC_phenoData_raw.csv')#/data_new/xiaolixing/pdac_revision/data/
#write.table(pdac1_expr,'CTC_exprData_raw.txt',sep = '\t',quote = T,row.names = T,col.names = T)#/data_new/xiaolixing/pdac_revision/data/
pdac1_expr_pro= cbind(rownames(pdac1_expr),pdac1_expr)
colnames(pdac1_expr_pro)[1]='GeneSymbol'
write.table(pdac1_expr_pro,'CTC_expr4est.txt',quote = F,row.names = F,sep = '\t')

#Tumor sample purity assessment####
file_dir='CTC_expr4est.txt'
filterCommonGenes(input.f = file_dir,output.f = "CTC_estimateGene.gct",id="GeneSymbol")
estimateScore("CTC_estimateGene.gct","CTC_estimat_score.gct",platform = "illumina")
estimat_score = read.table('CTC_estimat_score.gct',skip = 2,header = T,check.names = F,row.names = 1)
estimat_score = estimat_score[,-1]
estimat_score = as.data.frame(t(estimat_score))
estimat_score$TumorPurity = cos(0.6049872018+0.0001467884*estimat_score$ESTIMATEScore)
estimat_score$type = ifelse(estimat_score$patient=='healthy donor','Norm PBMC','PDAC CTC')
#mycompare = combn(names(table(estimat_score$patient)),2,simplify = F)
mycompare= list(c('Norm PBMC','PDAC CTC'))
pdf('../../Figures/CTC_Tumor_purity.pdf',width = 4,height = 6)
p1 =ggplot(estimat_score, aes(x = type, y = TumorPurity))+ 
  labs(y="Tumor purity",x= NULL,title = "PDAC CTC")+  
  geom_violin(aes(fill = type),position=position_dodge(0.2),width=0.5,size=0.4,
               outlier.alpha = 1, outlier.size = 0.5)+ 
  theme_bw() + mytheme+
  scale_fill_manual(values = c('Norm PBMC' = '#85c7e8','PDAC CTC' = '#8fbd7a'))+
  scale_y_continuous(labels = scales::percent)+
  stat_compare_means(label =  "p.format",
                     method = 't.test',comparisons = mycompare,
                     hide.ns = F)
print(p1)
dev.off()

estimat_score$patient = pdac1_phe$patient
write.csv(estimat_score,'CTC_tumor_purity_estimatResult_raw.csv')
estimat_score_pdac = subset(estimat_score,patient !='healthy donor')
estimat_score_hd = subset(estimat_score,patient =='healthy donor')

estimat_score_pdac = subset(estimat_score_pdac,TumorPurity >0.5)
estimat_score_hd = subset(estimat_score_hd,TumorPurity <0.5)
estimat_score_final = rbind(estimat_score_pdac,estimat_score_hd)
write.csv(estimat_score_final,'CTC_tumor_purity_estimatResult_final.csv')

#PCA####
estimat_score_final = read.csv('CTC_tumor_purity_estimatResult_final.csv',header = T,row.names = 1)
count_expr= read.table('CTC_exprData_raw.txt',header = T,row.names = 1,sep = '\t')
count_expr_filtered = count_expr[,match(rownames(estimat_score_final),colnames(count_expr))]

estimat_score_final$type = ifelse(estimat_score_final$patient=='healthy donor','HD',
                                  ifelse(estimat_score_final$patient=='localized PDAC','locPDAC','metPDAC'))
write.csv(estimat_score_final,'CTC_tumor_purity_estimatResult_final.csv')


#PCA
group <- interaction(estimat_score_final$type)
dge = DGEList(counts=count_expr_filtered,group=group)
plotMDS(dge, top=10,col = as.numeric(group))
mds_dat = plotMDS(dge, top=80,col = as.numeric(group),gene.selection ='common')
pca.data <- data.frame(sample =  rownames(mds_dat$distance.matrix.squared), 
                       Type=estimat_score_final$type, PC1=mds_dat$x,PC2=mds_dat$y)
pca.data$Type2 = ifelse(pca.data$Type=='HD','HD','PDAC CTC')

pdf('../../Figures/CTC_PCA_plot.pdf',width = 4,height = 4.5)  
p1= ggscatter(pca.data,x="PC1", y="PC2", color="Type2", ellipse=T, ellipse.type="convex",
              size=3, palette=c('#8fbd7a','#85c7e8', '#EB6000'), 
              label=NULL,repel=F, main="PCA plot")+
  ylab('Principal Component 2(19%)')+
  xlab('Principal Component 1(55%)')+
  theme_bw()+mytheme
print(p1)
dev.off()

pdf('../../Figures/CTC_PCA_plot2.pdf',width = 4,height = 4.5)  
p1= ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=F, ellipse.type="convex",
              size=3, palette=c('#8fbd7a','#85c7e8', '#EB6000'), 
              label=NULL,repel=F, main="PCA plot")+
  ylab('Principal Component 2(19%)')+
  xlab('Principal Component 1(55%)')+
  theme_bw()+mytheme
print(p1)
dev.off()

#Immune cell composition####
#counts2TPM
pdac1_phe = read.csv('CTC_phenoData_raw.csv',header = T,row.names = 1)
pdac1_expr = read.table('CTC_exprData_raw.txt',header = T,row.names = 1,sep = '\t')
exons_length1 = read.csv('/data_new/xiaolixing/pdac_revision/data/code/gencode.v38.annotation_v2.csv',header = T)
exons_length1 = read.csv('gencode.v38.annotation_v2.csv',header = T,row.names = 1)
epxr = cbind(rownames(pdac1_expr),pdac1_expr)
colnames(epxr)[1]='GeneSymbol'
epxr$GeneSymbol = str_split(epxr$GeneSymbol,'[.]',simplify = T)[,1]
epxr = epxr[!duplicated(epxr$GeneSymbol),]
rownames(epxr) = epxr$GeneSymbol
mergedataTP = epxr[,-1]
le = exons_length1$Length[match(rownames(mergedataTP),exons_length1$gene_name)]
mergedataTP = mergedataTP[-which(is.na(le)),]#33126 81
eff_length = le[-which(is.na(le))] #33125

table(is.na(le))
countToTpm <- function(counts, effLen)
{
  rate <- log(counts+1) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

tpms <- apply(mergedataTP,2,countToTpm,eff_length)
tpms[1:3,1:3]
dim(tpms)
write.table(tpms,'CTC_expr_tpm.txt',quote = F,row.names = T,sep = '\t')

tpms_ciber= cbind(rownames(tpms),tpms)
colnames(tpms_ciber)[1]='GeneSymbol'
write.table(tpms_ciber,'CTC_expr4ciber.txt',quote = F,row.names = F,sep = '\t')
tpms_ciber[1:3,1:5]

#cibersort calculate
source('../Cibersort.R')
pkgs <- c("matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "cowplot","ggpubr","bslib","ggthemes")
lapply(pkgs, library, character.only = T)
tpms_expr= read.table('CTC_expr_tpm.txt',header = T,row.names = 1,sep = '\t')
tpms_expr[1:4,1:4]
tpms_expr_filtered = tpms_expr[,match(rownames(estimat_score_final),colnames(tpms_expr))]
tpms_expr_ciber= cbind(rownames(tpms_expr_filtered),tpms_expr_filtered)
colnames(tpms_expr_ciber)[1]='GeneSymbol'
write.table(tpms_expr_ciber,'CTCFilterd_expr4ciber.txt',quote = F,row.names = F,sep = '\t')

result1 <- CIBERSORT('LM22.txt','CTCFilterd_expr4ciber.txt', perm = 1000, QN = F) 
#write.csv(result1,'ciber_Output_perm1000.csv')

cibersort_raw <- read.csv("ciber_Output_perm1000.csv",row.names = 1,header = T)
immCell_four_type <- read.table("Cibersort_four_types.txt", header = T, row.names = NULL, sep = "\t")
head(immCell_four_type)
all(rownames(cibersort_raw)==rownames(estimat_score_final))
colnames(cibersort_raw)
cibersort_raw = cibersort_raw[,1:22]
cibersort_raw$group = estimat_score_final$patient
cibersort_raw$sample= rownames(cibersort_raw)
TME_four_new = melt(cibersort_raw)
colnames(TME_four_new) = c("Group","Sample","Immune.cells","Composition")
names(table(TME_four_new$Immune.cells))

TME_four_new2 = left_join(TME_four_new, immCell_four_type, by = "Immune.cells") 
TME_four_final = aggregate(Composition ~ Sample+Group + Types, data = TME_four_new2, sum)

TME_four_final[1:4,]
#write.table(TME_four_final,'CTC_immCell4Type.txt',quote = F,row.names = T,col.names = T,sep ='\t')


TME_four_final$Group = as.factor(TME_four_final$Group)
comlist= combn(unique(TME_four_final$Group), 2, simplify =FALSE)
pdf('../../Figures/CTC_wbc_composition.pdf',width = 8,height = 4.5)
box_four_immtypes <- ggplot(TME_four_final, aes(x = Group, y = Composition))+ 
  labs(y="Cell composition",x= NULL,title = "PDAC CTC")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.2),width=0.3,size=0.4,
               outlier.alpha = 1, outlier.size = 0.5)+ 
  theme_bw() + mytheme+
  scale_fill_manual(values = c('healthy donor' = '#EB6000','localized PDAC' = '#8fbd7a','metastatic PDAC'='#85c7e8'))+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(~ Types,scales = "free",ncol = 4)+ 
  stat_compare_means(  method = "anova")+
  stat_compare_means(label =  "p.signif",
                     hide.ns = T)
print(box_four_immtypes) 
dev.off()


#pseudobulk EMT DEG analysis(epitheilial)####
sce.tumor_norm = readRDS('epi.tumor_norm.rds')
emt_signature = read.csv('emt_gene_signature.csv',header = T,row.names = 1)
pesudo_counts <- sce.tumor_norm@assays$RNA@
pesudo_metadata <- sce.tumor_norm@meta.data
pesudo_metadata$cluster_id <- factor(sce.tumor_norm$cell_type)
pesudo_metadata$sample_id <- factor(sce.tumor_norm$orig.ident)
sce_bulk <- SingleCellExperiment(assay = list(counts = pesudo_counts),
                                 colData = pesudo_metadata)
#saveRDS(sce_bulk,'epi.tumor_norm.sceBulk.rds')
#sce_bulk = readRDS('epi.tumor_norm.sceBulk.rds')
group <- colData(sce_bulk)[,c("cluster_id","sample_id",'DataSource')]
#pb <- aggregate(t(counts(sce_bulk)), list(group$sample_id), sum)
rownames(pb) <- pb$Group.1
#pb_t <- as.data.frame(t(pb[,-1]))
#write.table(pb_t,'epi.tumor_norm.pseudoBulk_count.txt',sep = '\t',quote = F,col.names = T,row.names = T)
pb_t = read.table('epi.tumor_norm.pseudoBulk_count.txt',sep = '\t',header = T,row.names = 1)
pb_metadat = as.data.frame(table(group$sample_id,group$cluster_id))
pb_metadat =subset(pb_metadat,Freq!=0)
colnames(pb_metadat) = c('Patient_id','type','NumberOfCells')
rownames(pb_metadat) = pb_metadat$Patient_id
table(group$DataSource)
pb_metadat$batch = "GSE156405"
pb_metadat$batch[grep('_TISSUE_',pb_metadat$Patient_id)]='GSE155698'
pb_metadat$batch[grep('SmartSeq2',pb_metadat$Patient_id)]='GSE81547'
table(pb_metadat$batch)
pb_metadat$batch = as.factor(pb_metadat$batch)
pb_metadat$type = as.factor(pb_metadat$type)

pb_t_final = pb_t[,match(rownames(pb_metadat),colnames(pb_t))]
all(colnames(pb_t_final)==rownames(pb_metadat))
#write.table(pb_metadat, 'combat_Expr_metadat.txt',quote = F,sep = '\t',col.names = T,row.names = T)
#write.table(pb_t_final,'epi.final.pseudoBulk_count.txt',sep = '\t',quote = F,col.names = T,row.names = T)
pb_t_final = read.table('epi.final.pseudoBulk_count.txt',sep = '\t',header = T,row.names = 1)
pb_metadat = read.table('combat_Expr_metadat.txt',sep = '\t',header = T,row.names = 1)
#add Bulk duct datasets(liver met)
#GSE193268
pancreatic_duct_cellline = read.table('GSE193268_total_count.txt.gz',header = T)
pancreatic_duct_cellline[1:4,]
Liver_metPDAC = pancreatic_duct_cellline[!duplicated(pancreatic_duct_cellline$name),1:4]
rownames(Liver_metPDAC) = Liver_metPDAC$name
Liver_metPDAC[1:4,]
Liver_metPDAC = Liver_metPDAC[,-1]
#write.table(Liver_metPDAC,'GSE193268_Liver_metPDAC_count.txt',quote = T,row.names = T,col.names = T)
Liver_metPDAC_meta = data.frame(Patient_id = colnames(Liver_metPDAC),
                            type = 'metPDAC',
                              batch = 'GSE193268Bulk',row.names =colnames(Liver_metPDAC) )
#add bulk GSE57973_duct_cell(normal) 
#featureCounts -T 4 -p -t exon -g gene_name  -a /data/xiaolixing/ref/gencode.v38.annotation.gtf -o featureCount/raw.featureCount *bam
# awk -F '\t' '{print $1,$7,$8,$9,$10}' OFS='\t' raw.featureCount >subread_matrix.out
norm_duct_cellline = read.table('GSE57973_duct_cell_featureCount.txt',sep = '\t',header = T)
norm_duct_cellline_filtered = norm_duct_cellline[rowSums(norm_duct_cellline[,2:5])!=0,]
norm_duct_cellline_filtered = norm_duct_cellline_filtered[!duplicated(norm_duct_cellline_filtered$Geneid),]
rownames(norm_duct_cellline_filtered) = norm_duct_cellline_filtered$Geneid
norm_duct_cellline_filtered = norm_duct_cellline_filtered[,-1]
#write.table(norm_duct_cellline_filtered,'GSE57973_norm_duct_count.txt',quote = T,row.names = T,col.names = T)
duct_norm_meta = data.frame(Patient_id = colnames(norm_duct_cellline_filtered),
                                type = 'Norm',
                                batch = 'GSE57973Bulk',
                            row.names =colnames(norm_duct_cellline_filtered) )

expr_list = list(epi_expr = pb_t_final,
                 norm_duct = norm_duct_cellline_filtered,
                 liver_metPDAC = Liver_metPDAC)
pdac_meta = list(epi_meta =pb_metadat[,-3],
                 norm_duct_meta = duct_norm_meta,
                 liver_metPDAC_meta = Liver_metPDAC_meta)
raw_pdac = list(expr_count_raw = expr_list,
                pdac_meta = pdac_meta)
#saveRDS(raw_pdac,'preudoCountRawData.rds')

raw_pdac = readRDS('preudoCountRawData.rds')
expr_list = raw_pdac$expr_count_raw
pb_t_final = expr_list$epi_expr
norm_duct_cellline_filtered = expr_list$norm_duct
Liver_metPDAC = expr_list$liver_metPDAC

meta_list = raw_pdac$pdac_meta
Liver_metPDAC_meta = meta_list$liver_metPDAC_meta
duct_norm_meta = meta_list$norm_duct_meta

#Merge emt gene expression matrix and prepare sample information table
pb_t_EMT = pb_t_final[na.omit(match(emt_signature$Gene,rownames(pb_t_final))),]
all(colnames(pb_t_EMT)==rownames(pb_metadat))
pb_metadat

Liver_metPDAC_EMT = Liver_metPDAC[na.omit(match(emt_signature$Gene,rownames(Liver_metPDAC))),]
duct_norm_EMT = norm_duct_cellline_filtered[na.omit(match(emt_signature$Gene,rownames(norm_duct_cellline_filtered))),]

emt_genes_inter = intersect(rownames(pb_t_EMT),rownames(Liver_metPDAC_EMT))
emt_genes_inter = intersect(emt_genes_inter,rownames(duct_norm_EMT))

pb_t_EMT_final = pb_t_EMT[emt_genes_inter,]
duct_norm_meta_EMT_final = duct_norm_EMT[emt_genes_inter,]
Liver_metPDAC_EMT_final = Liver_metPDAC_EMT[emt_genes_inter,]
  
epi_pdac = bind_cols(pb_t_EMT_final,Liver_metPDAC_EMT_final,duct_norm_meta_EMT_final)
epi_pdac_meta = bind_rows(pb_metadat[,-3],Liver_metPDAC_meta,duct_norm_meta)
  table(epi_pdac_meta$batch)

#remove batch effect
library(FactoMineR)
library(factoextra)
library(sva)
pdf('../../Figures/preuso_epi_rmBatch.pdf',width = 15,height = 4)
pre.pca <- PCA(t(epi_pdac),graph = FALSE)
P1 = fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = epi_pdac_meta$batch,
             addEllipses = TRUE,
             legend.title="Group"  )+mytheme+ggtitle('Before integration')
#ggsave('../../Figures/preuso_epi_beforeBatch.pdf')
epi_pdac = as.matrix(epi_pdac)
combat_Expr <- ComBat_seq(epi_pdac,batch = epi_pdac_meta$batch)
combat.pca <- PCA(t(combat_Expr),graph = FALSE)
P2 = fviz_pca_ind(combat.pca,
             geom= "point",
             col.ind = epi_pdac_meta$batch,
             addEllipses = TRUE,
             legend.title="Group"  )+mytheme+ggtitle('Remove batch effect')
print(P1+P2)
dev.off()
#ggsave('../../Figures/preuso_epi_afterBatch.pdf')
#write.table(combat_Expr,'combat_Expr_emtExpr.txt',quote = F,sep = '\t',col.names = T,row.names = T)
#write.table(epi_pdac_meta,'combat_meta_all.txt',quote = F,sep = '\t',col.names = T,row.names = T)
epi_pdac_meta = read.csv('combat_meta_all.txt',header = T,sep = '\t',row.names = 1)

#EMT gene DEG analysis(CTC)####
rm(list =  ls())
#Differential expression analysis with DEGseq2
#CTC EMT-DEG####
estimat_score_final = read.csv('CTC_tumor_purity_estimatResult_final.csv',header = T,row.names = 1)
emt_signature = read.csv('emt_gene_signature.csv',header = T,row.names = 1)
count_expr = read.table('CTC_exprData_raw.txt',sep = '\t',header = T,row.names = 1)
estimat_score_final$sampleSource = ifelse(estimat_score_final$sampleType=='NormPBMC',
                                          'NormPBMC','CTC')
count_expr_raw = count_expr[,match(row.names(estimat_score_final),colnames(count_expr))]
count_expr_emt = count_expr[na.omit(match(emt_signature$Gene,rownames(count_expr))),
                            match(row.names(estimat_score_final),colnames(count_expr))]
dim(count_expr_emt)#310 45
all(row.names(estimat_score_final)==colnames(count_expr_emt))

#DESeq2
countData=count_expr_emt
countData <- countData[rowMeans(countData)>1,] 
condition <- factor(estimat_score_final$sampleSource)
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, design = ~ condition)
dds1 <- DESeq(dds) 
res_CTC <- results(dds1,contrast = c('condition','CTC','NormPBMC'))

#sc emt deg####
combat_Expr = read.table('combat_Expr_emtExpr.txt',sep = '\t',header = T,row.names = 1)
pb_metadat = read.table('combat_meta_all.txt',sep = '\t',header = T,row.names = 1)
pb_metadat$sampleType = ifelse(pb_metadat$type=='Norm','Norm','EpiCancer')

countData=combat_Expr
phe_meta = pb_metadat
countData <- countData[rowMeans(countData)>1,] 
table(emt_signature$Category[match(rownames(countData),emt_signature$Gene)])
condition <- factor(phe_meta$sampleType)
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds1 <- DESeq(dds) 
res_epiCancer <- results(dds1,contrast = c('condition','EpiCancer','Norm'))

degseq2_res = list(ctc_deg = res_CTC,
                   duct_deg = res_epiCancer)
saveRDS(degseq2_res,'DEGseq2_diffRaw.rds')
#FindSig####

FindSig = function(resdat,foldChange_cutoff,padj_cutoff){
  res_raw <- data.frame(resdat, stringsAsFactors = FALSE, check.names = FALSE)
  resSig<- res_raw[which(abs(res_raw$log2FoldChange) > foldChange_cutoff & res_raw$padj < padj_cutoff),] 
  resSig$gene = rownames(resSig)
  resSig$change = ifelse(resSig$log2FoldChange>foldChange_cutoff,'Up-regulation','Down-regulation')
  resSig$type = emt_signature$Category[match(resSig$gene,emt_signature$Gene)]
  return(resSig)
}
res_CTCSig = FindSig(res_CTC,1,0.05)
res_EPISig = FindSig(res_epiCancer,1,0.05)
dim(res_CTCSig)#197   9
dim(res_EPISig)#28  9
degseq2_resSig = list(CTCSig_lfc1 = res_CTCSig,
                      DuctSig_lfc1 = res_EPISig)
saveRDS(degseq2_resSig,'DEGseq2_diffSig.rds')

#GO enrich####
#venn plot####
EMT_deg_Sig= readRDS('DEGseq2_diffSig.rds')
res_CTCSig = EMT_deg_Sig$CTCSig_lfc1
res_EPISig = EMT_deg_Sig$DuctSig_lfc1

res_CTCSig$cellSource = 'CTC'
res_EPISig$cellSource = 'Duct'

resSig_all = rbind(res_CTCSig,res_EPISig)


#compare_gene_type 
Gene_ID_epi <- bitr(resSig_all$gene, fromType="SYMBOL", 
                    toType="ENTREZID", 
                    OrgDb="org.Hs.eg.db")
head(Gene_ID_epi)
resSig_all$ENTREZID = Gene_ID_epi$ENTREZID[match(resSig_all$gene , Gene_ID_epi$SYMBOL)]
resSig_all = resSig_all[!is.na(resSig_all$ENTREZID),]
ctc_epi_up = resSig_all$ENTREZID[which(resSig_all$cellSource=='CTC'&resSig_all$type=='Epi'&resSig_all$change=='Up-regulation')]
ctc_epi_dn=  resSig_all$ENTREZID[which(resSig_all$cellSource=='CTC'&resSig_all$type=='Epi'&resSig_all$change=='Down-regulation')]
ctc_mes_up=  resSig_all$ENTREZID[which(resSig_all$cellSource=='CTC'&resSig_all$type=='Mes'&resSig_all$change=='Up-regulation')]
ctc_mes_dn=  resSig_all$ENTREZID[which(resSig_all$cellSource=='CTC'&resSig_all$type=='Mes'&resSig_all$change=='Down-regulation')]
duct_epi_up= resSig_all$ENTREZID[which(resSig_all$cellSource=='Duct'&resSig_all$type=='Epi'&resSig_all$change=='Up-regulation')]
duct_epi_dn= resSig_all$ENTREZID[which(resSig_all$cellSource=='Duct'&resSig_all$type=='Epi'&resSig_all$change=='Down-regulation')]
duct_mes_up= resSig_all$ENTREZID[which(resSig_all$cellSource=='Duct'&resSig_all$type=='Mes'&resSig_all$change=='Up-regulation')]
duct_mes_dn= resSig_all$ENTREZID[which(resSig_all$cellSource=='Duct'&resSig_all$type=='Mes'&resSig_all$change=='Down-regulation')]
Gene_ID_list = list(ctc_epi = c(ctc_epi_up,ctc_epi_dn),ctc_mes = c(ctc_mes_up,ctc_mes_dn),
                    duct_epi = c(duct_epi_up,duct_epi_dn),duct_mes =c(duct_mes_up,duct_mes_dn))
names(Gene_ID_list) = c('ctc_epi','ctc_mes',
                           'duct_epi','duct_mes')
ck <- compareCluster(geneCluster = Gene_ID_list,
                     pvalueCutoff=1,OrgDb = org.Hs.eg.db,
                     fun = enrichGO)#
edox = setReadable(ck, 'org.Hs.eg.db', 'ENTREZID')

pdf('../../Figures/emtSig_cnet_geneGroup_GO.pdf',width = 10,height = 10)
p1= cnetplot(edox,  colorEdge = TRUE) 
print(p1)
dev.off()
p1$data$name
write.csv(edox@compareClusterResult,'emtSig_cnet_geneGroup_resGO.csv')
write.csv(p1$data,'emtSig_cnet_geneGroup_resGO_cyto.csv')
saveRDS(edox,'emtSig_cnet_geneGroup_GO.rds')
