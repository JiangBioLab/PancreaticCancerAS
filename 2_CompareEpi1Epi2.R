rm(list = ls())
#Epi1 vs Epi2####
#!/usr/local/bin/R --vanilla
setwd('/home/lxxiao/xiaolixing/pdac/10x/')

.libPaths(c("/home/lxxiao/R/x86_64-pc-linux-gnu-library/4.0","/usr/local/lib64/R/library"))
.libPaths()
#multiprocess####
library(KEGG.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
library(ggplot2)
library(ggpie)
load('scefinal.cluster34.res4.Rdata')
sce = subset(sce_change,selfCelltype %in%c('Epithelial1','Epithelial2'))
Idents(sce) = factor(Idents(sce),levels = c('Epithelial1','Epithelial2'))
Epithelial1_names = names(sce$selfCelltype)[which(sce$selfCelltype=='Epithelial1')]
Epithelial2_names = names(sce$selfCelltype)[which(sce$selfCelltype=='Epithelial2')]
table(sce$selfCelltype)
length(Epithelial1_names)
length(Epithelial2_names)

markers_df <- FindMarkers(object = sce, 
                          logfc.threshold = 0.25,ident.1 = 'Epithelial1',ident.2 = 'Epithelial2',
                          slot = "counts",
                          min.pct = 0.1)
head(markers_df)
markers_df_sig = subset(markers_df,p_val_adj<0.05)
markers_df_sig = subset(markers_df_sig,abs(avg_log2FC)>2)
markers_df_sig$gene =rownames(markers_df_sig)
write.csv(markers_df_sig,'Epi1_vs_Epi2_DEG.csv')

#goenrich
markers_df_sig = read.csv('markers_df_sig',header = T,row.names = 1)
sig_gene <- data.frame(SYMBOL=markers_df_sig$gene,
                    group=ifelse(markers_df_sig$avg_log2FC>0,'Epithelial1','Epithelial2'))
table(group$group)
Gene_ID <- bitr(sig_gene$SYMBOL, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
sig_gene  <- merge(Gene_ID,sig_gene,by='SYMBOL')
head(sig_gene)
go_enrich = compareCluster(
  ENTREZID~group, 
  data=sig_gene, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
pdf('../../Figures/Epi1_vs_Epi2_GOenrich.pdf',width = 4,height = 6)
pgo = dotplot(data_GO_sim, showCategory=5,font.size = 8)
print(pgo)
dev.off()

gene1 = subset(sig_gene,group =='Epithelial1' )
gene2 = subset(sig_gene,group =='Epithelial2' )
gene_top10 = c(gene2$SYMBOL[1:5],gene1$SYMBOL[1:5])
pdf('../../Figures/Epi1_vs_Epi2_DEGtop10.pdf',width = 6,height = 6)
pdot = DotPlot(sce, features = unique(gene_top10))+ coord_flip()+
  scale_color_gradient(low = 'grey70',high = 'red')
print(pdot)
dev.off()
table(sce$selfCelltype,sce$cell_type2)

pie_data = sce@meta.data
pie_data$selfCelltype = factor(pie_data$selfCelltype,levels = c('Epithelial1','Epithelial2'))

epi1_piedata = subset(pie_data,selfCelltype=='Epithelial1')
epi2_piedata = subset(pie_data,selfCelltype=='Epithelial2')
colnames(pie_data)
pdf('../../Figures/Epi1_vs_Epi2_cellType_pie.pdf',width = 6,height = 4)
epi1_pie=ggpie(
  data = epi1_piedata, group_key = "cell_type", count_type = "full",
  label_info = "ratio",label_type = "circle",label_size = 4,label_pos = "out"
)
epi2_pie=ggpie(
  data = epi2_piedata, group_key = "cell_type", count_type = "full",
  label_info = "ratio",label_type = "circle",label_size = 4,label_pos = "out"
)
print(epi1_pie+epi2_pie)
dev.off()

epiNorm_piedata = subset(pie_data,cell_type=='Norm')
epi1_norm = subset(epi1_piedata,cell_type=='Norm')
epi2_norm = subset(epi2_piedata,cell_type=='Norm')

#ductal
table(epiNorm_piedata$selfCelltype,epiNorm_piedata$cell_type2)
pvalue = phyper(138,145,173,241,lower.tail = F)##3.83445e-16
qvalue <- p.adjust(pvalue,method='fdr')
#ductal
pvalue = phyper(89,126,192,241,lower.tail = F)
qvalue <- p.adjust(pvalue,method='fdr')
pdf('../../Figures/Epi1_vs_Epi2_norm_bar.pdf',width = 6,height = 4)

pbar = ggplot(epiNorm_piedata, aes(cell_type2))+
  geom_bar(aes(fill=selfCelltype), position = position_dodge(preserve = 'single'))+
  scale_fill_manual(values = c("Epithelial1"  =  '#0B68AB',"Epithelial2"  =  "#8fbd7a"))+
  theme_classic()
print(pbar)
dev.off()

table(epiNorm_piedata$selfCelltype=="Epithelial1",epiNorm_piedata$cell_type2)


pdf('../../Figures/Epi1_vs_Epi2_norm_pie.pdf',width = 6,height = 4)
epi1_pnorm=ggpie(
  data = epi1_norm, group_key = "cell_type2", count_type = "full",
  label_info = "ratio",label_type = "circle",label_size = 4,label_pos = "out"
)
epi2_pnorm=ggpie(
  data = epi2_norm, group_key = "cell_type2", count_type = "full",
  label_info = "ratio",label_type = "circle",label_size = 4,label_pos = "out"
)
print(epi1_pnorm+epi2_pnorm)
dev.off()

head(epi1_norm)
pdf('../../Figures/Epi1_vs_Epi2_Technology_pie.pdf',width = 6,height = 4)
epi1_pTechnology=ggpie(
  data = epi1_piedata, group_key = "Technology", count_type = "full",
  label_info = "ratio",label_type = "circle",label_size = 4,label_pos = "in"
)
epi2_pTechnology=ggpie(
  data = epi2_piedata, group_key = "Technology", count_type = "full",
  label_info = "ratio",label_type = "circle",label_size = 4,label_pos = "in"
)
print(epi1_pTechnology+epi2_pTechnology)
dev.off()

smart_epi = subset(pie_data,Technology=='SmartSeq2')
table(smart_epi$selfCelltype)
smart_data = as.data.frame(table(smart_epi$selfCelltype))
smart_data$Var1 = as.factor(smart_data$Var1)
pdf('../../Figures/Epi1_vs_Epi2_Smartseq2_bar.pdf')
psmr = ggplot(smart_data, aes( Var1, Freq))+geom_col()
  theme_classic()
print(psmr)
dev.off()