rm (list = ls())
setwd('/home/lxxiao/xiaolixing/ctc/data/')
library(ggplot2)
library(stringr)
library(Seurat)
library(dplyr)
library(reshape2)
library(pheatmap)
list.files()

#patient emt score####
ctc_sce = readRDS('pdac_1_GSE144561/pdac1_expr_filted.rds')
ctc_sce = ctc_sce[emt_signature$Gene,]
table(ctc_sce$patient)
ctc_met = subset(ctc_sce,patient == 'metastatic PDAC')
ctc_loc = subset(ctc_sce,patient == 'localized PDAC')
ctc_norm = subset(ctc_sce,patient == 'healthy donor')
#ctc
#ctc pdac_vs_hd EMTdeg####
markers_df_loc <- FindMarkers(object = sce_file, 
                              ident.1 = 'localized PDAC', 
                              ident.2 = 'healthy donor',
                              logfc.threshold = 1,  # logfc.threshold:类群中基因的平均表达量相对于所有其他类群的平均表达量的最小log2倍数
                              slot = "counts",
                              min.pct = 0.25)
markers_df_loc$gene = rownames(markers_df_loc)
markers_df_loc$cluster = paste0('loc_vs_hd')

markers_df_met <- FindMarkers(object = sce_file, 
                              ident.1 = 'metastatic PDAC', 
                              ident.2 = 'healthy donor',
                              logfc.threshold = 1,
                              slot = "counts",
                              min.pct = 0.25)
dim(markers_df_met)
print(x = head(markers_df_met))
markers_df_met$gene = rownames(markers_df_met)
markers_df_met$cluster = paste0('met_vs_hd')

markers_df = rbind(markers_df_loc,markers_df_met)#128

markers_df$gene_EMT_type = emt_signature$Category[match(markers_df$gene,emt_signature$Gene)]
print(x = head(markers_df))
table(markers_df$gene_EMT_type)
foldChange = 1
padj =0.05
ctc_degs$type <- ifelse( abs(ctc_degs$avg_logFC)<foldChange,'Non',
                          ifelse(ctc_degs$avg_logFC>foldChange,'Up','Down'))
table(ctc_degs$type)

ctc_degs= markers_df
ctc_degs = subset(markers_df,markers_df$type != 'Non')
#markers_df  = read.csv('pdac_1_GSE144561/fx_pdac_vs_hd_degs.csv')
#venn####
write.csv(ctc_degs,'pdac_1_GSE144561/fx_pdac_vs_hd_degs.csv')
emt_signature = read.csv('/home/lxxiao/xiaolixing/pdac/10x/emt_gene_signature.csv',header = T,row.names = 1)
table(ctc_degs$cluster)

#sc_emt_deg####
library(stringi)
sc_sce= readRDS('/home/lxxiao/xiaolixing/pdac/10x/10xsmart_cancer_norm_CNV_rmSMART_fdup.rds')
sce.epi.tumor.norm.fdup= readRDS('/home/lxxiao/xiaolixing/pdac/10x/10xsmart_cancer_norm_CNV_rmSMART_fdup.rds')
emt_signature = read.csv('/home/lxxiao/xiaolixing/pdac/10x/emt_gene_signature.csv',header = T,row.names = 1)

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
markers_df_loc$gene = rownames(markers_df_loc)
markers_df_loc$cluster = paste0('loc_vs_hd')

markers_df_met <- FindMarkers(object = sce_file, 
                              ident.1 = 'Epi_met_Tumor', 
                              ident.2 = 'Epi_Norm',
                              logfc.threshold = 1,
                              slot = "counts",
                              min.pct = 0.25)
print(x = head(markers_df_met))
markers_df_met$gene = rownames(markers_df_met)
markers_df_met$cluster = paste0('met_vs_hd')

markers_df = rbind(markers_df_loc,markers_df_met)#LOG2FC差异的是312个，这里画火山图logfc=0是366个
markers_df$gene_EMT_type = emt_signature$Category[match(markers_df$gene,emt_signature$Gene)]

foldChange = 1
padj =0.05
markers_df$type <- ifelse( abs(markers_df$avg_logFC)<foldChange,'Non',
                          ifelse(markers_df$avg_logFC>foldChange,'Up','Down'))
table(markers_df$type)

sc_degs= markers_df
sc_all_sig = subset(markers_df,markers_df$type != 'Non')
write.csv(sc_degs,'/home/lxxiao/xiaolixing/pdac/10x/fx_pdac_vs_hd_degs.csv')




all_new_inter = intersect(ctc_all_sig$gene,
                          sc_all_sig$gene)

#venn####
ctc_degs = read.csv('pdac_1_GSE144561/fx_pdac_vs_hd_degs.csv',header = T,row.names = 1)
sc_degs = read.csv('/home/lxxiao/xiaolixing/pdac/10x/fx_pdac_vs_hd_degs.csv',header = T,row.names = 1)
emt_signature = read.csv('/home/lxxiao/xiaolixing/pdac/10x/emt_gene_signature.csv',header = T,row.names = 1)
all_inter= intersect(ctc_degs$gene,sc_degs$gene)
emt_deg_count = readRDS('pdac_1_GSE144561/emt_deg_count.rds')
expr_sc_met  = emt_deg_count$expr_sc_met
expr_sc_loc  =emt_deg_count$expr_sc_loc
expr_ctc_met =emt_deg_count$expr_ctc_met
expr_ctc_loc =emt_deg_count$expr_ctc_loc

emt_cv = function(x,y){
  gene_type = emt_signature$Category[match(rownames(x),emt_signature$Gene)]
  y_n = which(gene_type == y)
  y_x = x[y_n,]
  #sc_x = scale(y_x)
  #sc_x[sc_x>0] = 1
  #sc_x[sc_x!=1] = 0
  #res_x = colMeans(sc_x)
  res_x = colMeans(y_x)
  return(res_x)
}

emt_pearson_data  = data.frame(type = c(rep('Epi_Tumor_met',ncol(expr_sc_met)),
                                        rep('Epi_Tumor_loc',ncol(expr_sc_loc)),
                                        rep('Epi_norm',ncol(expr_sc_norm)),
                                        rep('CTC_met',ncol(expr_ctc_met)),
                                        rep('CTC_loc',ncol(expr_ctc_loc)),
                                        rep('CTC_norm',ncol(expr_ctc_norm))),
                               Epi_cv = c(emt_cv(expr_sc_met[all_inter,],'Epi'),#67 intersect_gene
                                          emt_cv(expr_sc_loc[all_inter,],'Epi'),
                                          emt_cv(expr_sc_norm[all_inter,],'Epi'),
                                          
                                          emt_cv(expr_ctc_met[all_inter,],'Epi'),
                                          emt_cv(expr_ctc_loc[all_inter,],'Epi'),
                                          emt_cv(expr_ctc_norm[all_inter,],'Epi')) ,
                               
                               Mes_cv =c(emt_cv(expr_sc_met[all_inter,],'Mes'),
                                         emt_cv(expr_sc_loc[all_inter,],'Mes'),
                                         emt_cv(expr_sc_norm[all_inter,],'Mes'),
                                         
                                         emt_cv(expr_ctc_met[all_inter,],'Mes'),
                                         emt_cv(expr_ctc_loc[all_inter,],'Mes'),
                                         emt_cv(expr_ctc_norm[all_inter,],'Mes')))
head(emt_pearson_data)
emt_pearson_data$patient = rownames(emt_pearson_data)
emt_pearson_data$emt = log10(emt_pearson_data$Epi_cv/emt_pearson_data$Mes_cv)
emt_pearson_data$loc_met_type = ifelse(emt_pearson_data$type %in% c('CTC_met','Epi_Tumor_met'),'met',
                                       ifelse(emt_pearson_data$type %in% c('CTC_loc','Epi_Tumor_loc'), 'loc','Norm'))

emt_pearson_data$cell_type = str_split(emt_pearson_data$type,'[_]',simplify = T)[,1]
emt_pearson_data$emt_type = ifelse(emt_pearson_data$emt< -0.37,'Mesenchymal',
                                   ifelse(emt_pearson_data$emt> 0.37,'Epithelial',
                                          'pEMT'))
table(emt_pearson_data_reRun$type,emt_pearson_data_reRun$emt_type)
table(emt_pearson_data_old$type,emt_pearson_data_old$emt_type)

plot_data = emt_pearson_data[order(emt_pearson_data$emt,decreasing = T),]
plot_data$patient = factor(plot_data$patient,levels = plot_data$patient)
quantile(plot_data$emt[which(plot_data$cell_type=='CTC')])
quantile(emt_pearson_data$emt[which(emt_pearson_data$cell_type=='CTC')],c(0.25,0.75))
quantile(emt_pearson_data$emt)

#plot_data = emt_pearson_data
#table(plot_data$type)
#plot_data$type = factor(plot_data$type,levels = c('Epi_norm','Epi_Tumor_loc','CTC_norm','CTC_loc','CTC_met','Epi_Tumor_met' ))
#plot_data = plot_data[order(plot_data$type),]
#plot_data$patient = factor(plot_data$patient,levels = plot_data$patient)
temp_plot_data = plot_data
ggplot(plot_data,aes(x = patient,y=emt,fill =type))+
  geom_bar(position = "stack", stat = "identity", width = 1)+
  #coord_flip()+
  geom_hline(yintercept = c(0.37,-0.37),lty="dashed")+
  xlab('')+theme_classic()+
  scale_fill_manual(values = c(
    'Epi_Tumor_loc'="#F8C77D",
    'Epi_Tumor_met'='#DC5B64',
    'Epi_norm'='#755032',
    "CTC_met" ="#5CA4DA"  ,
    'CTC_loc' ='#87C46A' ,
    'CTC_norm' ='#8D709B') )+
  theme(   panel.background = element_blank(),
           panel.border = element_blank(),
           panel.grid = element_blank(), 
           axis.title = element_text(color='black', size=1),
           axis.ticks.length = unit(0.1, "lines"),
           #axis.ticks = element_blank(),
           axis.line = element_line(color = 'black',size = 0.2),
           axis.text.x = element_text(vjust = 0.5,hjust = 0.4),
           #legend.title = element_blank(),
           legend.text = element_text(size=10),
           legend.key = element_blank(),
           legend.key.size = unit(0.5, 'cm'))



#GO####
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)

#R.utils::setOption("clusterProfiler.download.method",'libcurl')
options(clusterProfiler.download.method = "libcurl")
R.utils::setOption("clusterProfiler.download.method","auto")

ctc_degs = read.csv('pdac_1_GSE144561/fx_pdac_vs_hd_degs.csv',header = T,row.names = 1)
sc_degs = read.csv('/home/lxxiao/xiaolixing/pdac/10x/fx_pdac_vs_hd_degs.csv',header = T,row.names = 1)
emt_signature = read.csv('/home/lxxiao/xiaolixing/pdac/10x/emt_gene_signature.csv',header = T,row.names = 1)
table(ctc_degs$cluster)
ctc_degs_unique = ctc_degs[!duplicated(ctc_degs$gene),]
sc_degs_unique = sc_degs[!duplicated(sc_degs$gene),]

inter_gene = intersect(ctc_degs_unique$gene,sc_degs_unique$gene)
ctc_unique = ctc_degs_unique$gene[-which(ctc_degs_unique$gene %in% inter_gene)]
sc_unique = sc_degs_unique$gene[-which(sc_degs_unique$gene %in% inter_gene)]

group <- data.frame(gene=c(inter_gene,ctc_unique,sc_unique),
                    group=c(rep('inter_gene',length(inter_gene)),
                            rep('ctc_unique',length(ctc_unique)),
                            rep('sc_unique',length(sc_unique))))
Gene_ID <- bitr(group$gene, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
#write.csv(data,'pdac_1_GSE144561/reRun_GO_data_input.csv')

