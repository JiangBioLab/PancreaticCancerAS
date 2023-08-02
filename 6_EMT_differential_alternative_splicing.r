rm (list = ls())
setwd('/home/lxxiao/xiaolixing/ctc/data/')
library(stringr)
library(ggplot2)

#pEMT_Mes####
source('/home/lxxiao/xiaolixing/pdac/result_dir/code/rmats_diff_func.R')
source('/home/lxxiao/xiaolixing/pdac/result_dir/code/rmats_allE_DaseFinder.R')

pEMT_pdacSC = read.csv('pdac_1_GSE144561/1.mapping/rmats_pEMTSC/SE.MATS.JC.txt',
                       header = T,row.names = 1,sep = '\t')
mes_pdacSC = read.csv('pdac_1_GSE144561/1.mapping/rmats_mesSC//SE.MATS.JC.txt',
                      header = T,row.names = 1,sep = '\t')

bam_pEMT = read.table('pdac_1_GSE144561/rmats_file/pdac1_pEMT.txt',sep = ',')
bam_mes = read.table('pdac_1_GSE144561/rmats_file/pdac1_mes.txt',sep = ',')
#bam_sc = read.table('pdac_1_GSE144561/rmats_file/sc_pdac_hd.txt',sep = ',')
bam_sc = read.table('pdac_1_GSE144561/rmats_file/sc_pdac_hd_reRun.txt',sep = ',')
bam_sc[1,] = str_split(bam_sc[1,],'[/]',simplify = T)[,9]

pEMT_SE= dase(pEMT_pdacSC,bam1 = bam_pEMT[1,],bam2 = bam_sc[1,],
              group1 = "pEMT",group2 = "HD",na_percent = 0.3,method = 'mean')
colnames(pEMT_SE$FilByPadj)
dim(pEMT_SE$FilByPadj)
data_plot = pEMT_SE$FilByPadj[,24:29]
pheatmap::pheatmap(data_plot,show_rownames = F,show_colnames = F)
mes_SE = dase(rmats_result_data = mes_pdacSC,bam1 = bam_mes[1,],bam2 = bam_sc[1,],
              group1 = "pEMT",group2 = "HD",na_percent = 0.3,method = 'mean')

source('/home/lxxiao/xiaolixing/pdac/result_dir/code/rmats_allE_DaseFinder.R')

dir_path = 'pdac_1_GSE144561/1.mapping/'
#lochd_all = as_type_result(dir_path,'rmats_mes_LocSC_v4/',group1 = 'Loc',group2 = 'HD',na_percent = 0.3,method = 'Mean')
#methd_all = as_type_result(dir_path,'rmats_mes_MetSC_v4/',group1 = 'Met',group2 = 'HD',na_percent = 0.3,method = 'Mean')

pEMThd_all = as_type_result(dir_path,dir_name = 'rmats_pEMTSC/',
                            bam1 = bam_pEMT[1,],bam2 = bam_sc[1,],
                            group1 = 'pEMT',group2 = 'HD',
                            na_percent = 0.3,method = 'Mean')
meshd_all = as_type_result(dir_path,'rmats_mesSC/',
                           bam1 = bam_mes[1,],bam2 = bam_sc[1,],
                           group1 = 'Mes',group2 = 'HD',
                           na_percent = 0.3,method = 'Mean')
#rmats_all_event####

#source('/home/lxxiao/xiaolixing/pdac/result_dir/code/rmats_allE_DaseFinder.R')
#dir_path = 'pdac_1_GSE144561/1.mapping/'
#files_list = list.files(paste0(dir_path,dir_name)) 
#lochd_all = as_type_result(dir_path,'rmats_mes_LocHD_v3/')
#methd_all = as_type_result(dir_path,'rmats_mes_MetHD_v3/')
#metLoc_all =as_type_result(dir_path,'rmats_mes_MetLoc/')
#
#
#lochd_all = as_type_result(dir_path,'rmats_mes_LocHD_v3/')
#methd_all = as_type_result(dir_path,'rmats_mes_MetSC/')
#metLoc_all =as_type_result(dir_path,'rmats_mes_MetLoc/')

#rmats_result_save####
rmats_mes_DifEventAll = list(pEMThd_all=pEMThd_all,
                             meshd_all=meshd_all)

rmats_mes_DifSEResult  = list(pEMTHD_dase=pEMThd_all$SE,
                              mesHD_dase = meshd_all$SE)

rmats_mes_SEPadjResult  = list(pEMThd_dase_Padj=pEMThd_all$SE$FilByPadj,
                               meshd_dase_Padj=meshd_all$SE$FilByPadj)
#saveRDS(rmats_mes_DifEventAll,'pdac_1_GSE144561/rmats_file/rmats_mes_DifEventAll.rds')
#saveRDS(rmats_mes_DifSEResult,'pdac_1_GSE144561/rmats_file/rmats_mes_DifSE.rds')
#saveRDS(rmats_mes_SEPadjResult,'pdac_1_GSE144561/rmats_file/rmats_mes_DifSEPadj.rds')
rmats_mes_DifEventAll= readRDS('pdac_1_GSE144561/rmats_file/rmats_mes_DifEventAll.rds')
rmats_mes_DifSEResult = readRDS('pdac_1_GSE144561/rmats_file/rmats_mes_DifSE.rds')
rmats_mes_SEPadjResult = readRDS('pdac_1_GSE144561/rmats_file/rmats_mes_DifSEPadj.rds')
GESP = read.csv('pdac_1_GSE144561/hunman_cell_surface_pro.csv',header = T)#TCSA Download
head(GESP)
emt_signature = read.csv('/home/lxxiao/xiaolixing/pdac/10x/emt_gene_signature.csv',header = T,row.names = 1)


#all_event_summarize###
sapply(rmats_mes_DifEventAll$pEMThd_all, function(x){dim(x[[3]])})[1,]
sapply(rmats_mes_DifEventAll$meshd_all, function(x){dim(x[[3]])})[1,]

all_event_summarize = data.frame(pEMT_vs_hd = sapply(rmats_mes_DifEventAll$pEMThd_all, function(x){dim(x[[3]])})[1,],
                                 mes_vs_hd = sapply(rmats_mes_DifEventAll$meshd_all, function(x){dim(x[[3]])})[1,],
                                 row.names = c('A3SS','A5SS','MXE','RI','SE'))




#event_inter####
temp_as_pEMT = sapply(rmats_mes_DifEventAll$pEMThd_all,function(x){x[[3]]})
temp_as_mes  = sapply(rmats_mes_DifEventAll$meshd_all,function(x){x[[3]]})

#table save
#colnames(temp_as_pEMT$SE)
#temp_as_pEMT_info = list(A3SS = temp_as_pEMT$A3SS[,c(1:23,42:50),],
#                         A5SS = temp_as_pEMT$A5SS[,c(1:23,43:51),],
#                         MXE = temp_as_pEMT$MXE[,c(1:25,42:50),],
#                         RI = temp_as_pEMT$RI[,c(1:24,43:51),],
#                         SE = temp_as_pEMT$SE[,c(1:23,40:48),]
#                         )
#write.csv(temp_as_pEMT_info$A3SS,'pdac_1_GSE144561/as_pEMT_A3SS.csv')
#write.csv(temp_as_pEMT_info$A5SS,'pdac_1_GSE144561/as_pEMT_A5SS.csv')
#write.csv(temp_as_pEMT_info$MXE,'pdac_1_GSE144561/as_pEMT_MXE.csv')
#write.csv(temp_as_pEMT_info$RI,'pdac_1_GSE144561/as_pEMT_RI.csv')
#write.csv(temp_as_pEMT_info$SE,'pdac_1_GSE144561/as_pEMT_SE.csv')
#
#
#colnames(temp_as_mes$SE)
#temp_as_mes_info = list(A3SS = temp_as_mes$A3SS[,c(1:23,62:70),],
#                        A5SS = temp_as_mes$A5SS[,c(1:23,63:71),],
#                        MXE =  temp_as_mes$MXE[,c(1:25,61:69),],
#                        RI =   temp_as_mes$RI[,c(1:24,69:77),],
#                        SE =   temp_as_mes$SE[,c(1:23,59:67),])
#
#write.csv(temp_as_mes_info$A3SS,'pdac_1_GSE144561/as_mes_A3SS.csv')
#write.csv(temp_as_mes_info$A5SS,'pdac_1_GSE144561/as_mes_A5SS.csv')
#write.csv(temp_as_mes_info$MXE,'pdac_1_GSE144561/as_mes_MXE.csv')
#write.csv(temp_as_mes_info$RI,'pdac_1_GSE144561/as_mes_RI.csv')
#write.csv(temp_as_mes_info$SE,'pdac_1_GSE144561/as_mes_SE.csv')


rmats_mes_diff = list(pEMTHD_dase_Padj = temp_as_pEMT,
                      MesHD_dase_Padj = temp_as_mes)



n_gene_pEMTHD = sapply(temp_as_pEMT,function(x){x[,'ENSEMBL.Gene.ID']})
n_gene_MesHD = sapply(temp_as_mes,function(x){x[,'ENSEMBL.Gene.ID']})

all_gene_summarize = data.frame(pEMT_vs_hd = sapply(n_gene_pEMTHD, function(x){length(unique(x))}),
                                mes_vs_hd = sapply(n_gene_MesHD, function(x){length(unique(x))}),
                                row.names = c('A3SS','A5SS','MXE','RI','SE'))
mut_as_gene = (all_event_summarize-all_gene_summarize)/all_gene_summarize
colnames(mut_as_gene) = c('pEMT_mult','mes_mult')
single_as_gene = 1-mut_as_gene
colnames(single_as_gene) = c('pEMT_single','mes_single')
as_gene = cbind(mut_as_gene,single_as_gene)
as_gene$as_type = rownames(as_gene)


library(reshape2)
colSums(all_event_summarize)
bar_data = all_event_summarize
bar_data$group = rownames(bar_data)
bar_data = melt(bar_data,'group')

library(RColorBrewer)
display.brewer.all()
col = brewer.pal(2,'YlOrBr')
mycolor = colorRampPalette(col)
mycolor(5)

ggplot(bar_data,aes(x = value,y=variable,fill =group))+
  geom_bar(position = "stack", stat = "identity", width = 0.7)+
  xlab('')+ylab('')+theme_classic()+
  #scale_fill_manual(values = mycolor(5))
  scale_fill_manual(values = c(
    'A3SS'="#dde5b6", 'A5SS'="#aaa1c8",'MXE'= "#81c3d7" ,
    'RI'="#967aa1", 'SE'="#ee6c4d" ))

all_event_summarize = apply(all_event_summarize, 2, function(x){round(100*x/sum(x),1)})

#DEG&DASE venn####
EMT_DEG = readRDS('pdac_1_GSE144561/emt_deg_result.rds')
temp_as_pEMT = sapply(rmats_mes_DifEventAll$pEMThd_all,function(x){x[[3]]})
temp_as_mes  = sapply(rmats_mes_DifEventAll$meshd_all,function(x){x[[3]]})

rmats_mes_diff = list(pEMTHD_dase_Padj = temp_as_pEMT,
                      MesHD_dase_Padj = temp_as_mes)

n_gene_pEMTHD = sapply(temp_as_pEMT,function(x){x[,'geneSymbol']})
n_gene_MesHD = sapply(temp_as_mes,function(x){x[,'geneSymbol']})

dase_gene_names = list(gene_pEMTHD = unique(unlist(n_gene_pEMTHD)),
                       gene_MesHD = unique(unlist(n_gene_MesHD)))
length(dase_gene_names$gene_pEMTHD)#598
length(dase_gene_names$gene_MesHD)#350

deg_gene_names = list(gene_pEMTHD = unique(rownames(EMT_DEG$pEMThd_deg)),
                      gene_MesHD = unique(rownames(EMT_DEG$Meshd_deg)))

library(ggvenn)
library(tidyverse)
venn_pEMT = list(DASE_pEMT=dase_gene_names$gene_pEMTHD,
                 DEG_pEMT =deg_gene_names$gene_pEMTHD)
ggvenn(venn_pEMT,show_percentage = T,show_elements = F,label_sep = ",",
       digits = 1,stroke_color = "white",
       fill_color = c("#E41A1C", "#1E90FF", "#FF8C00"),
       set_name_color = c("#E41A1C", "#1E90FF"))#,"#FF8C00","#984EA3", "#4DAF4A", "#984EA3"

venn_Mes = list(DASE_Mes=dase_gene_names$gene_MesHD,
                DEG_Mes =deg_gene_names$gene_MesHD)
ggvenn(venn_Mes,show_percentage = T,show_elements = F,label_sep = ",",
       digits = 1,stroke_color = "white",
       fill_color = c("#E41A1C", "#1E90FF", "#FF8C00"),
       set_name_color = c("#E41A1C", "#1E90FF"))#,"#FF8C00","#984EA3", "#4DAF4A", "#984EA3"
ggsave('pdac_1_GSE144561/figures/venn_AS_Exp_pEMT.pdf')
ggsave('pdac_1_GSE144561/figures/venn_AS_Exp_Mes.pdf')

venn_gene_Mes = intersect(venn_Mes$DASE_Mes,venn_Mes$DEG_Mes)
#write.table(t(venn_gene_Mes),quote = F,row.names = F,col.names = F,sep = ' ')

venn_gene_pEMT = intersect(venn_pEMT$DASE_pEMT,venn_pEMT$DEG_pEMT)
#write.table(t(venn_gene_pEMT),quote = F,row.names = F,col.names = F,sep = ' ')

mut_gene_summarize = list(pEMT_vs_hd =lapply(n_gene_pEMTHD, function(x){names(table(x))[which(table(x)!=1)]}),
                          mes_vs_hd = lapply(n_gene_MesHD, function(x){names(table(x))[which(table(x)!=1)]}))

sapply(mut_gene_summarize$pEMT_vs_hd, function(x){length(unique(x))})
mut_pEMT_data = data.frame(geneSymbol = unlist(mut_gene_summarize$pEMT_vs_hd),
                           asType =c(rep('A3SS',2),rep('A5SS',0),rep('MXE',27),rep('RI',8),rep('SE',84)),
                           emt_type = rep('pEMT',length(unlist(mut_gene_summarize$pEMT_vs_hd))))

sapply(mut_gene_summarize$mes_vs_hd, function(x){length(unique(x))})
mut_Mes_data = data.frame(geneSymbol = unlist(mut_gene_summarize$mes_vs_hd),
                          asType =c(rep('A3SS',3),rep('A5SS',2),rep('MXE',18),rep('RI',3),rep('SE',36)),
                          emt_type = rep('pEMT',length(unlist(mut_gene_summarize$mes_vs_hd))))

mut_data = rbind(mut_pEMT_data,mut_Mes_data) 

mut_as_gene = (all_event_summarize-all_gene_summarize)/all_gene_summarize
colnames(mut_as_gene) = c('pEMT_mult','mes_mult')
single_as_gene = 1-mut_as_gene
colnames(single_as_gene) = c('pEMT_single','mes_single')
as_gene = cbind(mut_as_gene,single_as_gene)
as_gene$as_type = rownames(as_gene)

#Manhattan Plot####
library(CMplot)
temp_as_pEMT = sapply(rmats_mes_DifEventAll$pEMThd_all,function(x){x[[3]]})
temp_as_mes  = sapply(rmats_mes_DifEventAll$meshd_all,function(x){x[[3]]})

rmats_mes_diff = list(pEMTHD_dase_Padj = temp_as_pEMT,
                      MesHD_dase_Padj = temp_as_mes)
cm_as_pEMT = lapply(temp_as_pEMT, function(x){x[,c(1:4,which(colnames(x) %in% c('delt','as_type')))]})
cm_as_pEMT = do.call(rbind,cm_as_pEMT)
head(cm_as_pEMT)
write.csv(cm_as_pEMT,'pdac_1_GSE144561/CMplot_as_pEMT.csv')

cm_as_mes  = lapply(temp_as_mes, function(x){x[,c(1:4,which(colnames(x) %in% c('delt','as_type')))]})
cm_as_mes  = do.call(rbind,cm_as_mes)
write.csv(cm_as_mes,'pdac_1_GSE144561/CMplot_as_Mes.csv')

#AS_event_inter####
temp_as_pEMT = sapply(rmats_mes_DifEventAll$pEMThd_all,function(x){x[[3]]})
temp_as_mes = sapply(rmats_mes_DifEventAll$meshd_all,function(x){x[[3]]})

temp_a3ss_pEMT = temp_as_pEMT$A3SS
temp_a3ss_mes  = temp_as_mes$A3SS

temp_a5ss_pEMT = temp_as_pEMT$A5SS
temp_a5ss_mes  = temp_as_mes$A5SS

temp_mxe_pEMT = temp_as_pEMT$MXE
temp_mxe_mes  = temp_as_mes$MXE

temp_ri_pEMT = temp_as_pEMT$RI
temp_ri_mes  = temp_as_mes$RI

temp_se_pEMT = temp_as_pEMT$SE
temp_se_mes = temp_as_mes$SE

#motif####

non_dif_pEMT = rmats_mes_DifEventAll$pEMThd_all$SE$result_all[rmats_mes_DifEventAll$pEMThd_all$SE$result_all$padj>0.05,3:10]
se_dif_pEMT_up   = rmats_mes_DifEventAll$pEMThd_all$SE$FilByPadj[rmats_mes_DifEventAll$pEMThd_all$SE$FilByPadj$delt>0, 3:10]
se_dif_pEMT_down = rmats_mes_DifEventAll$pEMThd_all$SE$FilByPadj[rmats_mes_DifEventAll$pEMThd_all$SE$FilByPadj$delt<0, 3:10]
colnames(non_dif_pEMT)
#write.table(non_dif_pEMT,col.names = T,row.names = F,sep = '\t',quote = F,'pdac_1_GSE144561/non_dif_pEMT.txt')
#write.table(se_dif_pEMT_up,col.names = T,row.names = F,sep = '\t',quote = F,'pdac_1_GSE144561/se_dif_pEMT_up.txt')
#write.table(se_dif_pEMT_down,col.names = T,row.names = F,sep = '\t',quote = F,'pdac_1_GSE144561/se_dif_pEMT_down.txt')


non_dif_Mes = rmats_mes_DifEventAll$meshd_all$SE$result_all[rmats_mes_DifEventAll$meshd_all$SE$result_all$padj>0.05,3:10]
se_dif_Mes_up    = rmats_mes_DifEventAll$meshd_all$SE$FilByPadj[rmats_mes_DifEventAll$meshd_all$SE$FilByPadj$delt>0,3:10]
se_dif_Mes_down = rmats_mes_DifEventAll$meshd_all$SE$FilByPadj[rmats_mes_DifEventAll$meshd_all$SE$FilByPadj$delt<0,3:10]
#write.table(non_dif_Mes,col.names = T,row.names = F,sep = '\t',quote = F,'pdac_1_GSE144561/non_dif_Mes.txt')
#write.table(se_dif_Mes_up,col.names = T,row.names = F,sep = '\t',quote = F,'pdac_1_GSE144561/se_dif_Mes_up.txt')
#write.table(se_dif_Mes_down,col.names = T,row.names = F,sep = '\t',quote = F,'pdac_1_GSE144561/se_dif_Mes_down.txt')

#rMAPS result####
EMT_DEG = readRDS('pdac_1_GSE144561/emt_deg_result.rds')
Mes_up_res = read.table('pdac_1_GSE144561/motif/Mes_resultMaps/pVal.up.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')
Mes_dn_res = read.table('pdac_1_GSE144561/motif/Mes_resultMaps/pVal.dn.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')

pEMT_up_res = read.table('pdac_1_GSE144561/motif/pEMT_resultMaps/pVal.up.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')
pEMT_dn_res = read.table('pdac_1_GSE144561/motif/pEMT_resultMaps/pVal.dn.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')

#plot
library(tidyverse)
library(gganimate)
library(ggplot2)
cols <- c(
  # "smallest_p_in_upstreamExon.3prime" = "#34a186",
  #"smallest_p_in_upstreamExonIntron" = "#f9cb45",
  
  "smallest_p_in_upstreamIntron" = "#b5182b",
  
  "smallest_p_in_targetExon.5prime" = "#4cb1c4",
  "smallest_p_in_targetExon.3prime" = "#ab96d2",
  
  "smallest_p_in_downstreamIntron" =  "#34a186"
  
)
colnames(Mes_up_res)
up_res_data =as.data.frame(ifelse(Mes_up_res<0.01,1,0)) 
rownames(up_res_data)
up_res_data = up_res_data[-grep('^motif.',rownames(up_res_data)),]
up_res_data = up_res_data[-which(rowSums(up_res_data[,3:6])==0),]
pheatmap::pheatmap(up_res_data[,3:6],show_colnames = F,treeheight_row = 0,treeheight_col = 0,
                   color = colorRampPalette(c("white", "#ee6055"))(50))

type = 'Mes'
dot_data = Mes_up_res[,3:6]
dot_data$name = rownames(dot_data)
dot_data$motif = str_split(rownames(dot_data),'[.]',simplify = T)[,2]
dot_data$RBP = str_split(rownames(dot_data),'[.]',simplify = T)[,1]

colnames( EMT_DEG$Meshd_tpm_fil)
dot_data$Mean_tpm=round(rowMeans(EMT_DEG$Meshd_tpm_fil[match(dot_data$RBP,rownames(EMT_DEG$Meshd_tpm_fil)),1:57]),1) 

dot_data = dot_data[-which(is.na(dot_data$Mean_tpm)),]
head(dot_data)

dot_data = melt(dot_data,c('motif','RBP','name','Mean_tpm'))
dot_data$logP = log10(dot_data$value)


ggplot(dot_data,aes(x = variable, 
                    y = motif))+
  geom_point(
    aes(size = logP, fill = variable), 
    pch = 21, color = "white", alpha = .9
  )+scale_fill_manual(values = cols)



#sum_delt_histogram####
delt_data_pEMT = lapply(temp_as_pEMT,
                        function(x){
                          df_delt = x[,c('delt','as_type')]
                          df_delt$group = rep('pEMT',nrow(df_delt))
                          df_delt})
delt_data_pEMT = do.call(rbind,delt_data_pEMT)
table(delt_data_pEMT$delt>0,delt_data_pEMT$as_type)
delt_data_mes = lapply(temp_as_mes,
                       function(x){
                         df_delt = x[,c('delt','as_type')]
                         df_delt$group = rep('Mes',nrow(df_delt))
                         df_delt})
delt_data_mes = do.call(rbind,delt_data_mes)
table(delt_data_mes$delt>0,delt_data_mes$as_type)

data_delt = rbind(delt_data_pEMT,delt_data_mes)
head(data_delt)
ggplot(data = data_delt,aes(x=delt))+
  geom_histogram(aes(fill=group),
                 bins = 50,
                 alpha=0.8)+
  facet_grid(~as_type)+
  scale_fill_manual(values = c("pEMT"="#a3cd5b",
                               "Mes"="#8ea0cc"),
                    labels=c("pEMT"="pEMT",
                             "Mes"="Mes"))+
  theme_bw()+
  theme(axis.title = element_blank())
#ggsave('pdac_1_GSE144561/figures/event_psi_histogram.pdf')

head(data_delt)
table(data_delt$delt>0,data_delt$group)
data_all = as.data.frame(as.matrix(table(data_delt$delt>0,data_delt$group)))
data_all$type = ifelse(data_all$Var1==FALSE,'LESS','OVER')
data_all$type = paste0()
data_all$value = c(-168,332,-285,592)
data_all$Var2 = factor(data_all$Var2,levels = c('pEMT','Mes'))
ggplot(data = data_all,aes(x=Var2,y = Freq,fill = type ))+
  geom_bar(stat="identity",position=position_dodge(0.75))+
  scale_fill_manual(values = c("LESS"="#a3cd5b",
                               "OVER"="#8ea0cc"),
                    labels=c("LESS"="ΔPSI <-0.1",
                             "OVER"="ΔPSI >0.1"))+
  theme_classic()+
  theme(axis.title = element_blank())
#ggsave('pdac_1_GSE144561/figures/all_delt_barplot.pdf')


colnames(temp_se_pEMT)
temp_se_pEMT_psi = temp_se_pEMT[,24:39]
colnames(temp_a3ss_mes)
temp_se_mes_psi = temp_se_mes[,24:61]
anno_col = data.frame(type =str_split(colnames(temp_se_pEMT_psi ),'_',simplify = T)[,1],
                      row.names = colnames(temp_se_pEMT_psi ))
pheatmap::pheatmap(temp_se_pEMT_psi,show_rownames = F,show_colnames = F,cluster_cols = T,
                   color = c('white',rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Reds")))(100))),
                   annotation_col = anno_col,
)

temp_se_mes_psi = temp_se_mes[,24:39]
colnames(temp_a3ss_mes)
temp_se_mes_psi = temp_se_mes[,24:61]
anno_col = data.frame(type =str_split(colnames(temp_se_pEMT_psi ),'_',simplify = T)[,1],
                      row.names = colnames(temp_se_pEMT_psi ))
pheatmap::pheatmap(temp_se_pEMT_psi,show_rownames = F,show_colnames = F,cluster_cols = T,
                   color = c('white',rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Reds")))(100))),
                   annotation_col = anno_col,
)

data_temp = temp_se_pEMT
colnames(data_temp)
data_temp_part1 = data_temp[,1:19]
data_temp_part1$PValue = data_temp$p_self
data_temp_part1$FDR = data_temp$padj


rmats_mes_diff = list(pEMTHD_dase_Padj = temp_as_pEMT,
                      MesHD_dase_Padj  = temp_as_mes)


tmp_event_pEMT =c()
tmp_event_mes = c()
for (i in 1:5) {
  print(i)
  i = 5#SE
  type = names(temp_as_pEMT)[i]
  gene_a = temp_as_pEMT[[i]]
  gene_b = temp_as_mes[[i]]
  inter_gene = intersect(gene_a$geneSymbol,gene_b$geneSymbol)
  length(inter_gene)
  inter_a  = gene_a[gene_a$geneSymbol %in% inter_gene,
                    c("ENSEMBL.Gene.ID","geneSymbol","chr","strand" ,"exonStart_0base" ,"exonEnd",
                      "upstreamES" ,"upstreamEE","downstreamES" , "downstreamEE" ,
                      'delt','as_type','event_num')]
  inter_b  = gene_b[gene_b$geneSymbol %in% inter_gene,
                    c("ENSEMBL.Gene.ID","geneSymbol","chr","strand" ,"exonStart_0base" ,"exonEnd",
                      "upstreamES" ,"upstreamEE","downstreamES" , "downstreamEE" ,
                      'delt','as_type','event_num')]
  dim(inter_a)
  dim(inter_b)
  trans_event = c()
  trans_event_pEMT = c()
  trans_event_mes = c()
  for (z in 1:nrow(inter_a)) {
    print(z)
    
    gsa = inter_a[z,]
    gs = inter_a[z,"geneSymbol"]
    gsb = inter_b[grep(gs,inter_b[,'geneSymbol']),]
    gsb = gsb[which(gsb[,'exonStart_0base']>=gsa[,'upstreamEE']),]
    gsb = gsb[which(gsb[,'exonEnd']<=gsa[,'downstreamES']),]
    #gsb = gsb[grep(gsa[,'exonEnd'],gsb[,'exonEnd']),]
    #gsb = gsb[grep(gsa[,'upstreamES'],gsb[,'upstreamES']),]
    #gsb = gsb[grep(gsa[,'downstreamES'],gsb[,'upstreamES']),]
    #gsb = gsb[grep(gsa[,'downstreamEE'],gsb[,'downstreamEE']),]
    if(dim(gsb)[1]==0){
      next
    }
    
    if (gsa[,'delt']< -0.25 & any(gsb[,'delt']> 0.25) ) {
      gsb2 = gsb[which(gsb[,'delt']>0),]
      event = data.frame(
        ENSEMBL.Gene.ID = gsa[,'ENSEMBL.Gene.ID'],
        geneSymbol = gsa[,'geneSymbol'],
        chr = gsa[,'chr'],
        strand_pEMT = gsa[,'strand'],
        upstreamEE = gsa[,'upstreamEE'],
        downstreamES = gsa[,'downstreamES'],
        
        exonStart_0base_pEMT = gsa[,'exonStart_0base'],
        exonEnd_pEMT = gsa[,'exonEnd'],
        event_id_pEMT = gsa[,'event_num'],
        
        strand_mes = gsb2[,'strand'],
        exonStart_0base_mes = gsb2[,'exonStart_0base'],
        exonEnd_mes = gsb2[,'exonEnd'],
        event_id_mes = gsb2[,'event_num'],
        
        delt_pEMT = gsa[,'delt'],
        delt_mes = gsb2[which(gsb[,'delt']>0),'delt'])
      
      trans_event = rbind(trans_event,event)
      trans_event_pEMT = rbind(trans_event_pEMT,gene_a[gsa[,'event_num'],])
      trans_event_mes = rbind(trans_event_mes,gene_b[gsb2[,'event_num'],])
    }else{
      if (gsa[,'delt']> 0.25 & any(gsb[,'delt']< -0.25) ) {
        gsb2 = gsb[which(gsb[,'delt']<0),]
        event = data.frame(
          ENSEMBL.Gene.ID = gsa[,'ENSEMBL.Gene.ID'],
          geneSymbol = gsa[,'geneSymbol'],
          chr = gsa[,'chr'],
          strand_pEMT = gsa[,'strand'],
          upstreamEE = gsa[,'upstreamEE'],
          downstreamES = gsa[,'downstreamES'],
          
          exonStart_0base_pEMT = gsa[,'exonStart_0base'],
          exonEnd_pEMT = gsa[,'exonEnd'],
          event_id_pEMT = gsa[,'delt'],
          
          strand_mes = gsb2[,'strand'],
          exonStart_0base_mes = gsb2[,'exonStart_0base'],
          exonEnd_mes = gsb2[,'exonEnd'],
          event_id_mes = gsb2[,'event_num'],
          
          delt_pEMT = gsa[,'delt'],
          delt_mes = gsb2[,'delt'])
        trans_event = rbind(trans_event,event)
        trans_event_pEMT = rbind(trans_event_pEMT,gene_a[gsa[,'event_num'],])
        trans_event_mes = rbind(trans_event_mes,gene_b[gsb2[,'event_num'],])
      }else{
        trans_event = trans_event
        trans_event_mes = trans_event_mes
        trans_event_pEMT= trans_event_pEMT
      }
    }
    
  }
  
  tmp_event_pEMT = rbind(tmp_event_pEMT,trans_event_pEMT)
  tmp_event_mes = rbind(tmp_event_mes,trans_event_mes)
}

#plot
i = 5#SE
type = names(temp_as_pEMT)[i]
gene_a = temp_as_pEMT[[i]]
gene_b = temp_as_mes[[i]]
inter_gene = intersect(gene_a$geneSymbol,gene_b$geneSymbol)
length(inter_gene)
inter_a  = gene_a[gene_a$geneSymbol %in% inter_gene,
                  c("ENSEMBL.Gene.ID","geneSymbol","chr","strand" ,"exonStart_0base" ,"exonEnd",
                    "upstreamES" ,"upstreamEE","downstreamES" , "downstreamEE" ,
                    'delt','as_type','event_num')]
inter_b  = gene_b[gene_b$geneSymbol %in% inter_gene,
                  c("ENSEMBL.Gene.ID","geneSymbol","chr","strand" ,"exonStart_0base" ,"exonEnd",
                    "upstreamES" ,"upstreamEE","downstreamES" , "downstreamEE" ,
                    'delt','as_type','event_num')]


dim(inter_a)
dim(inter_b)
trans_event = c()
trans_event_pEMT = c()
trans_event_mes = c()
for (z in 1:nrow(inter_a)) {
  print(z)
  gsa = inter_a[z,]
  gs = inter_a[z,"geneSymbol"]
  gsb = inter_b[grep(gs,inter_b[,'geneSymbol']),]
  #gsb = gsb[which(gsb[,'exonStart_0base']>=gsa[,'upstreamEE']),]
  #gsb = gsb[which(gsb[,'exonEnd']<=gsa[,'downstreamES']),]
  gsb = gsb[grep(gsa[,'exonStart_0base'],gsb[,'exonStart_0base']),]
  gsb = gsb[grep(gsa[,'exonEnd'],gsb[,'exonEnd']),]
  #gsb = gsb[grep(gsa[,'upstreamES'],gsb[,'upstreamES']),]
  #gsb = gsb[grep(gsa[,'downstreamES'],gsb[,'upstreamES']),]
  #gsb = gsb[grep(gsa[,'downstreamEE'],gsb[,'downstreamEE']),]
  
  pEMT = gsa
  mes = gsb
  colnames(pEMT) = paste0(colnames(pEMT),'.pEMT')
  colnames(mes) = paste0(colnames(mes),'.Mes')
  print(paste('dim:',dim(gsb)[1]))
  if(dim(gsb)[1]==0){
    next
  }else{
    if(dim(gsb)[1]==1){
      event = cbind(pEMT,mes[,c('delt.Mes','event_num.Mes',
                                'exonStart_0base.Mes','exonEnd.Mes')]) 
      trans_event = rbind(trans_event,event)
      trans_event_pEMT = rbind(trans_event_pEMT,gene_a[gsa[,'event_num'],])
      trans_event_mes = rbind(trans_event_mes,gene_b[gsb[,'event_num'],])
    }else{
      event = cbind(pEMT,mes[,c('delt.Mes','event_num.Mes',
                                'exonStart_0base.Mes','exonEnd.Mes')]) 
      trans_event = rbind(trans_event,event)
      trans_event_pEMT = rbind(trans_event_pEMT,gene_a[gsa[,'event_num'],])
      trans_event_mes = rbind(trans_event_mes,gene_b[gsb[,'event_num'],])
    }
    
  }
}



library(stringr)
#SE_res_inter####
SE_event_inter = trans_event
SE_event_mes   = trans_event_mes
SE_event_pEMT  = trans_event_pEMT
colnames(SE_event_mes)
colnames(SE_event_pEMT)


library(ggplot2)
library(stringr)
library(ggsignif)

colnames(SE_event_inter)
SE_event_inter[which(SE_event_inter$delt.pEMT>0&SE_event_inter$delt.Mes<0),c(2,3,5,6,13,15)]
SE_event_inter[which(SE_event_inter$delt.pEMT<0&SE_event_inter$delt.Mes>0),c(2,3,5,6)]
#dif_pemt_mes
colnames(SE_event_mes)
colnames(SE_event_pEMT)
dim(SE_event_pEMT)

SE_mes_mat = SE_event_mes[,c(1:6,24:49)]
SE_pEMT_mat = SE_event_pEMT[,c(1:6,24:39)]

SE_mes_mat$event = paste0(SE_mes_mat$geneSymbol,':',SE_mes_mat$chr,':',
                          SE_mes_mat$exonStart_0base,'-',SE_mes_mat$exonEnd,
                          ':',SE_mes_mat$strand)
SE_mes_mat = SE_mes_mat[!duplicated(SE_mes_mat$event),]
SE_pEMT_mat$event = paste0(SE_pEMT_mat$geneSymbol,':',SE_pEMT_mat$chr,':',
                           SE_pEMT_mat$exonStart_0base,'-',SE_pEMT_mat$exonEnd,
                           ':',SE_mes_mat$strand)
SE_pEMT_mat = SE_pEMT_mat[!duplicated(SE_pEMT_mat$event),]



#inter_se_event_Result####
SE_mat = merge(SE_mes_mat,SE_pEMT_mat[,-c(1:6)],by = 'event')
rownames(SE_mat) = SE_mat$event

SE_event_res = list(
  SE_merge = SE_mat,
  SE_event_inter=SE_event_inter,
  SE_event_mes  =SE_event_mes,
  SE_event_pEMT =SE_event_pEMT
)
saveRDS(SE_event_res,'pdac_1_GSE144561/SE_event_intersect.rds')
#SE_event_res = readRDS('pdac_1_GSE144561/SE_event_intersect.rds')
#dim(SE_event_res$SE_merge )
#dim(SE_mat)
#SE_event_res$SE_merge = SE_mat
SE_merge = SE_event_res$SE_merge
SE_event_inter = SE_event_res$SE_event_inter
SE_event_inter$event = paste0(
  SE_event_inter$chr.pEMT,":",
  SE_event_inter$exonStart_0base.pEMT,"-",
  SE_event_inter$exonEnd.Mes,":",
  SE_event_inter$strand.pEMT
)
table(duplicated(SE_event_inter$event))
SE_event_inter_dup = SE_event_inter[!duplicated(SE_event_inter$event),
                                    c(1:6,11,14,18)]
write.csv(SE_event_inter_dup,'pdac_1_GSE144561/SE_event_inter_114.csv')

library(ggtern)

inter_deg_dase_gene = intersect(venn_gene_Mes,venn_gene_pEMT)

deg_gene_names = list(gene_pEMTHD = unique(rownames(EMT_DEG$pEMThd_deg)),
                      gene_MesHD = unique(rownames(EMT_DEG$Meshd_deg)))

inter_deg_dase_gene = intersect(SE_event_inter$geneSymbol.pEMT,unique(unlist(deg_gene_names)))

SE_event_inter[which(SE_event_inter$geneSymbol.pEMT %in% 'MBNL1'),]
SE_event_mes[SE_event_mes$event_num =='event33983' , c('geneSymbol','chr','strand','upstreamES','downstreamEE')]
SE_event_pEMT[SE_event_pEMT$event_num =='event31945',c('geneSymbol','chr','strand','upstreamES','downstreamEE')]

#sashimi-plot.py -b input_bam.tsv -c chr11:35189834-35214861 -g CD44.gtf  
#-M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 
#-P palette.txt -o AT1G73660_2
#Get all position####
data_pos = data.frame(geneSymbol = SE_event_inter$geneSymbol.pEMT,
                      chr = SE_event_inter$chr.pEMT,
                      strand.pEMT = SE_event_inter$strand.pEMT,
                      upstreamES = SE_event_inter$upstreamES.pEMT,
                      downstreamEE = SE_event_inter$downstreamEE.pEMT
)
head(data_pos)
data_pos$pos = paste0(data_pos$geneSymbol,'.',data_pos$chr,':',data_pos$upstreamES,'-',data_pos$downstreamEE)
data_pos = data_pos[!duplicated(data_pos$pos),]
write.table(data_pos$pos,'pdac_1_GSE144561/rmats_file/sashimi/data_pos.tsv',
            sep = '\t',quote = F,col.names = F,row.names = F)

colnames(SE_mat)
mes_mat = SE_mat[,8:33]
pEMT_mat= SE_mat[,34:41]
hd_mat  = SE_mat[,42:49]

#bam_tsv_sashimiPlot####
colnames(temp_as_mes$SE)
meshd_bam = data.frame(row.names = colnames(temp_as_mes$SE)[24:58],
                       path = str_split(colnames(temp_as_mes$SE)[24:58],'[_]',simplify = T)[,2],
                       type = str_split(colnames(temp_as_mes$SE)[24:58],'[_]',simplify = T)[,1]
) 
colnames(temp_as_pEMT$SE)
pEMThd_bam = data.frame(row.names = colnames(temp_as_pEMT$SE)[24:39],
                        path = str_split(colnames(temp_as_pEMT$SE)[24:39],'[_]',simplify = T)[,2],
                        type = str_split(colnames(temp_as_pEMT$SE)[24:39],'[_]',simplify = T)[,1]
) 

all_samp =  c(colnames(mes_mat),colnames(pEMT_mat),colnames(hd_mat))
bam_all = data.frame(row.names = all_samp,
                     path = c(paste0('/home/lxxiao/xiaolixing/ctc/data/pdac_1_GSE144561/1.mapping/',str_split(colnames(mes_mat),'[_]',simplify = T)[,2]),
                              paste0('/home/lxxiao/xiaolixing/ctc/data/pdac_1_GSE144561/1.mapping/',str_split(colnames(pEMT_mat),'[_]',simplify = T)[,2]),
                              paste0(str_split(colnames(hd_mat),'[_]',simplify = T)[,2],'_merge.bam')
                     ),
                     type = str_split(all_samp,'[_]',simplify = T)[,1]
) 

write.table(meshd_bam,'pdac_1_GSE144561/rmats_file/meshd_bam.tsv',sep = '\t',row.names = T,col.names = F,quote = F)
write.table(pEMThd_bam,'pdac_1_GSE144561/rmats_file/pEMThd_bam.tsv',sep = '\t',row.names = T,col.names = F,quote = F)
write.table(bam_all,'pdac_1_GSE144561/rmats_file/bam_all.tsv',sep = '\t',row.names = T,col.names = F,quote = F)
path_ctc = c(paste0('/home/lxxiao/xiaolixing/ctc/data/pdac_1_GSE144561/1.mapping/',str_split(colnames(mes_mat),'[_]',simplify = T)[,2]),
             paste0('/home/lxxiao/xiaolixing/ctc/data/pdac_1_GSE144561/1.mapping/',str_split(colnames(pEMT_mat),'[_]',simplify = T)[,2])
             #paste0(str_split(colnames(hd_mat),'[_]',simplify = T)[,2],'_merge.bam')
)
write.table(path_ctc,'pdac_1_GSE144561/1.mapping//bam.tsv',
            sep = '\t',row.names = F,col.names = F,quote = F)

path_sc = c(#paste0('/home/lxxiao/xiaolixing/ctc/data/pdac_1_GSE144561/1.mapping/',str_split(colnames(mes_mat),'[_]',simplify = T)[,2]),
  #paste0('/home/lxxiao/xiaolixing/ctc/data/pdac_1_GSE144561/1.mapping/',str_split(colnames(pEMT_mat),'[_]',simplify = T)[,2])
  paste0(str_split(colnames(hd_mat),'[_]',simplify = T)[,2],'_merge.bam')
)
write.table(path_sc,'/home/lxxiao/xiaolixing/pdac/scnorm/data/merge//bam.tsv',
            sep = '\t',row.names = F,col.names = F,quote = F)

color_p = data.frame(row.names = 1:3,color=c("#F8C77D",'#DC5B64','#755032'))
write.table(color_p,'pdac_1_GSE144561/rmats_file/palette.txt',sep = '\t',row.names = F,col.names = F,quote = F)



#RBP_cor_heatmap####
EMT_DEG  = readRDS('pdac_1_GSE144561/emt_deg_result.rds')
expr_tpm = readRDS('pdac_1_GSE144561/expr_tpm4deg_EMT.rds')
expr_tpm_matrix = cbind(expr_tpm$ctc_Mes_tpm,expr_tpm$ctc_pEMT_tpm)
expr_tpm_matrix =as.data.frame(cbind(expr_tpm_matrix,expr_tpm$sc_ductal_tpm)) 

Mes_up_res = read.table('pdac_1_GSE144561/motif/Mes_resultMaps/pVal.up.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')
Mes_dn_res = read.table('pdac_1_GSE144561/motif/Mes_resultMaps/pVal.dn.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')

pEMT_up_res = read.table('pdac_1_GSE144561/motif/pEMT_resultMaps/pVal.up.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')
pEMT_dn_res = read.table('pdac_1_GSE144561/motif/pEMT_resultMaps/pVal.dn.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')

SE_event_res = readRDS('pdac_1_GSE144561/SE_event_intersect.rds')

mes_up_fill = ifelse(Mes_up_res>0.05,0,1)
mes_up_fill = mes_up_fill[rowSums(mes_up_fill)!=0,]
mes_dn_fill = ifelse(Mes_dn_res>0.05,0,1)
mes_dn_fill = mes_dn_fill[rowSums(mes_dn_fill)!=0,]
inter_motif = intersect(rownames(mes_up_fill),rownames(mes_dn_fill)) 

pEMT_up_fill = ifelse(pEMT_up_res>0.05,0,1)
pEMT_up_fill = pEMT_up_fill[rowSums(pEMT_up_fill)!=0,]
pEMT_dn_fill = ifelse(pEMT_dn_res>0.05,0,1)
pEMT_dn_fill = pEMT_dn_fill[rowSums(pEMT_dn_fill)!=0,]
inter_motif = intersect(rownames(pEMT_up_fill),rownames(pEMT_dn_fill))
#up_data = log(Mes_up_res[inter_motif,]+1)
#dn_data = log(Mes_dn_res[inter_motif,]+1) 
#identical(rownames(up_data), rownames(dn_data))
#identical(colnames(up_data), colnames(dn_data))
#cor_rbp


RBP_name = str_split(inter_motif,'[.]',simplify = T)[,1]
intersect(unique(unlist(venn_Mes)),RBP_name)
intersect(unique(venn_Mes$DEG_Mes),RBP_name)
intersect(unique(venn_pEMT$DEG_pEMT),RBP_name)

inter_motif_gene = str_split(inter_motif,'[.]',simplify = T)[,1]
expr_matrix = expr_tpm_matrix[intersect(inter_motif_gene,rownames(expr_tpm_matrix)),]

inter_motif_gene_final = rownames(expr_matrix)
SE_event_gene =SE_event_res$SE_event_inter$geneSymbol.pEMT
expr_deg = rownames(EMT_DEG$Meshd_deg) 

RBP_dase = intersect(inter_motif_gene_final,SE_event_gene)
RBP_deg = intersect(inter_motif_gene_final,expr_deg)
RBP_tpm = expr_matrix

#cor
head(RBP_tpm)
colnames(RBP_tpm)
psi_data = SE_event_res$SE_merge
colnames(psi_data)[1:41] = str_split(colnames(psi_data)[1:41],'Aligned.sorted',simplify = T)[,1]
colnames(psi_data)[42:49]= str_split(colnames(psi_data)[42:49],'/data/merge/',simplify = T)[,2]
colnames(psi_data)[42:49]= str_split(colnames(psi_data)[42:49],'_merge.bam',simplify = T)[,1]

sp_id = c(str_split(colnames(RBP_tpm)[1:57],'mes_',simplify = T)[,2],
          str_split(colnames(RBP_tpm)[58:76],'pEMT_',simplify = T)[,2])
SRR_ID = EMT_DEG$phe$SRRID[match(sp_id[1:76],EMT_DEG$phe$description)]
colnames(RBP_tpm) = c(paste0(rep('Mes',length(1:57)),'_',SRR_ID[1:57]),
                      paste0(rep('pEMT',length(58:76)),'_',SRR_ID[58:76]),
                      str_split(colnames(RBP_tpm)[77:84],'_',simplify = T)[,2])
RBP_tpm_final = RBP_tpm[,colnames(psi_data)[8:49]]

library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)
rbp = t(RBP_tpm_final)
as = (t(psi_data[,8:49]))
cor_res = cor(rbp,as)#method = 'pearson'
cor_res_fil = ifelse(abs(cor_res)<0.6,0,1)
cor_res_fil2 = cor_res[rowSums(cor_res_fil)!=0, colSums(cor_res_fil)!=0]

cor_res_fil3 = ifelse(abs(cor_res_fil2)<0.6,0,
                      ifelse(cor_res_fil2>0.6,1,-1))
rbp_cor_list = list(
  rbp = t(RBP_tpm_final),
  as = (t(psi_data[,8:49])),
  cor_res = cor_res,
  cor_plot_data =cor_res_fil3
)
#saveRDS(rbp_cor_list,'pdac_1_GSE144561/RBP_SE_cor.rds')#Mes
#saveRDS(rbp_cor_list,'pdac_1_GSE144561/RBP_SE_cor_pEMT.rds')
pheatmap(cor_res_fil3,#angle_col = '45',
         treeheight_col = 9,
         treeheight_row = 9)
#ggsave('pdac_1_GSE144561/figures/RBP_cor_SE.pdf')#5:9
#ggsave('pdac_1_GSE144561/figures/RBP_cor_SE_pEMT.pdf')
pheatmap(t(cor_res_fil3),#angle_col = '45',
         treeheight_col = 9,
         treeheight_row = 9)
#ggsave('pdac_1_GSE144561/figures/RBP_cor_SE_t_Mes.pdf')#9:5
#ggsave('pdac_1_GSE144561/figures/RBP_cor_SE_t_pEMT.pdf')#9:5

#triangle heatmap####
#c("pEMT"="#a3cd5b", "Mes"="#8ea0cc")
library(circlize)
library(ComplexHeatmap)
rbp_cor_list_Mes = readRDS('pdac_1_GSE144561/RBP_SE_cor.rds')
rbp_cor_list_pEMT = readRDS('pdac_1_GSE144561/RBP_SE_cor_pEMT.rds')

RBP_fill_final = rownames(rbp_cor_list_Mes$cor_plot_data)
up_gene = str_split(rownames(Mes_up_res),'[.]',simplify = T)[,1]
dn_gene = str_split(rownames(Mes_dn_res),'[.]',simplify = T)[,1]
up_motif = rownames(Mes_up_res)[match(RBP_fill_final,up_gene)]
dn_motif =rownames(Mes_dn_res)[match(RBP_fill_final,dn_gene)] 
inter_motif = intersect(up_motif,dn_motif)
up_data_mes = -log10(Mes_up_res[inter_motif,])
dn_data_mes = -log10(Mes_dn_res[inter_motif,]) 
#write.csv(up_data_mes,'pdac_1_GSE144561/Motif_up_mes.csv')
#write.csv(dn_data_mes,'pdac_1_GSE144561/Motif_down_mes.csv')


mes_rbp = str_split(rownames(up_data_mes),'[.]',simplify = T)[,1]


RBP_fill_final = rownames(rbp_cor_list_pEMT$cor_plot_data)
up_gene = str_split(rownames(pEMT_up_res),'[.]',simplify = T)[,1]
dn_gene = str_split(rownames(pEMT_dn_res),'[.]',simplify = T)[,1]
up_motif = rownames(pEMT_up_res)[match(RBP_fill_final,up_gene)]
dn_motif =rownames(pEMT_dn_res)[match(RBP_fill_final,dn_gene)] 
inter_motif = intersect(up_motif,dn_motif)
up_data_pEMT = -log10(pEMT_up_res[inter_motif,])
dn_data_pEMT = -log10(pEMT_dn_res[inter_motif,]) 
pemt_rbp = str_split(rownames(dn_data_pEMT),'[.]',simplify = T)[,1]
#write.csv(up_data_pEMT,'pdac_1_GSE144561/Motif_up_pEMT.csv')
#write.csv(dn_data_pEMT,'pdac_1_GSE144561/Motif_down_pEMT.csv')

setdiff(mes_rbp,pemt_rbp) #"PCBP1"
setdiff(pemt_rbp,mes_rbp) #"RBM42"
rbp_list = list(mes = mes_rbp,
                pEMT = pemt_rbp)
ggvenn(rbp_list,show_elements = F,label_sep = ",",
       digits = 1,stroke_color = "white",
       fill_color = c("#E41A1C", "#1E90FF", "#FF8C00"),
       set_name_color = c("#E41A1C", "#1E90FF"))
ggsave('pdac_1_GSE144561/figures/rbp_venn.pdf')


up_data_raw = up_data_mes
dn_data_raw = dn_data_mes

up_data_raw = up_data_pEMT
dn_data_raw = dn_data_pEMT

up_data = ifelse(up_data_raw> -log10(0.05),1,0)
dn_data = ifelse(dn_data_raw> -log10(0.05),1,0)

RightOrder <- rev(rownames(up_data))
up_data = up_data[RightOrder,]
dn_data = dn_data[RightOrder,]

#up_data = up_data_raw
#dn_data = dn_data_raw


identical(rownames(up_data), rownames(dn_data))
identical(colnames(up_data), colnames(dn_data))


UpColor <- colorRamp2(breaks = c(0,1), colors = c("white","#8ea0cc"))### #FFFADD
DnColor <- colorRamp2(breaks = c(0,1), colors = c("white","#a3cd5b"))


DiagFunc <- function(up, down){
  function(j, i, x, y, width, height, fill){
    grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                 unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                 gp = gpar(fill = DnColor(down[i, j]), col = "grey"))
    grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                 unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                 gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
  }
}


p1 <- Heatmap(dn_data, column_title = "Mes RBP motif Pval", #"Mes RBP motif Pval",
              rect_gp = gpar(type = "none"),
              show_heatmap_legend = F,
              show_column_names = F,
              cluster_rows = T,
              cluster_columns = F, 
              #left_annotation = row_an, 
              cell_fun = DiagFunc(up = dn_data, down =up_data ) 
)



lgd <- list(Legend(title ="Mes_upstream",#"Mes_upstream"
                   col_fun = DnColor, 
                   at = c(0, 1), 
                   direction = "horizontal" 
),
Legend(title = "Mes_downStream",#"Mes_downStream", 
       col_fun = UpColor, 
       at = c(0, 1),
       direction = "horizontal" 
) )

draw(p1, annotation_legend_list = lgd,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE)
up_data[1,]
mes_motif_row =
  
  p1



gene_pEMTHD = c(n_gene_pEMTHD$A3SS,n_gene_pEMTHD$A5SS,n_gene_pEMTHD$MXE,n_gene_pEMTHD$RI,n_gene_pEMTHD$SE)

gene_pEMTHD = unique(gene_pEMTHD)
length(gene_pEMTHD)#350

gene_MesHD = c(n_gene_MesHD$A3SS,n_gene_MesHD$A5SS,n_gene_MesHD$MXE,n_gene_MesHD$RI,n_gene_MesHD$SE)
gene_MesHD = unique(gene_MesHD)
length(gene_MesHD)#598


#dase_Padj####
temp_as_pEMT = sapply(rmats_mes_DifEventAll$pEMThd_all,function(x){x[[3]]})
temp_as_mes = sapply(rmats_mes_DifEventAll$meshd_all,function(x){x[[3]]})


rmats_mes_diff = list(pEMTHD_dase_Padj = temp_as_pEMT,
                      MesHD_dase_Padj = temp_as_mes)

n_gene_name_pEMTHD = sapply(temp_as_pEMT,function(x){x[,'geneSymbol']})
n_gene_name_MesHD = sapply(temp_as_mes,function(x){x[,'geneSymbol']})

n_gene_pEMTHD = sapply(temp_as_pEMT,function(x){x[,'ENSEMBL.Gene.ID']})
n_gene_MesHD = sapply(temp_as_mes,function(x){x[,'ENSEMBL.Gene.ID']})

gene_pEMTHD = c(n_gene_pEMTHD$A3SS,n_gene_pEMTHD$A5SS,n_gene_pEMTHD$MXE,n_gene_pEMTHD$RI,n_gene_pEMTHD$SE)
gene_pEMTHD = unique(gene_pEMTHD)
length(gene_pEMTHD)#350

gene_MesHD = c(n_gene_MesHD$A3SS,n_gene_MesHD$A5SS,n_gene_MesHD$MXE,n_gene_MesHD$RI,n_gene_MesHD$SE)
gene_MesHD = unique(gene_MesHD)
length(gene_MesHD)#598

#kegg####
library(clusterProfiler)
library(org.Hs.eg.db)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'libcurl')
options(clusterProfiler.download.method = "libcurl")

group <- data.frame(gene=c(gene_pEMTHD,gene_MesHD),
                    group=c(rep('pEMT',length(gene_pEMTHD)),rep('Mes',length(gene_MesHD))))
Gene_ID <- bitr(group$gene, fromType="ENSEMBL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
data  <- merge(Gene_ID,group,by.x='ENSEMBL',by.y='gene')
#write.csv(data,'pdac_1_GSE144561/GO_data_input.csv')
data_KEGG <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichKEGG",
  organism="hsa",
  pvalueCutoff = 0.05
)
data_GO_sim <- simplify(data_KEGG, 
                        cutoff=0.7, 
                        by="p.adjust", 
                        select_fun=min)
dotplot(data_GO_sim, showCategory=5,font.size = 8)

#go_Group####
library(clusterProfiler)
library(org.Hs.eg.db)
group <- data.frame(gene=c(gene_pEMTHD,gene_MesHD),
                    group=c(rep('pEMT',length(gene_pEMTHD)),rep('Mes',length(gene_MesHD))))
table(group$group)
Gene_ID <- bitr(group$gene, fromType="ENSEMBL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
data  <- merge(Gene_ID,group,by.x='ENSEMBL',by.y='gene')
data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
data_GO_sim <- simplify(data_GO, 
                        cutoff=0.7, 
                        by="p.adjust", 
                        select_fun=min)
dotplot(data_GO_sim, showCategory=5,font.size = 8)
