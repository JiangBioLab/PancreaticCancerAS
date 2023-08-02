rm (list = ls())
setwd('/home/lxxiao/xiaolixing/ctc/data/')
library(stringr)
library(ggplot2)
setwd('/home/lxxiao/xiaolixing/ctc/data/')
library(stringr)
library(ggplot2)
library(pheatmap)
#loading data#
rmats_mes_DifEventAll= readRDS('pdac_1_GSE144561/rmats_file/rmats_mes_DifEventAll.rds')
rmats_mes_DifSEResult = readRDS('pdac_1_GSE144561/rmats_file/rmats_mes_DifSE.rds')
rmats_mes_SEPadjResult = readRDS('pdac_1_GSE144561/rmats_file/rmats_mes_DifSEPadj.rds')
GESP = read.csv('pdac_1_GSE144561/hunman_cell_surface_pro.csv',header = T)#TCSA Download
#phe_ctc = read.csv('pdac_1_GSE144561/phenoData_v3.csv',header = T,row.names = 1)
phe_ctc = read.csv('pdac_1_GSE144561/phenoData_final.csv',header = T,row.names = 1)
head(GESP)
emt_signature = read.csv('/home/lxxiao/xiaolixing/pdac/10x/emt_gene_signature.csv',header = T,row.names = 1)
SE_inter_event_res = readRDS('pdac_1_GSE144561/SE_event_intersect.rds')

gene_annotation_gtf = read.table('~/xiaolixing/ref/gencode.v38.gtf.annotation_rewrite.txt',header = T,sep = '\t')

SE_inter_event_114 = SE_inter_event_res$SE_merge
SE_gene = unique(SE_inter_event_114$geneSymbol)#91

mes_cordata = readRDS('pdac_1_GSE144561/RBP_SE_cor.rds')#Mes
pEMT_cordata = readRDS('pdac_1_GSE144561/RBP_SE_cor_pEMT.rds')

mes_cordata_R =  as.data.frame(mes_cordata$cor_res)
pEMT_cordata_R =  as.data.frame(pEMT_cordata$cor_res)
rbp_inter = intersect(rownames(mes_cordata_R),rownames(pEMT_cordata_R))

RBP_inter_genename = intersect(rbp_inter,SE_gene)
RBP_inter_envent = SE_inter_event_114[which(SE_inter_event_114$geneSymbol %in% RBP_inter_genename),]




#mul_iso transition####
library(stringi)
pEMThd_result = rmats_mes_SEPadjResult$pEMThd_dase_Padj  #482
meshd_result = rmats_mes_SEPadjResult$meshd_dase_Padj    #245
colnames(SE_inter_event_114)
event_114 = SE_inter_event_114[,8:48]

length(unique(SE_inter_event_114$geneSymbol)) 
#sankey_plot####
# install.packages("remotes")
# remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(tidyverse)

get_mode = function(data,Sam_type){
  colnames(data) = str_split(colnames(data),'_',simplify = T)[,1]    
  temp = data[,which(colnames(data)==Sam_type)]                      
  colnames(temp) = str_split(colnames(temp),'[.]',simplify = T)[,1]  
  event_mean = rowMeans(temp)                                        
  mode_type = ifelse(event_mean>0.95,'included',ifelse(event_mean<0.05,'excluded','middle')) 
  return(mode_type)
}

mes_type = get_mode(event_114,'Mes')
pEMT_type = get_mode(event_114,'pEMT')
HD_type = get_mode(event_114,'HD')

mydata = data.frame(hd_type = HD_type,
                    pEMT_type = pEMT_type,
                    mes_type = mes_type,
                    row.names = rownames(event_114))

#write.csv(mydata,'pdac_1_GSE144561/SE_mode_table.csv')
names(table(mydata$mes_type))

df <- mydata %>%
  make_long(hd_type,
            pEMT_type,
            mes_type)

table(df$x)
#df$x = factor(df$x,levels = c("HD_type","pEMT_type","mes_type"))
#df$node = factor(df$node,levels = c("middle","included","excluded"))
#调整node的factor level 可以调整node的顺序（倒序）
color =list(
  #"bimodal"=       "#d00000",
  "excluded"=      "#F8B62C",
  "included"=      "#3E99CC",
  "middle"=        "#98CA6C"
  #"uncategorized"= "#136f63"
)

ggplot(df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey(flow.alpha = 0.50, node.color = 1) +
  scale_fill_manual( values =color )+
  theme_sankey(base_size = 16)

#ggsave('pdac_1_GSE144561/figures/sankey_114event.pdf')#3 7.5


#stack_barplot1####
library(reshape2)
hd_2_pEMT = ifelse(mydata$hd_type==mydata$pEMT_type,'No_mode_change','mode_change')
pEMT_2_mes = ifelse(mydata$pEMT_type==mydata$mes_type,'No_mode_change','mode_change')
table(hd_2_pEMT)
table(pEMT_2_mes)

data_bar = data.frame(hd_2_pEMT = c(64,50),
                      pEMT_2_mes = c(28,86),
                      type = c('mode_change','No_mode_change'),
                      row.names = c('mode_change','No_mode_change'))

library(reshape2)
library(ggplot2)
color =list(
  "No_mode_change"="#d00000",
  "mode_change"= "#ffba08"
)
data = melt(data_bar,'type')
head(data)

data$percent = data$value/114*100

col_fill = list(
  "No_mode_change"="#d00000",
  "mode_change"= "#ffba08"
)
ggplot(data, aes( x = variable,y=value,fill = type))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  scale_fill_manual(values = col_fill)+
  theme_bw()
#ggsave('pdac_1_GSE144561/figures/barplot_mode_event114.pdf')



#mode type length####
hd_mode_change = mydata[which(hd_2_pEMT=='mode_change'),1:2]
pEMT_mode_change = mydata[which(pEMT_2_mes=='mode_change'),2:3]
event_inter114 = SE_inter_event_res$SE_event_inter
event_inter114$event = paste(event_inter114$geneSymbol.pEMT,":",
                             event_inter114$chr.pEMT,":",
                             event_inter114$exonStart_0base.pEMT,"-",
                             event_inter114$exonEnd.pEMT,":",
                             event_inter114$strand.pEMT,
                             sep = '')

data_hd_mode_change_info = event_inter114[match(rownames(hd_mode_change),event_inter114$event),]
data_pEMT_mode_change_info = event_inter114[match(rownames(pEMT_mode_change),event_inter114$event),]

hd_mode_change_length = data_hd_mode_change_info[,c('delt.pEMT','delt.Mes','exonStart_0base.pEMT','exonEnd.pEMT','event')]
hd_mode_change_length$exon_length = hd_mode_change_length$exonEnd.pEMT-hd_mode_change_length$exonStart_0base.pEMT
rownames(hd_mode_change_length) = hd_mode_change_length$event
all(
  rownames(hd_mode_change) == hd_mode_change_length$event
)
hd_mode_change_length = cbind(hd_mode_change_length,
                              hd_mode_change)

dim(hd_mode_change_length)#64  5
table(hd_mode_change_length$exon_length>100 , hd_mode_change_length$pEMT_type)
head(hd_mode_change_length)
boxplot_data = hd_mode_change_length

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x <= (qnt[1] - H)] <- NA
  y[x >= (qnt[2] + H)] <- NA
  y
}

remove_outliers(boxplot_data$exon_length)
colnames(input)
library(dplyr)
library(ggpubr)
input <- boxplot_data %>% group_by(pEMT_type) %>%  mutate(value = remove_outliers(exon_length))
head(input)
input2 <- na.omit(input)
table(input2$pEMT_type)

box_comp_list = list(
  c('excluded','included'),
  c('included','middle'),
  c('excluded','middle')
)
ggplot(input2,aes(x = pEMT_type,y = value,fill =pEMT_type)) + 
  #geom_violin(trim =T,scale = 'width')+
  #geom_jitter(width = 0.1) +
  geom_boxplot(width = 0.4) +
  stat_boxplot(geom = "errorbar",width=0.15)+
  geom_signif(comparisons = box_comp_list,
              test = 't.test',
              map_signif_level = c("***"<0.001, "**"<0.01, "*"<0.05),
              y_position = c(220,220,245))+
  scale_fill_brewer()+
  ylab('Exon length')+
  theme_bw()
#ggsave('pdac_1_GSE144561/figures/exon_length_pemt.pdf') #4:5


#pemt2mes
pEMT_mode_change_length = data_pEMT_mode_change_info[,c('delt.pEMT','delt.Mes','exonStart_0base.pEMT','exonEnd.pEMT','event')]
pEMT_mode_change_length$exon_length = pEMT_mode_change_length$exonEnd.pEMT-pEMT_mode_change_length$exonStart_0base.pEMT
rownames(pEMT_mode_change_length) = pEMT_mode_change_length$event
all(
  rownames(pEMT_mode_change) == pEMT_mode_change_length$event
)
pEMT_mode_change_length = cbind(pEMT_mode_change_length,
                                pEMT_mode_change)

dim(pEMT_mode_change_length)#28  8
table(pEMT_mode_change_length$exon_length>100 , pEMT_mode_change_length$pEMT_type)

boxplot_data = pEMT_mode_change_length
input <- boxplot_data %>% group_by(mes_type) %>%  mutate(value = remove_outliers(exon_length))
head(input)
input2 <- na.omit(input)
table(input2$mes_type)

box_comp_list = list(
  c('excluded','included'),
  c('included','middle'),
  c('excluded','middle')
)
ggplot(input2,aes(x = mes_type,y = value,fill =mes_type)) + 
  #geom_violin(trim =T,scale = 'width')+
  #geom_jitter(width = 0.1) +
  geom_boxplot(width = 0.4) +
  stat_boxplot(geom = "errorbar",width=0.15)+
  geom_signif(comparisons = box_comp_list,
              test = 't.test',
              map_signif_level = c("***"<0.001, "**"<0.01, "*"<0.05),
              y_position = c(220,220,245))+
  scale_fill_brewer()+
  ylab('Exon length')+
  theme_bw()
#ggsave('pdac_1_GSE144561/figures/exon_length_pemt.pdf') #4:5



#nochange events####
data_nochange = mydata[which(hd_2_pEMT =='No_mode_change' & pEMT_2_mes=='No_mode_change'),] #全都是middle
data_middle_event = rownames(data_nochange) 
event_inter114 = SE_inter_event_res$SE_event_inter #122 18
event_inter114$event = paste(event_inter114$geneSymbol.pEMT,":",
                             event_inter114$chr.pEMT,":",
                             event_inter114$exonStart_0base.pEMT,"-",
                             event_inter114$exonEnd.pEMT,":",
                             event_inter114$strand.pEMT,
                             sep = '')
data_middle_event_info = event_inter114[match(data_middle_event,event_inter114$event),]


middle_event_length = data_middle_event_info[,c('delt.pEMT','delt.Mes','exonStart_0base.pEMT','exonEnd.pEMT')]
middle_event_length$exon_length = middle_event_length$exonEnd.pEMT-middle_event_length$exonStart_0base.pEMT
dim(middle_event_length)#45 5
table(middle_event_length$exon_length>100)
table(middle_event_length$delt.pEMT>0)
table(middle_event_length$delt.Mes>0)

data_change_event_info = event_inter114[-match(data_middle_event,event_inter114$event),]
data_change_event_info = data_change_event_info[,c('delt.pEMT','delt.Mes','exonStart_0base.pEMT','exonEnd.pEMT')]
data_change_event_info$exon_length = data_change_event_info$exonEnd.pEMT-data_change_event_info$exonStart_0base.pEMT
dim(data_change_event_info)#77 5
table(data_change_event_info$exon_length>100)
table(data_change_event_info$delt.pEMT>0,
      data_change_event_info$delt.Mes>0)




bar_data_nochange = as.data.frame(table(middle_event_length$exon_length>100))
bar_data_nochange$percent = bar_data_nochange$Freq/45
bar_data_nochange$type = 'nochange'

ggplot(bar_data_nochange,aes(x = type,y = percent,fill = Var1))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  #scale_fill_manual(values = col_fill)+
  theme_bw()

bar_data_change = as.data.frame(table(data_change_event_info$exon_length>100))
bar_data_change$percent = bar_data_change$Freq/77
bar_data_change$type = 'change'
ggplot(bar_data_change,aes(x = type,y = percent,fill = Var1))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  #scale_fill_manual(values = col_fill)+
  theme_bw()

#ggsave('pdac_1_GSE144561/figures/exon_length_nochange.pdf')
#ggsave('pdac_1_GSE144561/figures/exon_length_change.pdf')

#delta PSI####
table(mydata$mes_type)
data_mydata = mydata
data_mydata$mode = ifelse(mydata$pEMT_type==mydata$mes_type,'No_mode_change','mode_change')
data_mydata = subset(data_mydata,data_mydata$mode=='mode_change')
data_mydata_info = event_inter114[match(rownames(data_mydata),event_inter114$event),]
data_mydata_length = data_mydata_info[,c('delt.pEMT','delt.Mes',
                                         'exonStart_0base.pEMT','exonEnd.pEMT',
                                         'event')]
data_mydata_length$exon_length = data_mydata_length$exonEnd.pEMT-data_mydata_length$exonStart_0base.pEMT
all(data_mydata_length$event==rownames(data_mydata))
rownames(data_mydata_length) = data_mydata_length$event
data_mydata_length = cbind(data_mydata_length,
                           data_mydata)
data_mydata_length$length_type = ifelse(data_mydata_length$exon_length<100,'short','long')
colnames(data_mydata_length)
ggdata=data_mydata_length[,c("delt.pEMT","delt.Mes","length_type",'mes_type')]
table(ggdata$mes_type)
col_fill = list("short"='#c44536',
                "long"='#0077b6')
color =list(
  "excluded"=      "#F8B62C",
  "included"=      "#3E99CC",
  "middle"=        "#98CA6C"
)
ggplot(ggdata) +
  geom_point(aes(x=delt.pEMT, y=delt.Mes, col=mes_type), size=1) +
  scale_color_manual(values =color ) +
  theme_classic()

rownames(data_mydata_length)[which(data_mydata_length$delt.pEMT>0 & data_mydata_length$delt.Mes<0)]

#find ending event
data_ex = mydata[which(mydata$hd_type !='middle'),]
data_ex
#PTPN6:chr12:6951463-6951520:+####
#boxplot
library(ggplot2)
library(ggprism)
library(ggsignif)
box_data = event_114
box_data = box_data[grep('PTPN6',rownames(box_data)),]
box_data = data.frame(type = str_split(colnames(box_data),'[_]',simplify = T)[,1],
                      PSI =  as.numeric(box_data[1,]) ,
                      sample_id =  str_split(colnames(box_data),'[_]',simplify = T)[,2],
                      row.names = colnames(box_data))
table(box_data$type)
box_data$type = factor(box_data$type,levels = c('Mes','pEMT','HD'))
box_comp_list = list(
  c('HD','Mes'),
  c('HD','pEMT')
)
ggplot(box_data,aes(x = type,y = PSI))+
  geom_violin()
ggplot(box_data,aes(x =type,y = PSI ))+stat_boxplot(geom = "errorbar",width=0.2)+
  #geom_boxplot(outlier.size = 0.1,aes(fill=type))+
  geom_violin(trim =T,aes(fill=type),scale = 'width')+
  ylim(0,1.5)+
  geom_signif(comparisons = box_comp_list,
              test = 't.test',
              map_signif_level = c("***"<0.001, "**"<0.01, "*"<0.05),
              y_position = c(1.0,1.1,1.2))+
  scale_fill_manual(values = c(
    'Mes'='#DC5B64',
    'pEMT'="#F8C77D",
    'HD'='#755032') )+
  theme_prism(axis_text_angle = 0,
              base_line_size = 0.1,
              base_fontface = "plain",
              base_family = "serif")+
  #ggtitle(event)+
  xlab('')
#ggsave('pdac_1_GSE144561/figures/violin_PTPN6.pdf')


#cell surface event length####

rmats_mes_DifEventAll= readRDS('pdac_1_GSE144561/rmats_file/rmats_mes_DifEventAll.rds')
rmats_mes_DifSEResult = readRDS('pdac_1_GSE144561/rmats_file/rmats_mes_DifSE.rds')
rmats_mes_SEPadjResult = readRDS('pdac_1_GSE144561/rmats_file/rmats_mes_DifSEPadj.rds')
GESP = read.csv('pdac_1_GSE144561/hunman_cell_surface_pro.csv',header = T)#TCSA Download

pEMT_data = rmats_mes_DifSEResult$pEMTHD_dase$FilByPadj
mes_data = rmats_mes_DifSEResult$mesHD_dase$FilByPadj

dim(pEMT_data)#482
dim(mes_data)#245

pEMT_surface_result =  pEMT_data[which(pEMT_data$geneSymbol %in% GESP$HGNC.Symbol),]
mes_surface_result =  mes_data[which(mes_data$geneSymbol %in% GESP$HGNC.Symbol),]

dim(pEMT_surface_result)#53
dim(mes_surface_result)#22

colnames(pEMT_surface_result)
pEMT_surface_table=  pEMT_surface_result[,c("ENSEMBL.Gene.ID","geneSymbol","chr","strand","exonStart_0base","exonEnd","delt","p_self","padj")]
mes_surface_table =  pEMT_surface_result[,c("ENSEMBL.Gene.ID","geneSymbol","chr","strand","exonStart_0base","exonEnd","delt","p_self","padj")]

pEMT_surface_table$exon_length = pEMT_surface_table$exonEnd-pEMT_surface_table$exonStart_0base+1
mes_surface_table$exon_length = mes_surface_table$exonEnd-mes_surface_table$exonStart_0base+1

get_pEMT_surface_seq = paste0("samtools faidx /home/lxxiao/xiaolixing/ref/GRCh38.p13.genome.fa ",
                              pEMT_surface_table$chr,":",
                              pEMT_surface_table$exonStart_0base,"-",
                              pEMT_surface_table$exonEnd," ",
                 ifelse(pEMT_surface_table$strand=='-','--reverse-complement','')
)

write.table(get_pEMT_surface_seq,quote = F,row.names = F,col.names = F,
            'pdac_1_GSE144561/get_pEMT_surface_seq.sh')
#bash get_pEMT_surface_seq.sh >get_pEMT_surface_seq_result.txt
#python pemt_seq.py

get_mes_surface_seq = paste0("samtools faidx /home/lxxiao/xiaolixing/ref/GRCh38.p13.genome.fa ",
                              mes_surface_table$chr,":",
                              mes_surface_table$exonStart_0base,"-",
                              mes_surface_table$exonEnd," ",
                              ifelse(mes_surface_table$strand=='-','--reverse-complement','')
)

write.table(get_mes_surface_seq,quote = F,row.names = F,col.names = F,
            'pdac_1_GSE144561/get_mes_surface_seq.sh')
#bash get_mes_surface_seq.sh >get_mes_surface_seq_result.txt
#python mes_seq.py

get_pEMT_seq_result = read.table('pdac_1_GSE144561/get_pEMT_surface_seq_dataFrame',
                                  sep = '\t',header = F)
colnames(get_pEMT_seq_result) = c('event_position','sequence')
pEMT_surface_table$sequence = get_pEMT_seq_result$sequence

get_mes_seq_result = read.table('pdac_1_GSE144561/get_mes_surface_seq_dataFrame',
                                 sep = '\t',header = F)
colnames(get_mes_seq_result) = c('event_position','sequence')
mes_surface_table$sequence = get_mes_seq_result$sequence

write.csv(pEMT_surface_table,'pdac_1_GSE144561/pEMT_surface_protein.csv')
write.csv(mes_surface_table,'pdac_1_GSE144561/mes_surface_protein.csv')


#surface_barplot####
bar_plot_pEMT = pEMT_surface_result[,c("geneSymbol"  ,"delt")]#,"Level1_mean" ,"Level2_mean"
bar_plot_mes  = mes_surface_result[,c("geneSymbol" ,"delt")]#,"Level1_mean" ,"Level2_mean" 

#bar_plot_pEMT$hd_type = ifelse(bar_plot_pEMT$Level2_mean>0.95,1,ifelse(bar_plot_pEMT$Level2_mean<0.05,0,0.5))
#bar_plot_mes$hd_type = ifelse(bar_plot_mes$Level2_mean>0.95,1,ifelse(bar_plot_mes$Level2_mean<0.05,0,0.5))

bar_plot_pEMT$type = as.factor(ifelse(bar_plot_pEMT$delt<0,'down','up'))
bar_plot_mes$type = as.factor(ifelse(bar_plot_mes$delt<0,'down','up'))

#table(bar_plot_pEMT$type,bar_plot_pEMT$hd_type != 0.5)
#table(bar_plot_mes$type,bar_plot_mes$hd_type != 0.5)

bar_plot_pEMT_data = as.data.frame(table(bar_plot_pEMT$type))
bar_plot_mes_data = as.data.frame(table(bar_plot_mes$type))

bar_plot_pEMT_data$Var1 = as.factor(bar_plot_pEMT_data$Var1)
ggplot(bar_plot_pEMT_data,aes(x = Var1,y = Freq))+
  geom_bar(stat = 'identity',width = 0.6,position = 'dodge')+
  geom_text(aes(label = Freq), vjust = - 0.2)+
  scale_color_manual(values = c(
    'down'='#DC5B64',
    'up'="#F8C77D",
    'HD'='#755032') )+
  theme_bw()+
  #ggtitle(event)+
  xlab('')
#ggsave('pdac_1_GSE144561/figures/CART_pEMT.pdf')
bar_plot_mes_data$Var1 = as.factor(bar_plot_mes_data$Var1)
ggplot(bar_plot_mes_data,aes(x = Var1,y = Freq))+
  geom_bar(stat = 'identity',width = 0.6,position = 'dodge')+
  geom_text(aes(label = Freq), vjust = - 0.2)+
  scale_color_manual(values = c(
    'down'='#DC5B64',
    'up'="#F8C77D",
    'HD'='#755032') )+
  theme_bw()+
  #ggtitle(event)+
  xlab('')
#ggsave('pdac_1_GSE144561/figures/CART_Mes.pdf')



#get AS sequence####
#1515 Server
#conda actiavate R42
#BiocManager::install("BSgenome")
library(BSGenome)
library("BSgenome.Hsapiens.UCSC.hg38")
events = rownames(mydata)
chr = str_split(events,':',simplify = T)[,2]
exon_position = str_split(events,':',simplify = T)[,3]
star_position = as.numeric(str_split(exon_position,'-',simplify = T)[,1])
end_position = as.numeric(str_split(exon_position,'-',simplify = T)[,2])
event_position = data.frame(
  chr = chr,
  star_position = star_position,
  end_position = end_position,
  direction = str_split(events,':',simplify = T)[,4],
  exon_length = end_position-star_position+1
)
rownames(event_position) = paste0(event_position$chr,":",
                                  event_position$star_position,"-",
                                  event_position$end_position)
ifelse(event_position$direction=='-','--reverse-complement','')
get_seq = paste0("samtools faidx /home/lxxiao/xiaolixing/ref/GRCh38.p13.genome.fa ",
                 event_position$chr,":",
                 event_position$star_position,"-",
                 event_position$end_position," ",
                 ifelse(event_position$direction=='-','--reverse-complement','')
)


write.table(get_seq,
            quote = F,
            row.names = F,
            col.names = F,'pdac_1_GSE144561/get_event_sequence.sh')
#bash get_event_sequence.sh >get_event_sequence_result.txt
#python get_seq_dataframe.py

get_event_seq_result = read.table('pdac_1_GSE144561/get_sequence_dataFrame',
                                  sep = '\t',header = F)
colnames(get_event_seq_result) = c('event_position','sequence')
get_event_seq_result$event_position = str_split(get_event_seq_result$event_position,
                                                '/r',simplify = T)[,1]
all(rownames(event_position)==get_event_seq_result$event_position)

mydata_event_position  = cbind(mydata,event_position$direction,event_position$exon_length)
mydata_event_position  = cbind(mydata_event_position,get_event_seq_result)
colnames(mydata_event_position)[4] = 'strand'
colnames(mydata_event_position)[5] = 'exon_length'
mydata_event_position$event = rownames(mydata_event_position)
mydata_event_position = mydata_event_position[,c("event","exon_length","event_position",'strand',"sequence",
                                                 "hd_type","pEMT_type","mes_type")]
#write.table(mydata_event_position,'pdac_1_GSE144561/get_seq_length.txt',quote = F,sep = '\t' )


pct_TAG = as.numeric(str_detect(mydata_event_position$sequence,'TAG'))
pct_TAA = as.numeric(str_detect(mydata_event_position$sequence,'TAA'))
pct_TGA = as.numeric(str_detect(mydata_event_position$sequence,'TGA'))
pct_data = data.frame(TAG = pct_TAG,TAA= pct_TAA,TGA = pct_TGA)
rownames(pct_data) = rownames(mydata_event_position)
#write.table(pct_data,'pdac_1_GSE144561/PTC_result.txt',quote = F,row.names = T,sep = '\t')
pct_data = cbind(pct_data,mydata_event_position)
pct_data$chr = str_split(pct_data$event_position,'[:]',simplify = T)[,1]
pct_data$chr = factor(pct_data$chr,
                      levels = c(paste0('chr',1:20),'chrX'))
colnames(pct_data)
pct_data_plot= pct_data[,1:3]

pheatmap(pct_data_plot,show_rownames = T,cluster_rows = F,cluster_cols = F,
         show_colnames = T,
         #border_color = 'grey60',
         color = c("black", 'grey70',"white"),
)
pct_data['PTPN6:chr12:6951463-6951520:+',]
'GCTGCACCTCCTCATTCCCTGCGCCCCCTTCCTCTCCGGAAGCCCCCAGGATGGTGAG'

#circos.heatmap####
library('circlize')
library('ComplexHeatmap')
col_TAG=colorRamp2(c(0,1), c("#ffffff", "#a6cee3"))
col_TAA=colorRamp2(c(0,1), c("#ffffff", "#3690c0"))
col_TGA=colorRamp2(c(0,1), c("#ffffff", "#045a8d"))
col_length = colorRamp2(c(0, 200), c("white","#fb9a99"))
col_ptc = colorRamp2(c(0, 3), c("#F7EFE5","#0868ac"))
split = pct_data$chr
circos.clear()
circos.par(gap.after=c(90))
circos.heatmap(pct_data$exon_length,dend.track.height = 0.04,track.height = 0.05,
               col =col_length,cluster = F,bg.border = c('grey60'),
               bg.lwd = 1, bg.lty = 1)
circos.track(
  ylim = c(0,3), panel.fun = function(x, y) {
    m = pct_data_plot[CELL_META$subset, , drop = FALSE]
    m = m[CELL_META$row_order, , drop = FALSE]
    n = nrow(m)
    # circos.boxplot 应用于矩阵的列，所以需要转置
    value = rowSums(pct_data_plot)
    circos.barplot(
      value, pos = 1:n - 0.5, #pch = 16, cex = 0.3, 
      col = ifelse(value == 0, 'white', 
                   ifelse(value==1,'#b2df8a',
                          ifelse(value==2,'#fdbf6f','#cab2d6')))) #CELL_META$sector.numeric.index)
    circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = 'white')
  }, 
  cell.padding = c(0.02, 0, 0.02, 0)
) 
circos.heatmap(pct_data$TAG,dend.track.height = 0.2,track.height = 0.1,
               col =col_TAG,cluster = F,bg.border = c('grey60'),
               bg.lwd = 1, bg.lty = 1) 
circos.heatmap(pct_data$TAA,dend.track.height = 0.4,track.height = 0.1,
               col =col_TAA,cluster = F,bg.border = c('grey60'),
               bg.lwd = 1, bg.lty = 1) 
circos.heatmap(pct_data$TGA,dend.track.height = 0.4,track.height = 0.1,
               col =col_TGA,cluster = F,bg.border = c('grey60'),
               bg.lwd = 1, bg.lty = 1) 

circos.clear() 
#ggsave('pdac_1_GSE144561/figures/PTC_result_cirheatmap.pdf')

#pheatmap(pct_data$exon_length,color = col_length) #get annotation bar
#PTCbar####
pct_data$No.PCT = rowSums(pct_data_plot)
colnames(pct_data)
pct_data[,c("TAG","TAA","TGA","No.PCT")]
table(pct_data$No.PCT)
pct_bar = as.data.frame(table(pct_data$No.PCT))
value = pct_bar$Var1
ggplot(pct_bar,aes(x = Var1,y =Freq))+
  geom_bar(stat = 'identity', fill = ifelse(value == 0, '#fddaec', 
                                            ifelse(value==1,'#b3cde3',
                                                   ifelse(value==2,'#ccebc5','#decbe4'))))+
  xlab('Number of PTC')+
  theme_classic()
#ggsave('pdac_1_GSE144561/figures/PTC_Number_barplot.pdf')

#PTCdencity####
cols = c("#fddaec","#b3cde3","#ccebc5","#decbe4",
         "#e5ce81","#f47720","#459943","#bdc3d2")
pct_point = pct_data[,c("No.PCT","exon_length")]
pct_point$exon_length[pct_point$exon_length>500]=500
pct_point$No.PCT = as.factor(pct_point$No.PCT)
ggplot(pct_point, aes(x=exon_length,fill =No.PCT))+ 
  geom_density(alpha=0.8)+
  scale_fill_manual(values = cols)+
  theme_classic()
#ggsave('pdac_1_GSE144561/figures/PTC_length_density.pdf') 



#MEME####
#https://meme-suite.org/meme/
#INCLUSION
meme_data = pct_data
meme_data$position = str_split(meme_data$event_position,'[:]',simplify = T)[,2]
meme_data$star_position = str_split(meme_data$position,'[-]',simplify = T)[,1]
meme_data$end_position = str_split(meme_data$position,'[-]',simplify = T)[,2]
meme_data$chr = str_split(meme_data$event,'[:]',simplify = T)[,2]


pct = c("TAG","TAA","TGA","No.PCT")
library('Biostrings')
result_ptc = c()
for (i in 1:nrow(meme_data)) {
  x <- DNAString(meme_data[i,'sequence'])
  
  TAG <- matchPattern("TAG", x)
  TAA <- matchPattern("TAA", x)
  TGA <- matchPattern("TGA", x)
  if (all(c(length(TAG@ranges@start),
            length(TAA@ranges@start),
            length(TGA@ranges@start))==0) ) {
    next
  }else{
    start_ptc = c(TAG@ranges@start,
                  TAA@ranges@start,
                  TGA@ranges@start)
    
    strand_ptc = meme_data$strand[i]
    if (strand_ptc=='+') {
      #upstream_ptc_40nts = base01-40
    #downstream_ptc_60nts = base01+62
    base01 = start_ptc-1+ as.numeric(meme_data$star_position[i])
    upstream_ptc_40nts = base01-40
    downstream_ptc_60nts = base01+62
    
    up_position = paste0(upstream_ptc_40nts,'-',base01)
    down_position = paste0(base01,'-',downstream_ptc_60nts)
    }else{
      
      exon_length = meme_data$exon_length[i]
      #downstream_ptc_60nts = base01+40
      #upstream_ptc_40nts  = base01-62
      upstream_position  = exon_length-start_ptc+1
      downstream_position= exon_length-start_ptc+1+2
      
      base01 = upstream_position-1+ as.numeric(meme_data$star_position[i])
      upstream_ptc_40nts = base01-62
      downstream_ptc_60nts  = base01+40
      
      up_position = paste0(base01,'-',downstream_ptc_60nts)
      down_position = paste0(upstream_ptc_40nts,'-',base01)
    }
    
    chr_ptc = as.character(meme_data$chr)[i]
    print(i)
    print(chr_ptc)
    
    temp_data = data.frame(
      num = i,
      start_ptc = start_ptc,
      base01 =base01,
      upstream_ptc_40nts = upstream_ptc_40nts,
      downstream_ptc_60nts = downstream_ptc_60nts,
      chr_ptc = chr_ptc,
      strand_ptc = strand_ptc,
      up_position = up_position,
      down_position = down_position
      
    )
  }
    result_ptc =rbind(result_ptc,temp_data)
}

result_ptc$position = paste0(result_ptc$chr_ptc,':',
                             result_ptc$upstream_ptc_40nts,"-",
                             result_ptc$downstream_ptc_60nts)
result_ptc = result_ptc[!duplicated(result_ptc$position),]
get_seq_ptc = paste0("samtools faidx /home/lxxiao/xiaolixing/ref/GRCh38.p13.genome.fa ",
                     result_ptc$chr_ptc,":",
                     result_ptc$upstream_ptc_40nts,"-",
                     result_ptc$downstream_ptc_60nts," ",
                      ifelse(result_ptc$strand_ptc=='-','--reverse-complement',''))
write.table(get_seq_ptc,quote = F,row.names = F,col.names = F,
            'pdac_1_GSE144561/get_seq_PTC.sh')


result_ptc$upstream_position = paste0(result_ptc$chr_ptc,':',
                             result_ptc$up_position,':',
                             result_ptc$strand_ptc)
result_ptc = result_ptc[!duplicated(result_ptc$upstream_position),]
get_seq_up = paste0("samtools faidx /home/lxxiao/xiaolixing/ref/GRCh38.p13.genome.fa ",
                     result_ptc$chr_ptc,":",
                     result_ptc$up_position,
                     ifelse(result_ptc$strand_ptc=='-',' --reverse-complement',''))
write.table(get_seq_up,quote = F,row.names = F,col.names = F,
            'pdac_1_GSE144561/get_seq_up.sh')

#bash get_seq_up.sh >get_seq_up_result.fa
# python get_ptc_100nt.py #get_up40nt_dataFrame
#meme get_seq_PTC_result.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -objfun classic -revcomp -markov_order 0
#
result_ptc$downstream_position = paste0(result_ptc$chr_ptc,':',
                                      result_ptc$down_position,':',
                                      result_ptc$strand_ptc)
result_ptc = result_ptc[!duplicated(result_ptc$upstream_position),]
get_seq_down = paste0("samtools faidx /home/lxxiao/xiaolixing/ref/GRCh38.p13.genome.fa ",
                    result_ptc$chr_ptc,":",
                    result_ptc$down_position,
                    ifelse(result_ptc$strand_ptc=='-',' --reverse-complement',''))
write.table(get_seq_down,quote = F,row.names = F,col.names = F,
            'pdac_1_GSE144561/get_seq_down.sh')

#bash get_seq_down.sh >get_seq_down_result.fa
# python get_ptc_100nt.py #get_up40nt_dataFrame

ptc_100nt_set = read.table('pdac_1_GSE144561/get_PTC100nt_dataFrame',sep = '\t')
table(nchar(ptc_100nt_set$V2[1:10]) )
ggseqlogo(ptc_100nt_set$V2,method = 'prob')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('pdac_1_GSE144561/figures/ptc_100nt_prob.pdf')
ggseqlogo(ptc_100nt_set$V2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))


