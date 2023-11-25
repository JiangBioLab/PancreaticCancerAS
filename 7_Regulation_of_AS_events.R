setwd('/mywork/path')
ps <- c('stringr', 'stringi', 'stats', 'reshape2','ggpubr','dplyr',
        'DESeq2','clusterProfiler','org.Hs.eg.db')
lapply(ps, function(x){library(x, character.only = T)}) 
#Motif analysis of exon skipping events in each group was performed using 
#rMAPS (non-differential genes were used as background)
#Loading rMAPS result
Mes_up_res = read.table('AS_analysis/rMAPS2 Motif/mes.pVal.dn.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')
Mes_dn_res = read.table('AS_analysis/rMAPS2 Motif/mes.pVal.up.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')

pEMT_up_res = read.table('AS_analysis/rMAPS2 Motif/pEMT.pVal.up.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')
pEMT_dn_res = read.table('AS_analysis/rMAPS2 Motif/pEMT.pVal.up.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')

epi_up_res = read.table('AS_analysis/rMAPS2 Motif/epi.pVal.dn.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')
epi_dn_res = read.table('AS_analysis/rMAPS2 Motif/epi.pVal.up.vs.bg.RNAmap.txt',header = T,row.names = 1,sep = '\t')

mes_up_fill = ifelse(Mes_up_res>0.05,0,1)
mes_up_fill = mes_up_fill[rowSums(mes_up_fill)!=0,]
mes_up_gene = str_split(rownames(mes_up_fill),'[.]',simplify = T)[,1]
mes_dn_fill = ifelse(Mes_dn_res>0.05,0,1)
mes_dn_fill = mes_dn_fill[rowSums(mes_dn_fill)!=0,]
mes_dn_gene = str_split(rownames(mes_dn_fill),'[.]',simplify = T)[,1]
inter_motif_mes = intersect(rownames(mes_up_fill),rownames(mes_dn_fill))

pEMT_up_fill = ifelse(pEMT_up_res>0.05,0,1)
pEMT_up_fill = pEMT_up_fill[rowSums(pEMT_up_fill)!=0,]
pEMT_up_gene = str_split(rownames(pEMT_up_fill),'[.]',simplify = T)[,1]
pEMT_dn_fill = ifelse(pEMT_dn_res>0.05,0,1)
pEMT_dn_fill = pEMT_dn_fill[rowSums(pEMT_dn_fill)!=0,]
pEMT_dn_gene = str_split(rownames(pEMT_dn_fill),'[.]',simplify = T)[,1]
inter_motif_pemt = intersect(rownames(pEMT_up_fill),rownames(pEMT_dn_fill))

epi_up_fill = ifelse(epi_up_res>0.05,0,1)
epi_up_fill = epi_up_fill[rowSums(epi_up_fill)!=0,]
epi_up_gene = str_split(rownames(epi_up_fill),'[.]',simplify = T)[,1]
epi_dn_fill = ifelse(epi_dn_res>0.05,0,1)
epi_dn_fill = epi_dn_fill[rowSums(epi_dn_fill)!=0,]
epi_dn_gene = str_split(rownames(epi_dn_fill),'[.]',simplify = T)[,1]
inter_motif_epi = intersect(rownames(epi_up_fill),rownames(epi_dn_fill))

inter_motif_all_list = list( Epithelial = inter_motif_epi,
                             pEMT = inter_motif_pemt,
                             Mesenchymal=inter_motif_mes)

inter_motif_all = c(inter_motif_epi,inter_motif_mes)
inter_motif_all = c(inter_motif_all,inter_motif_pemt)
inter_motif_all_gene = unique(str_split(inter_motif_all,'[.]',simplify = T)[,1])
inter_motif_all_gene=inter_motif_all_gene[-1]
#saveRDS(inter_motif_all_list,'AS_analysis/inter_motif_all_list.RDS')


#We calculated the correlation between RBP expression and alternative splicing events, 
#and selected the top 10 strongly correlated RBPs.

EMT_result = read.csv('EMscore_tanh.csv',header = T,row.names = 1)
CTC_info =read.csv('CTC_tumor_purity_estimatResult_final.csv',header = T,row.names = 1)
EMT_result_ctc = subset(EMT_result,type2!='CTCHD'&sampleSource=='CTC')
EMT_result_ctc$SRRID=CTC_info$SRRID[match(EMT_result_ctc$sampleID,rownames(CTC_info))]
EMT_result_ductNorm = subset(EMT_result,dataSource=='GSE57973Bulk')
count_expr = read.table('CTC_exprData_raw.txt',sep = '\t',header = T,row.names = 1)

CTC_reGet = count_expr[,which(colnames(count_expr)%in% EMT_result_ctc$sampleID)]
colnames(CTC_reGet)=EMT_result_ctc$SRRID[match(colnames(CTC_reGet),EMT_result_ctc$sampleID)]
colnames(CTC_reGet)=paste0(colnames(CTC_reGet),"_",EMT_result_ctc$EMtype[match(colnames(CTC_reGet),EMT_result_ctc$SRRID)])

epi_count = CTC_reGet[,grep('Epithelial',colnames(CTC_reGet))]
pEMT_count = CTC_reGet[,grep('pEMT',colnames(CTC_reGet))]
Mes_count = CTC_reGet[,grep('Mesenchymal',colnames(CTC_reGet))]

duct_Expr = read.table('GSE57973_duct_cell_featureCount.txt',sep = '\t',header = T,row.names = 1)
Norm_DuctGet = duct_Expr
colnames(Norm_DuctGet) = paste0(str_split(colnames(Norm_DuctGet) ,'[A]',simplify = T)[,1],'_DuctNorm')

count_list = list(Epithelial = epi_count, pEMT = pEMT_count,
                  Mesenchymal = Mes_count, NormDuct = Norm_DuctGet)
#saveRDS(count_list,'CTC_Duct_exprCount.rds')

inter_expr_gene = intersect(rownames(CTC_reGet),rownames(Norm_DuctGet))
combat_Expr = bind_cols(epi_count[inter_expr_gene,],
                        pEMT_count[inter_expr_gene,],
                        Mes_count[inter_expr_gene,],
                        Norm_DuctGet[inter_expr_gene,])
phe_meta = data.frame(row.names = colnames(combat_Expr),
                      sampleType= str_split(colnames(combat_Expr),'[_]',simplify = T)[,2],
                      batch = c(rep('ctc',30),
                                rep('duct',4)))
combat_Expr_tpm = as.matrix(combat_Expr)
count_matrix<- ComBat_seq(combat_Expr_tpm,batch = phe_meta$batch)
#saveRDS(count_matrix,'AS_analysis/ctc_duct_ComBat_count_matrix.rds')

#load psiData
event_out_final = readRDS('AS_analysis/DASE_degOut_EMFinal.rds')
psi_SE = lapply(event_out_final,function(x){x[[5]]})
psi_SE = lapply(psi_SE,function(x){
  rownames(x) = x$event_name
  psi_x= x[,grep('bam',colnames(x))]
  colnames(psi_x) = paste0(str_split(colnames(psi_x),'Aligned.sortedByCoord.out.bam',simplify = T)[,1],
                           str_split(colnames(psi_x),'Aligned.sortedByCoord.out.bam',simplify = T)[,2])
  return(psi_x)
})

all(rownames(psi_SE$Epithelial)==rownames(psi_SE$pEMT))
psi_data = do.call(cbind,psi_SE)
psi_data = psi_data[,-grep('DuctNorm',colnames(psi_data))[1:8]]
colnames(psi_data)[1:30] = paste0(str_split(colnames(psi_data)[1:30],'[.]',simplify = T)[,2],"_",
                                  str_split(colnames(psi_data)[1:30],'[.]',simplify = T)[,1])
colnames(psi_data)[31:34] = str_split(colnames(psi_data)[31:34],'[.]',simplify = T)[,2]


inter_motifExp_tmp = count_matrix[rownames(count_matrix)%in% inter_motif_all_gene ,]
inter_motifExp_tmp = inter_motifExp_tmp[rowSums(inter_motifExp_tmp)!=0,]

psi_data_tmp = psi_data[,match(colnames(inter_motifExp_tmp),colnames(psi_data))]
all(colnames(psi_data_tmp)==colnames(inter_motifExp_tmp))

aimed_group = c('Epithelial','pEMT','Mesenchymal')
duct_expr = inter_motifExp_tmp[,grep('DuctNorm',colnames(inter_motifExp_tmp))]


motif_DEG_SE_cor=c()
RBP_final_save = c()
for (i in aimed_group) {
  print(i)
  exp_i =inter_motifExp_tmp[,grep(i,colnames(inter_motifExp_tmp))]
  expr_tmp = cbind(exp_i,duct_expr)
  
  x_tmp = t(expr_tmp)
  y_tmp = t(psi_data_tmp[,match(rownames(x_tmp),colnames(psi_data_tmp))])
  all(rownames(x_tmp)==rownames(y_tmp))
  motif_cor <- rcorr(x = x_tmp,  y =y_tmp,type = "pearson")
  
  cormat = round(motif_cor$r, 2)
  pmat = round(motif_cor$P, 2)
  
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  motif_cor_Res = flattenCorrMatrix(cormat, pmat)
  motif_cor_Res$rowtype = 'RBP'
  motif_cor_Res$rowtype[grep('ENSG',motif_cor_Res$row)]  ='event'  
  
  motif_cor_Res$coltype = 'RBP'
  motif_cor_Res$coltype[grep('ENSG',motif_cor_Res$column)]  ='event'
  
  motif_cor_Res$type = ifelse(motif_cor_Res$rowtype==motif_cor_Res$coltype,
                              'same','dif')
  motif_cor_ResFinal = subset(motif_cor_Res,type=='dif')
  motif_cor_ResFinal = motif_cor_ResFinal[which(abs(motif_cor_ResFinal$cor)>0.6 &
                                                  motif_cor_ResFinal$p<0.05),]
  
  motif_cor_ResFinal_top = motif_cor_ResFinal[order(motif_cor_ResFinal$p),]
  motif_cor_ResFinal_top = motif_cor_ResFinal_top[1:500,]
  RBP_final = unique(motif_cor_ResFinal_top$row)[1:10]
  event_final=unique(motif_cor_ResFinal_top$column)
  motif_cor_ResFinal_top = motif_cor_ResFinal_top[which(motif_cor_ResFinal_top$row %in% RBP_final),]

  motif_cor_ResFinal_cor =  dcast(motif_cor_ResFinal_top,column~row,value.var = 'cor')
  rownames(motif_cor_ResFinal_cor) = motif_cor_ResFinal_cor$column
  motif_cor_ResFinal_cor = motif_cor_ResFinal_cor[,-1]
  n_event = rownames(motif_cor_ResFinal_cor)
  n_RBP = colnames(motif_cor_ResFinal_cor)
  motif_cor_mat= as.data.frame(cormat[n_event,n_RBP])
  
  exon_positon = paste(str_split(rownames(motif_cor_mat),'[.]',simplify = T)[,4],':',
                       str_split(rownames(motif_cor_mat),'[.]',simplify = T)[,6],'-',
                       str_split(rownames(motif_cor_mat),'[.]',simplify = T)[,7],':',
                       str_split(rownames(motif_cor_mat),'[.]',simplify = T)[,5],'(',
                       str_split(rownames(motif_cor_mat),'[.]',simplify = T)[,3],')',
                       sep = "")
  gene_name=str_split(rownames(motif_cor_mat),'[.]',simplify = T)[,3]
  rownames(motif_cor_mat) = exon_positon
  
  pmt = data.frame(pmat[n_event,n_RBP])
  if (!is.null(pmt)){
    ssmt <- pmt< 0.01
    pmt[ssmt] <-'**'
    smt <- pmt >0.01& pmt <0.05
    pmt[smt] <- '*'
    pmt[!ssmt&!smt]<- ''
  } else {
    pmt <- F
  }
  
  p_dat = motif_cor_mat
  dim(p_dat)
  p_dat[abs(p_dat)<0.8]=0
  pmt[p_dat==0] = ''
  
  path_tmp = paste0('../../Figures/',i,'_RBP_event_correlation.pdf')
  name_tmp = paste0('Correlation between RBP expression and shared SE events',"(",i,")")
  pheatmap::pheatmap(p_dat,cluster_cols = T,cluster_rows = T,
                     show_colnames = T,angle_col = 45,
                     show_rownames = T,#annotation_row =ann_row,
                     treeheight_row = 0,treeheight_col = 0,
                     main = name_tmp,
                     na_col = "#DDDDDD",display_numbers = pmt,
                     legend_breaks=seq(-1,1,0.5),
                     color = colorRampPalette(colors = c("#c6def1",'white',"#c74787"))(100),
                     filename =path_tmp ,
                     width = 7,
                     height = 8)
  
  motif_DEG_SE_cor[[i]]=motif_cor
  RBP_final_save[[i]] = n_RBP
}

#saveRDS(motif_DEG_SE_cor,'AS_analysis/RBP_correlation.rds')
#saveRDS(RBP_final_save,'AS_analysis/RBP_correlation_choosed.rds')

#motifs that significantly regulate differential events
RBP_final_save = readRDS('AS_analysis/RBP_correlation_choosed.rds')
RBP_final_save$Epithelial

up_motif_res = list(Epithelial=epi_up_res,
                    pEMT=pEMT_up_res,
                    Mesenchymal=Mes_up_res)
dn_motif_res = list(Epithelial=epi_dn_res,
                    pEMT=pEMT_dn_res,
                    Mesenchymal=Mes_dn_res)

aimed_group = c("Epithelial","pEMT","Mesenchymal")
for (i in aimed_group) {
  RBP_fill_final = RBP_final_save[[i]]
  up_res = up_motif_res[[i]]
  dn_res = dn_motif_res[[i]]
  up_gene = str_split(rownames(up_res),'[.]',simplify = T)[,1]
  dn_gene = str_split(rownames(dn_res),'[.]',simplify = T)[,1]
  
  up_motif = rownames(up_res)[match(RBP_fill_final,up_gene)]
  dn_motif =rownames(dn_res)[match(RBP_fill_final,dn_gene)] 
  
  inter_motif = intersect(up_motif,dn_motif)
  up_data_tmp = -log10(up_res[inter_motif,])
  dn_data_tmp = -log10(dn_res[inter_motif,]) 
  rbp = str_split(rownames(up_res),'[.]',simplify = T)[,1]
  
  up_data_raw = up_data_tmp
  dn_data_raw = dn_data_tmp
  
  up_data = ifelse(up_data_raw> -log10(0.05),1,0)
  dn_data = ifelse(dn_data_raw> -log10(0.05),1,0)
  
  RightOrder <- rev(rownames(up_data))
  up_data = up_data[RightOrder,]
  dn_data = dn_data[RightOrder,]
  
  identical(rownames(up_data), rownames(dn_data))
  identical(colnames(up_data), colnames(dn_data))
  
  UpColor <- colorRamp2(breaks = c(0,1), colors = c("white","#8ea0cc"))
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
  col_title= paste0(i,' RBP motif -log10(Pval)')
  p1 <- Heatmap(dn_data, column_title = col_title,rect_gp = gpar(type = "none"),
                cell_fun = DiagFunc(up = dn_data, down =up_data ) )
  lgd <- list(Legend(title ='delta PSI>0',col_fun = DnColor, at = c(0, 1), direction = "horizontal" ),
  Legend(title = "delta PSI<0",col_fun = UpColor, at = c(0, 1), direction = "horizontal" ) )
  main_plotpath = paste0('../../Figures/Diagonal_heatmap_',i,'.pdf')
  print(main_plotpath)
  pdf(main_plotpath,height = 3,width = 4.3)
  draw(p1, annotation_legend_list = lgd,
       annotation_legend_side = "bottom",
       heatmap_legend_side = "bottom",
       merge_legend = TRUE)
  dev.off()
  
}

#MBNL1_PSI
#A sashimi plot was drawn based on the position of MBNL1 alternatively spliced exons.
#Next, we explored the expression of MBNL1 and DTU.
MBNL1_PSI_norm = SE_event_inter[,grep('DuctNorm',colnames(SE_event_inter))[1:4]]
MBNL1_PSI_ctc = SE_event_inter[,grep('sortedByCoord.out.bam$',colnames(SE_event_inter))]

MBNL1_PSI = cbind(MBNL1_PSI_ctc,MBNL1_PSI_norm)
MBNL1_PSI$event_name = rownames(MBNL1_PSI)

MBNL1_PSI_dat = melt(MBNL1_PSI,'event_name')

colnames(MBNL1_PSI_dat) = c('event_name','sampleID','PSI')
MBNL1_PSI_dat$type = str_split(MBNL1_PSI_dat$sampleID,'[.]',simplify = T)[,1]
MBNL1_PSI_dat$type[grep('_DuctNorm',MBNL1_PSI_dat$sampleID)]='Norm_duct'
MBNL1_PSI_dat$type = factor(MBNL1_PSI_dat$type,
                            levels = c('Epithelial','pEMT', 'Mesenchymal','Norm_duct'))
mycomp = list(c('Epithelial','Norm_duct'),c('pEMT','Norm_duct'),c('Mesenchymal','Norm_duct'))
pdf('../../Figures/MBNL1_psi_boxplot.pdf',height = 6,width = 4)
p1 = ggplot(MBNL1_PSI_dat,aes(x = type,y = PSI,color=type))+
  geom_boxplot()+
  geom_signif( comparisons = mycomp, 
    test = "wilcox.test",y_position= c(1.5,1.4,1.3),
    map_signif_level = TRUE, textsize = 4 )+
  facet_wrap('event_name')
print(p1)
dev.off()

#MBNL1 expr
count_list = readRDS('CTC_Duct_exprCount.rds')
intersect_gene= intersect(rownames(count_list$Epithelial),
                          rownames(count_list$NormDuct))
count_list_expr = lapply(count_list,function(x){
  x = x[match(intersect_gene,rownames(x)),]
  return(x)
})

count_expr = do.call(cbind,count_list_expr)
ids=annoGene(rownames(count_expr),'SYMBOL','human')
gene_ins_proteinCoding = unique(ids$SYMBOL[which(ids$biotypes=='protein_coding')])
count_expr_pro = count_expr[which(rownames(count_expr)%in%gene_ins_proteinCoding),]
phe_meta = data.frame(row.names = colnames(count_expr_pro),
                      sampleType= str_split(colnames(count_expr_pro),'[_]',simplify = T)[,2],
                      batch = c(rep('ctc',30), rep('duct',4)))
combat_Expr_tpm = as.matrix(count_expr_pro)
countData<- ComBat_seq(combat_Expr_tpm,batch = phe_meta$batch)

#Calculate TPM
countToTpm <- function(counts, effLen)
{
  rate <- log(counts+1) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
mergedataTP = countData

exons_length1 = read.csv('gencode.v38.annotation_v2.csv',header = T,row.names = 1)
le = exons_length1$Length[match(rownames(mergedataTP),exons_length1$gene_name)]

mergedataTP = mergedataT
eff_length = le
tpms <- apply(mergedataTP,2,countToTpm,eff_length)
tpms_scale = log2(tpms+1)
mbnl1_expr = tpms_scale[which(rownames(tpms_scale)=='MBNL1'),]
mbnl1_expr_dat = data.frame(sampleID = names(mbnl1_expr),
                            TPM = mbnl1_expr,
                            type = str_split(names(mbnl1_expr),'[_]',simplify = T)[,2])

mbnl1_expr_dat$type = factor(mbnl1_expr_dat$type,
                             levels = c('Epithelial','pEMT', 
                                        'Mesenchymal','DuctNorm'))
mbnl1_expr_dat = mbnl1_expr_dat[order(mbnl1_expr_dat$type),]
mbnl1_expr_dat$sampleID = factor(mbnl1_expr_dat$sampleID,
                                 levels = mbnl1_expr_dat$sampleID)
pdf('../../Figures/MBNL1_exprTPM_barplot.pdf',width = 9,height = 2.7)
p1 =ggplot(mbnl1_expr_dat,aes(x =sampleID,y=TPM,fill =type  ))+
  geom_bar(stat = 'identity',width = 0.6,position = 'dodge')
print(p1)
dev.off()


#DTU
mySwitchList = readRDS('AS_analysis/isoformSwitchResult/mySwitchList_importRdata.rds')

mbnl1_switch <- preFilter(
  switchAnalyzeRlist = mySwitchList,
  geneExpressionCutoff = 10,
  isoformExpressionCutoff = 3,
  removeSingleIsoformGenes = TRUE
)

all_switchAnalyzed <- isoformSwitchTestSatuRn( mbnl1_switch, dIFcutoff = 0.1,
                                                 reduceToSwitchingGenes=T )

### exampleSwitchListAnalyzedSubset is created above

pdf('../../Figures/IsoUsage_MBNL1.pdf',width = 6,height = 9)
p1 = switchPlotIsoUsage(all_switchAnalyzed, gene = 'MBNL1',
                        condition1 = 'Epithelial',condition2 = 'DuctNorm')
p2 =switchPlotIsoUsage(all_switchAnalyzed, gene = 'MBNL1',localTheme =theme_bw(),
                       condition1 = 'pEMT',condition2 = 'DuctNorm')
p3 = switchPlotIsoUsage(all_switchAnalyzed, gene = 'MBNL1',
                        condition1 = 'Mesenchymal',condition2 = 'DuctNorm')
print(p1/p2/p3)
dev.off()

#saveRDS(all_switchAnalyzed,'AS_analysis/all_switchAnalyzed')

#Analyze exon skipping conditions and classify them into completely retained (included), 
#completely skipped (excluded) and partially retained (middle)
event_out_final = readRDS('AS_analysis/DASE_degOut_EMFinal.rds')
psi_SE = lapply(event_out_final,function(x){x[[5]]})
psi_SE_mode = lapply(psi_SE,function(x){
  ctc_mean  = x$ctc_mean
  duct_mean =x$duct_mean
  mode_ctc = ifelse(ctc_mean>0.95,'included',ifelse(ctc_mean<0.05,'excluded','middle')) 
  mode_duct = ifelse(duct_mean>0.95,'included',ifelse(duct_mean<0.05,'excluded','middle')) 
  mode_dat = x[,1:10]
  mode_dat$mode_ctc = mode_ctc
  mode_dat$mode_duct= mode_duct
  return(mode_dat)
})
psi_SE_mode_rbind = do.call(rbind,psi_SE_mode)
psi_SE_mode_rbind$emtype = str_split(rownames(psi_SE_mode_rbind),'[.]',simplify = T)[,1]
#write.csv(psi_SE_mode_rbind,'AS_analysis/psi_se_mode.csv')
#saveRDS(psi_SE_mode,'AS_analysis/psi_se_mode.rds')

all(psi_SE_mode$Epithelial$mode_duct==psi_SE_mode$pEMT$mode_duct)
psi_SE_mode_dat = data.frame(
  Norm = psi_SE_mode$Epithelial$mode_duct,
  Epithelial = psi_SE_mode$Epithelial$mode_ctc,
  pEMT = psi_SE_mode$pEMT$mode_ctc,
  Mesenchymal = psi_SE_mode$Mesenchymal$mode_ctc)

mydata = psi_SE_mode_dat
hd_2_epi = ifelse(mydata$Norm==mydata$Epithelial,'No_mode_change','mode_change')
epi_2_pEMT = ifelse(mydata$Epithelial==mydata$pEMT,'No_mode_change','mode_change')
pEMT_2_mes = ifelse(mydata$pEMT==mydata$Mesenchymal,'No_mode_change','mode_change')
data_bar = data.frame(hd_2_epi = as.numeric(table(hd_2_epi)),
                      epi_2_pEMT=as.numeric(table(epi_2_pEMT)),
                      pEMT_2_mes = as.numeric(table(pEMT_2_mes)),
                      type = c('mode_change','No_mode_change'),
                      row.names = c('mode_change','No_mode_change'))
color =list(
  "No_mode_change"="grey70",
  "mode_change"= "#8B9CC5"
)
data = melt(data_bar,'type')
data <- data %>%
  group_by(variable) %>%
  mutate(pct=prop.table(value))
pdf('../../Figures/se_barplot_mode.pdf',width = 4,height = 4)
p1 = ggplot(data, aes( x = variable,y=pct,fill = type))+
  geom_bar(position="fill", stat = "identity", width = 0.4) +
  facet_wrap('variable',scales = 'free_x')
print(p1)
dev.off()

hd_2_epi= ifelse(mydata$Norm==mydata$Epithelial,'No_mode_change','mode_change')
epi_2_pEMT = ifelse(mydata$Epithelial==mydata$pEMT,'No_mode_change','mode_change')
pEMT_2_mes = ifelse(mydata$pEMT==mydata$Mesenchymal,'No_mode_change','mode_change')
table(hd_2_epi)
table(pEMT_2_mes)
epi_mode_change = mydata[which(hd_2_epi=='mode_change'),1:2]
pEMT_mode_change = mydata[which(epi_2_pEMT=='mode_change'),2:3]
mes_mode_change = mydata[which(pEMT_2_mes=='mode_change'),3:4]

epi_mode_change_result = as.data.frame(table(epi_mode_change$Epithelial))
pEMT_mode_change_result =  as.data.frame(table(pEMT_mode_change$pEMT))
mes_mode_change_result = as.data.frame(table(mes_mode_change$Mesenchymal))

epi_mode_change_result$type = 'Epithelial'
pEMT_mode_change_result$type = 'pEMT'
mes_mode_change_result$type = 'Mesenchymal'

bar_data_mod  = bind_rows(epi_mode_change_result,
                          pEMT_mode_change_result,
                          mes_mode_change_result)
bar_data_mod <- bar_data_mod %>%
  group_by(type) %>%
  mutate(pct=prop.table(Freq))

pdf('../../Figures/se_barplot_mode_change.pdf',width = 4,height = 4)
p1 = ggplot(bar_data_mod, aes( x = type,y=pct,fill = Var1))+
  geom_bar(position="fill", stat = "identity", width = 0.4) +
  facet_wrap('type',scales = 'free_x')
print(p1)
dev.off()

#find mode-changed event
psi_SE_mode = readRDS('AS_analysis/psi_se_mode.rds')
data_ex =lapply(psi_SE_mode,function(x){x[which(x$mode_duct !='middle'),]}) 
data_ex_tmp = data.frame(
  event  = rownames(data_ex$Epithelial),
  norm =data_ex$Epithelial$mode_duct,
  epi = data_ex$Epithelial$mode_ctc,
  pEMT= data_ex$pEMT$mode_ctc,
  mes = data_ex$Mesenchymal$mode_ctc)
data_ex$Epithelial[2,]
#'CYB5R3.chr22:42646440-42647701'
all_switchAnalyzed = readRDS('AS_analysis/all_switchAnalyzed')

CYB5R3_switchAnalyzed = subsetSwitchAnalyzeRlist(all_switchAnalyzed,
                                                 all_switchAnalyzed$isoformFeatures$gene_name=='CYB5R3')


CYB5R3_switchAnalyzed_ORF <- extractSequence(
  CYB5R3_switchAnalyzed, 
  pathToOutput = 'AS_analysis/isoformSwitchResult/cyb5r3/',
  removeShortAAseq =F,removeORFwithStop=F,
  outputPrefix='cellSurface_isoform',
  writeToFile=T
)
CYB5R3_switchAnalyzed_ORF <- analyzeDeepTMHMM(
  switchAnalyzeRlist   = CYB5R3_switchAnalyzed,
  pathToDeepTMHMMresultFile = "AS_analysis/isoformSwitchResult/cyb5r3/TMRs.gff3",
  showProgress=FALSE
)
pdf('../../Figures/CYB5R3_switchAnalyzedEpi.pdf',width = 7,height = 4)
p1 = switchPlot(CYB5R3_switchAnalyzed_ORF, gene = 'CYB5R3',
                condition1 = 'Epithelial',condition2 = 'DuctNorm')
print(p1)
dev.off()

pdf('../../Figures/CYB5R3_switchAnalyzedpEMT.pdf',width = 7,height = 4)
p2 =switchPlot(CYB5R3_switchAnalyzed_ORF, gene = 'CYB5R3',
               condition2 = 'pEMT',condition1 = 'DuctNorm')
print(p2)
dev.off()

pdf('../../Figures/CYB5R3_switchAnalyzedMes.pdf',width = 7,height = 4)
p3 =switchPlot(CYB5R3_switchAnalyzed_ORF, gene = 'CYB5R3',
               condition2 = 'Mesenchymal',condition1 = 'DuctNorm')
print(p3)
dev.off()

#Obtain sequence information of events on chromosomes based on their location
get_event_seq_result = read.table('AS_analysis/get_sequence_dataFrame',
                                  sep = '\t',header = F)
colnames(get_event_seq_result) = c('exon_position','sequence')
get_event_seq_result$event_position = str_split(get_event_seq_result$exon_position,
                                                '/r',simplify = T)[,1]
all(rownames(event_position)==get_event_seq_result$event_position)
get_event_seq_result$event_name = events$event_name

mydata_event_position  = cbind(event_position$direction,event_position$exon_length)
mydata_event_position  = cbind(mydata_event_position,get_event_seq_result)

colnames(mydata_event_position)[1] = 'strand'
colnames(mydata_event_position)[2] = 'exon_length'

mydata_event_position = mydata_event_position[,c("event_name","exon_length","event_position",'strand',"sequence")]
mydata_event_position$hd_type = psi_SE$Epithelial$duct_mean
mydata_event_position$epi_type= psi_SE$Epithelial$ctc_mean
mydata_event_position$pEMT_type= psi_SE$pEMT$ctc_mean
mydata_event_position$mes_type= psi_SE$Mesenchymal$ctc_mean

pct_TAG = as.numeric(str_detect(mydata_event_position$sequence,'TAG'))
pct_TAA = as.numeric(str_detect(mydata_event_position$sequence,'TAA'))
pct_TGA = as.numeric(str_detect(mydata_event_position$sequence,'TGA'))
pct_data = data.frame(TAG = pct_TAG,TAA= pct_TAA,TGA = pct_TGA)
rownames(pct_data) = mydata_event_position$event_name

pct_data = cbind(pct_data,mydata_event_position)
pct_data$chr = str_split(pct_data$event_position,'[:]',simplify = T)[,1]
pct_data$chr = factor(pct_data$chr,levels = c(paste0('chr',1:22),'chrX'))
colnames(pct_data)
#write.table(pct_data,'AS_analysis/PTC_result.txt',quote = F,row.names = T,sep = '\t')

#pheatmap(pct_data$exon_length,color = col_length) #get annotation bar
#PTC Statistics
pct_data$No.PCT = rowSums(pct_data_plot)

pct_data[,c("TAG","TAA","TGA","No.PCT")]

pct_bar = as.data.frame(table(pct_data$No.PCT))
value = pct_bar$Var1
pdf('../../Figures//PTC_Number_barplot.pdf',width = 4,height = 4)
p1 = ggplot(pct_bar,aes(x = Var1,y =Freq))+
  geom_bar(stat = 'identity')+
  xlab('Number of PTC')+
  theme_classic()
print(p1)
dev.off()
#dencity
pdf('../../Figures//PTC_length_density.pdf',width = 4,height = 4) 
pct_point = pct_data[,c("No.PCT","exon_length")]
pct_point$exon_length[pct_point$exon_length>500]=500
pct_point$No.PCT = as.factor(pct_point$No.PCT)
p1 = ggplot(pct_point, aes(x=exon_length,fill =No.PCT))+ 
  geom_density(alpha=0.8)+
  theme_classic()
print(p1)
dev.off()

#Obtain PTC upstream 40nt and downstream 60nt sequences
#bash get_seq_PTC.sh >get_100ntPTC_result.fa
#python get_ptc_100nt.py #get_100ntPTC_dataFrame
ptc_100nt_set = read.table('AS_analysis//get_100ntPTC_dataFrame',sep = '\t')
pdf('../../Figures/seqlogo.pdf',width = 12,height = 3)
p1 = ggseqlogo(ptc_100nt_set$V2,method = 'prob')+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = .5))
print(p1)
dev.off()