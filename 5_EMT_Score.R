rm(list = ls())
library(limma)
library(edgeR)
library(stringr)
library(stringi)
library(DESeq2)
library(singscore)

#singleScore analysis
#load data
degseq2_res = readRDS('DEGseq2_diffRaw.rds')
res_CTC = degseq2_res$ctc_deg
res_epiCancer = degseq2_res$duct_deg

degseq2_resSig = readRDS('DEGseq2_diffSig.rds')
res_CTCSig = degseq2_resSig$CTCSig_lfc1
res_EPISig = degseq2_resSig$DuctSig_lfc1

resEM = function(resSig,resRaw,expr_m){
  epi_up = resSig$gene[which(resSig$type=='Epi' &resSig$change=='Up-regulation')]
  epi_dn = resSig$gene[which(resSig$type=='Epi' &resSig$change=='Down-regulation')]
  mes_up = resSig$gene[which(resSig$type=='Mes' &resSig$change=='Up-regulation')]
  mes_dn = resSig$gene[which(resSig$type=='Mes' &resSig$change=='Down-regulation')]
  
  stable_genes = setdiff(rownames(resRaw),resSig$gene)[1:5]

  EMT_score = function(expr,upGene,dnGene,stable_genes){#
    measured <-  unique(c(stable_genes,upGene, dnGene))# 
    small_expr  <-  expr[measured, ]
    rankData_st <-  rankGenes(small_expr,stableGenes = stable_genes)# 
    scoredf_st <- simpleScore(rankData_st, upSet = upGene, downSet = dnGene)
    return(scoredf_st)
  }
  epi_score = EMT_score(expr_m,epi_up,epi_dn,stable_genes)
  mes_score = EMT_score(expr_m,mes_up,mes_dn,stable_genes)
  
  epi_final = epi_score$TotalScore
  mes_final = mes_score$TotalScore
  final_EM = data.frame( epi_score = epi_final,mes_score = mes_final,
                      EM_score = epi_final/mes_final,
                      sampleID = rownames(epi_score))
  resRet = list(epi_score,mes_score,final_EM)
  names(resRet)=c('epi_score','mes_score','final_EM')
  return(resRet)
}
#singleScore Result####
ctc_EM = resEM(resRaw = res_CTC,
               resSig = res_CTCSig,
               expr_m = count_expr_raw)
epi_EM = resEM(resRaw = res_epiCancer,
               resSig = res_EPISig,
               expr_m = combat_Expr)
#hexbin plot####
pdf('../../Figures/hexbin_plot_major.pdf',width = 8,height = 4)
p1 = plotScoreLandscape(ctc_EM$epi_score, isInteractive = F,
                   ctc_EM$mes_score, 
                   scorenames = c('CTC-EPI','CTC-MES'),hexMin = 10)
p2 =plotScoreLandscape(epi_EM$epi_score, epi_EM$mes_score, 
                   scorenames = c('Duct-EPI','Duct-MES'),hexMin = 10)
print(p1+p2)
dev.off()
#plt CTC result####
plot_dat_ctc = ctc_EM$final_EM
plot_dat_ctc$sampleSource = 'CTC'
estimat_score_final$sampleID = rownames(estimat_score_final)
plot_dat_ctc$type = estimat_score_final$type[match(plot_dat_ctc$sampleID,rownames(estimat_score_final))]
plot_dat_ctc$type  = factor(plot_dat_ctc$type,levels = c('locPDAC','metPDAC','HD')) 


#plt all result####
all_EM = rbind(plot_dat_sc,plot_dat_ctc)
plot_dat_sc$EMtanh = tanh(plot_dat_sc$EM_score/2)
plot_dat_ctc$EMtanh = tanh(plot_dat_ctc$EM_score/2)
plot_data = all_EM[order(all_EM$EM_score,decreasing = T),]
plot_data$sampleID = factor(plot_data$sampleID,levels = plot_data$sampleID)
plot_data$type2= as.factor(paste0(plot_data$sampleSource,plot_data$type))
plot_data$EMtanh = tanh(plot_data$EM_score/2)
#save result####
write.csv(plot_data,'EMscore_tanh2.csv')
write.csv(plot_data,'../../Tables/EMscore_tanh2.csv')

plot_data = read.csv('EMscore_tanh2.csv',header = T,row.names = 1)

p = ggplot(plot_data, aes(x=EMtanh))+ 
  geom_density(aes(fill = sampleSource), alpha=0.4)+theme_bw() +mytheme
a = plot_dat_sc$EMtanh[-24]
b = plot_dat_ctc$EMtanh
scdensity = as.data.frame(density(a,from=-1,to = 1)[c("x", "y")])
ctcdensity =as.data.frame(density(b,from=-1,to = 1)[c("x", "y")])
comp = ctcdensity$y>scdensity$y
diff(comp)
cross =c(NA, diff(comp))
round(ctcdensity[which(cross != 0), ],2)
x1 = -0.05
y1 = 0.13
x2 = 0.76
y2 = 0.54
dev.new()
#plt density####
pdf('../../Figures/EMscore_density.pdf',width = 3,height = 3)
p1 = ggplot(plot_data, aes(x=EMtanh))+ 
  geom_density(aes(fill = sampleSource), alpha=0.4)+
  geom_vline(xintercept = c(-0.05,0.76),lty="dashed")+theme_bw() +mytheme
print(p1)
dev.off()

plot_data$type2 = factor(plot_data$type2,
                         levels = c('CTCHD','CTClocPDAC','CTCmetPDAC',
                                    'DuctNorm','DuctlocPDAC','DuctmetPDAC'))
#plt corelation####
pdf('../../Figures/EMscore_summary_CorLmplot.pdf',width = 5,height = 5)
p1 = ggplot(plot_data, aes(y=mes_score,x=epi_score,
                           color =type2) )+ 
  geom_point(size=1.2) +
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson",
           label.x = -0.75, label.y = 0.27,color='black')+
  scale_color_manual(values = c(
    'DuctlocPDAC'="#F8C77D",
    'DuctmetPDAC'='#DC5B64',
    'DuctNorm'='#755032',
    "CTCmetPDAC" ="#5CA4DA"  ,
    'CTClocPDAC' ='#87C46A' ,
    'CTCHD' ='#8D709B'))+
  facet_wrap(~type2)+theme_bw()+mytheme
print(p1)
dev.off()

#plt data_source_plot####
plot_data$EMtype = ifelse(plot_data$EMtanh<x1,'Mesenchymal',
                          ifelse(plot_data$EMtanh>x2,'Epithelial','pEMT'))
plot_data$EMtype = factor(plot_data$EMtype,levels = c('Epithelial','pEMT','Mesenchymal'))
plot_data$dataSource =phe_meta$batch[match(plot_data$sampleID,phe_meta$Patient_id)] 
plot_data$dataSource[grep('CTC',plot_data$sampleSource)] ='GSE144561'#CTC dataset
#write.csv(plot_data,'EMscore_tanh2.csv')
#write.csv(plot_data,'../../Tables/EMscore_tanh2.csv')

#plt corelation_emtgene####
count_expr = read.table('CTC_exprData_raw.txt',sep = '\t',header = T,row.names = 1)
estimat_score_final = read.csv('CTC_tumor_purity_estimatResult_final.csv',header = T,row.names = 1)
emt_signature = read.csv('emt_gene_signature.csv',header = T,row.names = 1)
estimat_score_final$sampleSource = ifelse(estimat_score_final$sampleType=='NormPBMC','NormPBMC','CTC')

combat_Expr = read.table('combat_Expr_emtExpr.txt',sep = '\t',header = T,row.names = 1)
pb_metadat = read.table('combat_meta_all.txt',sep = '\t',header = T,row.names = 1)
head(pb_metadat)
pb_metadat$sampleType = ifelse(pb_metadat$type=='Norm','Norm','EpiCancer')

EMT_res = read.csv('../../Tables/EMscore_tanh2.csv',header = T,row.names =1 )
EMT_deg_Sig= readRDS('DEGseq2_diffSig.rds')
res_CTCSig = EMT_deg_Sig$CTCSig_lfc1
res_EPISig = EMT_deg_Sig$DuctSig_lfc1
res_CTCSig$cellSource = 'CTC'
res_EPISig$cellSource = 'Duct'

resSig_all = rbind(res_CTCSig,res_EPISig)

epi_ctc = resSig_all$gene[which(resSig_all$type=='Epi' &resSig_all$cellSource=='CTC')]
epi_duct = resSig_all$gene[which(resSig_all$type=='Epi' &resSig_all$cellSource=='Duct')]
mes_ctc = resSig_all$gene[which(resSig_all$type=='Mes' &resSig_all$cellSource=='CTC')]
mes_duct = resSig_all$gene[which(resSig_all$type=='Mes' &resSig_all$cellSource=='Duct')]

epi_shared = intersect(epi_ctc,epi_duct)
mes_shared = intersect(mes_ctc,mes_duct)

epi_shared = c("TACSTD2", "EHF","OR7E14P")
mes_shared =  c("VCAM1","DDR2","SFRP1" ,  "ECM2")

resSig_all_epiShared = subset(resSig_all,gene %in% epi_shared)
resSig_all_mesShared = subset(resSig_all,gene %in% mes_shared)

#test_gene = epi_shared
test_gene = mes_shared

ctc_expr_epiShared =count_expr[match(test_gene,rownames(count_expr)),]
sc_expr_epiShared= combat_Expr[match(test_gene,rownames(combat_Expr)),]

expr_epiShared = cbind(ctc_expr_epiShared,sc_expr_epiShared)
colnames(expr_epiShared)
expr_epiShared = expr_epiShared[,match(EMT_res$sampleID,colnames(expr_epiShared))]


expr_epiShared_t = as.data.frame(t(expr_epiShared))
all(rownames(expr_epiShared_t)==EMT_res$sampleID)
expr_epiShared_t$EMT_tanh = EMT_res$EMtanh
expr_epiShared_t$sampleType = EMT_res$sampleSource
expr_epiShared_t = na.omit(expr_epiShared_t)

expr_epiShared_ctc = subset(expr_epiShared_t,sampleType=='CTC')
expr_epiShared_Duct = subset(expr_epiShared_t,sampleType=='Duct')



library(Hmisc)

Dct_cor <- rcorr(as.matrix(scale(expr_epiShared_Duct[,-ncol(expr_epiShared_Duct)])), 
                 type = "pearson")
Dct_res = round(Dct_cor$r, 2)

CTC_cor <- rcorr(as.matrix(scale(expr_epiShared_ctc[,-ncol(expr_epiShared_ctc)])), 
                 type = "pearson")
CTC_res = round(CTC_cor$r, 2)


res = data.frame(Duct = Dct_res[-nrow(Dct_res),ncol(Dct_res)],
                 CTC = CTC_res[-nrow(CTC_res),ncol(CTC_res)],
                 row.names = rownames(CTC_res)[-nrow(CTC_res)])
res[res>0.3]=0.3
res[res<0.1]=0.1


pmt = data.frame(Duct = Dct_cor$P[-nrow(Dct_res),ncol(Dct_res)],
                 CTC = CTC_cor$P[-nrow(CTC_res),ncol(CTC_res)],
                 row.names = rownames(CTC_res)[-nrow(CTC_res)])
if (!is.null(pmt)){
  ssmt <- pmt< 0.01
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}
pdf('../../Figures/EMscore_cor_Expr_epithelial.pdf',height = 1.5,width = 2)
p1  = pheatmap::pheatmap(res,cluster_cols = F,cluster_rows = F,angle_col = 0,
                         display_numbers = pmt,main = 'Epithelial gene',
                         legend_breaks=seq(-0.1,0.3,0.1),
                         color = colorRampPalette(colors = c("#c6def1", 'white',"#c74787"))(100))
print(p1)
dev.off()

pdf('../../Figures/EMscore_cor_Expr_mesenchymal.pdf',height = 1.5,width = 2)
p1  = pheatmap::pheatmap(res,cluster_cols = F,cluster_rows = F,angle_col = 0,
                         display_numbers = pmt,main = 'Mesenchymal gene',
                         legend_breaks=seq(-0.1,0.3,0.1),
                         color = colorRampPalette(colors = c("#c6def1", 'white',"#c74787"))(100))
print(p1)
dev.off()

#survival####
#tcga####
library(estimate)
tcga_path_ex = 'E:/307/data/TCGA/rna expr/rna_sd/'
tcga_path_cli = 'E:/307/data/TCGA/clinical1/'
paad_ex=read.table(paste0(tcga_path_ex,'PAAD_HiSeqV2_PANCAN/HiSeqV2_PANCAN'),header = T,row.names =1,sep = '\t')
paad_ex[1:4,1:5]
paad_cli = read.table(paste0(tcga_path_cli,'paad_clinicalMatrix'),header = T,sep = '\t')
rownames(paad_cli) = gsub('[-]','.',paad_cli$sampleID)
paad_cli[1:2,]
paad_cli_need = paad_cli[,match(c('sample_type_id','bcr_patient_barcode',
                            "vital_status",
                            'days_to_birth',
                            "days_to_last_followup",
                            "days_to_death",
                            "gender",
                            "sample_type",
                            'alcoholic_exposure_category','tobacco_smoking_history'),colnames(paad_cli))]
#paad_cli_need = paad_cli_need[!is.na(paad_cli_need$days_to_last_followup),]

write.table(paad_ex,'tcag_paad.txt',sep = '\t',col.names = T,row.names = T,quote = F)
write.table(paad_cli_need,'tcag_paad_clinical.txt',sep = '\t',col.names = T,row.names = T,quote = F)

paad_cli_need = read.table('tcag_paad_clinical.txt',
                           sep = '\t',header = T,row.names = 1)#196
paad_cli_need$alcoholic_exposure_category[paad_cli_need$alcoholic_exposure_category=='']='empty'
paad_cli_need[is.na(paad_cli_need)]='empty'

paad_cli_rm = paad_cli_need[which(paad_cli_need$alcoholic_exposure_category != 'Daily Drinker'),]
paad_cli_rm = paad_cli_rm[which(paad_cli_rm$tobacco_smoking_history %in% c('empty','1','2','3')),]
paad_ex_rm = paad_ex[,na.omit(match(rownames(paad_cli_rm),colnames(paad_ex)))]
dim(paad_ex_rm)#20530   134
paad_cli_rm =paad_cli_need[na.omit(match(colnames(paad_ex_rm),rownames(paad_cli_need))),]
write.table(paad_ex_rm,'tcag_paad_filtered.txt',sep = '\t',col.names = T,row.names = T,quote = F)
write.table(paad_cli_rm,'tcag_paad_filtered__clinical.txt',sep = '\t',col.names = T,row.names = T,quote = F)

vld_tcga= resEM(resRaw = res_epiCancer,
                resSig = res_EPISig,
                expr_m = paad_ex_rm)
plotScoreLandscape(vld_tcga$epi_score, isInteractive =F,
                   vld_tcga$mes_score, scorenames = c('TCGA-EPI','TCGA-MES'),hexMin = 10)
plot_dat_vld = vld_tcga$final_EM
plot_dat_vld =plot_dat_vld[plot_dat_vld$EM_score!=Inf,]
plot_dat_vld$sampleSource = 'TCGA'
row.names(plot_dat_vld) = gsub('[.]','-',plot_dat_vld$sampleID)
plot_dat_vld$type = as.numeric(str_split(rownames(plot_dat_vld),'[-]',simplify = T)[,4])
table(plot_dat_vld$type)
plot_dat_vld$type = ifelse(plot_dat_vld$type<= 9,'tumor','norm')
plot_dat_vld$days_to_birth = paad_cli_rm$days_to_birth[match(rownames(plot_dat_vld),paad_cli_rm$sampleID)]
plot_dat_vld$EMtanh = round(tanh(plot_dat_vld$EM_score/2),2)
plot_dat_vld_norm = subset(plot_dat_vld,type=='norm')

outlier_values <- boxplot.stats(plot_dat_vld_norm$EMtanh)$out
ggplot(plot_dat_vld_norm,mapping = aes(x= type,y = EMtanh,fill=type))+
  geom_boxplot()+
  geom_hline(yintercept = c(x1,x2),lty="dashed")+
  theme_bw()+mytheme

plot_dat_vld$EMtype= ifelse(plot_dat_vld$EMtanh<x1,'Mesenchymal',
                            ifelse(plot_dat_vld$EMtanh>x2,'Epithelial','pEMT'))
plot_dat_vld$EMtype =factor(plot_dat_vld$EMtype,levels = c('Epithelial','pEMT','Mesenchymal'))

pdf('../../Figures/EMscore_Validate_TCGA.pdf',width = 8,height = 5)
p2 = ggplot(plot_dat_vld,mapping = aes(x= EMtype,y = EMtanh,fill=type))+
  geom_boxplot( position = position_dodge2(preserve = "single") )+
  geom_hline(yintercept = c(x1,x2),lty="dashed")+
  theme_bw()+mytheme
p1 = plotScoreLandscape(vld_tcga$epi_score, isInteractive =F,
                        vld_tcga$mes_score, scorenames = c('TCGA-EPI','TCGA-MES'),hexMin = 10)
print(p1+p2)
dev.off()
#plot_dat_vld$`Tumor purity` = estimat_score$TumorPurity[match(plot_dat_vld$sampleID,rownames(estimat_score))]
write.csv(plot_dat_vld,'EMscore_validate_TCGA.csv')
write.csv(plot_dat_vld,'../../Tables/EMscore_validate_TCGA.csv')
dim(plot_dat_vld)
#os tcga####
paad_cli_rm = read.table('tcag_paad_filtered__clinical.txt',sep = '\t',
                         header = T,row.names = 1)
paad_cli_rm$EMtype = plot_dat_vld$EMtype[match(rownames(paad_cli_rm),plot_dat_vld$sampleID)]
paad_cli_rm$EMtanh = plot_dat_vld$EMtanh[match(rownames(paad_cli_rm),plot_dat_vld$sampleID)]
paad_cli_em = paad_cli_rm[!is.na(paad_cli_rm$EMtype),]
paad_cli_em = na.omit(paad_cli_em)

paad_cli_em$OS <- as.integer(
  ifelse( is.na(paad_cli_em$days_to_death),
          paad_cli_em$days_to_last_followup,
          paad_cli_em$days_to_death))
paad_cli_em$event=ifelse(paad_cli_em$vital_status=='LIVING',0,1)

#median(paad_cli_em$EMtanh)
paad_cli_em$EMtype2 = ifelse(paad_cli_em$EMtanh>median(paad_cli_em$EMtanh),
                             'High(EM_score)','Low(EM score)')
sfit <- survfit(Surv(OS, event)~EMtype2, data=paad_cli_em)
print(sfit)
surv_pvalue(sfit,method = 'survdiff')
pdf('../../Figures/EMscore_Validate_TCGA_survival.pdf',height = 5,width = 5)
p1 = ggsurvplot(sfit, data = paad_cli_em,surv.median.line = "hv",pval = TRUE, 
           ggtheme = theme_bw()+mytheme, 
           palette = c("#E7B800", "#2E9FDF",'red'))
print(p1)
dev.off()
#validate GSE40174####
emt_signature = read.csv('/home/lxxiao/xiaolixing/pdac/10x/emt_gene_signature.csv',header = T,row.names = 1)
validate_pdac3 = read.table('GSE40174_human_processed_data.txt',
                              header = T,row.names = 1,sep = '\t')
gene_name = str_split(rownames(validate_pdac3),'[.]',simplify = T)[,1]
validate_pdac3 = cbind(gene_name,validate_pdac3) 
validate_pdac3 = validate_pdac3[!duplicated(validate_pdac3$gene_name),]
validate_pdac3 = validate_pdac3[,-c(1:2)]  
validate_pdac3 = validate_pdac3[which(str_split(validate_pdac3$patient,'[.]',simplify = T)[,2]=="EpCAM"),]
validate_pdac3_phe = data.frame(sampleID =colnames(validate_pdac3),
                                type = str_split(colnames(validate_pdac3),'[_]',simplify = T)[,1],
                                sampleSource =str_split(colnames(validate_pdac3),'[.]',simplify = T)[,2] )


emM_PDAC3 = validate_pdac3[,which(validate_pdac3_phe$sampleSource=="EpCAM")]
vld_pdac3= resEM(resRaw = res_CTC,
       resSig = res_CTCSig,
      expr_m = emM_PDAC3)
pdf('../../Figures/EMscore_Validate_pdac3.pdf',width = 10,height = 10)
p1 = plotScoreLandscape(vld_pdac3$epi_score, isInteractive = F,
                   vld_pdac3$mes_score, 
                   scorenames = c('CTC-EPI','CTC-MES'),hexMin = 10)
plot_dat_vld = vld_pdac3$final_EM
plot_dat_vld$sampleSource = 'CTC'
plot_dat_vld$type = validate_pdac3_phe$type[match(plot_dat_vld$sampleID,
                                                  validate_pdac3_phe$sampleID)]
table(plot_dat_vld$type)
plot_dat_vld$EMtanh = tanh(plot_dat_vld$EM_score/2)
OutVals = boxplot(plot_dat_vld$EM_score)$out
plot_dat_vld_rmOt = plot_dat_vld[!plot_dat_vld$EM_score %in% OutVals,]
p2 = ggplot(plot_dat_vld_rmOt,mapping = aes(x= type,y = EMtanh,fill=type))+
  geom_boxplot( position = position_dodge2(preserve = "single") )+
  geom_hline(yintercept = c(x1,x2),lty="dashed")+
  theme_bw()+mytheme
print(p1+p2)
dev.off()

#validate GSE114704#### 
validate_pdac4 = read.table('GSE114704_allSamples_htseqcount_readcouts.tsv.gz',
                            header = T,row.names = 1,sep = '\t')
validate_pdac4_phe = data.frame(sampleID = colnames(validate_pdac4),
                                type = 'CTC')
validate_pdac4_phe$type[grep('CTC',colnames(validate_pdac4))]='CTC'
validate_pdac4_phe$type[grep('LM',colnames(validate_pdac4))]='Liver_metastasis'
validate_pdac4_phe$type[grep('PT',colnames(validate_pdac4))]='Primary_Tumor'

validate_pdac4_ctc = validate_pdac4[,validate_pdac4_phe$sampleID[which(validate_pdac4_phe$type=='CTC')]]
validate_pdac4_tumor = validate_pdac4[,validate_pdac4_phe$sampleID[which(validate_pdac4_phe$type!='CTC')]]

vld_pdac4_ctc= resEM(resRaw = res_CTC,
                 resSig = res_CTCSig,
                 expr_m = validate_pdac4_ctc)
vld_pdac4_tumor= resEM(resRaw = res_epiCancer,
                     resSig = res_EPISig,
                     expr_m = validate_pdac4_tumor)
pdf('../../Figures/EMscore_Validate_pdac4.pdf',width = 12,height = 10)

p1=plotScoreLandscape(vld_pdac4_ctc$epi_score, isInteractive = F,
                   vld_pdac4_ctc$mes_score, 
                   scorenames = c('CTC4-EPI','CTC4-MES'),hexMin = 10)
p2=plotScoreLandscape(vld_pdac4_tumor$epi_score, isInteractive = F,
                   vld_pdac4_tumor$mes_score, 
                   scorenames = c('TUMOR4-EPI','TUMOR4-MES'),hexMin = 10)
vld_pdac4_ctc$final_EM$sampleSource= 'CTC'
vld_pdac4_tumor$final_EM$sampleSource= 'Tumor'

plot_dat_vld = rbind(vld_pdac4_ctc$final_EM,vld_pdac4_tumor$final_EM)
plot_dat_vld$type = validate_pdac4_phe$type[match(plot_dat_vld$sampleID,
                                                  validate_pdac4_phe$sampleID)]
table(plot_dat_vld$type,plot_dat_vld$sampleSource)
OutVals = boxplot(plot_dat_vld$EM_score)$out
plot_dat_vld_rmOt = plot_dat_vld[!plot_dat_vld$EM_score %in% OutVals,]
plot_dat_vld_rmOt$EMtanh = tanh(plot_dat_vld_rmOt$EM_score/2)
plot_dat_vld_rmOt$EMtype= ifelse(plot_dat_vld_rmOt$EMtanh<x1,'Mesenchymal',
                            ifelse(plot_dat_vld_rmOt$EMtanh>x2,'Epithelial','pEMT'))
plot_dat_vld_rmOt$EMtype =factor(plot_dat_vld_rmOt$EMtype,levels = c('Epithelial','pEMT','Mesenchymal'))

p3 = ggplot(plot_dat_vld_rmOt,mapping = aes(x= EMtype,y = EMtanh,fill=type))+
  geom_boxplot( position = position_dodge2(preserve = "single") )+
  geom_hline(yintercept = c(x1,x2),lty="dashed")+
  theme_bw()+mytheme
p4 = ggplot(plot_dat_vld_rmOt,mapping = aes(x= sampleID,y = EM_score,fill=type))+
  geom_bar(position = "stack", stat = "identity", width = 0.8)+theme_bw()+mytheme
p5 = ggplot(plot_dat_vld_rmOt,mapping = aes(x= epi_score,y = mes_score,color=type))+
  geom_point()+theme_bw()+mytheme
print(p1+p2+p3)
dev.off()
