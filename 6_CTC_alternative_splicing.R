setwd('/mywork/path')
ps <- c('stringr', 'stringi', 'stats', 'reshape2','ggpubr',
        'dplyr','survival','survminer',
        'DESeq2','clusterProfiler','org.Hs.eg.db')
lapply(ps, function(x){library(x, character.only = T)}) 
#Load sample information and prepare files required by rMATS
EMT_result = read.csv('EMscore_tanh.csv',header = T,row.names = 1)
CTC_info =read.csv('CTC_tumor_purity_estimatResult_final.csv',header = T,row.names = 1)
EMT_result_ctc = subset(EMT_result,type2!='CTCHD'&sampleSource=='CTC')
EMT_result_ctc$SRRID=CTC_info$SRRID[match(EMT_result_ctc$sampleID,rownames(CTC_info))]
EMT_result_ductNorm = subset(EMT_result,dataSource=='GSE57973Bulk')
EMT_result_ctc$bamid = paste0(EMT_result_ctc$SRRID,'Aligned.sortedByCoord.out.bam')
write.table(paste0(EMT_result_ctc$SRRID,'Aligned.sortedByCoord.out.bam',collapse = ','),
            '/fastqc/bam4rmats',quote = F,row.names = F,col.names = F)
write.csv(EMT_result_ctc,'AS_analysis/EMT_result_ctc_Info.csv')

#rMATS was used to obtain alternative splicing information for each sample.
#Here we integrate the AS results (including A3SS, A5SS, MXE, SE, RI) of 
#each group (Epithelial, pEMT and Mesenchymal)
rmats_resultGet= function(dir_path,dir_name,dataset,per_row,bamfile){
  files_list = list.files(paste0(dir_path,dir_name)) 
  as_files = grep('.MATS.JC.txt',files_list,value = T)
  as_path = paste0(dir_path,dir_name,as_files)
  
  as_list = c()
  for (i in 1:5) {
    print(i)
    as_res = read.csv(as_path[i],header = T,row.names = 1,sep = '\t')
    bam4rmats = read.table(paste0(dir_path,bamfile),sep = ',')
    rmats_data = as_res
    rmats_data$event_num = paste0('event',1:nrow(rmats_data))
    rownames(rmats_data) = rmats_data$event_num 
    
    Level1_data = apply(str_split(rmats_data$IncLevel1,'[,]',simplify = T),2,as.numeric)
    Level1_data = as.data.frame(Level1_data)
    rownames(Level1_data)=paste0('event',1:nrow(Level1_data))
    colnames(Level1_data)=bam4rmats[1,]
    
    na_psi_v1 = as.data.frame(is.na(Level1_data))
    print(table(rowSums(na_psi_v1)<4))
    if(dataset == 'CTC'){
      ctc_type = EMT_result_ctc$EMtype[match(bam4rmats[1,],EMT_result_ctc$bamid)]
      ctc_type = as.factor(ctc_type)
      
      print(levels(ctc_type))
      event_merge = c()
      as_type_list = c()
      for (n in 1:length(levels(ctc_type)) ){
        z = levels(ctc_type)[n]
        print(z)
        
        Level1_data_tmp = Level1_data[,which(ctc_type==z)]
        na_psi_v1_tmp = as.data.frame(is.na(Level1_data_tmp))
        
        na_per_row = rowSums(na_psi_v1_tmp)/ncol(na_psi_v1_tmp)
        event_merge = Level1_data_tmp[na_per_row<per_row,]
        
        
        print(paste0(nrow(event_merge),' events were saved'))
        
        event_final = rownames(event_merge)
        event_info = rmats_data[event_final,1:which(colnames(rmats_data)=='IncLevel1')-1]
        
        
        event_info$event_name = apply(event_info[,1:which(colnames(event_info)=="ID.1")-1], 1, 
                                      function(x){
                                        y = paste0(x,collapse = '.')
                                        return(y)})
        
        rmats_data_final = cbind(event_info,event_merge)    
        
        as_type_list[[z]]=rmats_data_final
      }
      as_type = str_split(str_split(as_path[i],dir_name,simplify = T)[2],'[.]',simplify = T)[,1]
      print(as_type)
      as_list[[i]]=as_type_list
      names(as_list)[i]=as_type
      
    }
    else{
      na_per_row = rowSums(na_psi_v1)/ncol(na_psi_v1)
      event_merge = Level1_data[na_per_row<per_row,]
      
      event_final = rownames(event_merge)
      event_info = rmats_data[event_final,1:which(colnames(rmats_data)=='IncLevel1')-1]
      
      event_info$event_name = apply(event_info[,1:which(colnames(event_info)=="ID.1")-1], 1, 
                                    function(x){
                                      y = paste0(x,collapse = '.')
                                      return(y)})
      
      rmats_data_final = cbind(event_info,event_merge)    
      
      as_type = str_split(str_split(as_path[i],dir_name,simplify = T)[2],'[.]',simplify = T)[,1]
      print(as_type)
      
      as_list[[i]]=rmats_data_final
      names(as_list)[i]=as_type
    }
    
  }
  return(as_list)
}

CTC_path= '/data_new/xiaolixing/pdac_revision/ctc_GSE144561/fastqc/'
CTC_rmats_name='rmats_PDAC30_ctc/'
CTC_rmats_resultGet = rmats_resultGet(dir_path=CTC_path,
                                      dir_name=CTC_rmats_name,
                                      dataset='CTC',
                                      bamfile = 'bam4rmats',
                                      per_row=0.85)

DuctNorm_path= '/data_new/xiaolixing/pdac_revision/validateNorm_GSE57973/fastqc/'
DuctNorm_name='rmats_Norm/'
DuctNorm_resultGet = rmats_resultGet(dir_path=DuctNorm_path,
                                     dir_name=DuctNorm_name,
                                     dataset='Duct',
                                     bamfile = 'bam4rmats',
                                     per_row=0.85)

CTC_Duct = list(ctc = CTC_rmats_resultGet,
                duct= DuctNorm_resultGet)

#saveRDS(CTC_Duct,'CTC_Duct_rmats_resultRaw.rds')

CTC_Duct = readRDS('CTC_Duct_rmats_resultRaw.rds')
#Intersect the normal samples with alternating events 
#of different groups to obtain a matrix of each type of AS, 
#which is used for subsequent difference analysis.
event_final_list = c()
for (z in 1:5) {
  CTC_type = CTC_rmats_resultGet[[z]]
  duct_tmp=DuctNorm_resultGet[[z]]
  n_sp= (which(colnames(duct_tmp)=='event_name')+1):ncol(duct_tmp)
  colnames(duct_tmp)[n_sp]=paste0(colnames(duct_tmp)[n_sp],'_DuctNorm')
  iner_list = c() 
  for (i in 1:3) {
    ctc_tmp=CTC_type[[i]]
    
    event_intersect = intersect(ctc_tmp$event_name,duct_tmp$event_name)
    print(paste(length(event_intersect),' event intersect find'))
    
    event_dat= cbind(ctc_tmp[match(event_intersect,ctc_tmp$event_name),],
                     duct_tmp[match(event_intersect,duct_tmp$event_name),n_sp])
    print(paste(dim(event_dat)[1],' event intersect'))
    EMTYPE = names(CTC_type)[i]
    iner_list[[EMTYPE]]=event_dat
  }
  astype = names(CTC_rmats_resultGet)[z]
  print(astype)
  event_final_list[[astype]]=iner_list
}
#saveRDS(event_final_list,'CTC_Duct_resultIntersectRaw.rds')

#Perform differential alternative splicing analysis
diffevent_list =  lapply(event_final_list, function(y){
  as_inner=lapply(y, function(z){
    event_info = z[,1:which(colnames(z)=='event_name')]
    event_duct = z[,(ncol(z)-3):ncol(z)]
    event_ctc = z[,(which(colnames(z)=='event_name')+1):(ncol(z)-4)]
    
    event_ctc_tmp = event_ctc
    
    print(paste0(dim(event_ctc_tmp)[2],' samples included'))
    
    event_ctc_impute = apply(event_ctc_tmp, 1, function(x){
      x[is.na(x)]=mean(x,na.rm = T)
      return(x)})
    event_ctc_impute = as.data.frame(t(event_ctc_impute))
    
    event_duct_impute = apply(event_duct, 1, function(x){
      x[is.na(x)]=mean(x,na.rm = T)
      return(x)})
    event_duct_impute = as.data.frame(t(event_duct_impute))
    
    ctc_mean_tmp = apply(event_ctc_impute, 1, function(x) mean(x,na.rm =T))
    duct_mean_tmp = apply(event_duct_impute, 1, function(x) mean(x,na.rm =T))
    
    
    event_tmp = cbind(event_ctc_impute,event_duct_impute)
    
    event_W.test = apply(event_tmp, 1, function(x){
      
      ctc_tmp = na.omit(as.numeric(x[1:(length(x)-4)]))
      duct_tmp = na.omit(as.numeric(x[(length(x)-3):length(x)]))
      
      if (any(c(length(ctc_tmp),length(duct_tmp))<2) ){
        w_p=NA
      }else{
        w_t = wilcox.test(ctc_tmp,duct_tmp,paired = F)
        w_p = w_t$p.value
      }})
    event_padj = round(p.adjust(event_W.test,method = 'BH'),3)
    
    event_tmp$ctc_mean = ctc_mean_tmp
    event_tmp$duct_mean = duct_mean_tmp
    event_tmp$delt = event_tmp$ctc_mean-event_tmp$duct_mean 
    event_tmp$p.value = event_W.test
    event_tmp$padj = event_padj
    
    event_dif_tmp = subset(event_tmp,padj<0.05 & abs(delt)>0.1)
    
    event_dif_tmp = cbind(event_info[rownames(event_dif_tmp),],
                          event_dif_tmp)
    return(event_dif_tmp)
    
  })
  return(as_inner)
  
})

#saveRDS(diffevent_list,'CTC_Duct_DASEresult.rds')

#Differential alternative splicing analysis results processing
#Calculate the differential alternative splicing landscape in each group
EMT_result = read.csv('EMscore_tanh.csv',header = T,row.names = 1)
CTC_info =read.csv('CTC_tumor_purity_estimatResult_final.csv',header = T,row.names = 1)
EMT_result_ctc = subset(EMT_result,type2!='CTCHD'&sampleSource=='CTC')
EMT_result_ctc$SRRID=CTC_info$SRRID[match(EMT_result_ctc$sampleID,rownames(CTC_info))]
EMT_result_ductNorm = subset(EMT_result,dataSource=='GSE57973Bulk')

CTC_Duct = readRDS('AS_analysis/CTC_Duct__rmats_resultRaw.rds')
event_final_list=readRDS('AS_analysis/CTC_Duct_resultIntersectRaw.rds')
diffevent_list = readRDS('AS_analysis/CTC_Duct_DASEresult.rds')

all_ctcRawevent_summarize =  sapply(CTC_Duct$ctc, function(x){lapply(x,function(y){dim(y)[1]})})
all_ductRawevent_summarize = sapply(CTC_Duct$duct, function(x){dim(x)[1]})
all_Rawevent_summarize = rbind(all_ctcRawevent_summarize,all_ductRawevent_summarize)
rownames(all_Rawevent_summarize)[4]='Duct'

all_Intevent_summarize = as.data.frame(sapply(event_final_list, function(x){lapply(x,function(y){dim(y)[1]})}))
all_Diffevent_summarize=do.call(cbind,lapply(diffevent_list,function(x){lapply(x,function(y){dim(y)[1]})}))

dif_se= diffevent_list$SE
dif_event = intersect(dif_se$Epithelial$event_name,dif_se$Mesenchymal$event_name)
dif_event = intersect(dif_event,dif_se$pEMT$event_name)

epiRaw  = lapply(event_final_list, function(x){x[[1]]})
mesRaw  = lapply(event_final_list, function(x){x[[2]]})
pEMTRaw = lapply(event_final_list, function(x){x[[3]]})
EMT_event_raw=list(epithelial = epiRaw, mesenchymal =mesRaw,pEMT =pEMTRaw)
#saveRDS(EMT_event_raw,'CTC_Duct_EM_RawEresult.rds')

epidiff = lapply(diffevent_list, function(x){x[[1]]})
mesdiff = lapply(diffevent_list, function(x){x[[2]]})
pEMTdiff= lapply(diffevent_list, function(x){x[[3]]})

diffevent_list_v2 = list(epithelial = epidiff,mesenchymal =mesdiff,pEMT =pEMTdiff)
#saveRDS(diffevent_list_v2,'CTC_Duct_EM_DASEresult.rds')

#EMT_event_raw = readRDS('CTC_Duct_EM_RawEresult.rds')
#diffevent_list_v2 = readRDS('CTC_Duct_EM_DASEresult.rds')

n_gene_epiHD = sapply(epidiff,function(x){x[,'geneSymbol']})
n_gene_pEMTHD = sapply(pEMTdiff,function(x){x[,'geneSymbol']})
n_gene_MesHD = sapply(mesdiff,function(x){x[,'geneSymbol']})

all_event_summarize = data.frame(epi_vs_hd  = sapply(epidiff, function(x){dim(x)[1]}),
                                 pEMT_vs_hd = sapply(pEMTdiff,function(x){dim(x)[1]}),
                                 mes_vs_hd  = sapply(mesdiff, function(x){dim(x)[1]}),
                                 row.names  = c('A3SS','A5SS','MXE','RI','SE'))
all_gene_summarize = data.frame(epi_vs_hd =sapply(n_gene_epiHD, function(x){length(unique(x))}),
                                pEMT_vs_hd = sapply(n_gene_pEMTHD, function(x){length(unique(x))}),
                                mes_vs_hd = sapply(n_gene_MesHD, function(x){length(unique(x))}),
                                row.names = c('A3SS','A5SS','MXE','RI','SE'))
gene_summarize = list(epi_vs_hd =sapply(n_gene_epiHD, function(x){unique(x)}),
                            pEMT_vs_hd = sapply(n_gene_pEMTHD, function(x){unique(x)}),
                            mes_vs_hd = sapply(n_gene_MesHD, function(x){unique(x)}))

mut_as_gene = list(epi_vs_hd =sapply(n_gene_epiHD, function(x){unique(x[duplicated(x)])}),
                         pEMT_vs_hd = sapply(n_gene_pEMTHD, function(x){unique(x[duplicated(x)])}),
                         mes_vs_hd = sapply(n_gene_MesHD, function(x){unique(x[duplicated(x)])}))


colSums(all_event_summarize)
bar_data = all_event_summarize
bar_data$group = factor(rownames(bar_data),levels = c('A3SS','A5SS','MXE','RI','SE'))
bar_data = reshape2::melt(bar_data,'group')
epi = subset(bar_data,variable=='epi_vs_hd')
mes = subset(bar_data,variable=='pEMT_vs_hd')
pemt = subset(bar_data,variable=='mes_vs_hd')

pdf('../../Figures/DASE_sumpie.pdf',width = 9,height = 3)
p1 = ggpie(epi,"value",lab.pos = "in",  color = "group",
      palette = c("#E01516", "#228B22", "#FDB462" ,"#8B658B", "#4876FF" ))
p2 = ggpie(pemt,"value", lab.pos = "in", color = "group",
           palette = c("#E01516", "#228B22", "#FDB462" ,"#8B658B", "#4876FF" ))
p3 = ggpie(mes,"value",lab.pos = "in",color = "group",
           palette = c("#E01516", "#228B22", "#FDB462" ,"#8B658B", "#4876FF" ))
print(p1+p2+p3)
dev.off()

#Find genes associated with multiple differential alternative splicing events

all_InterEvent_summarize = data.frame(epi_vs_hd  = sapply(epiRaw, function(x){dim(x)[1]}),
                                 pEMT_vs_hd = sapply(pEMTRaw,function(x){dim(x)[1]}),
                                 mes_vs_hd  = sapply(mesRaw, function(x){dim(x)[1]}),
                                 row.names  = c('A3SS','A5SS','MXE','RI','SE'))
bar_data_interRaw = all_InterEvent_summarize
bar_data_interRaw$group = factor(rownames(bar_data_interRaw),
                                 levels = c('A3SS','A5SS','MXE','RI','SE'))
bar_data_interRaw = reshape2::melt(bar_data_interRaw,'group')
epi_inter = subset(bar_data_interRaw,variable=='epi_vs_hd')
mes_inter = subset(bar_data_interRaw,variable=='pEMT_vs_hd')
pemt_inter= subset(bar_data_interRaw,variable=='mes_vs_hd')

pdf('../../Figures/eventRaw_inter_sumpie.pdf',width = 9,height = 3)
p1 = ggpie(epi_inter,"value", lab.pos = "in", color = "group",
           palette = c("#E01516", "#228B22", "#FDB462" ,"#8B658B", "#4876FF" ))
p2 = ggpie(pemt_inter,"value",lab.pos = "in", color = "group",
           palette = c("#E01516", "#228B22", "#FDB462" ,"#8B658B", "#4876FF" ))
p3 = ggpie(mes_inter,"value", lab.pos = "in", color = "group",
           palette = c("#E01516", "#228B22", "#FDB462" ,"#8B658B", "#4876FF" ))
print(p1+p2+p3)
dev.off()

n_gene_epiHD = sapply(epidiff,function(x){x[,'geneSymbol']})
n_gene_pEMTHD = sapply(pEMTdiff,function(x){x[,'geneSymbol']})
n_gene_MesHD = sapply(mesdiff,function(x){x[,'geneSymbol']})

mut_as_gene = list(epi_vs_hd =sapply(n_gene_epiHD, function(x){length(unique(x[duplicated(x)]))}),
                   pEMT_vs_hd = sapply(n_gene_pEMTHD, function(x){length(unique(x[duplicated(x)]))}),
                   mes_vs_hd = sapply(n_gene_MesHD, function(x){length(unique(x[duplicated(x)]))}))
num_mut_gene =do.call(rbind,mut_as_gene)
all_as_gene = list(epi_vs_hd =sapply(n_gene_epiHD, function(x){length(unique(x))}),
                   pEMT_vs_hd = sapply(n_gene_pEMTHD, function(x){length(unique(x))}),
                   mes_vs_hd = sapply(n_gene_MesHD, function(x){length(unique(x))}))
num_all_gene =do.call(rbind,all_as_gene)
num_sinle_gene=num_all_gene-num_mut_gene

num_sinle_gene = as.data.frame(num_sinle_gene)
num_sinle_gene$gene_type = 'sinle_gene'
num_sinle_gene$group_type = rownames(num_sinle_gene)
num_mut_gene = as.data.frame(num_mut_gene)
num_mut_gene$gene_type = 'multi_gene'
num_mut_gene$group_type = rownames(num_mut_gene)

num_gene = rbind(num_sinle_gene,num_mut_gene)
num_gene_melt= reshape2::melt(num_gene, id.vars =c('group_type','gene_type'))
num_gene_melt$group_type = factor(num_gene_melt$group_type,
                                   levels=c("epi_vs_hd", "pEMT_vs_hd","mes_vs_hd"))
num_gene_melt$gene_type = factor(num_gene_melt$gene_type,
                                  levels=c("sinle_gene","multi_gene"))

pdf('../../Figures/DASE_gene_type.pdf',width = 6,height = 2.5)
p1= ggplot(num_gene_melt,aes(x = variable,y=value,fill =gene_type))+
  geom_bar(position = "fill",stat = 'identity',width = 0.5,colour='black')+
  facet_wrap(~group_type)+
  theme(axis.text.x = element_text( hjust = 1,angle = 45 ))+
  scale_fill_manual(values = c('sinle_gene'="#dde5b6", multi_gene= "#81c3d7" ))
print(p1)
dev.off()

#Statistical distribution density of delta PSI
temp_as_epi = lapply(epidiff,function(x){x[,c('delt','event_name')]})
temp_as_epi = do.call(rbind,temp_as_epi)
temp_as_epi$type = 'Epithelial'
temp_as_pemt = lapply(pEMTdiff,function(x){x[,c('delt','event_name')]})
temp_as_pemt = do.call(rbind,temp_as_pemt)
temp_as_pemt$type = 'pEMT'
temp_as_mes = lapply(mesdiff,function(x){x[,c('delt','event_name')]})
temp_as_mes = do.call(rbind,temp_as_mes)
temp_as_mes$type = 'Mesenchymal'

data_delt = bind_rows(temp_as_epi,temp_as_pemt,temp_as_mes)
data_delt$as_type = str_split(rownames(data_delt),'[.]',simplify = T)[,1]
head(data_delt)
pdf('../../Figures/DASE_delta_histogram.pdf',width = 6,height = 2)
p1 = ggplot(data = data_delt,aes(x=delt))+
  geom_histogram(aes(fill=type), bins = 50, alpha=0.8)+
  facet_grid(~as_type)+theme_bw()+theme(axis.title = element_blank(),
        axis.text.x = element_text( hjust = 1,angle = 90 ))
print(p1)
dev.off()

#The intersection of differential alternative splicing events 
#between different groups is taken for subsequent comparisons between groups.
Diffevent_info_list = c()
for (z in 1:5) {
  as_type = names(diffevent_list)[z]
  print(as_type)
  tmp_type = diffevent_list[[z]]
  intersect_event = intersect(tmp_type$Epithelial$event_name,
                              tmp_type$Mesenchymal$event_name)
  intersect_event = intersect(intersect_event,
                              tmp_type$pEMT$event_name)
  tmp_info = tmp_type$Epithelial[,1:which(colnames(tmp_type$Epithelial)=='event_name')]
  tmp_info = tmp_info[which(tmp_info$event_name %in%intersect_event ),]
  
  Diffevent_info_list[[as_type]]=tmp_info
  }
saveRDS(Diffevent_info_list,'AS_analysis/Diffevent_info_list.rds')

Diffevent_info_list = readRDS('AS_analysis/Diffevent_info_list.rds')
bar_plot = data.frame(event_type = names(Diffevent_info_list),
                      event_num = sapply(Diffevent_info_list,function(x){dim(x)[1]}))
bar_plot = bar_plot[order(bar_plot$event_num),]
bar_plot$event_type = factor(bar_plot$event_type,levels = bar_plot$event_type)
pdf('../../Figures/DASE_intersectEvent_barplot.pdf',width = 6,height = 4)
p1 = ggplot(bar_plot,aes(x = event_num,y =event_type ,fill=event_num))+
  geom_bar(stat = 'identity',width = 0.6,position = 'dodge')+
  geom_text(aes(label=event_num),size=4,vjust=0)+
  scale_fill_gradientn(colours = c("#C7E9B4", "#1D91C0","#225EA8"))+
  theme_classic()+theme(legend.position = 'top')
print(p1)
dev.off()

#Count the number of genes
total_gene_name = sapply(Diffevent_info_list,function(x){unique(x$geneSymbol)})
library(venn)
names(total_gene_name)
pdf('../../Figures/total_gene_venn.pdf',width = 6,height = 6)
p1= venn(total_gene_name,zcolor=2:7,box=F)
print(p1)
dev.off()

#Count the number of alternative splices
pdf('../../Figures/DASE_intersectGene_barplot.pdf',width = 6,height = 4)
p1 = ggplot(bar_plot,aes(x = event_num,y =event_type ,fill=event_num))+
  geom_bar(stat = 'identity',width = 0.6,position = 'dodge')+
  geom_text(aes(label=event_num),size=4,vjust=0)+
  scale_fill_gradientn(colours = c("#C7E9B4", "#1D91C0","#225EA8"))+
  theme_classic()+theme(legend.position = 'top')
print(p1)
dev.off()

event_func = function(x){
  inter_res = c()
  for (i in 1:5) {
    event_type = names(x)[i]
    dif_info  = Diffevent_info_list[[i]]
    inter_x = x[[i]]
    inter_x = inter_x[match(dif_info$event_name,inter_x$event_name),]
    inter_res[[event_type]]=inter_x
  }
  return(inter_res)
}
epidiff_inter =event_func(epidiff)
pEMTdiff_inter=event_func(pEMTdiff)
mesdiff_inter =event_func(mesdiff)

event_inter = list(epidiff_inter=epidiff_inter,
                   pEMTdiff_inter = pEMTdiff_inter,
                   mesdiff_inter =mesdiff_inter)
saveRDS(event_inter,'AS_analysis/Diffevent_infoEvent_list.rds')

#Obtain genes with significantly different top200 events for GO enrichment analysis
event_inter = readRDS('AS_analysis/Diffevent_infoEvent_list.rds')
top_dif_event = lapply(event_inter,function(x){
  x_ase = lapply(x,function(y){
    y_ase = y[,c("event_name","GeneID","delt","padj")]
    return(y_ase)
  })
  x_ase= do.call(rbind,x_ase)
  return(x_ase)
})
top_dif_event = do.call(rbind,top_dif_event)
top_dif_event =top_dif_event[order(top_dif_event$padj),]
top_dif_event = top_dif_event[1:300,]
top_dif_event = top_dif_event[order(abs(top_dif_event$delt),
                                            decreasing = T),]
top_dif_event =top_dif_event[1:200,]
top_dif_eventGene =unique(top_dif_event$GeneID)#93

write.table(top_dif_eventGene,row.names = F,col.names = F,quote = F)
#Enrichment analysis using Metascape website

#Analysis of differentially expressed genes between each group (epithelial cells, pEMT and mesenchymal cells) 
#and normal samples (normal duct samples)
estimat_score_final = read.csv('CTC_tumor_purity_estimatResult_final.csv',header = T,row.names = 1)
emt_signature = read.csv('emt_gene_signature.csv',header = T,row.names = 1)
count_expr = read.table('CTC_exprData_raw.txt',sep = '\t',header = T,row.names = 1)
count_expr[1:4,1:4]

CTC_reGet = count_expr[,which(colnames(count_expr)%in% EMT_result_ctc$sampleID)]
colnames(CTC_reGet)=EMT_result_ctc$SRRID[match(colnames(CTC_reGet),EMT_result_ctc$sampleID)]
colnames(CTC_reGet)=paste0(colnames(CTC_reGet),"_",EMT_result_ctc$EMtype[match(colnames(CTC_reGet),EMT_result_ctc$SRRID)])

epi_count = CTC_reGet[,grep('Epithelial',colnames(CTC_reGet))]
pEMT_count = CTC_reGet[,grep('pEMT',colnames(CTC_reGet))]
Mes_count = CTC_reGet[,grep('Mesenchymal',colnames(CTC_reGet))]

duct_Expr = read.table('GSE57973_duct_cell_featureCount.txt',sep = '\t',header = T,row.names = 1)
Norm_DuctGet = duct_Expr
colnames(Norm_DuctGet) = paste0(str_split(colnames(Norm_DuctGet) ,'[A]',simplify = T)[,1],'_Norm')

gene_ins = intersect(rownames(CTC_reGet),rownames(Norm_DuctGet))
ids=annoGene(gene_ins,'SYMBOL','human')
gene_ins_proteinCoding = unique(ids$SYMBOL[which(ids$biotypes=='protein_coding')])

Mes_count_procoding = Mes_count[na.omit(match(gene_ins_proteinCoding,rownames(Mes_count))),]
epi_count_procoding = epi_count[na.omit(match(gene_ins_proteinCoding,rownames(epi_count))),]
pEMT_count_procoding = pEMT_count[na.omit(match(gene_ins_proteinCoding,rownames(pEMT_count))),]
Norm_count_procoding = Norm_DuctGet[na.omit(match(gene_ins_proteinCoding,rownames(Norm_DuctGet))),]

diff_geneExp = function(x,y,type){
  combat_Expr= cbind(x,y)
  phe_meta = data.frame(row.names = colnames(combat_Expr),
                        sampleType= str_split(colnames(combat_Expr),'[_]',simplify = T)[,2],
                        batch = c(rep('ctc',ncol(x)),
                                  rep('duct',ncol(y))))
  combat_Expr_tpm = as.matrix(combat_Expr)
  countData<- ComBat_seq(combat_Expr_tpm,batch = phe_meta$batch)
  
  condition <- factor(phe_meta$sampleType)
  colData <- data.frame(row.names=colnames(countData), condition)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
  dds1 <- DESeq(dds) 
  res_deg <- results(dds1,contrast = c('condition',type,'Norm'))
  return(res_deg)
}

epi_deg = diff_geneExp(x=epi_count_procoding,y=Norm_count_procoding,type='Epithelial')
pEMT_deg= diff_geneExp(x=pEMT_count_procoding,y=Norm_count_procoding,type='pEMT')
Mes_deg = diff_geneExp(x=Mes_count_procoding,y=Norm_count_procoding,type='Mesenchymal')

CTC_All_deg =list(Epithelial = epi_deg,pEMT=pEMT_deg,Mesenchymal=Mes_deg)

FindSig = function(resdat,foldChange_cutoff=1,padj_cutoff=0.05){
  res_raw <- data.frame(resdat, stringsAsFactors = FALSE, check.names = FALSE)
  resSig<- res_raw[which(abs(res_raw$log2FoldChange) > foldChange_cutoff & res_raw$padj < padj_cutoff),] 
  resSig$gene = rownames(resSig)
  resSig$change = ifelse(resSig$log2FoldChange>foldChange_cutoff,'Up-regulation','Down-regulation')
  #resSig$type = emt_signature$Category[match(resSig$gene,emt_signature$Gene)]
  return(resSig)
}
CTC_All_degSig =lapply(CTC_All_deg,FindSig)

pdf('../../Figures/DASE_DEG_venn.pdf',width = 5,height = 5)
p1 = ggvenn::ggvenn(CTC_All_degSig_gene,show_percentage = F,
                    fill_color = c("#FDB462","#a3cd5b","#8ea0cc", "red"))
print(p1)
dev.off()

#saveRDS(CTC_All_degSig,'AS_analysis/DEG_ctc_vs_duct.rds')

#Remove events from differentially expressed genes
event_inter = readRDS('AS_analysis/Diffevent_infoEvent_list.rds')

event_out_final=c()
for (i in 1:3) {
  em_type = names(CTC_All_degSig_gene)[i]
  print(em_type)
  deg = CTC_All_degSig_gene[[i]]
  dase = event_inter[[i]]
  event_out = lapply(dase, function(x){
    event_out_num = which(x$geneSymbol %in% deg)
    print(length(event_out_num))
    if (length(event_out_num)==0) {
      event_out=x
    }else{
      event_out=x[-event_out_num,]
    }
    return(event_out)
  })
  event_out_final[[em_type]]=event_out
}

#saveRDS(event_out_final,'AS_analysis/CTC_Duct_DASE_degOut.rds')
#event_out_final = readRDS('AS_analysis/CTC_Duct_DASE_degOut.rds')
epi_out_event = event_out_final$Epithelial
pemt_out_event = event_out_final$pEMT
mes_out_event = event_out_final$Mesenchymal
event_out_deg_final = c()
event_out_deg_info = c()
event_out_deg_info_gene = c()
for (i in 1:5) {
  as_type = names(epi_out_event)[i]
  print(as_type)
  e = epi_out_event[[i]]
  p = pemt_out_event[[i]]
  m = mes_out_event[[i]]
  
  event_out_inter =intersect(e$event_name,p$event_name)
  event_out_inter =intersect(event_out_inter,m$event_name)
  
  e_in = e[which(e$event_name %in% event_out_inter),]
  p_in = p[which(p$event_name %in% event_out_inter),]
  m_in = m[which(m$event_name %in% event_out_inter),]
  
  as_in = list(Epithelial =e_in,
               pEMT = p_in,
               Mesenchymal = m_in)
  event_out_deg_final[[as_type]] =  as_in
  
  event_out_deg_info[[as_type]] = event_out_inter
  event_out_deg_info_gene[[as_type]] =unique(str_split(event_out_inter,'[.]',simplify = T)[,1]) 
  
}

epi_event_outDeg =lapply(event_out_deg_final,function(x){x[[1]]})
pEMT_event_outDeg =lapply(event_out_deg_final,function(x){x[[2]]})
Mes_event_outDeg =lapply(event_out_deg_final,function(x){x[[3]]})
event_EM_final = list(Epithelial = epi_event_outDeg,
                      pEMT = pEMT_event_outDeg,
                      Mesenchymal =Mes_event_outDeg)
#saveRDS(event_EM_final,'AS_analysis/DASE_degOut_EMFinal.rds')
#saveRDS(event_out_deg_final,'AS_analysis/CTC_Duct_DASE_degOut_interFinal.rds')
#saveRDS(event_out_deg_info_gene,'AS_analysis/CTC_Duct_DASE_degOut_interInfo_gene.rds')

event_out_deg_final = readRDS('AS_analysis/CTC_Duct_DASE_degOut_interFinal.rds')
event_out_deg_info_gene =  readRDS('AS_analysis/CTC_Duct_DASE_degOut_interInfo_gene.rds')
bar_plot = data.frame(event_type = names(event_out_deg_info),
                      event_num = sapply(event_out_deg_info,length))
bar_plot = bar_plot[order(bar_plot$event_num),]
bar_plot$event_type = factor(bar_plot$event_type,levels = bar_plot$event_type)
pdf('../../Figures/DASEoutDEG_intersectEvent_barplot.pdf',width = 6,height = 4)
p1 = ggplot(bar_plot,aes(x = event_num,y =event_type ,fill=event_num))+
  geom_bar(stat = 'identity',width = 0.6,position = 'dodge')+
  geom_text(aes(label=event_num),size=4,vjust=0)+
  scale_fill_gradientn(colours = c("#C7E9B4", "#1D91C0","#225EA8"))+
  theme_classic()+theme(legend.position = 'top')
print(p1)
dev.off()

#Analyze the distribution of differential splicing exons on chromosomes
event_out_deg_final = readRDS('AS_analysis/CTC_Duct_DASE_degOut_interFinal.rds')
event_out_deg_info_gene =  readRDS('AS_analysis/CTC_Duct_DASE_degOut_interInfo_gene.rds')

chr_size = read.table('AS_analysis//hg38.chrom.sizes.txt',sep = '\t',header = F,row.names = 1)
chr_size$chr = rownames(chr_size)
colnames(chr_size)[1]='Total_length_bp'

event_chr_info=c()
for (i in 1:5) {
  as_type = names(event_out_deg_final)[i]
  as_tmp = event_out_deg_final[[i]]
  chr_info_tmp = lapply(as_tmp,function(x){
    x_tmp = x[,c("event_name","chr","delt")]
    return(x_tmp)
  })
  chr_info_tmp= do.call(rbind,chr_info_tmp)
  chr_info_tmp$emtype = str_split(rownames(chr_info_tmp),'[.]',simplify = T)[,1]
  chr_info_tmp$astype = as_type
  event_chr_info = rbind(event_chr_info,chr_info_tmp)
  
}
saveRDS(event_chr_info,'AS_analysis/event_chr_info.rds')

#Analyzing the correlation between chromosome length and 
#the number of alternative splices

dat_length_cor =subset(event_chr_info,emtype=='Epithelial')
dat_length_cor$chr_length = dfcol$length[match(dat_length_cor$chr,dfcol$x)]
dat_cor = as.data.frame(table(dat_length_cor$chr,dat_length_cor$emtype))
dat_length_cor$event_total_count = dat_cor$Freq[match(dat_length_cor$chr,dat_cor$Var1)]

dat_coras = as.data.frame(table(dat_length_cor$chr,dat_length_cor$astype))
dat_coras$length = log10(dfcol$length[match(dat_coras$Var1,dfcol$x)])
dat_coras$event_total_count = dat_length_cor$event_total_count[match(dat_coras$Var1,dat_length_cor$chr)]
dat_coras = dat_coras[order(dat_coras$length),]
dat_coras$Var1 = as.factor(dat_coras$Var1)

pdf('../../Figures/DASE_Chr_correlation.pdf',width = 8,height = 3)
g <- ggplot(dat_coras, aes(length, Freq))
p1= g + geom_jitter(aes(col=Var2, size=Freq)) + 
  scale_size(range = c(1,3))+
  facet_wrap('Var2',scales = 'free_y')+
  geom_smooth(aes(col=Var2), method="lm", se=T)+
  stat_cor(method = 'pearson',alternative = 'greater')+
  scale_color_manual(values  = c("#E01516", "#228B22", "#FDB462" ,
                                 "#8B658B", "#4876FF"))+
  theme_bw()
print(p1)
dev.off()

#Analysis of correlations between 
#differential alternative splicing events and EMT stage

event_out_deg_final = readRDS('AS_analysis/CTC_Duct_DASE_degOut_interFinal.rds')
event_out_deg_info_gene =  readRDS('AS_analysis/CTC_Duct_DASE_degOut_interInfo_gene.rds')

event_out_comb = c()
for (i in 1:5) {
  astype = names(event_out_deg_final)[i]
  dat = event_out_deg_final[[i]]
  #get all info from each AS type
  dat_info = lapply(dat,function(x){
    x = x[,c("GeneID","geneSymbol","event_name","ctc_mean","delt","padj")]
    return(x)
  })

  dat_info_bind = do.call(rbind,dat_info)
  dat_info_bind$EM = str_split(rownames(dat_info_bind),'[.]',simplify = T)[,1]
  dat_info_bind$as_type = astype
  event_out_comb[[astype]]=dat_info_bind
}

event_out_comb_final = do.call(rbind,event_out_comb)
event_out_tmp_delt = dcast(event_out_comb_final,event_name+as_type~EM,value.var = 'delt')
rownames(event_out_tmp_delt) = paste(event_out_tmp_delt$as_type,'_',event_out_tmp_delt$event_name)
event_out_tmp_delt = event_out_tmp_delt[,c('Epithelial','Mesenchymal','pEMT')]
event_out_tmp_delt$sd = apply(event_out_tmp_delt,1,sd)
sd_event = rownames(event_out_tmp_delt)
event_out_tmp_T = dcast(event_out_comb_final,event_name+as_type~EM,value.var = 'ctc_mean')
rownames(event_out_tmp_T) = paste(event_out_tmp_T$as_type,'_',event_out_tmp_T$event_name)
event_out_tmp_T =event_out_tmp_T[which(rownames(event_out_tmp_T) %in%sd_event ),]
event_out_tmp_T$gene_ID = str_split(event_out_tmp_T$event_name,'[.]',simplify = T)[,1]
event_out_tmp_T$gene_name =  str_split(event_out_tmp_T$event_name,'[.]',simplify = T)[,3]


dupgene = na.omit(event_out_tmp_T$gene_name)[duplicated(na.omit(event_out_tmp_T$gene_name))]
event_out_tmp_T$dupGene =ifelse(event_out_tmp_T$gene_name %in%dupgene,'mutiEvent',
                                'singleEvent')
event_out_tmp_T$genelab = event_out_tmp_T$gene_name
event_out_tmp_T$genelab[which(event_out_tmp_T$dupGene!='mutiEvent')]=NA

pdf('../../Figures/DASE_eventTernary.pdf',width = 5,height = 5)
p1=ggtern(event_out_tmp_T, aes(x = Mesenchymal, y = pEMT,z = Epithelial, color= as_type))+ 
  geom_mask()+geom_point(aes(size = ifelse(dupGene=='mutiEvent',1,0.5)),alpha=0.8, shape = 16)+
  scale_size(range = c(2,5))
print(p1)
dev.off()


#Analyze events from cell surface protein genes
event_out_deg_final = readRDS('AS_analysis/CTC_Duct_DASE_degOut_interFinal.rds')
event_out_deg_info_gene =  readRDS('AS_analysis/CTC_Duct_DASE_degOut_interInfo_gene.rds')
GESP = read.csv('AS_analysis//hunman_cell_surface_pro.csv',header = T)

sur_all = sapply(event_out_deg_info_gene, function(z){
  e_id = z
  s_gene = intersect(GESP$ENSEMBL.Gene.ID,e_id)
  m_data = length(s_gene)
  return(m_data)
})
sur_all = as.data.frame(sur_all)
sur_all$type = rownames(sur_all)
sur_all = sur_all[order(sur_all$sur_all,decreasing = T),]
sur_all$type = factor(sur_all$type,levels = sur_all$type)
sum(sur_all$sur_all)#77 cell-surface gene
pdf('../../Figures/DASE_surfaceAllGene_barplot.pdf',width = 2,height = 2.5)
p1 = ggplot(sur_all,aes(y =sur_all,x =  type,fill=type))+
  geom_bar(stat = 'identity',width = 0.6)+
  theme_bw()
print(p1)
dev.off()


gene_info_epiHD = lapply(epi_event_outDeg,function(x){x[,c('event_name',"GeneID",'geneSymbol','ctc_mean','delt')]})
gene_info_pEMTHD = lapply(pEMT_event_outDeg,function(x){x[,c('event_name',"GeneID",'geneSymbol','ctc_mean','delt')]})
gene_info_MesHD = lapply(Mes_event_outDeg,function(x){x[,c('event_name',"GeneID",'geneSymbol','ctc_mean','delt')]})

surface_pro = function(x,type){
  y  = lapply(x, function(z){
    e_id = str_split(z$GeneID,'[.]',simplify = T)[,1]
    s_gene = intersect(GESP$ENSEMBL.Gene.ID,e_id)
    m_data = z[which(e_id %in% s_gene),]
    
    m_data$gene_sum = length(s_gene)
    m_data$event_sum= nrow(m_data)
    return(m_data)
  })
  res = do.call(rbind ,y)
  res$as_type = str_split(rownames(res),'[.]',simplify = T)[,1]
  res$deltaType= ifelse(res$delt>0,'dPSI>0','dPSI<0')
  res$dateType=type
  return(res)
}

epi_surface = surface_pro(gene_info_epiHD,'Epithelial')
pEMT_surface = surface_pro(gene_info_pEMTHD,'pEMT')
mes_surface = surface_pro(gene_info_MesHD,'Mesenchymal')
surface_raw= list(epithelial =epi_surface, pEMT =pEMT_surface,Mesenchymal =mes_surface )
#saveRDS(surface_raw,'AS_analysis/surface_rawIntersect_info.rds')

bar_func = function(x,type){
  bar_data=as.data.frame(table(x$deltaType,x$as_type))
  bar_data$type = type
  return(bar_data)
}
epi_bar = bar_func(epi_surface,'Epithelial')
pEMT_bar = bar_func(pEMT_surface,'pEMT')
mes_bar = bar_func(mes_surface,'Mesenchymal')
bar_data = bind_rows(epi_bar,pEMT_bar,mes_bar)
pdf('../../Figures/DASE_surfacePro_barplot.pdf',width = 8,height = 4)
p1 =ggplot(bar_data,aes(x =Var1 ,y = Freq,fill=type))+
  geom_bar(stat = 'identity',width = 0.6,position = 'dodge')+
  facet_wrap(~Var2, scales = "free_y")+
  ggtitle('Neojunction events in cell-surface protein')+
  ylab('Number of events')+
  xlab("")
print(p1)
dev.off()

box_dat = subset(event_out_tmp_T,dupGene=='mutiEvent')
rownames(box_dat)=box_dat$event_name
box_dat$gene_id = str_split(rownames(box_dat),'[.]',simplify = T)[,1]

box_dat=box_dat[,c('Epithelial','pEMT','Mesenchymal','gene_name','as_type','gene_id')]
box_dat_melt = melt(box_dat,id.vars = c('gene_name','as_type','gene_id'))

colnames(box_dat_melt)=c('Gene','AS_type','GeneID','Group','PSI')

inter_gene = unique(box_dat_melt$GeneID)
inter_gene_id = GESP$ENSEMBL.Gene.ID[which(GESP$ENSEMBL.Gene.ID %in% inter_gene)]

surGeneEvent_tmp = c()
gene_surGet=c()
for (i in 1:5) {
  astype = names(event_out_deg_final)[i]
  dat = event_out_deg_final[[i]]
  
  x = dat[[1]]
  emtype= names(dat)[1]
  rownames(x)= x$event_name
  x_gene_id = str_split(x$GeneID,'[.]',simplify = T)[,1]
  x_data = as.data.frame(x[which(x_gene_id %in% inter_gene_id),
                           c('ctc_mean','event_name')])
  x_em = x_data
  colnames(x_em)=paste(colnames(x_em),emtype,sep = '_')
  dat_em = x_em
  for (i in 2:3) {
    x = dat[[i]]
    emtype= names(dat)[i]
    rownames(x)= x$event_name
    x_gene_id = str_split(x$GeneID,'[.]',simplify = T)[,1]
    x_data = as.data.frame(x[which(x_gene_id %in% inter_gene_id),
                             c('ctc_mean','event_name')])
    x_em = x_data
    colnames(x_em)=paste(colnames(x_em),emtype,sep = '_')
    x_Norm=  x[which(x_gene_id %in% inter_gene_id),c('duct_mean','event_name')]
    colnames(x_Norm)=paste(colnames(x_Norm),'DuctNorm',sep = '_')
    dat_em=cbind(dat_em,x_em) 
  }
  dat_em  = cbind(dat_em,x_Norm)
  dat_em = dat_em[,-grep('event_nam',colnames(dat_em))]
  rownames(dat_em) = paste(astype,rownames(dat_em),sep = '_')
  
  dat_em$eventName = rownames(dat_em) 
  dat_em = melt(dat_em,variable.name =c('eventName'))
  colnames(dat_em)=c("eventName", "SampleName", "PSI")
  surGeneEvent_tmp[[astype]]=dat_em
  
  gene_surGet[[astype]] =ifelse(is.na(dat_em$eventName),
       NULL,str_split(dat_em$eventName,'[.]',simplify = T)[,3]) 
}
surface_aimedGene = list(surGenePSI=surGeneEvent_tmp,
                         surGeneName=gene_surGet)

surGeneEvent =do.call(rbind,surGeneEvent_tmp)
surGeneEvent$astype = str_split(surGeneEvent$eventName,'[_]',simplify = T)[,1]
surGeneEvent$gene = str_split(surGeneEvent$eventName,'[.]',simplify = T)[,3]
surGeneEvent$group = str_split(surGeneEvent$SampleName,'[_]',simplify = T)[,3]
surGeneEvent$sample = str_split(surGeneEvent$SampleName,'[_]',simplify = T)[,1]
surGeneEvent$group  = factor(surGeneEvent$group ,
                             levels = c('Epithelial','pEMT','Mesenchymal','DuctNorm'))

pdf('../../Figures/surPro_mutiEvent.pdf',width = 5,height = 4)
p1 = ggplot(surGeneEvent,aes(y=PSI,x=group))+
  geom_point(aes(color =eventName))+
  scale_color_manual(values  = c(mycol))+
  geom_line(aes(group=eventName,color =eventName ), linetype="dashed")+
  facet_wrap('gene',nrow = 3)+theme_bw()
print(p1)
dev.off()
#saveRDS(surface_aimedGene,'AS_analysis/surface_aimedGeneAll_info.rds')
#saveRDS(surGeneEvent,'AS_analysis/surface_aimedGene_info.rds')

#Surface DTU analysis
library(IsoformSwitchAnalyzeR)
library(patchwork)
know2ens = read.table('AS_analysis/genecode.v38.cds.information.txt',sep = '\t',header = T)
know2ens$transcript_id = str_split(know2ens$transcript_id ,'[.]',simplify = T)[,1]
EMT_result = read.csv('../../Tables/EMscore_tanh2.csv',header = T,row.names = 1)
CTC_info =read.csv('CTC_tumor_purity_estimatResult_final2.csv',header = T,row.names = 1)
EMT_result_ctc = subset(EMT_result,type2!='CTCHD'&sampleSource=='CTC')
EMT_result_ctc$SRRID=CTC_info$SRRID[match(EMT_result_ctc$sampleID,rownames(CTC_info))]

all_salmon = importIsoformExpression(
  parentDir ='AS_analysis/salmon_all/',
  addIsofomIdAsColumn = TRUE)

sampleID = colnames(all_salmon$abundance)[2:ncol(all_salmon$abundance)]
sampleID = str_split(sampleID,'[_]',simplify = T)[,1]
sampleID[31:34]=paste0(sampleID[31:34],'_DuctNorm')
ctcidtype = EMT_result_ctc$EMtype[match(sampleID[1:30],EMT_result_ctc$SRRID)]
sampleID[1:30] = paste(sampleID[1:30] ,ctcidtype,sep='_')

myDesign <- data.frame(
  sampleID = colnames(all_salmon$abundance)[2:ncol(all_salmon$abundance)],
  condition = str_split(sampleID,'[_]',simplify = T)[,2]
  
)

mySwitchList <- importRdata(
  isoformCountMatrix   = all_salmon$counts,
  isoformRepExpression = all_salmon$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = 'AS_analysis/gencode.v44.annotation.gtf.gz',
  isoformNtFasta       = 'AS_analysis/gencode.v44.transcripts.fa.gz',
  fixStringTieAnnotationProblem = TRUE)
#saveRDS(mySwitchList,'AS_analysis/isoformSwitchResult/mySwitchList_importRdata.rds')
mySwitchList = readRDS('AS_analysis/isoformSwitchResult/mySwitchList_importRdata.rds')
mySwitchListFiltered <- preFilter(switchAnalyzeRlist = mySwitchList)
"orfAnalysis" %in% names( mySwitchListAnalyzed )
mySwitchListAnalyzed <- isoformSwitchTestSatuRn( mySwitchListFiltered, 
                                                 reduceToSwitchingGenes=T )

#Analyze cell-surface protein_ORF
#Here we first extract ORF information based on cell surface gene IDs, and then retain ORFs 
#that intersect with exon positions based on the positions of relevant alternatively spliced exons.
surface_rawEvent = readRDS('AS_analysis/surface_rawIntersect_info.rds')
all(surface_rawEvent$epithelial$event_name==surface_rawEvent$pEMT$event_name)
surface_event = surface_rawEvent$epithelial
surface_event$exonStart0base = str_split(surface_event$event_name,'[.]',simplify = T)[,6]
surface_event$exonEnd = str_split(surface_event$event_name,'[.]',simplify = T)[,7]#89
surface_event$strand = str_split(surface_event$event_name,'[.]',simplify = T)[,5]#89

isoformFeatures_dat = mySwitchListAnalyzed$isoformFeatures
orfAnalysis_dat = mySwitchListAnalyzed$orfAnalysis
orfAnalysis_dat$strand = ifelse(orfAnalysis_dat$stopDistanceToLastJunction<0,'-','+')

cell_surface_isodat = c()
for (i in 1:nrow(surface_event)) {
  print(i)
  exonStart0base = as.numeric(surface_event$exonStart0base[i])
  exonEnd = as.numeric(surface_event$exonEnd[i])
  exonStrand= surface_event$strand[i]
  cell_surface_geneID= surface_event$geneSymbol[i]
  cell_surface_geneID2= surface_event$GeneID[i]
  cell_surface_isoID = isoformFeatures_dat$isoform_id[grep(cell_surface_geneID,
                                                           isoformFeatures_dat$gene_name)]
  if(length(cell_surface_isoID)==0){
    next
  }
  
  cell_surface_isoFinal = c()
  orfStart_tmp= c()
  orfEnd_tmp= c()
  orfStrand_tmp =c()
  for (x in unique(cell_surface_isoID)) {
    print(x)
    orfStart =orfAnalysis_dat$orfStartGenomic[grep(x,orfAnalysis_dat$isoform_id)]
    orfEnd   =orfAnalysis_dat$orfEndGenomic[grep(x,orfAnalysis_dat$isoform_id)]
    orfStrand = orfAnalysis_dat$strand[grep(x,orfAnalysis_dat$isoform_id)]
    
    if(is.na(orfStart)|is.na(orfEnd)){
      next
    }else{
      if(orfStrand>0){
        orf_small=orfStart
        orf_big  =orfEnd
      }else{
        orf_small=orfEnd
        orf_big  =orfStart
      }
    }
    
    if(is.na(orf_big)|is.na(orf_small)){
      next
    }else{
      if(orf_small> exonEnd|exonStart0base>orf_big){
        next
      }else{
        cell_surface_isoFinal = x
        orfStart_tmp =c(orfStart_tmp,orfStart)
        orfEnd_tmp   =c(orfEnd_tmp,orfEnd)
        orfStrand_tmp = c(orfStrand_tmp,orfStrand)
        print(paste0('final:',cell_surface_isoFinal))
      }
    }
  }
  if(is.null(cell_surface_isoFinal)){
    next
  }
  cell_surface_isodat_tmp = data.frame(genename = cell_surface_geneID,
                                       geneID = cell_surface_geneID2,
                                       isoID =cell_surface_isoFinal ,
                                       exonStart0base = exonStart0base,
                                       exonEnd = exonEnd,
                                       exonStrand= exonStrand,
                                       orfStart = orfStart_tmp,
                                       orfEnd =orfEnd_tmp,
                                       orfStrand =orfStrand_tmp)
  cell_surface_isodat = rbind(cell_surface_isodat,cell_surface_isodat_tmp)
}

cell_surface_isodat$transcript_name= know2ens$transcript_name[match(cell_surface_isodat$isoID,know2ens$transcript_id)]

write.csv(cell_surface_isodat,'AS_analysis/cell_surface_isoORFdat.csv')
write.csv(cell_surface_isodat,'../../Tables//cell_surface_isoORFdat.csv')

mySwitchListAnalyzed_cellSurface = subsetSwitchAnalyzeRlist(mySwitchListAnalyzed,
                                                            mySwitchListAnalyzed$isoformFeatures$isoform_id %in% unique(cell_surface_isodat$isoID) )
#saveRDS(mySwitchListAnalyzed_cellSurface,'AS_analysis/isoformSwitchResult/mySwitchList_cellSurface_exonIn.rds')


mySwitchListAnalyzed_cellSurface_ORF <- extractSequence(
  mySwitchListAnalyzed_cellSurface, 
  pathToOutput = 'AS_analysis/isoformSwitchResult/my_surface_orf/',
  removeShortAAseq =F,removeORFwithStop=F,
  outputPrefix='cellSurface_isoform',
  writeToFile=T 
)

mySwitchListAnalyzed_cellSurface_ORF = mySwitchListAnalyzed_cellSurface
### Add CPC2 analysis
mySwitchListAnalyzed_cellSurface_ORF <- analyzeCPC2(
  switchAnalyzeRlist   = mySwitchListAnalyzed_cellSurface_ORF,
  pathToCPC2resultFile = "AS_analysis/isoformSwitchResult/my_surface_orf/result_cpc2_mycellSurface.txt",
  removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)
### Add PFAM analysis
mySwitchListAnalyzed_cellSurface_ORF <- analyzePFAM(
  switchAnalyzeRlist   = mySwitchListAnalyzed_cellSurface_ORF,
  pathToPFAMresultFile = "AS_analysis/isoformSwitchResult/my_surface_orf/result_pFam_mycellSurface.txt",
  showProgress=T
)
### Add DeepTMHMM analysis
mySwitchListAnalyzed_cellSurface_ORF <- analyzeDeepTMHMM(
  switchAnalyzeRlist   = mySwitchListAnalyzed_cellSurface_ORF,
  pathToDeepTMHMMresultFile = "AS_analysis/isoformSwitchResult/my_surface_orf/TMRs.gff3",
  showProgress=FALSE
)
summary(mySwitchListAnalyzed_cellSurface_ORF)# 21 isoforms from 21 genes
### Add SignalP analysis
mySwitchListAnalyzed_cellSurface_ORF <- analyzeSignalP(
  switchAnalyzeRlist   = mySwitchListAnalyzed_cellSurface_ORF,
  pathToSignalPresultFile = "AS_analysis/isoformSwitchResult/my_surface_orf/signalp5output_protein_type.txt"
)
summary(mySwitchListAnalyzed_cellSurface_ORF)# 22isoforms from 21 genes
### Add DeepLoc2 analysis
mySwitchListAnalyzed_cellSurface_ORF <- analyzeDeepLoc2(
  switchAnalyzeRlist = mySwitchListAnalyzed_cellSurface_ORF,
  pathToDeepLoc2resultFile = "AS_analysis/isoformSwitchResult/my_surface_orf/results_DeepLoc2.csv",
  quiet = FALSE
)

unique(mySwitchListAnalyzed_cellSurface_ORF$topologyAnalysis$region_type)
topologyAnalysis_res = mySwitchListAnalyzed_cellSurface_ORF$topologyAnalysis

pdf('../../Figures/ORF_pie_surface.pdf',height = 4,width = 4)
p1 =ggpie::ggpie(topologyAnalysis_res,group_key = 'region_type', count_type = "full")
print(p1)
dev.off()

out_side_iso =mySwitchListAnalyzed_cellSurface_ORF$topologyAnalysis$isoform_id[which(mySwitchListAnalyzed_cellSurface_ORF$topologyAnalysis$region_type=="outside")]
pdf('../../Figures/cell_surface_ORF_MCAM-201.pdf',width =9,height = 2.5)
mySwitchListAnalyzed_cellSurface_ORF$signalPeptideAnalysis$isoform_id
cell_surface_isodat$genename[grep("ENST00000264036.6",cell_surface_isodat$isoID)]
cell_surface_isodat$transcript_name[grep("ENST00000264036.6",cell_surface_isodat$isoID)]
p1 =switchPlotTranscript(mySwitchListAnalyzed_cellSurface_ORF, gene = 'MCAM')+
  theme(legend.position = "bottom",legend.direction= 'horizontal')+
  ggtitle("MCAM-201")
print(p1)
dev.off()

pdf('../../Figures/cell_surface_ORF_LMBRD1-229.pdf',width =9,height = 2.5)
cell_surface_isodat$genename[grep("ENST00000649934.3",cell_surface_isodat$isoID)]
cell_surface_isodat$transcript_name[grep("ENST00000649934.3",cell_surface_isodat$isoID)]

p2 =switchPlotTranscript(mySwitchListAnalyzed_cellSurface_ORF, gene = 'LMBRD1')+
  theme(legend.position = "bottom",legend.direction= 'horizontal')+
  ggtitle("LMBRD1-229")
print(p2)
dev.off()

pdf('../../Figures/cell_surface_ORF_HLA-B-249.pdf',width =9,height = 2.5)
cell_surface_isodat$genename[grep("ENST00000399213.2",cell_surface_isodat$isoID)]
cell_surface_isodat$transcript_name[grep("ENST00000399213.2",cell_surface_isodat$isoID)]

p3 =switchPlotTranscript(mySwitchListAnalyzed_cellSurface_ORF, gene = 'PI4KA')+
  theme(legend.position = "bottom",legend.direction= 'horizontal')+
  ggtitle("PI4KA-202")
print(p3)
dev.off()

#Alternative splicing events and survival
OS.final = read.csv('AS_analysis/OS.final.Filtered.csv',header = T,row.names = 1)
event_out_deg_final = readRDS('AS_analysis/CTC_Duct_DASE_degOut_interFinal.rds')

eventpsi_data_tmp = lapply(event_out_deg_final,function(x){
  x_tmp_Epithelial = x[['Epithelial']]
  x_tmp_pEMT = x[['pEMT']]
  x_tmp_Mesenchymal = x[['Mesenchymal']]
  
  x_tmp_all=bind_cols(x_tmp_Epithelial[,grep('sortedByCoord.out.bam$',colnames(x_tmp_Epithelial))],
                      x_tmp_pEMT[,grep('sortedByCoord.out.bam$',colnames(x_tmp_pEMT))],
                      x_tmp_Mesenchymal[,grep('sortedByCoord.out.bam$',colnames(x_tmp_Mesenchymal))])
  rownames(x_tmp_all) = x_tmp_Epithelial$event_name
  duct_mean = x_tmp_epi$duct_mean
  i=1
  row_psi = x_tmp_all[i,]
  type_tmp=ifelse(as.numeric(row_psi[1,])<median(as.numeric(row_psi[1,])),'low','high')
  type_tmp_dat= data.frame(type = type_tmp,row.names = colnames(row_psi))
  colnames(type_tmp_dat)=rownames(x_tmp_all)[i]
  x_type = type_tmp_dat
  
  for (i in 2:nrow(x_tmp_all)) {
    print(i)
    row_psi = x_tmp_all[i,]
    type_tmp=ifelse(as.numeric(row_psi[1,])<median(as.numeric(row_psi[1,])),'low','high')
    type_tmp_dat= data.frame(type = type_tmp,row.names = colnames(row_psi))
    colnames(type_tmp_dat) = rownames(x_tmp_all)[i]
    x_type = cbind(x_type,type_tmp_dat)
  }
  rownames(x_type) = str_split(rownames(x_type),'[Aligned]',simplify = T)[,1]
  x_type = x_type[match(OS.final$SRRID,rownames(x_type)),]
  return(x_type)
})
saveRDS(eventpsi_data_tmp,'AS_analysis/eventpsi_data_type.rds')

all(rownames(eventpsi_data_tmp$A3SS)==OS.final$SRRID)

phe = OS.final
mySurv=with(phe,Surv(OS.days., status))
log_rank_p <- lapply(eventpsi_data_tmp , function(x){
  event = apply(x,2,function(event){
    phe$group=event
    if(length(table(phe$group))==1){
      p.val = NA
    }else{
      data.survdiff=survdiff(mySurv~group,data=phe)
      p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    }
    return(p.val)
  })
  return(event)
})
log_rank_p_cutoff = lapply(log_rank_p,function(x){x[which(x<0.05)]})

for (i in 1:5) {
  event = eventpsi_data_tmp[[i]]
  pval_sig= log_rank_p_cutoff[[i]]
  event = event[,names(pval_sig)]
  as_type = names(log_rank_p_cutoff)[i]
  for (z in colnames(event)) {
    type = event[,z]
    exon_loc = paste0(as_type,' ',
                      str_split(z,'[.]',simplify = T)[,3],'.',
                      str_split(z,'[.]',simplify = T)[,4],':',
                      str_split(z,'[.]',simplify = T)[,6],'-',
                      str_split(z,'[.]',simplify = T)[,7],':',
                      str_split(z,'[.]',simplify = T)[,5])
    phe$group = type
    fit = survfit(Surv(OS.days.,status)~group ,data = phe)
    surPlot = ggsurvplot(fit, data = phe,surv.median.line = "hv",pval = TRUE,
                         palette = c('#d7191c','#2b83ba'),
                         legend.title = exon_loc)
    path_to_save =gsub(':','_',exon_loc)
    path_to_save =gsub('[+]','.',path_to_save)
    pdf(file = paste0('../../Figures/',path_to_save,'.pdf'),width = 3,height = 3)
    print(surPlot)
    dev.off()
  }
}

