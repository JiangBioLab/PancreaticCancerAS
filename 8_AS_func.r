dase = function(rmats_result_data,bam1,bam2,
                group1,group2,na_percent,method){
  library(stringr)
  library(limma)
  library(fdrtool)
  rmats_data = rmats_result_data
  rmats_data$event_num = paste0('event',1:nrow(rmats_data))
  rownames(rmats_data) = rmats_data$event_num 
  
  Level1_data = as.data.frame(apply(strsplit2(rmats_data$IncLevel1,','),2,as.numeric))
  rownames(Level1_data)=paste0('event',1:nrow(Level1_data))
  colnames(Level1_data)=paste0(group1,'_',bam1)
  
  Level2_data = as.data.frame(apply(strsplit2(rmats_data$IncLevel2,','),2,as.numeric))
  rownames(Level2_data)=paste0('event',1:nrow(Level2_data))
  colnames(Level2_data)=paste0(group2,'_',bam2)
  
  rmats_psi_raw = cbind(Level1_data,Level2_data)
  na_psi = as.data.frame(is.na(rmats_psi_raw))
  na_per_col = colSums(na_psi)/nrow(na_psi)
  
  rmats_psi = rmats_psi_raw[,na_per_col<0.9]
  print(paste('Samples Name:',paste(colnames(rmats_psi),collapse = ' '),sep = ' '))
  
  Level1_data_per = Level1_data[,intersect(colnames(rmats_psi),colnames(Level1_data))]
  na_psi = as.data.frame(is.na(Level1_data_per))
  na_per_row = rowSums(na_psi)/ncol(na_psi)
  Level1_data_per = Level1_data_per[na_per_row<na_percent,]
  
  Level2_data_per = Level2_data[,intersect(colnames(rmats_psi),colnames(Level2_data))]
  na_psi = as.data.frame(is.na(Level2_data_per))
  na_per_row = rowSums(na_psi)/ncol(na_psi)
  Level2_data_per = Level2_data_per[na_per_row<na_percent,]
  
  event_final = intersect(rownames(Level1_data_per),rownames(Level2_data_per))
  Level1_data_final = Level1_data_per[event_final,]
  Level2_data_final = Level2_data_per[event_final,]
  rmats_psi_na = cbind(Level1_data_final,Level2_data_final)
  
  
  if(method =='knn'){
    library(impute)
  set.seed(2022)
  #rmats_psi_na  = rmats_psi[na_per_row<na_percent,]
  #Level1_data_na = Level1_data[rownames(rmats_psi_na),intersect(colnames(rmats_psi_na),colnames(Level1_data))]
  #Level2_data_na = Level2_data[rownames(rmats_psi_na),intersect(colnames(rmats_psi_na),colnames(Level2_data))]
  
  data_imp = impute.knn(as.matrix(rmats_psi_na),colmax =1)
  data_imp = as.data.frame(data_imp$data)
  rmats_merge_psi = data_imp
  
  }else{
    
    #rmats_psi_na  = rmats_psi[na_per_row<na_percent,]
    
    Level1_data_na = Level1_data_final[rownames(rmats_psi_na),intersect(colnames(rmats_psi_na),colnames(Level1_data_final))]
    Level2_data_na = Level2_data_final[rownames(rmats_psi_na),intersect(colnames(rmats_psi_na),colnames(Level2_data_final))]
    
  Level1_data_imp = apply(Level1_data_na, 1, function(x){
    x[is.na(x)]=mean(x,na.rm = T)
    return(x)})
  Level1_data_imp = as.data.frame(t(Level1_data_imp))
  dim(Level1_data_imp)
  
  #Level2_data$mean = rowMeans(Level2_data_na,na.rm = T)
  Level2_data_imp = apply(Level2_data_na, 1, function(x){
    x[is.na(x)]=mean(x,na.rm = T)
    return(x)})
  Level2_data_imp  = as.data.frame(t(Level2_data_imp))
  dim(Level2_data_imp)
  rmats_merge_psi = cbind(Level1_data_imp,Level2_data_imp)
  }
  
  
  #dim(rmats_merge_psi) #746
  #table(is.na(rmats_merge_psi))
  
  psi_obse = rmats_merge_psi[,colnames(Level1_data_imp)]
  psi_othe = rmats_merge_psi[,colnames(Level2_data_imp)]
  
  p = c()
  delt = c()
  logFC = c()
  Level1_mean = c()
  Level2_mean = c()
  print("-----p_value calculate------")
  for (i in 1:nrow(rmats_merge_psi)) {
    psi_obs = as.numeric(psi_obse[i,])
    psi_oth = as.numeric(psi_othe[i,])
    if(all(psi_obs==psi_oth)){
      w_p = NA
    }else{
      w_t = wilcox.test(psi_obs,psi_oth)
      w_p = w_t$p.value
      #w_t = t.test(psi_1,psi_0.1,paired = F)
      #w_p = w_t$p.value
    }
    p = c(p,w_p)
    delt = c(delt,round(mean(psi_obs)-mean(psi_oth),3))
    logFC = c(logFC,log10(mean(psi_obs+0.001)/mean(psi_oth+0.001)))
    Level1_mean = c(Level1_mean,mean(psi_obs))
    Level2_mean = c(Level2_mean,mean(psi_oth))
  }
  print(paste0("-----finished(events N= ",nrow(psi_obse),")------"))
  #apply(psi_fna, 1, function(x){ w_t = wilcox.test(as.numeric(x[1:42]),as.numeric(x[43:63])) w_p = w_t$p.value return(w_p)})
  p_data = data.frame(event_num = rownames(psi_obse),
                      p_self = p,
                      delt = delt,
                      logFC = logFC,
                      Level1_mean = Level1_mean,
                      Level2_mean = Level2_mean)
  p_narm = p_data[!is.na(p_data$p_self),]
  #table(p_narm$p_self<0.05)
  #padj = round(p.adjust(p_narm$p_self,method = 'BH'),3)
  
  fdr=fdrtool(p_narm$p_self,statistic="pvalue",plot = F)
  #table(fdr$qval<0.05)
  padj = round(fdr$qval,3)
  #table(padj<0.05)
  p_narm$padj = padj
  rmats_dase = rmats_merge_psi[p_narm$event_num,]
  #all(rownames(rmats_dase) ==p_narm$event_num )
  rmats_result = cbind(rmats_dase,p_narm)
  rmats_data_result = rmats_data[p_narm$event_num,]
  
  #all(rownames(rmats_data_result)==rownames(rmats_result))
  res_final = cbind(rmats_data_result,rmats_result)
  res_final_p    = res_final[which(res_final$p_self<0.05 & abs(res_final$logFC )>1),]
  
  res_final_padj = res_final[which(res_final$padj<0.05 & abs(res_final$logFC)>1 ),]
  
  print(paste0(nrow(res_final_padj),' diff-AS were found(FDR)'))
  result = list(result_all = res_final, 
                FilByPvalue= res_final_p,
                FilByPadj  = res_final_padj)
  return(result)
}
