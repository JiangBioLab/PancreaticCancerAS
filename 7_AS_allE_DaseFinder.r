as_type_result = function(dir_path,dir_name,bam1,bam2,group1,group2,na_percent,method){
  library(stringr)
  library(limma)
  files_list = list.files(paste0(dir_path,dir_name)) 
  as_files = grep('.MATS.JC.txt',files_list,value = T)
  as_path = paste0(dir_path,dir_name,as_files)
  dase_list = c()
  source('/home/lxxiao/xiaolixing/pdac/result_dir/code/rmats_diff_func.R')
  for (i in 1:5) {
    as_res = read.csv(as_path[i],header = T,row.names = 1,sep = '\t')
    dase_final = dase(rmats_result_data = as_res,
                      group1 =group1,
                      group2 = group2,
                      bam1 = bam1,
                      bam2 = bam2,
                      na_percent = na_percent,
                      method = method)
    dase_final$result_all$ENSEMBL.Gene.ID =str_split(dase_final$result_all$GeneID,'[.]',simplify = T)[,1] 
    dase_final$FilByPvalue$ENSEMBL.Gene.ID =str_split(dase_final$FilByPvalue$GeneID,'[.]',simplify = T)[,1] 
    dase_final$FilByPadj$ENSEMBL.Gene.ID =str_split(dase_final$FilByPadj$GeneID,'[.]',simplify = T)[,1] 
    
    dase_final$result_all$as_type  =str_split(as_files,'[.]',simplify = T)[i,1]
    dase_final$FilByPvalue$as_type =str_split(as_files,'[.]',simplify = T)[i,1]
    dase_final$FilByPadj$as_type   =str_split(as_files,'[.]',simplify = T)[i,1]
    dase_list[[i]] = dase_final
    names(dase_list)[i] = str_split(as_files,'[.]',simplify = T)[i,1]
  }
  return(dase_list)
}