#下载、安装coloc
if(!require("remotes"))
  install.packages("remotes")
install.packages("dplyr")
library(remotes)
#install_github("chr1swallace/coloc",build_vignettes=TRUE) 这个下载方式不太行
install.packages("coloc")
library("coloc")
library(dplyr)
library(bigreadr)

#导入表型1（GWAS）数据
#格式（如果表型是二分类变量（case和control），输入文件二选一）
#rs编号rs_id、P值pval_nominal、SNP的效应值beta、效应值方差varbeta(即se)（用这个）
#rs编号rs_id、P值pval_nominal、case在所有样本中的比例s
gwas<-fread2(file="/public/home/biostat06/TWAS/Preeclampsia_gwas_summary/GCST90301704.tsv",header=T)
gwas_bind<-cbind(gwas$SNPID,gwas$p_value,gwas$beta,gwas$standard_error)
colnames(gwas_bind)<-c("rs_id","pval_nominal","beta","varbeta")

#导入表型2（eQTL）数据
#格式（如果表型是连续型变量，输入文件三选一）
#rs编号rs_id、P值pval_nominal、表型的标准差sdY
#rs编号rs_id、P值pval_nominal、效应值beta,效应值方差 varbeta, 样本量N,次等位基因频率 MAF
#rs编号rs_id、P值pval_nominal、次等位基因频率 MAF（用这个）
#相同基因提出
main_folder <- "/public/home/biostat06/TWAS/49_different_tissues"
output_main_folder <- "/public/home/biostat06/TWAS/49coloccccc/list_gene"
subfolders <- list.files(path = main_folder,pattern = "tsv.gz", fall.names = TRUE)
for (subfolder in subfolders) {
 dir_name <- basename(subfolder)
 dir_out <- file.path(output_main_folder,die_name)
 if (!dir.exists(dir_out)){
   dir.create(dir_out)
   }
   df <- fread2(subfolder,stringsAsFactors = FALSE,header = T)
   df_unique <- unique(df$molecular_trait_i)
   for (i in df_unique){
     df_sub <- subset(df,molecular_trait_id==i)
     write.csv(df_sub, file = paste0(dir_out, i, ".csv"), row.names = FALSE)
   }
}
#单个
#eqtl <- fread2(file="/public/home/biostat06/TWAS/49_different_tissues/Whole_Blood.tsv.gz", header=T)
#eqtl_unique<-unique(eqtl$molecular_trait_id)
#for (i in eqtl_unique){
#eqtl_sub<-subset(eqtl,molecular_trait_id==i)
#filename <- paste0("/public/home/biostat06/TWAS/49coloccccc/Whole Blood/", i, ".csv")
#write.csv(eqtl_sub, file = filename, row.names = FALSE)
#}
#找出每个基因所需要的列
main_folder <-"/public/home/biostat06/TWAS/49coloccccc/list_gene"
output_main_folder <-"/public/home/biostat06/TWAS/49coloccccc/sub"
subfolders <- list.dirs(main_folder, full.names = TRUE)
for (subfolder in subfolders) {
  subfolder_name <- basename(subfolder)
  output_subfolder <- file.path(output_main_folder, subfolder_name)
  if (!dir.exists(output_subfolder)) {
    dir.create(output_subfolder)
  }
file_list <- list.files(path = subfolder, pattern = "*.csv", full.names = TRUE)
 for (file_path in file_list) {
  df <- fread2(file_path, stringsAsFactors = FALSE, header = T)
  eqtl_bind <- cbind(df$rsid,df$pvalue,df$maf)
  colnames(eqtl_bind) <- c("rs_id","pval_nominal","MAF")
  file_name <- basename(file_path)
  output_file_path <- file.path(output_subfolder, file_name)
  write.csv(eqtl_bind, file = output_file_path, row.names = FALSE)
 }
}

#合并GWAS和eQTL数据
main_folder<-"/public/home/biostat06/TWAS/49coloccccc/sub"
output_main_folder <-"/public/home/biostat06/TWAS/coloc_result/merge"
subfolders <- list.dirs(main_folder, full.names = TRUE)
for (subfolder in subfolders) {
  subfolder_name <- basename(subfolder)
  output_subfolder <- file.path(output_main_folder, subfolder_name)
  if (!dir.exists(output_subfolder)) {
    dir.create(output_subfolder)
  }
  file_list <- list.files(path = subfolder, pattern = "*.csv", full.names = TRUE)
  for (file_path in file_list) {
    df <- fread2(file_path, stringsAsFactors = FALSE, header = T)
    input <- merge(df, gwas_bind, by="rs_id", all=FALSE, suffixes=c("_eqtl","_gwas"))
    file_name <- basename(file_path)
    output_file_path <- file.path(output_subfolder, file_name)
    write.csv(eqtl_bind, file = output_file_path, row.names = FALSE)
  }
}
#input <- merge(eqtl_bind, gwas_bind, by="rs_id", all=FALSE, suffixes=c("_eqtl","_gwas"))


#共定位分析

result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_gwas, type="cc", s=0.33, N=50000),
                    dataset2=list(pvalues=input$pval_nominal_eqtl, type="quant", N=10000), 
                    MAF=input$maf)
#dataset1的type="cc"指的是GWAS的表型是二分类（case和control）；dataset2的type="quant"指的是eQTL的表型（基因表达量）是连续型；N指样本量
#s=case/总数，N=总数

# 筛选共定位的位点
need_result=result$results %>% filter(SNP.PP.H4 > 0.95) #很多文献PPA > 0.95的位点是共定位位点，有些会放松到0.75。这里假定后验概率大于0.95为共定位位点

