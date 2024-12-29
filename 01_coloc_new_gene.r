# Perform coloc analysis
# Load packages
library(bigreadr)
library(dplyr)
library(coloc)
library(optparse)
library(ggplot2)
library(optparse)
# Input parameters
args_list = list(
  make_option("--tissue", type="character", default=NULL,
              help="INPUT: tissue information",
              metavar="character")
)
opt_parser <- OptionParser(option_list = args_list)
opt <- parse_args(opt_parser)

# Fix parameter
COLOC_P1 <- 1E-04 
COLOC_P2 <- 1E-04 
COLOC_P3 <- 1E-05
#WIN_PLOT <- 2E05

# Function 1: process summary_info.txt
summary_file <- list(summ1 = '/public/home/biostat06/TWAS/Preeclampsia_gwas_summary/outcome.txt', 
                     type1 = 'cc',
                     summ2 = paste0('/public/home/biostat06/TWAS/49_different_tissues/',opt$tissue,'.tsv.gz'), 
                     type2 = 'quant'
                     #use_model1 = summary_info[5, 2],
                     #use_model2 = summary_info[6, 2],
                     #PPH4 = as.numeric(summary_info[7, 2])
                     )
gwas_file <- fread2(summary_file$summ1)
eqtl_file <- fread2(summary_file$summ2)
save(eqtl_file,file=paste0('/public/home/biostat06/TWAS/49_different_tissues/',opt$tissue,'.RData'))
gene_list <- unique(eqtl_file$molecular_trait_id)

system(paste0('mkdir /public/home/biostat06/TWAS/COLOC_RESULT/coloc_result/',opt$tissue))

coloc_result <- plyr::alply(gene_list,1,function(gene){
    # intersect snp set
    eqtl_file_sub <- eqtl_file %>% filter(molecular_trait_id == gene)
    sub_snp <- eqtl_file_sub$rsid
    gwas_file_sub <- gwas_file %>% filter(SNP %in% sub_snp)
    eqtl_file_sub <- eqtl_file_sub %>% filter(rsid %in% gwas_file_sub$SNP)
    # 
    eqtl_file_sub <- eqtl_file_sub[!duplicated(eqtl_file_sub$rsid),]
    gwas_file_sub <- gwas_file_sub[!duplicated(gwas_file_sub$SNP),]
    
    # define sample size
    N_eqtl <- 1000
    N_gwas <- 296824 # 16743 case
    # define dataset
    MAF1 <- ifelse(gwas_file_sub$FRQ<0.5, gwas_file_sub$FRQ, 1-gwas_file_sub$FRQ)
    gwas_dat <- list(snp=gwas_file_sub$SNP,
                     beta=gwas_file_sub$BETA,
                     varbeta=gwas_file_sub$SE^2,
                     N=N_gwas,
                     MAF=MAF1,
                     type=summary_file$type1)

    MAF2 <- ifelse(eqtl_file_sub$maf<0.5, eqtl_file_sub$maf, 1-eqtl_file_sub$maf)
    eqtl_dat <- list(snp=eqtl_file_sub$rsid,
                     beta=eqtl_file_sub$beta,
                     varbeta=eqtl_file_sub$se^2,
                     N=N_eqtl,
                     MAF=MAF2,
                     type=summary_file$type2)
    coloc_result <- coloc.abf(gwas_dat,eqtl_dat,p1=COLOC_P1,p2=COLOC_P2,p12=COLOC_P3)
    PPH <- as.vector(coloc_result$summary)
    res <- data.frame(coloc_result$results[, c("snp", "SNP.PP.H4")],nsnps=PPH[1],
                      PP.H0.abf=PPH[2],PP.H1.abf=PPH[3],
                      PP.H2.abf=PPH[4],PP.H3.abf=PPH[5],
                      PP.H4.abf=PPH[6],gene=gene)
    fwrite2(res,file=paste0('/public/home/biostat06/TWAS/COLOC_RESULT/coloc_result/',opt$tissue,'/',gene,'.txt'),
            row.names=F,col.names=T,quote=F,sep='\t')
})


