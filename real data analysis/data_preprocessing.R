library(data.table)
#################################### merge GTEX and GCTA ######################################
merge_gtex_gcta_target <- function(){
  setwd("D:/gbm/real data/target data")
  gene_expre <- fread("TCGA-GBM.star_counts.tsv", header = T, sep = '\t', data.table = F)

  filtered_exprs <- gene_expre#[rowMeans((2^(gene_expre[,-1])-1) > 1) > 0.1, ]
  dim(filtered_exprs)

  exp_map <- fread("gencode.v36.annotation.gtf.gene.probemap", header = T, sep = '\t', data.table = F)
  head(exp_map)[1:5, 1:5]
  exp_map <- exp_map[ , c(1, 2)]
  gene_expre <- merge(gene_expre, exp_map, by.x  = "Ensembl_ID", by.y = "id" )
  # 
  dlbc.phe <- fread("TCGA-GBM.clinical.tsv", header = T, sep = '\t', data.table = F)
  dlbc.phe$sample[1:5]

  # 
  rownames(dlbc.phe) <- dlbc.phe$sample

  # 
  table(dlbc.phe$sample_type.samples)

  # Primary Tumor
  dlbc.phe.t <- filter(dlbc.phe, sample_type.samples == "Primary Tumor")
  # 
  merge_phe_fpkm <- intersect(rownames(dlbc.phe.t), colnames(gene_expre))
  merge_phe_fpkm = c('gene',merge_phe_fpkm )
  # 
  filtered_exprs <- gene_expre[ , merge_phe_fpkm]
  dim(filtered_exprs)
  #saveRDS(filtered_exprs, file='dlbc.exp.rds')
  
  
  
  #################################### GTEx ######################################
  setwd("D:/gbm/real data/GTEX")
  # 
  gtex.exp <- fread("gtex_RSEM_Hugo_norm_count", header = T, sep = '\t', data.table = F)
  # gtex.exp[1:5, 1:4]
  # # saveRDS(gtex.exp, file = 'gtex.exp.rds')
  gtex.pro <- fread("probeMap_hugo_gencode_good_hg38_v23comp.probemap", header = T, sep = '\t', data.table = F)
  head(gtex.pro)
  dim(gtex.pro) #
  
  # 
  gtex.pro <- gtex.pro[, c(1,2)]
  # 
  # 
  gtex.fpkm.pro <- merge(gtex.pro, gtex.exp, by.y ="sample", by.x = "id" )
  # save(gtex.fpkm.pro,file = 'gtex.fpkm.pro.Rdata')
  #load('gtex.fpkm.pro.Rdata')
  gtex.fpkm.pro <- get("gtex.fpkm.pro")
  gtex.phe <- fread("GTEX_phenotype", header = T, sep = '\t', data.table = F)
  rownames(gtex.phe) <- gtex.phe$Sample
  
  # 
  colnames(gtex.phe) <- c("Sample", "body_site_detail (SMTSD)", "primary_site", "gender", "patient", "cohort")
  table(gtex.phe$primary_site)
  gtex.phe.s <- filter(gtex.phe, primary_site == "Brain")
  
  # 
  merge_phe_fpkm_gtex <- intersect(rownames(gtex.phe.s), colnames(gtex.fpkm.pro)) # 444
  gtex.s <- gtex.fpkm.pro[ , c("gene", merge_phe_fpkm_gtex)]
  
  # 
  gtex.s <- distinct(gtex.s, gene, .keep_all = T)
  rownames(gtex.s) <- gtex.s$gene
  gtex.s <- gtex.s[ , -1]
  dim(gtex.s)
  all.data <- merge(gtex.s,t(filtered_exprs), by = 0)
  head(gtex.s)[1:5, 1:4]
  head(filtered_exprs)[1:5, 1:4]
  all.data <- column_to_rownames(all.data, "Row.names")
  #saveRDS(all.data, file = 'target_all.data.rds')
  all.data
  list(all.data = filtered_exprs, gene_expre = filtered_exprs)
}
merge_gtex_gcta_source1 <- function(){
  setwd("D:/gbm/real data/source 1")
  
  gene_expre <- fread("HiSeqV2", header = T, sep = '\t', data.table = F)
  
  filtered_exprs <- gene_expre#[rowMeans((2^(gene_expre[,-1])-1) > 1) > 0.1, ]
  dim(filtered_exprs)
  head(filtered_exprs)[1:5, 1:5]
  b = filtered_exprs$sample
  filtered_exprs<-as.data.frame(t(filtered_exprs))
  
  colnames(filtered_exprs) = b
  filtered_exprs$sample<-rownames(filtered_exprs)
  head(filtered_exprs)[1:5, 1:5]
  
  # 
  dlbc.phe <- fread("TCGA.GBM.sampleMap_GBM_clinicalMatrix", header = T, sep = '\t', data.table = F)
  dlbc.phe$sampleID[1:5]
  
  # 
  rownames(dlbc.phe) <- dlbc.phe$sampleID
  
  #
  table(dlbc.phe$sample_type)
  # extract brain tumor data
  dlbc.phe.t <- filter(dlbc.phe, sample_type == "Primary Tumor")
  
  #
  merge_phe_fpkm <- intersect(rownames(dlbc.phe.t), filtered_exprs$sample)
  
  #
  filtered_exprs <- filtered_exprs[ , filtered_exprs$sample %in% merge_phe_fpkm]
  
  dim(filtered_exprs)
  
  #saveRDS(filtered_exprs, file='dlbc.exp.rds')
  
  
  
  #################################### GTEx ######################################
  setwd("D:/gbm/real data/GTEX")
  # 
  gtex.exp <- fread("gtex_RSEM_Hugo_norm_count", header = T, sep = '\t', data.table = F)
  # gtex.exp[1:5, 1:4]
  # # saveRDS(gtex.exp, file = 'gtex.exp.rds')
  gtex.pro <- fread("probeMap_hugo_gencode_good_hg38_v23comp.probemap", header = T, sep = '\t', data.table = F)
  head(gtex.pro)
  dim(gtex.pro) #
  
  # 
  gtex.pro <- gtex.pro[, c(1,2)]
  # 
  # 
  gtex.fpkm.pro <- merge(gtex.pro, gtex.exp, by.y ="sample", by.x = "id" )
  # save(gtex.fpkm.pro,file = 'gtex.fpkm.pro.Rdata')
  #load('gtex.fpkm.pro.Rdata')
  gtex.fpkm.pro <- get("gtex.fpkm.pro")
  gtex.phe <- fread("GTEX_phenotype", header = T, sep = '\t', data.table = F)
  rownames(gtex.phe) <- gtex.phe$Sample
  
  # 
  colnames(gtex.phe) <- c("Sample", "body_site_detail (SMTSD)", "primary_site", "gender", "patient", "cohort")
  table(gtex.phe$primary_site)
  gtex.phe.s <- filter(gtex.phe, primary_site == "Brain")
  
  # 
  merge_phe_fpkm_gtex <- intersect(rownames(gtex.phe.s), colnames(gtex.fpkm.pro)) # 444
  gtex.s <- gtex.fpkm.pro[ , c("gene", merge_phe_fpkm_gtex)]
  
  # 
  gtex.s <- distinct(gtex.s, gene, .keep_all = T)
  rownames(gtex.s) <- gtex.s$gene
  gtex.s <- gtex.s[ , -1]
  dim(gtex.s)
  all.data <- merge(gtex.s,t(filtered_exprs), by = 0)
  head(gtex.s)[1:5, 1:4]
  head(filtered_exprs)[1:5, 1:4]
  all.data <- column_to_rownames(all.data, "Row.names")
  #saveRDS(all.data, file = 'source1_all.data.rds')
  all.data
  list(all.data = filtered_exprs, gene_expre = filtered_exprs)
}
merge_gtex_gcta_source2 <- function(){
  setwd("D:/gbm/real data/source 2")
  
  gene_expre <- fread("TCGA-HNSC.star_counts.tsv", header = T, sep = '\t', data.table = F)
  
  filtered_exprs <- gene_expre#[rowMeans((2^(gene_expre[,-1])-1) > 1) > 0.1, ]
  dim(filtered_exprs)
  
  exp_map <- fread("gencode.v36.annotation.gtf.gene.probemap", header = T, sep = '\t', data.table = F)
  head(exp_map)[1:5, 1:5]
  exp_map <- exp_map[ , c(1, 2)]
  gene_expre <- merge(gene_expre, exp_map, by.x  = "Ensembl_ID", by.y = "id" )
  # 
  dlbc.phe <- fread("TCGA-HNSC.clinical.tsv", header = T, sep = '\t', data.table = F)
  dlbc.phe$sample[1:5]
  
  # 
  rownames(dlbc.phe) <- dlbc.phe$sample
  
  # 
  table(dlbc.phe$sample_type.samples)
  
  # Primary Tumor
  dlbc.phe.t <- filter(dlbc.phe, sample_type.samples == "Primary Tumor")
  # 
  merge_phe_fpkm <- intersect(rownames(dlbc.phe.t), colnames(gene_expre))
  merge_phe_fpkm = c('gene',merge_phe_fpkm )
  # 
  filtered_exprs <- gene_expre[ , merge_phe_fpkm]
  dim(filtered_exprs)
  
  #saveRDS(filtered_exprs, file='dlbc.exp.rds')
  
  
  
  #################################### GTEx ######################################
  setwd("D:/gbm/real data/GTEX")
  # 
  gtex.exp <- fread("gtex_RSEM_Hugo_norm_count", header = T, sep = '\t', data.table = F)
  # gtex.exp[1:5, 1:4]
  # # saveRDS(gtex.exp, file = 'gtex.exp.rds')
  gtex.pro <- fread("probeMap_hugo_gencode_good_hg38_v23comp.probemap", header = T, sep = '\t', data.table = F)
  head(gtex.pro)
  dim(gtex.pro) #
  
  # 
  gtex.pro <- gtex.pro[, c(1,2)]
  # 
  # 
  gtex.fpkm.pro <- merge(gtex.pro, gtex.exp, by.y ="sample", by.x = "id" )
  # save(gtex.fpkm.pro,file = 'gtex.fpkm.pro.Rdata')
  #load('gtex.fpkm.pro.Rdata')
  gtex.fpkm.pro <- get("gtex.fpkm.pro")
  gtex.phe <- fread("GTEX_phenotype", header = T, sep = '\t', data.table = F)
  rownames(gtex.phe) <- gtex.phe$Sample
  
  # 
  colnames(gtex.phe) <- c("Sample", "body_site_detail (SMTSD)", "primary_site", "gender", "patient", "cohort")
  table(gtex.phe$primary_site)
  gtex.phe.s <- filter(gtex.phe, primary_site == "Brain")
  
  # 
  merge_phe_fpkm_gtex <- intersect(rownames(gtex.phe.s), colnames(gtex.fpkm.pro)) # 444
  gtex.s <- gtex.fpkm.pro[ , c("gene", merge_phe_fpkm_gtex)]
  
  # 
  gtex.s <- distinct(gtex.s, gene, .keep_all = T)
  rownames(gtex.s) <- gtex.s$gene
  gtex.s <- gtex.s[ , -1]
  dim(gtex.s)
  all.data <- merge(gtex.s,t(filtered_exprs), by = 0)
  head(gtex.s)[1:5, 1:4]
  head(filtered_exprs)[1:5, 1:4]
  all.data <- column_to_rownames(all.data, "Row.names")
  #saveRDS(all.data, file = 'source2_all.data.rds')
  all.data
  list(all.data = filtered_exprs, gene_expre = filtered_exprs)
}
result = merge_gtex_gcta_source1()
source1_all.data = result$all.data
source1_gene_expre = result$gene_expre

result = merge_gtex_gcta_source2()
source2_all.data = result$all.data
source2_gene_expre = result$gene_expre

result = merge_gtex_gcta_target()
target_all.data = result$all.data
target_gene_expre = result$gene_expre
#################################### significant differential expression analysis ###################
DE_fun <- function(all.data){
  data_tumor <- all.data[,grep('TCGA',colnames(all.data))]
  data_normal <- all.data[,grep('GTEX',colnames(all.data))]
  
  #
  exp_gbm <- cbind(data_tumor,data_normal)
  
  group <- c(rep('tumor',ncol(data_tumor)),rep('normal',ncol(data_normal)))
  group <- factor(group,levels=c("normal","tumor"))
  
  colData<-data.frame(row.names=colnames(exp_gbm),
                      group=group)
  colData$group<-factor(colData$group,levels=c("normal","tumor"))
  head(colData)
  exp_gbm_int <- 2^(exp_gbm) - 1
  exp_gbm_int <- apply(exp_gbm_int,2,as.integer)
  rownames(exp_gbm_int)<-rownames(exp_gbm)
  ############
  
  library(DESeq2)
  
  dds<-DESeqDataSetFromMatrix(countData=exp_gbm_int,
                              colData=colData,
                              design=~group)
  #DESeq2
  dds <- DESeq(dds)
  
  res<-results(dds,contrast=c("group",rev(levels(group))))
  resOrdered<-res[order(res$padj),]
  
  
  DEG<-as.data.frame(resOrdered)
  DEG_deseq2<-na.omit(DEG)
  DEG_deseq2
  #
  head(DEG_deseq2)
  
  #save(DEG_deseq2,file='DEG_deseq2.Rdata')
}
source1_DEG_deseq2 = DE_fun(source1_all.data)
source2_DEG_deseq2 = DE_fun(source2_all.data)
target_DEG_deseq2 = DE_fun(target_all.data)

#################################### merge protein data ######################################

analysis1_data <- function(RPPA_RBN, target_DEG_deseq2, target_gene_expre){

  RPPA_RBN = RPPA_RBN[which(RPPA_RBN$peptide_target=='EGFR'),]
  t(RPPA_RBN )
  RPPA_RBN = data.frame(sample = rownames(t(RPPA_RBN ))[-1], EGFR_pro = t(RPPA_RBN )[-1])
  dim(RPPA_RBN)
  #load("target_DEG_deseq2.Rdata")
  logFC=2.5
  P.Value=0.01
  k1<-(target_DEG_deseq2$pvalue<P.Value)&(target_DEG_deseq2$log2FoldChange<-logFC)
  k2<-(target_DEG_deseq2$pvalue<P.Value)&(target_DEG_deseq2$log2FoldChange>logFC)
  target_DEG_deseq2<-mutate(target_DEG_deseq2,change=ifelse(k1,"down",ifelse(k2,"up","stable")))
  table(target_DEG_deseq2$change)
  #
  deg_filter<-function(df){
    rownames(df)[df$change!="stable"]
  }
  
  #
  all_degs<-deg_filter(target_DEG_deseq2)
  #
  
  target_gene_expre<-target_gene_expre[rownames(target_gene_expre)%in%all_degs,]
  dim(target_gene_expre)
  target_gene_expre = t(target_gene_expre)
  target_gene_expre = cbind(c(rownames(target_gene_expre)),target_gene_expre )
  colnames(target_gene_expre)[1] = 'sample'
  target_gene_expre = as.data.frame(target_gene_expre)
  dim(target_gene_expre)
  target_gene_expre = target_gene_expre[rownames(target_gene_expre) %in% RPPA_RBN$sample ,]
  dim(target_gene_expre)
  data <-left_join(RPPA_RBN ,target_gene_expre,by="sample")
  data
}
#################################### target data
setwd("D:/gbm/GDC TCGA Glioblastoma")
RPPA_RBN <- fread("TCGA-GBM.protein.tsv", header = T, sep = '\t', data.table = F)
head(RPPA_RBN)[1:5, 1:5]
data0 = analysis1_data(RPPA_RBN, target_DEG_deseq2, target_gene_expre)
#################################### source 1
setwd("D:/gbm/TCGA Glioblastoma (GBM) (25 datasets)")
RPPA_RBN <- fread("RPPA_RBN", header = T, sep = '\t', data.table = F)
head(RPPA_RBN)[1:5, 1:5]
data1 = analysis1_data(RPPA_RBN, source1_DEG_deseq2, source1_gene_expre)
#################################### source 2
setwd("D:/gbm/GDC TCGA Head and Neck Cancer (HNSC)")
RPPA_RBN <- fread("TCGA-HNSC.protein.tsv", header = T, sep = '\t', data.table = F)
head(RPPA_RBN)[1:5, 1:5]
data2 = analysis1_data(RPPA_RBN, source2_DEG_deseq2, source2_gene_expre)
#################################### Retained genes within EGFR pathway ######################################

library(clusterProfiler)
library(org.Hs.eg.db)

# 
gene_list <- colnames(data0)[-(1:2)]
gene_list_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

ego <- enrichGO(gene = gene_list_entrez,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE) 
egfr_related <- ego[grep("epidermal growth factor receptor signaling pathway", ego$Description), ]
egfr_related$geneID
egfr_genes <- unlist(strsplit(egfr_related$geneID, "/"))
#################################### final analysis data ######################################
g_0 = colnames(data0)[!colnames(data0) %in% c('sample','EGFR_pro')]
g_1 = colnames(data1)[!colnames(data1) %in% c('sample','EGFR_pro')]
g_2 = colnames(data2)[!colnames(data2) %in% c('sample','EGFR_pro')]

id = intersect(intersect(g_1,g_2),g_0)

analysis2_data <- function(data, id, egfr_genes){
  g = colnames(data)[!colnames(data) %in% c('sample','EGFR_pro')]
  g = intersect(g,id)
  posvec = NULL
  for(j in 1:length(egfr_genes)){
    pos = which(g==egfr_genes[j])
    posvec = c(posvec,pos)
  }
  g = g[posvec]
  
  x_s1 <- as.matrix(data[,!colnames(data) %in% c('sample','EGFR_pro')])
  y_s1 <- as.numeric(data$EGFR_pro)
  x_s1 = apply(x_s1[,posvec], c(1, 2), as.numeric) 
  data.frame(X = x_s1, Y = y_s1)
}
#################################### target data
target_data = analysis2_data(data0, egfr_genes)
#################################### source 1
source1_data = analysis2_data(data1, egfr_genes)
#################################### source 2
source2_data = analysis2_data(data2, egfr_genes)