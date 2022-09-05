############## 1. Library preparation ################
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library("ggpubr")
library(ggsci)
library(factoextra)
library(RColorBrewer)
############## 2. Function preparation ###########
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}
############## 3. Data preparation ####################
count <- read.table("rawdata/vlidation_gene_count.csv",header = T,sep = ",",row.names = 1)
gene_list <- count[,25]
gene_length <- count[,26]
count <- count[,1:24]

### clinical data
metadata <- read.table("rawdata/validation_metadata.csv",header = T,sep = ",")
metadata <- metadata[,c(1:5)]
rownames(metadata) <- metadata$ori_name
#### count to tpm ###
tpm <- read.csv("rawdata/validation_tpm.csv",row.names = 1)
rownames_tpm <- rownames(tpm)
tpm <- apply(tpm,2,function(x) as.integer(x))
rownames(tpm) <- rownames_tpm

human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

### deci
Deci_150_validation_DEGs <- read_excel("results/validation/Deci_150_validation_DEGs.xls")

Deci_30_validation_DEGs <- read_excel("results/validation/Deci_30_validation_DEGs.xls")

### filter
summary(abs(Deci_30_validation_DEGs$log2FoldChange))
#logFC_cutoff <- with(Deci_30_validation_DEGs,mean(abs(log2FoldChange)) + 1*sd(abs(log2FoldChange)))
Deci_30_validation_diffMat = Deci_30_validation_DEGs [Deci_30_validation_DEGs$pvalue<=0.05 & abs(Deci_30_validation_DEGs$log2FoldChange)>=0.50,]

summary(abs(Deci_150_validation_DEGs$log2FoldChange))
logFC_cutoff <- with(Deci_150_validation_DEGs,mean(abs(log2FoldChange)) + 1*sd(abs(log2FoldChange)))

Deci_150_validation_diffMat = Deci_150_validation_DEGs [Deci_150_validation_DEGs$pvalue<=0.05 & abs(Deci_150_validation_DEGs$log2FoldChange)>=0.12,]


dim(Deci_150_validation_diffMat)
dim(Deci_30_validation_diffMat)
### 11 67

## convert
## deci 30
deci_30_gene_symbol <-HGNC[HGNC$Ensembl.gene.ID %in% Deci_30_validation_diffMat$...1,]
deci_30_gene_symbol <- deci_30_gene_symbol$NCBI.gene.ID
length(deci_30_gene_symbol)

reactome.total_30 <- enrichPathway(gene       = deci_30_gene_symbol ,
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.03,
                                qvalueCutoff  = 1)
write.table(data.frame(reactome.total_30),paste0("results/validation/deci_30_validation_total_gene_reactome_test.xls"),sep = "\t",quote = F,row.names = F)


## deci 150
deci_150_gene_symbol <-HGNC[HGNC$Ensembl.gene.ID %in% Deci_150_validation_diffMat$...1,]
deci_150_gene_symbol <- deci_150_gene_symbol$NCBI.gene.ID
length(deci_150_gene_symbol)

reactome.total_150 <- enrichPathway(gene       = deci_150_gene_symbol ,
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.03,
                                qvalueCutoff  = 1)
write.table(data.frame(reactome.total_150),paste0("results/validation/deci_150_validation_total_gene_reactome_test.xls"),sep = "\t",quote = F,row.names = F)

pathway <- c("Cyclin E associated events during G1/S transition", 
             "Cyclin A:Cdk2-associated events at S phase entry",
             "Mitotic G1 phase and G1/S transition",
             "Transcriptional Regulation by TP53",
             "G2/M Transition",
             "Cdc20:Phospho-APC/C mediated degradation of Cyclin A",
             "Mitotic G2-G2/M phases")

reactome.total_30@result <- reactome.total_30@result[reactome.total_30@result$Description %in% pathway, ]

pdf(paste0("results/validation/new_deci_30_validation_total_gene_Reactome_dotplot.pdf"),width = 8,height = 10)
dotplot(reactome.total_30)#点状图
dev.off()  
  
  
  
