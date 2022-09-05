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
### count
count <- read.table("rawdata/baiduII_gene_count_.txt",header = T, row.names = 1)
nrow(count)
nrow(na.omit(count))
row_list <- c()

for (rowname in rownames(count)) {
  row_list <- c(row_list , strsplit(rowname, "_")[[1]][1])
}
rownames(count) <- row_list
### T vs N
DEG_Tumor_vs_Normal_sig <- as.data.frame(read_excel("rawdata/DEG_Tumor_vs_Normal_sig.xlsx"))
summary(abs(DEG_Tumor_vs_Normal_sig$log2FC))
row_list <- c()

for (rowname in DEG_Tumor_vs_Normal_sig$EnsemblGene_GeneSymbol) {
  row_list <- c(row_list , strsplit(rowname, "_")[[1]][1])
}
dim(DEG_Tumor_vs_Normal_sig)
DEG_Tumor_vs_Normal_sig$ensemble <- row_list
## only check tvsm sig gene
count <- count[rownames(count) %in% DEG_Tumor_vs_Normal_sig$ensemble,]
dim(count)
## clinical data 
metadata <- read.csv(file = "rawdata/baiduII_clinical_data.csv", header = TRUE)
metadata_gender <- metadata[,c("sample.ID............BDESCC2..", "Gender",
                               "Location.1","Smoking.history","TNM.stage.the.Eighth.Edition.",
                               "Grade","Age")]
colnames(metadata_gender)[1] <- "Sample_ID"
metadata_gender = metadata_gender[order(metadata_gender$Sample_ID),]

## change to integer
data = round(as.matrix(count))

### tpm
tpm <- read.csv("rawdata/baiduII_tpm.csv",row.names = 1,header = T)


############## 4. DEG list preparation (DESeq2) ##################
sorted_colname <- sort(colnames(data))
group_list <- gsub(".*T","T",sorted_colname)
group_list <- gsub(".*N","N",group_list)
condition = group_list

subject <- gsub("T","",sorted_colname)
subject <- gsub("N","",subject)
subject <- gsub("X","",subject)

coldata <- data.frame(row.names = sorted_colname, group_list)
coldata <- cbind(coldata, subject)
colnames(coldata) <- c("condition","subject")
coldata$ori_id <- rownames(coldata)
coldata$subject <- as.numeric(coldata$subject)
coldata <- merge(x =metadata_gender, y =coldata, by.x = "Sample_ID", by.y = "subject")
coldata$condition <- ifelse(coldata$condition == "T", "Tumor", "Normal")

### visualization tpm prep
tpm <- tpm[,coldata$ori_id]
tpm_T <- tpm[,colnames(tpm) %in% coldata[coldata$condition == "Tumor","ori_id"]]
tpm_T <- tpm_T[,sort(colnames(tpm_T))]
dim(tpm_T)
tpm_T <- round(as.matrix(tpm_T))

data_T <- data[,colnames(data) %in% coldata[coldata$condition == "Tumor","ori_id"]]
dim(data_T)
coldata_T <- coldata[coldata$condition == "Tumor",]
coldata_T$Sample_ID <- as.integer(coldata_T$Sample_ID)
dim(coldata_T)

## add age group
coldata_T$age_group <- ifelse(coldata_T$Age >= 60,"OV60","BL60")
coldata_T_ov60 <- coldata_T[coldata_T$age_group == "OV60",]
coldata_T_bl60 <- coldata_T[coldata_T$age_group == "BL60",]

data_T_ov60 <- data_T[,colnames(data_T) %in% coldata_T_ov60[coldata_T_ov60$condition == "Tumor","ori_id"]]
dim(data_T_ov60)
data_T_bl60 <- data_T[,colnames(data_T) %in% coldata_T_bl60[coldata_T_bl60$condition == "Tumor","ori_id"]]
dim(data_T_bl60) 

data_ov60 <- data[,colnames(data) %in% coldata[coldata$Age >=60,"ori_id"]]
data_bl60 <- data[,colnames(data) %in% coldata[coldata$Age <60,"ori_id"]]

coldata_ov60 <- coldata[coldata$Age >=60,]
coldata_bl60 <- coldata[coldata$Age < 60,]
############## 5. DEseq OV60 ########
dds_gender_condition <- DESeqDataSetFromMatrix(countData = data_T_ov60,
                                               colData = coldata_T_ov60,
                                               design = ~ Gender )
dds.temp_ <- dds_gender_condition
nrow(dds.temp_)
summary(rowMedians(counts(dds_gender_condition)))
summary(rowSums(counts(dds_gender_condition)))
dds <- dds_gender_condition[rowMedians(counts(dds_gender_condition)) > 403 ,]
dds.temp <- DESeq(dds)
res <- results(dds.temp)
res <- res[order(res$padj),]
res = data.frame(res)
res <- na.omit(res)
res$fdr = p.adjust(res$pvalue,method="BH")

res = res[order(res$padj), ]
dim(res)
summary(abs(res$log2FoldChange))
logFC_cutoff <- with(res,mean(abs(log2FoldChange)) + 1*sd(abs(log2FoldChange)))
diffMat = res[res$pvalue<=0.05 & abs(res$log2FoldChange) >= 0.2,]
diffMat <- na.omit(diffMat)
dim(diffMat)

write.table(res,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/1. DEGs/BaiduII_FvsM_Degs.xls",sep = "\t",quote = F)
write.table(diffMat,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/1. DEGs/BaiduII_FvsM_Sig_Degs.xls",sep = "\t",quote = F)
diffMat <- read.table("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/1. DEGs/BaiduII_FvsM_Sig_Degs.xls")
dim(diffMat)

############## 6. DEseq BL60 ########
dds_gender_condition <- DESeqDataSetFromMatrix(countData = data_T_bl60,
                                               colData = coldata_T_bl60,
                                               design = ~ Gender )
dds.temp_ <- dds_gender_condition
nrow(dds.temp_)
summary(rowMedians(counts(dds_gender_condition)))
summary(rowSums(counts(dds_gender_condition)))

dds <- dds_gender_condition[rowMedians(counts(dds_gender_condition)) > 380.8 ,]
dds.temp <- DESeq(dds)
res <- results(dds.temp)
res <- res[order(res$padj),]
res = data.frame(res)
res <- na.omit(res)
res$fdr = p.adjust(res$pvalue,method="BH")

res = res[order(res$padj), ]
dim(res)
summary(abs(res$log2FoldChange))

logFC_cutoff <- with(res,mean(abs(log2FoldChange)) + 1*sd(abs(log2FoldChange)))
diffMat = res[res$pvalue<=0.05 & abs(res$log2FoldChange) >= 1,]
diffMat <- na.omit(diffMat)
dim(diffMat)

write.table(res,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_BL60/1. DEGs/BaiduII_FvsM_Degs.xls",sep = "\t",quote = F)
write.table(diffMat,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_BL60/1. DEGs/BaiduII_FvsM_Sig_Degs.xls",sep = "\t",quote = F)

dim(diffMat)
############## 7. Heatmap ####################
### ov60 ####
coldata_ov60$TNM.stage.the.Eighth.Edition. <- ifelse(coldata_ov60$TNM.stage.the.Eighth.Edition. %in% c("IIIA","IIIB"),3,
                                                ifelse(coldata_ov60$TNM.stage.the.Eighth.Edition. %in% c("IIA","ⅡA","IIB"),2,
                                                       ifelse(coldata_ov60$TNM.stage.the.Eighth.Edition. %in% c("IA","IB"),1,
                                                              ifelse(coldata_ov60$TNM.stage.the.Eighth.Edition. == "ⅣA",4,"stop"))))

color_TNM <- brewer.pal(4,'YlOrRd')
names(color_TNM) <- unique(coldata_ov60$TNM.stage.the.Eighth.Edition.)
color_grade <- brewer.pal(3,'YlGnBu')
names(color_grade) <- c("G1","G2","G3")
color_smoke <- brewer.pal(4,'RdPu')
names(color_smoke) <- unique(coldata_ov60$Smoking.history)

col_ha = HeatmapAnnotation(
  df=coldata_ov60[,c(
    "Gender" ,                
    "Smoking.history"       ,
    "TNM.stage.the.Eighth.Edition.", 
    "Grade" ,                       
    "condition" )],
  col =  list(
    condition = c("Tumor" =  pal_jco()(5)[2], "Normal" =  pal_jco()(5)[5]),
    TNM.stage.the.Eighth.Edition. = color_TNM,
    Gender = c("Female"= pal_jco()(5)[4],"Male" =  pal_jco()(5)[1]),
    Grade = color_grade , Smoking.history = color_smoke),
  annotation_legend_param = list(direction = "horizontal"))

diffmat_genes <- HGNC[HGNC$Ensembl.gene.ID %in% rownames(diffMat),]
diffmat_genes <- na.omit(diffmat_genes)
diffmat_genes <- diffmat_genes[!duplicated(diffmat_genes$Ensembl.gene.ID),]
dim(diffmat_genes)

int <- intersect(rownames(data_T_ov60), diffmat_genes$Ensembl.gene.ID)
plot.mat = scale_mat(data_ov60[int,],"row")
dim(plot.mat)
length(diffmat_genes[diffmat_genes$Ensembl.gene.ID %in% int,"Approved.symbol"])
rownames(plot.mat) <- diffmat_genes[diffmat_genes$Ensembl.gene.ID %in% int,"Approved.symbol"]

p3 = Heatmap(plot.mat,cluster_rows = T,cluster_columns = T,row_dend_side = "right",row_names_side = "left",
             show_row_names = F,top_annotation = col_ha,show_column_names = F,
             heatmap_legend_param = list(title = "expression levels",legend_direction = "horizontal"),
             #row_split = vocano[row.names(plot.mat),"change"],
             column_split = coldata_ov60[,c("condition","Gender")],
             row_names_gp = gpar(fontsize = 10)   ,                  
             col = colorRamp2(c(-2,0,2), c("purple", "black", "yellow")))

pdf("results/heatmap_ov60_all.pdf",w =12 ,h=10)
draw(p3)
dev.off()
### bl60 #####
coldata_bl60$TNM.stage.the.Eighth.Edition. <- ifelse(coldata_bl60$TNM.stage.the.Eighth.Edition. %in% c("IIIA","IIIB"),3,
                                                       ifelse(coldata_bl60$TNM.stage.the.Eighth.Edition. %in% c("IIA","ⅡA","IIB"),2,
                                                              ifelse(coldata_bl60$TNM.stage.the.Eighth.Edition. %in% c("IA","IB"),1,
                                                                     ifelse(coldata_bl60$TNM.stage.the.Eighth.Edition. == "ⅣA",4,"stop"))))

color_TNM <- brewer.pal(4,'YlOrRd')
names(color_TNM) <- unique(coldata_bl60$TNM.stage.the.Eighth.Edition.)
color_grade <- brewer.pal(3,'YlGnBu')
names(color_grade) <- c("G1","G2","G3")
color_smoke <- brewer.pal(4,'RdPu')
names(color_smoke) <- unique(coldata_bl60$Smoking.history)

# get plot matrix
diffmat_genes <- HGNC[HGNC$Ensembl.gene.ID %in% rownames(diffMat),]
diffmat_genes <- na.omit(diffmat_genes)
diffmat_genes <- diffmat_genes[!duplicated(diffmat_genes$Ensembl.gene.ID),]
dim(diffmat_genes)

int <- intersect(rownames(data_T_bl60), diffmat_genes$Ensembl.gene.ID)
plot.mat = scale_mat(data_bl60[int,],"row")
dim(plot.mat)
length(diffmat_genes[diffmat_genes$Ensembl.gene.ID %in% int,"Approved.symbol"])
rownames(plot.mat) <- diffmat_genes[diffmat_genes$Ensembl.gene.ID %in% int,"Approved.symbol"]

col_ha = HeatmapAnnotation(
  df=coldata_bl60[,c(
    "Gender" ,                
    "Smoking.history"       ,
    "TNM.stage.the.Eighth.Edition.", 
    "Grade" ,                       
    "condition" )],
  col =  list(
    condition = c("Tumor" =  pal_jco()(5)[2], "Normal" =  pal_jco()(5)[5]),
    TNM.stage.the.Eighth.Edition. = color_TNM,
    Gender = c("Female"= pal_jco()(5)[4],"Male" =  pal_jco()(5)[1]),
    Grade = color_grade , Smoking.history = color_smoke),
  annotation_legend_param = list(direction = "horizontal"))


p3 = Heatmap(plot.mat,cluster_rows = T,cluster_columns = T,row_dend_side = "right",row_names_side = "left",
             show_row_names = F,top_annotation = col_ha,show_column_names = F,
             heatmap_legend_param = list(title = "expression levels",legend_direction = "horizontal"),
             #row_split = vocano[row.names(plot.mat),"change"],
             column_split = coldata_bl60[,c("condition","Gender")],
             row_names_gp = gpar(fontsize = 10)   ,                  
             col = colorRamp2(c(-2,0,2), c("purple", "black", "yellow")))

pdf("results/heatmap_bl60_all.pdf",w =12 ,h=10)
draw(p3)
dev.off()


######## 8. Reactome ########
## OV60
diffmat_ov60 <- diffMat
diffmat_genes <- HGNC[HGNC$Ensembl.gene.ID %in% rownames(diffmat_ov60),]
diffmat_genes <- na.omit(diffmat_genes)
diffmat_genes_ov60 <- diffmat_genes[!duplicated(diffmat_genes$Ensembl.gene.ID),]
dim(diffmat_genes_ov60)
diffmat_ov60 <- diffmat_ov60[diffmat_genes_ov60$Ensembl.gene.ID,]
rownames(diffmat_ov60) <- diffmat_genes_ov60$Approved.symbol
reactome.ov60 <- enrichPathway(gene       = diffmat_genes$NCBI.gene.ID ,
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 1)
## OV60
diffmat_bl60 <- diffMat
diffmat_genes <- HGNC[HGNC$Ensembl.gene.ID %in% rownames(diffmat_bl60),]
diffmat_genes <- na.omit(diffmat_genes)
diffmat_genes_bl60 <- diffmat_genes[!duplicated(diffmat_genes$Ensembl.gene.ID),]
dim(diffmat_genes_bl60)
reactome.bl60 <- enrichPathway(gene       = diffmat_genes_bl60$NCBI.gene.ID ,
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 1)
reactome.ov60@result <- reactome.ov60@result[reactome.ov60@result$Description != "HDR through Homologous Recombination (HRR) or Single Strand Annealing (SSA)",]
######
edox_ov60 <- setReadable(reactome.ov60, 'org.Hs.eg.db', 'ENTREZID')
FC <- setNames(diffmat_ov60$log2FoldChange,rownames(diffmat_ov60))
pdf("results/reactome_ov60_cnetplot.pdf",width = 8, height = 6)
cnetplot(edox_ov60, foldChange = FC,node_label = "all",layout = "gem") + scale_color_viridis(option="C")
dev.off()

pdf("results/reactome_ov60_upsetplot.pdf")
upsetplot(reactome.ov60)
dev.off()

###### save gene set
length(diffmat_genes[diffmat_genes$Ensembl.gene.ID %in% int,"Approved.symbol"])
rownames(diffmat_new ) <- diffmat_genes[diffmat_genes$Ensembl.gene.ID %in% int,"Approved.symbol"]
up_genes <- rownames(diffmat_new)[diffmat_new$log2FoldChange>= 0]
down_genes <- rownames(diffmat_new)[diffmat_new$log2FoldChange<0]
write.table(up_genes,"results/baiduII_ov60_new_rna_deg_up.txt",row.names = F, quote = F)
write.table(down_genes,"results/baiduII_ov60_new_rna_deg_down.txt",row.names = F, quote = F)


