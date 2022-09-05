############## 0.library and function #########
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(annotatr)
library(ggsci)
library(ggrepel)
library(methylKit)
library(dplyr)
library(GenomicRanges)
library(plyranges)
library(IRanges)
library("ggpubr")
library("qpcR") 
library(RColorBrewer)
library(ComplexHeatmap)
library(colorRamp)
library(circlize)

############## 1.data import ####
#### 1.1 clinical data #####
metadata <- read.csv(file = "~/onedrive/Work/phD/phd_project/TME_gender/rawdata/baiduII_clinical_data.csv", header = TRUE)
metadata <- metadata[metadata$Age >= 60,]
metadata <- metadata[metadata$sample.ID............BDESCC2..!= 178,]
anno_data <- metadata[,c("sample.ID............BDESCC2..", "Gender",
                         "Smoking.history","TNM.stage.the.Eighth.Edition.",
                         "Grade")]
colnames(anno_data)[1] <- "Sample_ID"

metadata_gender <- metadata[,c("sample.ID............BDESCC2..", "Gender")]
colnames(metadata_gender)<- c("Sample_ID", "Gender")
metadata_gender = metadata_gender[order(metadata_gender$Sample_ID),]
dim(metadata_gender)
#### 1.2 female #####
female.file.list=list.files("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/10.Methylation/paired_data/female")
female_file_mydiff_final_list <- NA
setwd("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/10.Methylation/paired_data/female")
# get myDiff for each pair and combine to one
for (file in female.file.list){
  patient_id <- sapply(strsplit(sapply(strsplit(file,"new_female"),"[[",2),"TvsN"),"[[",1)
  if (patient_id %in% metadata_gender$Sample_ID){
    female_DMR_myDiff <- readRDS(file)
    female_DMR_myDiff=getMethylDiff(female_DMR_myDiff,difference=0,qvalue=0.05)
    female_DMR_myDiff <- as(female_DMR_myDiff,"GRanges")
    
    female_DMR_myDiff <- as.data.frame(female_DMR_myDiff)
    female_DMR_myDiff <- female_DMR_myDiff[,c(1,2,3,8)]
    colnames(female_DMR_myDiff)[4] <- patient_id
    if(length(female_file_mydiff_final_list) == 1){
      female_file_mydiff_final_list <- as.data.frame(female_DMR_myDiff)
    } else {
      female_file_mydiff_final_list <- merge(female_file_mydiff_final_list,as.data.frame(female_DMR_myDiff),
                                             by.x = c("seqnames",  "start",     "end"),by.y = c("seqnames",  "start",     "end"),all.x=TRUE,all.y=TRUE)
    }}
}
female_file_mydiff_final_list <- female_file_mydiff_final_list[female_file_mydiff_final_list$seqnames %in% c(paste0("chr",rep(1:22,1)),"chrX","chrY"),]
colnames(female_file_mydiff_final_list)[4:35] <- gsub("X","",colnames(female_file_mydiff_final_list[4:35]))
nrow(female_file_mydiff_final_list)
## Filter with mean > 0 and #NA < 5
female_file_mydiff_final_list_ <- female_file_mydiff_final_list[rowSums(is.na(female_file_mydiff_final_list)) <= 5,]
female_file_mydiff_final_list_ <- female_file_mydiff_final_list_[abs(rowMeans(female_file_mydiff_final_list_[,4:23],na.rm=T)) >= 0,]
nrow(female_file_mydiff_final_list_)
# 182511
saveRDS(female_file_mydiff_final_list,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/10.Methylation/female_intersect_mydiff_final_list.rds")
female_file_mydiff_final_list <- readRDS("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/10.Methylation/female_intersect_mydiff_final_list.rds")
nrow(female_file_mydiff_final_list_)
## convert back to GRanges
female_mydiff_GR <- as(female_file_mydiff_final_list_,"GRanges")
female_mydiff_df <- female_file_mydiff_final_list_
####  1.3 male  ######
male.file.list=list.files("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/10.Methylation/paired_data/male")
male_file_mydiff_final_list <- NA
setwd("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/10.Methylation/paired_data/male")
for (file in male.file.list){
  patient_id <- sapply(strsplit(sapply(strsplit(file,"new_male"),"[[",2),"TvsN"),"[[",1)
  if (patient_id %in% metadata_gender$Sample_ID){
    male_DMR_myDiff <- readRDS(file)
    male_DMR_myDiff=getMethylDiff(male_DMR_myDiff,difference=0,qvalue=0.05)
    male_DMR_myDiff <- as(male_DMR_myDiff,"GRanges")
    
    male_DMR_myDiff <- as.data.frame(male_DMR_myDiff)
    male_DMR_myDiff <- male_DMR_myDiff[,c(1,2,3,8)]
    colnames(male_DMR_myDiff)[4] <- patient_id
    if(length(male_file_mydiff_final_list) == 1){
      male_file_mydiff_final_list <- as.data.frame(male_DMR_myDiff)
    } else {
      male_file_mydiff_final_list <- merge(male_file_mydiff_final_list,as.data.frame(male_DMR_myDiff),
                                           by.x = c("seqnames",  "start",     "end"),by.y = c("seqnames",  "start",     "end"),all.x=TRUE,all.y=TRUE)
    }}
}
nrow(male_file_mydiff_final_list)
male_file_mydiff_final_list <- male_file_mydiff_final_list[male_file_mydiff_final_list$seqnames %in% c(paste0("chr",rep(1:22,1)),"chrX","chrY"),]
colnames(male_file_mydiff_final_list)[4:69] <- gsub("X","",colnames(male_file_mydiff_final_list[4:69]))
## Filter with mean > 0 and #NA < 18
male_file_mydiff_final_list_ <- male_file_mydiff_final_list[rowSums(is.na(male_file_mydiff_final_list)) <= 9,]
nrow(male_file_mydiff_final_list_)
# 127625
saveRDS(male_file_mydiff_final_list,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/10.Methylation/male_intersect_mydiff_final_list.rds")
#male_file_mydiff_final_list_ <- readRDS("~/onedrive/Work/phD/phd_project/TME_gender/results/baiduII_OV60/10. Methylation/male_intersect_mydiff_final_list.rds")
nrow(male_file_mydiff_final_list_)
## convert back to GRanges
male_mydiff_GR <- as(male_file_mydiff_final_list_,"GRanges")
male_mydiff_df <- male_file_mydiff_final_list_
####  1.4 fvsm #####
fvsm_myDiff <- readRDS("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/10.Methylation/all_99_DMR_treatment_gender_myDiff.rds")
fvsm_myDiff=getMethylDiff(fvsm_myDiff ,difference=0,qvalue=0.05)
fvsm_myDiff_gr <- as(fvsm_myDiff,"GRanges")
fvsm_myDiff_df <- as.data.frame(fvsm_myDiff_gr)


##############  2. Compare TvsN DMR in genders ####
#### 2.1 intersection FandM gr vs FvsM DMR #####
## get intersection of TvsN
int_female_gr <-  subsetByOverlaps(female_mydiff_GR,male_mydiff_GR)
int_male_gr <-  subsetByOverlaps(male_mydiff_GR,female_mydiff_GR)
mcols(int_female_gr) <- cbind(mcols(int_female_gr),mcols(int_male_gr))
int_tvsn_gr <- int_female_gr
#### 2.2 only female #####
only_f_gr <- GenomicRanges::setdiff(female_mydiff_GR,int_female_gr)
only_f_gr <- subsetByOverlaps(female_mydiff_GR,only_f_gr)

only_f_tvsn_gr <- subsetByOverlaps(only_f_gr,fvsm_myDiff_gr)
only_f_fvsm_gr <- subsetByOverlaps(fvsm_myDiff_gr,only_f_gr)

mcols(only_f_tvsn_gr) <- cbind(mcols(only_f_tvsn_gr),mcols(only_f_fvsm_gr)[,"meth.diff"])
only_f_gr <- only_f_tvsn_gr
colnames(mcols(only_f_gr))[32] <- "meth.diff.gender"

only_f_gr <- readRDS("~/onedrive/Work/phD/phd_project/TME_gender/results/baiduII_OV60/10. Methylation/only_f_gr.rds")

#### 2.3 only male #####
only_m_gr <- GenomicRanges::setdiff(male_mydiff_GR,int_male_gr)
only_m_gr <- subsetByOverlaps(male_mydiff_GR,only_m_gr)

## get gr for FvsM methylation
only_m_tvsn_gr <- subsetByOverlaps(only_m_gr,fvsm_myDiff_gr)
only_m_fvsm_gr <- subsetByOverlaps(fvsm_myDiff_gr,only_m_gr)

mcols(only_m_tvsn_gr) <- cbind(mcols(only_m_tvsn_gr),mcols(only_m_fvsm_gr)[,"meth.diff"])
only_m_gr <- only_m_tvsn_gr
colnames(mcols(only_m_gr))[66] <- "meth.diff.gender"
saveRDS(only_m_gr,"~/onedrive/Work/phD/phd_project/TME_gender/results/baiduII_OV60/10. Methylation/only_m_gr.rds")

############## 3. Annotation ####################
### 3.1 Prep #######
# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
           'hg19_genes_intronexonboundaries'
)
# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg19', annotations = annots)

#### 3.2 Female only #####
only_f_gr <- only_f_gr[only_f_gr@seqnames %in% c(paste0("chr", rep(1:22)),"chrX", "chrY"),]
only_f_gr_annotated = annotate_regions(
  regions = only_f_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
#### 3.3  Male only #######
only_m_gr <- only_m_gr[only_m_gr@seqnames %in% c(paste0("chr", rep(1:22)),"chrX", "chrY"),]
only_m_gr_annotated = annotate_regions(
  regions = only_m_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

# Coerce to a data.frame

only_f_annotated_df = data.frame(only_f_gr_annotated)
only_m_annotated_df = data.frame(only_m_gr_annotated )

############## 4. Gender meth diff vs chr vs gene plot #####
###### 4.1 Label gender type #######
only_f_annotated_df$type <- "Female only"
#fm_int_sig_annotated_df$type <- "Intersection"
only_m_annotated_df$type <- "Male only"

df_meth_diff_gene <- rbind(only_f_annotated_df[,c("annot.seqnames","type","meth.diff.gender","annot.symbol")],
                           only_m_annotated_df[,c("annot.seqnames","type","meth.diff.gender","annot.symbol")])
#fm_int_sig_annotated_df[,c("annot.seqnames","type","meth.diff.gender","annot.symbol")])
df_meth_diff_gene <- na.omit(df_meth_diff_gene)


######  check promoter average tvsn meth diff for each gene ####
#### OV60####
only_f_promoter <- only_f_annotated_df[only_f_annotated_df$seqnames %in% c(paste0("chr",rep(1:22,1)),"chrX","chrY"),]
only_f_promoter <- only_f_promoter[only_f_promoter$annot.type == 'hg19_genes_promoters', ]
unique(only_f_promoter$annot.symbol)

only_m_promoter <- only_m_annotated_df[only_m_annotated_df$seqnames %in% c(paste0("chr",rep(1:22,1)),"chrX","chrY"),]
only_m_promoter <- only_m_promoter[only_m_promoter$annot.type == 'hg19_genes_promoters', ]
unique(only_m_promoter$annot.symbol)

only_f_promoter$meth.diff.tn.mean <- rowMeans(only_f_promoter[,c(6:37)],na.rm=T)
only_m_promoter$meth.diff.tn.mean <- rowMeans(only_m_promoter[,c(6:71)],na.rm=T)

only_f_promoter$type <- "Female only"
only_m_promoter$type <- "Male only"

df_meth_diff_gene_promoter_df <- rbind(only_f_promoter[,c("annot.seqnames","type","meth.diff.gender","annot.symbol","meth.diff.tn.mean")],
                                       only_m_promoter[,c("annot.seqnames","type","meth.diff.gender","annot.symbol","meth.diff.tn.mean")])
df_meth_diff_gene_promoter_df <- na.omit(df_meth_diff_gene_promoter_df)

df_meth_diff_gene_promoter_df_ <- df_meth_diff_gene_promoter_df %>% 
  group_by(annot.symbol) %>% 
  summarise(type = type,chr = annot.seqnames,avg_gender = mean(meth.diff.gender),avg_tn = mean(meth.diff.tn.mean))
df_meth_diff_gene_promoter_df_$abs_meth_diff <- abs(df_meth_diff_gene_promoter_df_$avg_gender)
df_meth_diff_gene_promoter_df_ <- distinct(df_meth_diff_gene_promoter_df_)
df_meth_diff_gene_promoter_df_ <- df_meth_diff_gene_promoter_df_[order(df_meth_diff_gene_promoter_df_$abs_meth_diff,decreasing = T ),]
nrow(df_meth_diff_gene_promoter_df_)
## 2525
## labels
df_meth_diff_gene_promoter_df_ <- as.data.frame(df_meth_diff_gene_promoter_df_ )
labeled_genes <- df_meth_diff_gene_promoter_df_[abs(df_meth_diff_gene_promoter_df_$avg_gender) > 50,]
p_test <- ggplot(df_meth_diff_gene_promoter_df_) +
  geom_point(aes(x = `avg_tn`, y = `avg_gender`,color = chr)) + labs(x = "Methylation difference between tumor and normal", y = "Methylation difference between gender")
p17 <- p_test + geom_label_repel(data = labeled_genes, max.overlaps = 90,
                                 aes(x = `avg_tn`, y = `avg_gender`,label = `annot.symbol`),
                                 size = 4) + facet_grid(~type)
pdf("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/10.Methylation/methylation_plot_methdiff_tvsn_fvsm_promoter_version_3_no178.pdf",width = 10,height = 5)
p17
dev.off()

### check import genes
df_meth_diff_gene_promoter_df_[df_meth_diff_gene_promoter_df_$annot.symbol == "CDK1",]

