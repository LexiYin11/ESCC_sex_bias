############# 0. Library ##################
library(circlize)
library(maftools)
library(viridis)
library(ggsci)
library(GenomicRanges)
library(splicejam)
library(ComplexHeatmap)

cytoband_bed <- read.cytoband(species = "hg19")$df

#exp_set <- read.table("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/gene_lists/expression_geneset_106.txt")
#exp_set <- exp_set$V1
############# 1. Data filter and export to gistic  ##############
cnv_file <- read.delim("~/onedrive/Work/phD/phd_project/TME_gender/rawdata/ESCC_Biseq_CNV_seq.txt")
cnv_file  <- cnv_file[grep("T",cnv_file$ID),]
cnv_file$ID <- gsub("X","",cnv_file$ID)
cnv_file$ID <- gsub("T","",cnv_file$ID)
cnv_file$cn <- 2^(cnv_file$seg.mean)*2
cnv_file$ID <- as.numeric(cnv_file$ID)

## clinical data
metadata <- read.csv(file = "~/onedrive/Work/phD/phd_project/TME_gender/rawdata/BaiduII_clinical_data.csv", header = TRUE)
metadata <- metadata[metadata$Age >= 60,]
metadata_gender <- metadata[,c("sample.ID............BDESCC2..", "Gender")]
colnames(metadata_gender)<- c("Sample_ID", "Gender")
metadata_gender = metadata_gender[order(metadata_gender$Sample_ID),]

cnv_file <- as.data.frame(cnv_file[cnv_file$ID %in% metadata_gender$Sample_ID,])
## female and male data
metadata_female <- metadata_gender[metadata_gender$Gender == 'Female',]
metadata_male <- metadata_gender[metadata_gender$Gender == 'Male',]

cnv_female <- as.data.frame(cnv_file[cnv_file$ID %in% metadata_female$Sample_ID,])
cnv_male <- as.data.frame(cnv_file[cnv_file$ID %in% metadata_male$Sample_ID,])

write.table(cnv_file,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/ESCC_Biseq_CNV_seq_ov60_all.txt", quote = F,row.names = F,sep = "\t")
write.table(cnv_male,"results/BaiduII_OV60/7. CNV/ESCC_Biseq_CNV_seq_ov60_male.txt", quote = F,row.names = F,sep = "\t")
write.table(cnv_female,"results/BaiduII_OV60/7. CNV/ESCC_Biseq_CNV_seq_ov60_female.txt", quote = F,row.names = F,sep = "\t")
############# 2. Import data from GISTIC2  ###############
#### gistic score total and female and male 
scores_total <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/all/scores.gistic")
scores_female <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/female/scores.gistic")
scores_male <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/male/scores.gistic")

#### amp and del resgions female and male 
regions_male <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/female/regions_track.conf_99.bed")
regions_female <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/male/regions_track.conf_99.bed")

regions_female$V4 <- sapply(strsplit(regions_female$V4,"-"),"[[",2)
regions_female$V4 <- sapply(strsplit(regions_female$V4,"P"),"[[",1)

regions_male$V4 <- sapply(strsplit(regions_male$V4,"-"),"[[",2)
regions_male$V4 <- sapply(strsplit(regions_male$V4,"P"),"[[",1)

#### total all thresholded by genes 
ori_total_cnv_thresholded_genes <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/all/all_thresholded.by_genes.txt")
total_cnv_thresholded_genes <- ori_total_cnv_thresholded_genes[,4:102]
total_cnv_thresholded_genes <- as.data.frame(total_cnv_thresholded_genes)
rownames(total_cnv_thresholded_genes) <- ori_total_cnv_thresholded_genes$`Gene Symbol`
## female
ori_female_thresholded_gene <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/female/all_thresholded.by_genes.txt")
female_thresholded_gene <- ori_female_thresholded_gene[,4:35]
female_thresholded_gene <- as.data.frame(female_thresholded_gene)
rownames(female_thresholded_gene) <- ori_female_thresholded_gene$`Gene Symbol`

## male
ori_male_thresholded_gene <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/male/all_thresholded.by_genes.txt")
male_thresholded_gene <- ori_male_thresholded_gene[,4:70]
male_thresholded_gene <- as.data.frame(male_thresholded_gene)
rownames(male_thresholded_gene) <- ori_male_thresholded_gene$`Gene Symbol`



############## 3. Create female circos ################
## Split Amp scores
scores_female_amp <-  scores_female[scores_female$Type == "Amp",] 
scores_female_amp$Chromosome <- paste0("chr",scores_female_amp$Chromosome)
scores_female_amp <- scores_female_amp[,c(2,3,4,6)]
## Split Del scores
scores_female_del <-  scores_female[scores_female$Type == "Del",] 
scores_female_del$Chromosome <- paste0("chr",scores_female_del$Chromosome)
scores_female_del <- scores_female_del[,c(2,3,4,6)]

## Split Amp track
regions_female_amp <- as.data.frame(regions_female[regions_female$V4 =="A",])
## Split Del track
regions_female_del <- as.data.frame(regions_female[regions_female$V4 =="D",])

## create amp score Grange object
female_amp_grange <- GRanges(seqnames = scores_female_amp$Chromosome, ranges = IRanges(scores_female_amp$Start,scores_female_amp$End))
values(female_amp_grange) <- scores_female_amp$`G-score`
colnames(values(female_amp_grange)) <- "G-score"

## create del score Grange object
female_del_grange <- GRanges(seqnames = scores_female_del$Chromosome, ranges = IRanges(scores_female_del$Start,scores_female_del$End))
values(female_del_grange) <- scores_female_del$`G-score`
colnames(values(female_del_grange)) <- "G-score"

## create amp track Grange object
regions_female_amp_grange <- GRanges(seqnames = regions_female_amp[,1], ranges = IRanges(as.numeric(regions_female_amp[,2]),as.numeric(regions_female_amp[,3])))
values(regions_female_amp_grange) <- regions_female_amp[,4]
colnames(values(regions_female_amp_grange)) <- "Region"

## create del track Grange object
regions_female_del_grange <- GRanges(seqnames = regions_female_del[,1], ranges = IRanges(as.numeric(regions_female_del[,2]),as.numeric(regions_female_del[,3])))
values(regions_female_del_grange) <- regions_female_del[,4]
colnames(values(regions_female_del_grange)) <- "Region"

############## 4. Create male circos ################
## Split Amp scores
scores_male_amp <-  scores_male[scores_male$Type == "Amp",] 
scores_male_amp$Chromosome <- paste0("chr",scores_male_amp$Chromosome)
scores_male_amp <- scores_male_amp[,c(2,3,4,6)]
## Split Del scores
scores_male_del <-  scores_male[scores_male$Type == "Del",] 
scores_male_del$Chromosome <- paste0("chr",scores_male_del$Chromosome)
scores_male_del <- scores_male_del[,c(2,3,4,6)]

## Split Amp track
regions_male_amp <- as.data.frame(regions_male[regions_male$V4 =="A",])
## Split Del track
regions_male_del <- as.data.frame(regions_male[regions_male$V4 =="D",])

## create amp score Grange object
male_amp_grange <- GRanges(seqnames = scores_male_amp$Chromosome, ranges = IRanges(scores_male_amp$Start,scores_male_amp$End))
values(male_amp_grange) <- scores_male_amp$`G-score`
colnames(values(male_amp_grange)) <- "G-score"

## create del score Grange object
male_del_grange <- GRanges(seqnames = scores_male_del$Chromosome, ranges = IRanges(scores_male_del$Start,scores_male_del$End))
values(male_del_grange) <- scores_male_del$`G-score`
colnames(values(male_del_grange)) <- "G-score"

## create amp track Grange object
regions_male_amp_grange <- GRanges(seqnames = regions_male_amp[,1], ranges = IRanges(as.numeric(regions_male_amp[,2]),as.numeric(regions_male_amp[,3])))
values(regions_male_amp_grange) <- regions_male_amp[,4]
colnames(values(regions_male_amp_grange)) <- "Region"

## create del track Grange object
regions_male_del_grange <- GRanges(seqnames = regions_male_del[,1], ranges = IRanges(as.numeric(regions_male_del[,2]),as.numeric(regions_male_del[,3])))
values(regions_male_del_grange) <- regions_male_del[,4]
colnames(values(regions_male_del_grange)) <- "Region"


############## 5. Extract gender only and intersection regions ########
#### 5.1 Amp ####
## female prepare
ori_female_amp_cytobands <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/female/amp_genes.conf_99.txt")
ori_female_amp_cytobands <- as.data.frame(t(ori_female_amp_cytobands))
female_amp_cytobands <- ori_female_amp_cytobands[,1:5]
colnames(female_amp_cytobands) <- female_amp_cytobands[1,]
female_amp_cytobands <- na.omit(female_amp_cytobands)
female_amp_cytobands <- female_amp_cytobands$cytoband[2:nrow(female_amp_cytobands)]

## male prepare
ori_male_amp_cytobands <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/male/amp_genes.conf_99.txt")
ori_male_amp_cytobands <- as.data.frame(t(ori_male_amp_cytobands))
male_amp_cytobands <- ori_male_amp_cytobands[,1:5]
colnames(male_amp_cytobands) <- male_amp_cytobands[1,]
male_amp_cytobands <- na.omit(male_amp_cytobands)
male_amp_cytobands <- male_amp_cytobands$cytoband[2:nrow(male_amp_cytobands)]

## female only regions
female_only_amp <- setdiff(female_amp_cytobands,male_amp_cytobands)
## Male only regions
male_only_amp <- setdiff(male_amp_cytobands,female_amp_cytobands)
## intersect regions
intersect_amp <- intersect(female_amp_cytobands,male_amp_cytobands)

#### 5.2 Del ####
## female
ori_female_del_cytobands <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/female/del_genes.conf_99.txt")
ori_female_del_cytobands <- as.data.frame(t(ori_female_del_cytobands))
female_del_cytobands <- ori_female_del_cytobands[,1:5]
colnames(female_del_cytobands) <- female_del_cytobands[1,]
female_del_cytobands <- na.omit(female_del_cytobands)
female_del_cytobands <- female_del_cytobands$cytoband[2:nrow(female_del_cytobands)]

## male
ori_male_del_cytobands <- data.table::fread("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/male/del_genes.conf_99.txt")
ori_male_del_cytobands <- as.data.frame(t(ori_male_del_cytobands))
male_del_cytobands <- ori_male_del_cytobands[,1:5]

colnames(male_del_cytobands) <- male_del_cytobands[1,]
male_del_cytobands <- na.omit(male_del_cytobands)
male_del_cytobands <- male_del_cytobands$cytoband[2:nrow(male_del_cytobands)]

## female only regions
female_only_del <- setdiff(female_del_cytobands,male_del_cytobands)
## Male only regions
male_only_del <- setdiff(male_del_cytobands,female_del_cytobands)
## intersect resgions
intersect_del <- intersect(female_del_cytobands,male_del_cytobands)

############## 6. Extract gender only genes ########
##### 6.1 Amp #####
## female 
female_only_amp_genes <- ori_female_amp_cytobands[ori_female_amp_cytobands$V1 %in% female_only_amp, ]
female_only_amp_cytoband_gene_df <-  data.frame(matrix(ncol = 2,nrow = 1))
for (i in 1:nrow(female_only_amp_genes)){
  cur_row <- female_only_amp_genes[i,]
  filter_cur_row <- cur_row[,colSums(is.na(cur_row) | cur_row == "") != nrow(cur_row)]
  ### 先只取三个
  for (gene in filter_cur_row[5:ncol(filter_cur_row)]){
    # if (ncol(filter_cur_row) < 7){
    #   rownum <- ncol(filter_cur_row)
    # } else{
    #   rownum <- 7
    # }
    # 
    # for (gene in filter_cur_row[5:rownum]){
    female_only_amp_cytoband_gene_df <- rbind(female_only_amp_cytoband_gene_df,c(filter_cur_row$V1[1],gene))
  }
}
female_only_amp_cytoband_gene_df <- female_only_amp_cytoband_gene_df[2:nrow(female_only_amp_cytoband_gene_df),]
colnames(female_only_amp_cytoband_gene_df) <- c("cytoband","gene")
dim(female_only_amp_cytoband_gene_df)
saveRDS(female_only_amp_cytoband_gene_df,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/female_only_amp_cytoband_gene_df.rds")
## Male
male_only_amp_genes <- ori_male_amp_cytobands[ori_male_amp_cytobands$V1 %in% male_only_amp, ]
male_only_amp_cytoband_gene_df <-  data.frame(matrix(ncol = 2,nrow = 1))
for (i in 1:nrow(male_only_amp_genes)){
  cur_row <- male_only_amp_genes[i,]
  filter_cur_row <- cur_row[,colSums(is.na(cur_row) | cur_row == "") != nrow(cur_row)]
  # if (ncol(filter_cur_row) < 7){
  #   rownum <- ncol(filter_cur_row)
  # } else{
  #   rownum <- 7
  # }
  # 
  # for (gene in filter_cur_row[5:rownum]){
  for (gene in filter_cur_row[5:ncol(filter_cur_row)]){
    male_only_amp_cytoband_gene_df <- rbind(male_only_amp_cytoband_gene_df,c(filter_cur_row$V1[1],gene))
  }
}
male_only_amp_cytoband_gene_df <- male_only_amp_cytoband_gene_df[2:nrow(male_only_amp_cytoband_gene_df),]
colnames(male_only_amp_cytoband_gene_df) <- c("cytoband","gene")
dim(male_only_amp_cytoband_gene_df)
saveRDS(male_only_amp_cytoband_gene_df,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/male_only_amp_cytoband_gene_df.rds")

##### 6.2 Del ####
## female 
female_only_del_genes <- ori_female_del_cytobands[ori_female_del_cytobands$V1 %in% female_only_del, ]
female_only_del_cytoband_gene_df <-  data.frame(matrix(ncol = 2,nrow = 1))
for (i in 1:nrow(female_only_del_genes)){
  cur_row <- female_only_del_genes[i,]
  filter_cur_row <- cur_row[,colSums(is.na(cur_row) | cur_row == "") != nrow(cur_row)]
  # if (ncol(filter_cur_row) < 7){
  #   rownum <- ncol(filter_cur_row)
  # } else{
  #   rownum <- 7
  # }
  # 
  # for (gene in filter_cur_row[5:rownum]){
  for (gene in filter_cur_row[5:ncol(filter_cur_row)]){
    female_only_del_cytoband_gene_df <- rbind(female_only_del_cytoband_gene_df,c(filter_cur_row$V1[1],gene))
  }
}
female_only_del_cytoband_gene_df <- female_only_del_cytoband_gene_df[2:nrow(female_only_del_cytoband_gene_df),]
colnames(female_only_del_cytoband_gene_df) <- c("cytoband","gene")
dim(female_only_del_cytoband_gene_df)
saveRDS(female_only_del_cytoband_gene_df,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/female_only_del_cytoband_gene_df.rds")

## Male
male_only_del_genes <- ori_male_del_cytobands[ori_male_del_cytobands$V1 %in% male_only_del, ]
male_only_del_cytoband_gene_df <-  data.frame(matrix(ncol = 2,nrow = 1))
for (i in 1:nrow(male_only_del_genes)){
  cur_row <- male_only_del_genes[i,]
  filter_cur_row <- cur_row[,colSums(is.na(cur_row) | cur_row == "") != nrow(cur_row)]
  # if (ncol(filter_cur_row) < 7){
  #   rownum <- ncol(filter_cur_row)
  # } else{
  #   rownum <- 7
  # }
  # 
  # for (gene in filter_cur_row[5:rownum]){
  for (gene in filter_cur_row[5:ncol(filter_cur_row)]){
    male_only_del_cytoband_gene_df <- rbind(male_only_del_cytoband_gene_df,c(filter_cur_row$V1[1],gene))
  }
}
male_only_del_cytoband_gene_df <- male_only_del_cytoband_gene_df[2:nrow(male_only_del_cytoband_gene_df),]
colnames(male_only_del_cytoband_gene_df) <- c("cytoband","gene")

dim(male_only_del_cytoband_gene_df)
saveRDS(male_only_del_cytoband_gene_df,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/male_only_del_cytoband_gene_df.rds")

############## 7. Extract Significant gene in intersect regions ###########
##### 7.1 Amp ####
## use annotation from female
intersection_amp_genes <- ori_female_amp_cytobands[ori_female_amp_cytobands$V1 %in% intersect_amp, ]
intersection_amp_cytoband_gene_df <-  data.frame(matrix(ncol = 2,nrow = 1))
for (i in 1:nrow(intersection_amp_genes)){
  cur_row <- intersection_amp_genes[i,]
  filter_cur_row <- cur_row[,colSums(is.na(cur_row) | cur_row == "") != nrow(cur_row)]
  # if (ncol(filter_cur_row) < 7){
  #   rownum <- ncol(filter_cur_row)
  # } else{
  #   rownum <- 7
  # }
  # 
  # for (gene in filter_cur_row[5:rownum]){
  for (gene in filter_cur_row[5:ncol(filter_cur_row)]){
    intersection_amp_cytoband_gene_df <- rbind(intersection_amp_cytoband_gene_df,c(filter_cur_row$V1[1],gene))
  }
}
intersection_amp_cytoband_gene_df <- intersection_amp_cytoband_gene_df[2:nrow(intersection_amp_cytoband_gene_df),]
colnames(intersection_amp_cytoband_gene_df) <- c("cytoband","gene")
dim(intersection_amp_cytoband_gene_df)

## Split female and male
total_amp_intersection_cnv_thresholded_genes <- as.data.frame(total_cnv_thresholded_genes)[intersection_amp_cytoband_gene_df$gene,]
female_ID <- metadata_gender[metadata_gender$Gender =="Female","Sample_ID"]
female_amp_thresholded_gene <- total_amp_intersection_cnv_thresholded_genes[,as.character(female_ID)]

male_ID <- metadata_gender[metadata_gender$Gender =="Male","Sample_ID"]
male_amp_thresholded_gene <- total_amp_intersection_cnv_thresholded_genes[,as.character(male_ID)]

## female CNV ratio
female_ratio <- data.frame(matrix(ncol = 2,nrow = nrow(female_amp_thresholded_gene)))
rownames(female_ratio) <- rownames(female_amp_thresholded_gene)
female_num <- nrow(metadata_female)
for (i in 1:nrow(female_amp_thresholded_gene)){
  cut_row <- female_amp_thresholded_gene[i,]
  ratio <- c(length(which(cut_row ==0)), length(which(cut_row !=0)))
  female_ratio[i,] <- ratio
}

## male CNV ratio
male_ratio <- data.frame(matrix(ncol = 2,nrow = nrow(male_amp_thresholded_gene)))
rownames(male_ratio) <- rownames(male_amp_thresholded_gene)
male_num <- nrow(metadata_male)
for (i in 1:nrow(male_amp_thresholded_gene)){
  cut_row <- male_amp_thresholded_gene[i,]
  ratio <- c(length(which(cut_row ==0)), length(which(cut_row !=0)))
  male_ratio[i,] <- ratio
}
## calculate chisq
total_chisq_result <- data.frame(matrix(ncol = 2,nrow = nrow(male_ratio)))
rownames(total_chisq_result) <- rownames(male_ratio)
for (i in 1:nrow(male_ratio)){
  cut_matrix <- rbind(female_ratio[i,],male_ratio[i,])
  cut_chisq_result <- fisher.test(cut_matrix)
  total_chisq_result[i,] <- c(cut_chisq_result$statistic,cut_chisq_result$p.value)
}
colnames(total_chisq_result) <- c("statistics", "pvalue")
total_chisq_result <- total_chisq_result[order(total_chisq_result$pvalue),]

amp_total_chisq_result <- total_chisq_result
## only look at the pvalue 
amp_sig_chisq_result <- amp_total_chisq_result[amp_total_chisq_result$pvalue <= 0.05,]
amp_sig_chisq_result$gene <- rownames(amp_sig_chisq_result)

## 
cytoband_genes <- as.data.frame(ori_total_cnv_thresholded_genes)
rownames(cytoband_genes) <- cytoband_genes$`Gene Symbol`
cytoband_genes <- cytoband_genes[amp_sig_chisq_result$gene,]
amp_sig_chisq_result$cytoband <- cytoband_genes$Cytoband
write.table(amp_sig_chisq_result,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/amp_sig_chisq_result.csv",row.names = F,quote = F)
##### 7.2 Del ####
## use annotation from female
intersection_del_genes <- ori_female_del_cytobands[ori_female_del_cytobands$V1 %in% intersect_del, ]
intersection_del_cytoband_gene_df <-  data.frame(matrix(ncol = 2,nrow = 1))
for (i in 1:nrow(intersection_del_genes)){
  cur_row <- intersection_del_genes[i,]
  filter_cur_row <- cur_row[,colSums(is.na(cur_row) | cur_row == "") != nrow(cur_row)]
  
  for (gene in filter_cur_row[5:ncol(filter_cur_row)]){
    intersection_del_cytoband_gene_df <- rbind(intersection_del_cytoband_gene_df,c(filter_cur_row$V1[1],gene))
  }
}
intersection_del_cytoband_gene_df <- intersection_del_cytoband_gene_df[2:nrow(intersection_del_cytoband_gene_df),]
colnames(intersection_del_cytoband_gene_df) <- c("cytoband","gene")
dim(intersection_del_cytoband_gene_df)

## Split female and male
total_del_intersection_cnv_thresholded_genes <- as.data.frame(total_cnv_thresholded_genes)[intersection_del_cytoband_gene_df$gene,]
female_ID <- metadata_gender[metadata_gender$Gender =="Female","Sample_ID"]
female_del_thresholded_gene <- total_del_intersection_cnv_thresholded_genes[,as.character(female_ID)]

male_ID <- metadata_gender[metadata_gender$Gender =="Male","Sample_ID"]
male_del_thresholded_gene <- total_del_intersection_cnv_thresholded_genes[,as.character(male_ID)]

## female CNV ratio
female_ratio <- data.frame(matrix(ncol = 2,nrow = nrow(female_del_thresholded_gene)))
rownames(female_ratio) <- rownames(female_del_thresholded_gene)
female_num <- nrow(metadata_female)
for (i in 1:nrow(female_del_thresholded_gene)){
  cut_row <- female_del_thresholded_gene[i,]
  ratio <- c(length(which(cut_row ==0)), length(which(cut_row !=0)))
  female_ratio[i,] <- ratio
}
#female_ratio <- na.omit(female_ratio)
## male CNV ratio
male_ratio <- data.frame(matrix(ncol = 2,nrow = nrow(male_del_thresholded_gene)))
rownames(male_ratio) <- rownames(male_del_thresholded_gene)
male_num <- nrow(metadata_male)
for (i in 1:nrow(male_del_thresholded_gene)){
  cut_row <- male_del_thresholded_gene[i,]
  ratio <- c(length(which(cut_row ==0)), length(which(cut_row !=0)))
  male_ratio[i,] <- ratio
}

## calculate chisq
del_total_chisq_result <- data.frame(matrix(ncol = 2,nrow = nrow(male_ratio)))
rownames(del_total_chisq_result) <- rownames(male_ratio)

for (i in 1:nrow(male_ratio)){
  cut_matrix <- rbind(female_ratio[i,],male_ratio[i,])
  cut_chisq_result <- fisher.test(cut_matrix)
  del_total_chisq_result[i,] <- c(cut_chisq_result$statistic,cut_chisq_result$p.value)
}
colnames(del_total_chisq_result) <- c("statistics", "pvalue")
del_total_chisq_result <- del_total_chisq_result[order(del_total_chisq_result$pvalue),]

## only look at the pvalue 
del_sig_chisq_result <- del_total_chisq_result[del_total_chisq_result$pvalue <= 0.05,]
del_sig_chisq_result$gene <- rownames(del_sig_chisq_result)

## 
cytoband_genes <- as.data.frame(ori_total_cnv_thresholded_genes)
rownames(cytoband_genes) <- cytoband_genes$`Gene Symbol`
cytoband_genes <- cytoband_genes[del_sig_chisq_result$gene,]
del_sig_chisq_result$cytoband <- cytoband_genes$Cytoband

## no significant

############## 8. Prep for gene annotation ############
#### 8.1 Female only #######
### Amp
amp_locations <- female_only_amp_cytoband_gene_df$cytoband
chr <- paste0("chr",unique(sapply(strsplit(amp_locations,"q"),"[[",1)))
chr <- unique(sapply(strsplit(chr,"p"),"[[",1))

amp_cytoband_bed <- cytoband_bed[cytoband_bed$V1 %in% chr,]
amp_cytoband_bed$new_cytoband <- paste0(sapply(strsplit(amp_cytoband_bed$V1,"chr"),"[[",2),amp_cytoband_bed$V4)
amp_cytoband_bed <- amp_cytoband_bed[amp_cytoband_bed$new_cytoband %in% amp_locations,]

amp_female_bed <- merge(female_only_amp_cytoband_gene_df,amp_cytoband_bed,by.x = "cytoband", by.y = "new_cytoband")
amp_female_bed <- amp_female_bed [,c(2,3,4,5)]
colnames(amp_female_bed) <- c("gene", "chr","start","end")
amp_female_bed <- amp_female_bed[,c("chr","start","end","gene")]
amp_female_bed$chr <- paste0("Amp_",amp_female_bed$chr)
### Del
del_locations <- female_only_del_cytoband_gene_df$cytoband
chr <- paste0("chr",unique(sapply(strsplit(del_locations,"q"),"[[",1)))
chr <- unique(sapply(strsplit(chr,"p"),"[[",1))

del_cytoband_bed <- cytoband_bed[cytoband_bed$V1 %in% chr,]
del_cytoband_bed$new_cytoband <- paste0(sapply(strsplit(del_cytoband_bed$V1,"chr"),"[[",2),del_cytoband_bed$V4)
del_cytoband_bed <- del_cytoband_bed[del_cytoband_bed$new_cytoband %in% del_locations,]

del_female_bed <- merge(female_only_del_cytoband_gene_df,del_cytoband_bed,by.x = "cytoband", by.y = "new_cytoband")
del_female_bed <- del_female_bed [,c(2,3,4,5)]
colnames(del_female_bed) <- c("gene", "chr","start","end")
del_female_bed <- del_female_bed[,c("chr","start","end","gene")]
del_female_bed$chr <- paste0("Del_",del_female_bed$chr)
### combine
female_gene_circos <- rbind(amp_female_bed,del_female_bed)
nrow(female_gene_circos)
#### 8.2 Male only #####
### Amp
amp_locations <- male_only_amp_cytoband_gene_df$cytoband
chr <- paste0("chr",unique(sapply(strsplit(amp_locations,"q"),"[[",1)))
chr <- unique(sapply(strsplit(chr,"p"),"[[",1))

amp_cytoband_bed <- cytoband_bed[cytoband_bed$V1 %in% chr,]
amp_cytoband_bed$new_cytoband <- paste0(sapply(strsplit(amp_cytoband_bed$V1,"chr"),"[[",2),amp_cytoband_bed$V4)
amp_cytoband_bed <- amp_cytoband_bed[amp_cytoband_bed$new_cytoband %in% amp_locations,]

amp_male_bed <- merge(male_only_amp_cytoband_gene_df,amp_cytoband_bed,by.x = "cytoband", by.y = "new_cytoband")
amp_male_bed <- amp_male_bed [,c(2,3,4,5)]
colnames(amp_male_bed) <- c("gene", "chr","start","end")
amp_male_bed <- amp_male_bed[,c("chr","start","end","gene")]
amp_male_bed$chr <- paste0("Amp_",amp_male_bed$chr)
### Del
del_locations <- male_only_del_cytoband_gene_df$cytoband
chr <- paste0("chr",unique(sapply(strsplit(del_locations,"q"),"[[",1)))
chr <- unique(sapply(strsplit(chr,"p"),"[[",1))

del_cytoband_bed <- cytoband_bed[cytoband_bed$V1 %in% chr,]
del_cytoband_bed$new_cytoband <- paste0(sapply(strsplit(del_cytoband_bed$V1,"chr"),"[[",2),del_cytoband_bed$V4)
del_cytoband_bed <- del_cytoband_bed[del_cytoband_bed$new_cytoband %in% del_locations,]

del_male_bed <- merge(male_only_del_cytoband_gene_df,del_cytoband_bed,by.x = "cytoband", by.y = "new_cytoband")
del_male_bed <- del_male_bed [,c(2,3,4,5)]
colnames(del_male_bed) <- c("gene", "chr","start","end")
del_male_bed <- del_male_bed[,c("chr","start","end","gene")]
del_male_bed$chr <- paste0("Del_",del_male_bed$chr)
### combine
male_gene_circos <- rbind(amp_male_bed,del_male_bed)
nrow(male_gene_circos)
#### 8.3 Intersection ####
### no del
### Amp
amp_locations <- amp_sig_chisq_result$cytoband
chr <- paste0("chr",unique(sapply(strsplit(amp_locations,"q"),"[[",1)))
chr <- unique(sapply(strsplit(chr,"p"),"[[",1))

amp_cytoband_bed <- cytoband_bed[cytoband_bed$V1 %in% chr,]
amp_cytoband_bed$new_cytoband <- paste0(sapply(strsplit(amp_cytoband_bed$V1,"chr"),"[[",2),amp_cytoband_bed$V4)
amp_cytoband_bed <- amp_cytoband_bed[amp_cytoband_bed$new_cytoband %in% amp_locations,]

amp_intersection_bed <- merge(amp_sig_chisq_result,amp_cytoband_bed,by.x = "cytoband", by.y = "new_cytoband")
amp_intersection_bed <- amp_intersection_bed [,c(4,5,6,7)]
colnames(amp_intersection_bed) <- c("gene", "chr","start","end")
amp_intersection_bed <- amp_intersection_bed[,c("chr","start","end","gene")]
amp_intersection_bed$chr <- paste0("Amp_",amp_intersection_bed$chr)
### combine
intersection_gene_circos <- amp_intersection_bed


###### 8.4 check intersection
intersect(female_gene_circos$gene,exp_set)
intersect(male_gene_circos$gene,exp_set)
intersect(intersection_gene_circos$gene,exp_set)
############## 6. Initialize circos ##############
#pdf("results/BaiduII_OV60/cnv_circos.pdf")
human1_cytoband = read.cytoband(species = "hg19")$df
human2_cytoband = read.cytoband(species = "hg19")$df

human1_cytoband[ ,1] = paste0("Del_", human1_cytoband[, 1])
human2_cytoband[ ,1] = paste0("Amp_", human2_cytoband[, 1])

cytoband = rbind(human1_cytoband, human2_cytoband)
head(cytoband)

chromosome.index = c(paste0("Del_chr", c(1:22)), 
                     rev(paste0("Amp_chr", c(1:22))))
##### Initialize  #####
circos.par(gap.after = c(rep(1, 21), 5, rep(1, 21), 5))
circos.initializeWithIdeogram(cytoband, plotType = NULL, 
                              chromosome.index = chromosome.index)
circos.track(ylim = c(0, 0.8), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = TRUE)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)

highlight.chromosome(paste0("Amp_chr", c(1:22)), 
                     col =  pal_jco()(5)[2] , track.index = 1)
highlight.chromosome(paste0("Del_chr", c(1:22)), 
                     col = pal_jco()(5)[1], track.index = 1)

circos.genomicIdeogram(cytoband)

############# 7. Add  female track ############
#colnames(scores_female_amp) <- c("chr",'start',"end","value1")
#colnames(scores_male_amp) <- c("chr",'start',"end","value1")
track_female_amp <- scores_female_amp
track_female_amp[,1] <- paste0("Amp_",scores_female_amp$Chromosome)

track_female_del <- scores_female_del
track_female_del[,1] <- paste0("Del_",scores_female_del$Chromosome)

final_female_track <- rbind(track_female_amp, track_female_del)
## Color design
color_assign_f <- colorRamp2(breaks = seq(0,0.5,0.5/1000), 
                             col = viridis(option = "B",2001)[800:1800])
### readRDS
final_female_track <- readRDS("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/final_female_track.rds")
circos.genomicTrack(final_female_track, ylim=c(-0.5, 2.5),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h",col = color_assign_f(value[[1]]),...)
                      
                    })

############# 8. Add male track ############
track_male_amp <- scores_male_amp
track_male_amp[,1] <- paste0("Amp_",scores_male_amp$Chromosome)

track_male_del <- scores_male_del
track_male_del[,1] <- paste0("Del_",scores_male_del$Chromosome)

final_male_track <- rbind(track_male_amp, track_male_del)

## Color design
color_assign_m <- colorRamp2(breaks = seq(0,0.5,0.5/1000), 
                             col = viridis(option = "D",1500)[200:1200])
### readRDS
final_male_track <- readRDS("~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/final_male_track.rds")
circos.genomicTrack(final_male_track, ylim=c(-0.5, 2.5),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h",col = color_assign_m(value[[1]]),...)
                      #circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "black")
                    })

############### 9.Add gene annotation ##########
## add female only track
# human_choose_female_only_genes <- c("IFNA1","AKT3","IFNA2","AGT","RPS6",
#                                     "ATM","CTLA4","CCL19","EIF2AK4","HMGB1")
# final_female_gene_circos <- female_gene_circos[female_gene_circos$gene %in% human_choose_female_only_genes,]
# final_female_gene_circos <- final_female_gene_circos[c(1:7,9,10,14),]
# circos.genomicLabels(final_female_gene_circos, labels.column = 2, side = "inside",
#                      cex = 0.2,col = "red",niceFacing = T)
# 
# ## add male only track
# human_choose_male_only_genes <- c("FGF6","IL1A","IL1B","FGF23","IL1F10",
#                                   "UGT1A1","UGT1A3","FCER1G","GSTM1","FCER1A")
# final_male_gene_circos <- male_gene_circos[male_gene_circos$gene %in% human_choose_male_only_genes,]
# 
# circos.genomicLabels(final_male_gene_circos, labels.column = 4, side = "inside",cex = 0.2,col = "blue")
## add intersection track
intersection_gene_circos <- readRDS("results/BaiduII_OV60//7. CNV/intersection_gene_circos.rds")
circos.genomicLabels(intersection_gene_circos, labels.column = 4, side = "inside",cex = 0.7,col = "blue")

## test together
final_gene_circos <- rbind(final_male_gene_circos,final_female_gene_circos,intersection_gene_circos)
final_gene_circos <- readRDS("results/BaiduII_BL60/final_gene_circos.rds")

circos.genomicLabels(final_gene_circos, labels.column = 7, side = "inside",cex = 0.7,col = "blue")

############### 10. Add legends #############
text(0, -1.1, "Del")
text(0, 1.1, "Amp")

lgd_female= Legend(at = c(0, 1.5, 3), col_fun = color_assign_f, 
                   title_position = "topleft", title = "Female",direction = "horizontal")

lgd_male= Legend(at = c(0, 1.5, 3), col_fun = color_assign_m, 
                 title_position = "topleft", title = "Female",direction = "horizontal")
lgd_list = packLegend( lgd_female,lgd_male)

draw(lgd_list)
############### 9. Clear ##########
circos.clear()
dev.off()
############## 10. Save result ########
saveRDS(sig_chisq_result,"results/BaiduII_OV60/CNV_chisq_results.rds")
saveRDS(intersection_gene_circos,"results/BaiduII_OV60//7. CNV/intersection_gene_circos.rds")
female_thresholded_sig_fvsm <- female_thresholded_gene[sig_fvsm_genes,]
male_thresholded_sig_fvsm <- male_thresholded_gene[sig_fvsm_genes,]
write.csv(male_thresholded_sig_fvsm,"results/BaiduII_OV60/male_thresholded_sig_fvsm.csv",row.names = T)

saveRDS(final_male_track,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/final_male_track.rds")
saveRDS(final_female_track,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/7. CNV/final_female_track.rds")
