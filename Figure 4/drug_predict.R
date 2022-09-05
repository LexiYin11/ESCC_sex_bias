########### 1. Library import  #########
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(tibble)
library(data.table)
library(drugbankR)
library(ggsci)

########### 2. Training expression data and response data prep ############
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='./rawdata/Training Data'
dir="/thinker/tb269store/liuLab/Yin_yin/training_data"

GDSC2_Expr = readRDS(file = file.path(dir,"GDSC2_Expr (RMA Normalized and Log Transformed).rds"))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res)

CTRP2_Expr = readRDS(file = file.path(dir,"CTRP2_Expr (TPM, log2(x+1) Transformed).rds"))
CTRP2_Res = readRDS(file = file.path(dir,"CTRP2_Res.rds"))

##  our data
tpm <- read.csv("rawdata/baiduII_tpm_t.csv",row.names = 1,header = T)
tpm <- read.csv("/thinker/tb269store/liuLab/Yin_yin/rawData/BaiduII/baiduII_tpm_t_for_CTRP_predict.csv",row.names = 1,header = T)
colnames(tpm) <- gsub("X","",colnames(tpm))
colnames(tpm) <- as.integer(colnames(tpm))
tpm <- as.matrix(tpm)


########### 3. Gene extraction ############
## extract gene ID list
gene_list_fvsm <- read.delim("results/BaiduII_OV60/gene_lists/all_gene_symbol.txt",header = F)
gene_list_fvsm <- gene_list_fvsm[!duplicated(gene_list_fvsm),]

## Merge two gene lists
gene_list_total <- union(gene_list_fvsm, gene_list_tvsn)
length(gene_list_total)

## tpm gene ID convert
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id","chromosome_name"), filters="ensembl_gene_id", values=rownames(tpm), mart=human)

tpm_gene_conversion <- data.frame(rownames(tpm), gene_coords$hgnc_symbol)
tpm_gene_conversion <- tpm_gene_conversion[tpm_gene_conversion[,2] %in% gene_list_total,]
tpm_gene_conversion <- tpm_gene_conversion[!duplicated(tpm_gene_conversion[,1]),]
tpm_gene_conversion <- tpm_gene_conversion[!duplicated(tpm_gene_conversion[,2]),]

tpm <- tpm[tpm_gene_conversion[,1],]
rownames(tpm) <- tpm_gene_conversion[,2]

###########（Archive） Other data symbol conversion  ###########
### convert symbols
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=rownames(tpm),mart= mart)
rownames(tpm) <- G_list$hgnc_symbol
length(G_list$hgnc_symbol)
### if there is duplication
#undup <- unique(G_list$hgnc_symbol)
#length(undup)
#rownames(tpm) <- toupper(rownames(tpm))

#length(intersect(rownames(tpm), rownames(CTRP2_Expr)))
G_list_CTRP2 <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=rownames(CTRP2_Expr),mart= mart)
nrow(G_list_CTRP2)
nrow(CTRP2_Expr)

### if there is duplication
unq <- G_list_CTRP2[!duplicated(G_list_CTRP2$hgnc_symbol),]
dup <- G_list_CTRP2$hgnc_symbol[duplicated(G_list_CTRP2$hgnc_symbol)]
#dup_hgnc <- G_list_CTRP2[!duplicated(G_list_CTRP2$hgnc_symbol),]
length(dup)
dup_G_list <- c()
for (gene_name in dup){
  list_1 <- G_list_CTRP2[G_list_CTRP2$hgnc_symbol == gene_name,1]
  chosen_one <- c()
  for (id in list_1){
    if (id %in% rownames(tpm) && length(chosen_one) == 0){
      chosen_one = id
    }
  }
  if (length(chosen_one) != 0){
      if (length(dup_G_list) == 0){
        dup_G_list <- G_list_CTRP2[G_list_CTRP2$ensembl_gene_id == chosen_one,]
      } else {
        dup_G_list <- rbind(dup_G_list,G_list_CTRP2[G_list_CTRP2$ensembl_gene_id == chosen_one,])
      }
  }
}
dup_G_list <- dup_G_list[!duplicated(dup_G_list),]
unq <- rbind(unq,dup_G_list)
CTRP2_Expr <- CTRP2_Expr[unq$hgnc_symbol,]
rownames(CTRP2_Expr) <- unq$ensembl_gene_id

########### 4. OncoPredict! ############

calcPhenotype(trainingExprData = CTRP2_Expr,
              trainingPtype = CTRP2_Res,
              testExprData = as.matrix(tpm),
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
########### 5. Result analysis  ##########
final_result <- read.csv('results/BaiduII_OV60/5.Drug_prediction/GDSC/GDSC2_DrugPredictions.csv',row.name = 1)
dim(final_result)
rownames(final_result) <- as.numeric(gsub("X","",rownames(final_result)))
## filted by age
metadata_survival <- metadata[,c("sample.ID............BDESCC2..", "overall.survival...............1..dead..0.alive.",
                                 "overall.survival........time","Gender","Age", "Location","Family.History.of.ESCC","Family.History.of.other.cancers","recurrence.or.metastasis.....1..dead.recurrence.metastasis..0..free.","TNM.stage.the.Eighth.Edition.","TNM.stage..the.7th.Edition.","Grade","Smoking.history","Drinking.history","T","N","M","TNM.status")]
colnames(metadata_survival)<- c("Sample_ID","Survival_status","Survival_time",  "Gender", "Age", "Location", "ESCC_family_history", "Other_cancer_family_history", "Recurrence_status","TNM_stage_8","TNM_stage_7","Grade","smoking","drinking","T","N","M","TNM.status")
metadata_survival = metadata_survival[order(metadata_survival$Sample_ID),]
age_metadata <- mutate(metadata_survival, AG = ifelse((Age < 60), "LT60", "OV60"),
                       AG = factor(AG)
)
new_metadata <- age_metadata[age_metadata$AG == "OV60",]
metadata_gender <- new_metadata[,c("Sample_ID","Gender")]
final_result <- final_result[as.character(metadata_gender$Sample_ID),]
##########  6. Top 10 drugs in total, Female, Male #############

###### total
colmean <- sort(colMeans(final_result, na.rm = T))
colmean <- as.data.frame(colmean)
colmean$med <- rownames(colmean)
colnames(colmean) <- c("IC50","med")
# round to 0.01
colmean$IC50 <- round(colmean$IC50,3)
# Keep drug list
drug_list <- colmean[order(colmean$IC50,decreasing = F),]
# Extract top10 drugs
colmean <- head(colmean[order(colmean$IC50,decreasing = F),],10)
p1 <- ggplot(colmean) + geom_bar(aes(x = reorder(med,
    IC50), y = IC50,fill = factor(IC50)), stat = "identity",
    show.legend = F) + scale_fill_brewer(palette="RdYlGn") + geom_text(aes(x = reorder(med, -IC50),
    y = IC50,label = IC50),position = position_dodge(width=0.9),
    vjust=-0.25) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs( x = "medicines")
p1
ggsave("results/BaiduII_OV60/5.Drug_prediction/GDSC/GDSC_total_med_predict.pdf")
colmean_total <- colmean
## save drug list
write.csv(as.data.frame(drug_list),'results/BaiduII_OV60/5.Drug_prediction/GDSC/GDSC_total_med_predict_drug_list.csv',row.names = F)
## make a table with annotation from DrugBank
drugbank_dataframe <- dbxml2df(xmlfile="rawdata/drugbank.xml", version="5.1.3") 

########### Female
## Baidu
metadata <- read.csv(file = "rawdata/baiduII_clinical_data.csv", header = TRUE)
metadata_gender <- metadata[,c("sample.ID............BDESCC2..", "Gender")]
colnames(metadata_gender)<- c("Sample_ID", "Gender")
metadata_gender = metadata_gender[order(metadata_gender$Sample_ID),]

female_result <- final_result[as.character(metadata_gender[metadata_gender$Gender == "Female", "Sample_ID"]),]
colmean <- colMeans(female_result, na.rm = T)
colmean <- as.data.frame(colmean)
colmean$med <- rownames(colmean)
colnames(colmean) <- c("IC50","med")
# round to 0.01
colmean$IC50 <- round(colmean$IC50,3)

colmean <- head(colmean[order(colmean$IC50,decreasing = F),],10)
p2 <- ggplot(colmean) + geom_bar(aes(x = reorder(med,
     IC50), y = IC50,fill = factor(IC50)), stat = "identity",
    show.legend = F) + scale_fill_brewer(palette="RdYlGn") + geom_text(aes(x = reorder(med, IC50),
     y = IC50,label = IC50),position = position_dodge(width=0.9),vjust=-0.25) + theme(axis.text.x = element_text(angle = 90, 
    vjust = 0.5, hjust=1)) + labs( x = "medicines") 

p2
ggsave("results/BaiduII_CTRP_Female_med_predict.pdf")
colmean_female <- colmean
#### Male
male_result <- final_result[as.character(metadata_gender[metadata_gender$Gender == "Male", "Sample_ID"]),]
dim(male_result)
rownames(male_result)
colmean <- colMeans(male_result, na.rm = T)
colmean <- as.data.frame(colmean)
colmean$med <- rownames(colmean)
colnames(colmean) <- c("IC50","med")
# get log and absolute value
#colmean$IC50 <- abs(log10(colmean$IC50 ))
# round to 0.01
colmean$IC50 <- round(colmean$IC50,3)
# filter value with IC50 beyond 3
colmean <- head(colmean[order(colmean$IC50,decreasing = F),],10)
p3 <- ggplot(colmean) + geom_bar(aes(x = reorder(med,IC50), 
                                     y = IC50,fill = factor(IC50)), 
                                 stat = "identity", show.legend = F) + scale_fill_brewer(palette="RdYlGn") + geom_text(aes(x = reorder(med, 
                                 IC50), y = IC50,label = IC50),
                                 position = position_dodge(width=0.9),
                                 vjust=-0.25) + theme(axis.text.x = element_text(angle = 90, 
                                  vjust = 0.5, hjust=1)) + labs( x = "medicines") 
p3
ggsave("results/BaiduII_CTRP_Male_med_predict.pdf")

colmean_male <- colmean
##########  7. Sig drugs in F vs M ##############
## pre-filter - test
final_result_ttest <- final_result[as.character(metadata_gender[, "Sample_ID"]),]
metadata_ttest <- cbind(metadata_gender,final_result_ttest)

ttest_result <- c()
for (x in 1:198){
  new_table <- metadata_ttest[,c(2,x+2)]
  colnames(new_table) <- c("Gender", "Drug")
  ttest_result <- c(ttest_result, compare_means(Drug ~ Gender, new_table)$p.adj)
}
result_after_filter <- data.frame(colnames(final_result), ttest_result)
result_after_filter <- result_after_filter[order(result_after_filter[,2]),]
## save the fvsm result
write.csv(result_after_filter,"results/BaiduII_OV60/5.Drug_prediction/GDSC/GDSC_FvsM_drugs_with_padj.csv", row.names = F)

## extract only significant drugs
#Baiduii
sig_med_fvsm <- result_after_filter[result_after_filter$ttest_result < 0.05,1]
## Save
write_delim(as.data.frame(sig_med_fvsm),"results/BaiduII_OV60/5.Drug_prediction/GDSC/GDSC_sig_drugs_fvsm")
## use ggplot detect diff
pdf("results/BaiduII_OV60/5.Drug_prediction/GDSC/GDSC_FvsM_sig_drug_boxplot.pdf",w=6,h=7)
for (medType in sig_med_fvsm){
  p <- ggboxplot(metadata_ttest, x = "Gender", y = medType
                 ,color = 'Gender', ylim = c(0,50),
                 ,legend="none",
                 add.params = list(size=1),xlab = "",
                 add = "jitter") +theme_bw() + stat_compare_means(aes(label = paste0("P = ", ..p.format..)),
                                                                  label.x = 1.5) 
  print(change_palette(p, "simpsons"))
}
dev.off()



## all in one
dat <- final_result_ttest[,sig_med_fvsm] %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Med_type,value = IC50,-Sample)
dat <- na.omit(dat)
rownames(metadata_gender) <- metadata_gender$Sample_ID
dat$Gender <- metadata_gender[dat$Sample,2]

pdf("results//BaiduII_OV60/5.Drug_prediction/GDSC/GDSC_fvsm_sig_allinone_plot.pdf",height = 6, width = 15)
p <- ggplot(dat,aes(reorder(Med_type,IC50),IC50,fill = Gender)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  #ylim(0,50) +
  labs(x = "Drugs", y = "IC50") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
   stat_compare_means(aes(label = ..p.signif..),
                                                                        label.x = 1.5) 
p + scale_fill_jco()
dev.off()

##########  8. Check FvsM in TOP 30 drugs ######
# get top30 drug names
colmean <- sort(colMeans(final_result, na.rm = T))
colmean <- as.data.frame(colmean)
colmean$med <- rownames(colmean)
colnames(colmean) <- c("IC50","med")
colmean$IC50 <- round(colmean$IC50,3)
colmean <- head(colmean[order(colmean$IC50,decreasing = F),],30)
####
## all in one
allinone_final_result <- final_result[,colmean$med]
dat <- allinone_final_result %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Med_type,value = IC50,-Sample)
rownames(metadata_gender) <- metadata_gender$Sample_ID
dat$Gender <- metadata_gender[dat$Sample,2]

pdf("results/BaiduII_OV60/5.Drug_prediction/GDSC/GDSC_top30_fvsm_allinone_plot.pdf",height = 6, width = 15)

p <- ggplot(dat,aes(reorder(Med_type,IC50),IC50,fill = Gender)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Drugs", y = "IC50") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  stat_compare_means(aes(label = ..p.signif..),
                                                                         label.x = 1.5) 
p + scale_fill_jco()
dev.off()

