#### 0. Library import ###############
library(maftools)
library("ggpubr")
library("grid")
library("ggplotify")
color_order <- c("Missense_Mutation" = "#0073C2FF", "Nonsense_Mutation" ="#EFC000FF", 
                 "Frame_Shift_Del" = "#868686FF", "Splice_Site" = "#CD534CFF",
                 "Frame_Shift_Ins" = "#7AA6DCFF","In_Frame_Del" = "#003C67FF",
                 "In_Frame_Ins" = "#8F7700FF", "Nonstop_Mutation" = "#3B3B3BFF")


#### 1. Import Data ########
#### 1.1  Baidu phase 1 mutation   #####
mut_baiduI <- read.csv2('rawdata/baiduI_snv_cnv.csv',header = T,sep = ',')
# clinical data
metadata_baiduI <- read.csv(file = "rawdata/ESCC_pheno_508_20180911.csv", header = TRUE)
metadata_baiduI_ov60 <- metadata_baiduI[metadata_baiduI$age >= 60,]
metadata_baiduI_ov60_gender <- metadata_baiduI_ov60[,c("tumor", "gender")]

## split female and male
mut_baiduI_over60 <- mut_baiduI[mut_baiduI$Tumor_Sample_Barcode %in% metadata_baiduI_ov60_gender$tumor,]
## Split female and male
pre_mut_baiduI_female <- mut_baiduI_over60[mut_baiduI_over60$Tumor_Sample_Barcode %in% metadata_baiduI_ov60_gender[metadata_baiduI_ov60_gender$gender == 2,1],]
pre_mut_baiduI_male <- mut_baiduI_over60[mut_baiduI_over60$Tumor_Sample_Barcode %in% metadata_baiduI_ov60_gender[metadata_baiduI_ov60_gender$gender == 1,1],]
## Read to maftools 

mut_baiduI_female <- read.maf(maf=pre_mut_baiduI_female)
mut_baiduI_male <- read.maf(maf=pre_mut_baiduI_male)

#### 1.2  Baidu phase 2 mutation ####
mut <- read.csv2('rawdata/155ESCC.maf',header = T,sep = '\t')
## Clinical data
metadata <- read.csv(file = "rawdata/baiduII_clinical_data.csv", header = TRUE)
## Filter age
metadata_ov60 <- metadata[metadata$Age >= 60,]
metadata_gender_ov60 <- metadata_ov60[,c("sample.ID............BDESCC2..", "Gender")]
colnames(metadata_gender_ov60)<- c("Sample_ID", "Gender")
metadata_gender_ov60 = metadata_gender_ov60[order(metadata_gender_ov60$Sample_ID),]

mut$patient <- as.numeric(sapply(strsplit(mut$Tumor_Sample_Barcode,"T"),"[[",1))

## modify mutation maf
mut_baiduII_over60 <- mut[mut$patient %in% metadata_gender_ov60$Sample_ID,]
mut_baiduII_over60$Tumor_Sample_Barcode <- as.numeric(gsub("T","",mut_baiduII_over60$Tumor_Sample_Barcode))
## Split female and male
pre_mut_baiduII_female <- mut_baiduII_over60[mut_baiduII_over60$patient %in% metadata_gender_ov60[metadata_gender_ov60$Gender == "Female",1],]
pre_mut_baiduII_male <- mut_baiduII_over60[mut_baiduII_over60$patient %in% metadata_gender_ov60[metadata_gender_ov60$Gender == "Male",1],]

## Read to maftools 
mut_baiduII_female <- read.maf(maf=pre_mut_baiduII_female )
mut_baiduII_male <- read.maf(maf=pre_mut_baiduII_male)

## Add clinical data
#mut_155@clinical.data$Gender <- metadata_gender$Gender
#### 2. Combine baiduI and baiduII #######
## female
pre_mut_female_combine <- rbind(pre_mut_baiduI_female,pre_mut_baiduII_female[,colnames(pre_mut_baiduI_female)])
mut_female_combine <- read.maf(maf=pre_mut_female_combine)
## male
pre_mut_male_combine <- rbind(pre_mut_baiduI_male,pre_mut_baiduII_male[,colnames(pre_mut_baiduI_male)])
mut_male_combine <- read.maf(maf=pre_mut_male_combine)
## total
pre_mut_baidu_IandII_combine <- rbind(pre_mut_female_combine ,pre_mut_male_combine )
mut_baidu_IandII_combine <- read.maf(maf=pre_mut_baidu_IandII_combine)
saveRDS(mut_baidu_IandII_combine,"results/BaiduII_OV60/8. Mutation/mut_baidu_IandII_combine.rds")
mut_baidu_IandII_combine <- readRDS("results/BaiduII_OV60/8. Mutation/mut_baidu_IandII_combine.rds")

## combine clinical info 
metadata_baiduI_ov60 <- metadata_baiduI_ov60[c("tumor","gender","Smoking.Status","Tumor.Grade","Stage.TNM")]
metadata_baiduI_ov60$gender <- ifelse(metadata_baiduI_ov60$gender == "2","Female","Male")
metadata_baiduI_ov60$Smoking.Status <- ifelse(metadata_baiduI_ov60$Smoking.Status =="0","never",
                                              ifelse(metadata_baiduI_ov60$Smoking.Status == "1","light",
                                                     ifelse(metadata_baiduI_ov60$Smoking.Status == "2","moderate",
                                                            ifelse(metadata_baiduI_ov60$Smoking.Status == "3","heavy","stop"))))
metadata_baiduI_ov60$Tumor.Grade <- ifelse(metadata_baiduI_ov60$Tumor.Grade == "1","G1",
                                           ifelse(metadata_baiduI_ov60$Tumor.Grade == "2","G2",
                                                  ifelse(metadata_baiduI_ov60$Tumor.Grade == "3","G3","stop")))
colnames(metadata_baiduI_ov60) <- c("Tumor_Sample_Barcode","Gender","Smoking status","Grade","TNM stage")
metadata_BaiduII_OV60 <- metadata_ov60[,c("sample.ID............BDESCC2..", "Gender",
                                          "Smoking.history",  "Grade","TNM.stage.the.Eighth.Edition.")]
colnames(metadata_BaiduII_OV60) <- c("Tumor_Sample_Barcode","Gender","Smoking status","Grade","TNM stage")
metadata_BaiduII_OV60$`TNM stage` <- ifelse(metadata_BaiduII_OV60$`TNM stage` %in% c("IIIA","IIIB"),3,
                                            ifelse(metadata_BaiduII_OV60$`TNM stage` %in% c("IIA","ⅡA","IIB"),2,
                                                   ifelse(metadata_BaiduII_OV60$`TNM stage` == "IB",1,
                                                          ifelse(metadata_BaiduII_OV60$`TNM stage` == "ⅣA",4,"stop"))))
metadata_baiduI_ov60$`data resource` <- "Baidu phase 1"
metadata_BaiduII_OV60$`data resource` <- "Baidu phase 2"
metadata_baiduI_BaiduII_OV60 <- rbind(metadata_baiduI_ov60,metadata_BaiduII_OV60)
write.csv(metadata_baiduI_BaiduII_OV60,"results/BaiduII_OV60/8. Mutation/metadata_baiduI_BaiduII_OV60.csv",col.names = F,row.names = F,quote = F)




#### 3. Compare Female vs Male ###########
##  baidu phase I phase II combine 
fvsm_compare_combine <- mafCompare(m1 = mut_female_combine, m2 = mut_male_combine, m1Name = 'Female', m2Name = 'Male', minMut = 9,useCNV = F)
saveRDS(fvsm_compare_combine,"results/BaiduII_OV60/8. Mutation/fvsm_compare_combine.rds")
fvsm_compare_combine <- readRDS("results/BaiduII_OV60/8. Mutation/fvsm_compare_combine.rds")
print(fvsm_compare_combine)
fvsm_compare_combine$results <- fvsm_compare_combine$results[order(fvsm_compare_combine$results$Female,fvsm_compare_combine$results$Male,decreasing = T),]
## sig genes (both cohort > 5%)
mut_sig_genes_combine <- fvsm_compare_combine$results[fvsm_compare_combine$results$pval <= 0.05 ,]

mut_sig_genes_combine <- mut_sig_genes_combine$Hugo_Symbol
mut_sig_genes_combine <- c("EPPK1", "KMT2C" ,"DNAH8", "TENM1" ,"KALRN",
                           "APOB" ,  "CACNA1E",  "COL3A1" , "ASPM","OTOGL")

fvsm_compare_combine$results <- fvsm_compare_combine$results[fvsm_compare_combine$results$Hugo_Symbol %in% mut_sig_genes_combine,]
## save 
write.csv(mut_sig_genes_combine,"~/onedrive/Work/phD/phd_project/TME_gender/results/BaiduII_OV60/gene_lists/mutation_baiduI_baiduII_geneset.txt",col.names = F,row.names = F,quote = F)


#### 4. Forest plot ########
pdf("results/BaiduII_OV60/8. Mutation/mutation_forest_fvsm_baidu_combined.pdf",width = 8, height = 8)
p3_forestplot <- forestPlot(mafCompareRes = fvsm_compare_combine, pVal = 0.05 ,color = c('royalblue', 'maroon'), geneFontSize = 0.8)
dev.off()


#### 5. Co-oncoplot ####

pdf("results/BaiduII_OV60/8. Mutation/mutation_oncoplot_fvsm_combine.pdf",width = 5, height = 20,useDingbats = FALSE)
oncoplot(maf  = mut_female_combine,minMut = 0.05,colors = color_order,fontSize = 0.4)
oncoplot(maf = mut_male_combine,minMut = 0.05,colors = color_order,fontSize = 0.4)
dev.off()

### female combine
metadata_baiduI_BaiduII_OV60_female <- metadata_baiduI_BaiduII_OV60[metadata_baiduI_BaiduII_OV60$Gender == "Female",]
mut_female_combine@clinical.data$smoking_status <- metadata_baiduI_BaiduII_OV60_female$`Smoking status`
mut_female_combine@clinical.data$grade <- metadata_baiduI_BaiduII_OV60_female$`Grade`
mut_female_combine@clinical.data$`TNM_stage` <- metadata_baiduI_BaiduII_OV60_female$`TNM stage`
mut_female_combine@clinical.data$`data` <- metadata_baiduI_BaiduII_OV60_female$`data resource`

pdf("results/BaiduII_OV60/8. Mutation/mutation_oncoplot_female_combine.pdf",width = 5, height = 20)
oncoplot(maf  = mut_female_combine,top = 50,colors = color_order,fontSize = 0.8,
         clinicalFeatures = c("data","smoking_status","grade","TNM_stage"),
         annotationColor = list(data = c("Baidu phase 1" = pal_jco()(5)[2],"Baidu phase 2" = pal_jco()(5)[5]),
                                smoking_status = color_smoke,grade = color_grade, TNM_stage = color_stage
         )
)
dev.off()
### female combine
metadata_baiduI_BaiduII_OV60_male <- metadata_baiduI_BaiduII_OV60[metadata_baiduI_BaiduII_OV60$Gender == "Male",]
mut_male_combine@clinical.data$smoking_status <- metadata_baiduI_BaiduII_OV60_male$`Smoking status`
mut_male_combine@clinical.data$grade <- metadata_baiduI_BaiduII_OV60_male$`Grade`
mut_male_combine@clinical.data$`TNM_stage` <- metadata_baiduI_BaiduII_OV60_male$`TNM stage`
mut_male_combine@clinical.data$`data` <- metadata_baiduI_BaiduII_OV60_male$`data resource`

pdf("results/BaiduII_OV60/8. Mutation/mutation_oncoplot_male_combine.pdf",width = 5, height = 20,useDingbats = FALSE)
oncoplot(maf  = mut_male_combine,top = 50,colors = color_order,fontSize = 0.8,
         clinicalFeatures = c("data","smoking_status","grade","TNM_stage"),
         annotationColor = list(data = c("Baidu phase 1" = pal_jco()(5)[2],"Baidu phase 2" = pal_jco()(5)[5]),
                                smoking_status = color_smoke,grade = color_grade, TNM_stage = color_stage
         ))
dev.off()

color_grade <- brewer.pal(3,'YlGnBu')
names(color_grade) <- c("G1","G2","G3")
color_smoke <- brewer.pal(4,'RdPu')
names(color_smoke) <- unique(metadata_baiduI_BaiduII_OV60$`Smoking status`)
color_stage <- brewer.pal(4,"OrRd")
names(color_stage) <- unique(metadata_baiduI_BaiduII_OV60$`TNM stage`)

mut_baidu_IandII_combine@clinical.data$gender <- metadata_baiduI_BaiduII_OV60$Gender
mut_baidu_IandII_combine@clinical.data$smoking_status <- metadata_baiduI_BaiduII_OV60$`Smoking status`
mut_baidu_IandII_combine@clinical.data$grade <- metadata_baiduI_BaiduII_OV60$`Grade`

pdf("results/BaiduII_OV60/8. Mutation/mutation_oncoplot_baiduI_baiduII_combine.pdf",width = 5, height = 10,useDingbats = FALSE)
oncoplot(maf  = mut_baidu_IandII_combine,minMut = 0.05,colors = color_order,fontSize = 0.3,
         sampleOrder = metadata_baiduI_BaiduII_OV60$Tumor_Sample_Barcode,clinicalFeatures = c("gender","smoking_status","grade"),
)
dev.off()

### mutation fvsm significant genes

pdf("results/BaiduII_OV60/8. Mutation/mutation_coOncoplot_fvsm.pdf",width = 18, height = 10)
coOncoplot(m1 = mut_female, m2 = mut_male, m1Name = 'Female', m2Name = 'Male',
           genes = mut_sig_genes,colors = color_order)
dev.off()

pdf("results/BaiduII_OV60/8. Mutation/mutation_coOncoplot_fvsm_combine.pdf",width = 10, height = 8)
coOncoplot(m1 = mut_male_combine, m2 = mut_female_combine, m1Name = 'Female', m2Name = 'Male',
           genes = mut_sig_genes_combine,colors = color_order,keepGeneOrder = T,
           )
dev.off()

### gene expression significant genes
pdf("results/BaiduII_OV60/8. Mutation/mutation_coOncoplot_rna_fvsm.pdf",width = 18, height = 40)
coOncoplot(m1 = mut_female, m2 = mut_male, m1Name = 'Female', m2Name = 'Male',
           genes = degs_fvsm)
dev.off()
### CNV significant genes
pdf("results/BaiduII_OV60/8. Mutation/mutation_coOncoplot_cnv_fvsm.pdf")
coOncoplot(m1 = mut_female, m2 = mut_male, m1Name = 'Female', m2Name = 'Male',
           genes = head(cnv_fvsm,20))
dev.off()

### CNV Exp overlapping significant genes
# Only NOVA1 exists and error
pdf("results/BaiduII_OV60/8. Mutation/mutation_coOncoplot_cnv_exp_fvsm.pdf")
coOncoplot(m1 = mut_female_combine, m2 = mut_male_combine, m1Name = 'Female', m2Name = 'Male',
           genes = "NOVA1")
dev.off()
#### pathway sig genes
mutation_cell_cycle <- c("TP53","CDKN2A","CCND1","MDM2","ATM","CDKN1A","CCNE1","FBXW7")
mutation_MYC <- c("MYC","MYCL","MXD1","MXI1","MXD4","MNT")
mutation_PI3K <- c("ERBB3","KIT","RET","IGF1R","ERBB4","PTEN","IRS2","ARAF",
                   "PIK3CB","PIK3CA","MAPK1","AKT")
pdf("results/BaiduII_OV60/8. Mutation/BaiduII_ov60_mutation_coOncoplot_pathway.pdf")
coOncoplot(m1 = mut_female_combine, m2 = mut_male_combine, m1Name = 'Female', m2Name = 'Male',
           genes = mutation_cell_cycle,colors = color_order,keepGeneOrder = F)
dev.off()

#### 6. CoBarplot ###########
pdf("results/BaiduII_OV60/8. Mutation/mutation_coBarplot_fvsm.pdf",height = 8)
coBarplot(m1 = mut_female_combine, m2 = mut_male_combine, m1Name = 'Female', m2Name = 'Male',
          genes = mut_sig_genes_combine,colors = color_order,orderBy = "m1",geneSize = 0.7)

dev.off()


#### 7. somatic interaction ###########

pdf("results/BaiduII_OV60/8. Mutation/mutation_somatic_int_fandm_combine.pdf",width = 10, height = 10)
somaticInteractions(maf = mut_female_combine, top = 50, pvalue = c(0.01, 0.05),colPal = "RdYlBu",fontSize = 0.5)
somaticInteractions(maf = mut_male_combine, top = 50, pvalue = c(0.01, 0.05),colPal = "RdYlBu",fontSize = 0.5)
dev.off()

sig_fvsm_compare_combine <- fvsm_compare_combine$results[fvsm_compare_combine$results$pval <= 0.05,]
sig_fvsm_compare_combine <- sig_fvsm_compare_combine[(sig_fvsm_compare_combine$Female + sig_fvsm_compare_combine$Female) > 10,]
pdf("results/BaiduII_OV60/8. Mutation/mutation_somatic_int_fvsm_combine.pdf",width = 8, height = 8)
somaticInteractions(maf = mut_baidu_IandII_combine,genes =mut_sig_genes_combine,  
                    pvalue = c(0.01, 0.05),colPal = "RdYlBu",fontSize = 0.7)
dev.off()
mut_sig_genes_combine




