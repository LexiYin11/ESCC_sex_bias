#### 0. Library import and functions ####
library(maftools)
library(gtools)
library(viridis)
get_gene_mutation_rate <- function(mut){
    gs = getGeneSummary(mut)
    nsamps = as.numeric(mut@summary[ID %in% "Samples", summary])
    gs.load = gs[,.(Hugo_Symbol, AlteredSamples)]
    gs.load[,AlteredSamples := round(AlteredSamples/nsamps, digits = 2) * 100]
    data.table::setDF(x = gs.load, rownames = gs.load$Hugo_Symbol)
    gs = gs[,colnames(gs)[!colnames(x = gs) %in% c('total', 'Amp', 'Del', 'CNV_total', 'MutatedSamples', 'AlteredSamples')], with = FALSE]
    
    selected_gs.load <- gs.load[gs.load$AlteredSamples >= 0,]
    
    gs <- as.data.frame(gs)
    gs <- gs[gs$Hugo_Symbol %in% selected_gs.load$Hugo_Symbol,]
    gs <- gs[order(match(gs$Hugo_Symbol, selected_gs.load$Hugo_Symbol)),]
    gs$load <- selected_gs.load$AlteredSamples
    edited_gs_df <- NA
    for (variant in  colnames(gs)[2:9]){
      cur_df <- as.data.frame(gs)[,c("Hugo_Symbol",variant,"load")]
      cur_df$type <- variant
      colnames(cur_df)[2] <- "value"
      if(length(edited_gs_df) == 1){
        edited_gs_df <- cur_df
      } else {
        edited_gs_df <- rbind( edited_gs_df, cur_df)
      }
    }
    final_df <- edited_gs_df[,c("Hugo_Symbol","load")]
    return(final_df)
}
color_order <- viridis(14,option ="C")
names(color_order) <- c(new_variant_types,"Multi_Hit")
color_order <- c("Missense_Mutation" = "#F0F921FF", "Nonsense_Mutation" ="#FBD424FF", 
                 "Frame_Shift_Del" = "#FDB32FFF", "Splice_Site" = "#F89441FF",
                 "Frame_Shift_Ins" = "#ED7953FF","In_Frame_Del" = "#DE5F65FF",
                 "In_Frame_Ins" = "#CC4678FF", "Nonstop_Mutation" = "#B52F8CFF",
                 "3'UTR" =   "#9C179EFF", "3'Flank" = "#7E03A8FF",
                 "5'UTR" =   "#5D01A6FF" , "5'Flank" = "#3B049AFF","Silent" = "#0D0887FF")

#### 1. Get Pathway gene sets #####
pathway_gene_list <- c()
for ( i in 1:10){
  pathway_genes <- na.omit(xlsx::read.xlsx("rawdata/pathway_genes.xlsx",sheetIndex = i,colIndex = 1))
  pathway_gene_list <- c(pathway_gene_list,list(pathway_genes[,1]))
}
names(pathway_gene_list) <- c("Cell_cycle","HIPPO","MYC","NOTCH","NRF2","PI3K",
                              "TGF-Beta","RTK-RAS","TP53","WNT")
#### 2. Import OV60 data  #######

pre_mut_female_combine_ov60 <- read.table("rawdata/pre_mut_female_combine_ov60.xls",header = T)
pre_mut_male_combine_ov60 <- read.table("rawdata/pre_mut_male_combine_ov60.xls",header = T)

new_variant_types <- unique(pre_mut_female_combine_ov60$Variant_Classification)
new_variant_types <- new_variant_types[!new_variant_types %in% c("IGR","Intron")]

mut_female_combine_ov60 <- read.maf(maf=pre_mut_female_combine_ov60,vc_nonSyn = new_variant_types)
mut_male_combine_ov60 <- read.maf(maf=pre_mut_male_combine_ov60,vc_nonSyn = new_variant_types)

### get gene mutation rate
female_ov60 <- get_gene_mutation_rate(mut_female_combine_ov60)
male_ov60 <- get_gene_mutation_rate(mut_male_combine_ov60)

dim(female_ov60)
dim(male_ov60)

male_ov60$load <- -male_ov60$load
merged_ov60 <- merge(female_ov60,male_ov60,by = "Hugo_Symbol",all = T)
merged_ov60 <- merged_ov60[!duplicated(merged_ov60),]
colnames(merged_ov60) <- c("Gene","Female","Male")
merged_ov60 <- na.replace(merged_ov60,0)
dim(merged_ov60)

#write.table(merged_ov60,"results/new_variant_type_mutation_rate_fandm_ov60.txt",sep = "\t",quote = F,row.names = F)
### compare FvsM
fvsm_compare_ov60 <- mafCompare(m1 = mut_female_combine_ov60, m2 = mut_male_combine_ov60,
                                m1Name = 'Female', m2Name = 'Male', minMut = 1,
                                useCNV = F)
fvsm_compare_ov60_pathway <- mafCompare(m1 = mut_female_combine_ov60, m2 = mut_male_combine_ov60,
                                m1Name = 'Female', m2Name = 'Male', minMut = 1,
                                useCNV = F,pathways = T)
fvsm_compare_ov60_sig_ori <- fvsm_compare_ov60$results[fvsm_compare_ov60$results$pval <= 0.05,]
fvsm_compare_ov60_sig  <- fvsm_compare_ov60_sig_ori$Hugo_Symbol

fvsm_compare_ov60_pathway

## Forest plot
p3_forestplot <- forestPlot(mafCompareRes = fvsm_compare_ov60, pVal = 0.05 ,color = c('royalblue', 'maroon'), geneFontSize = 0.8)
## save
#write.csv(fvsm_compare_ov60_sig_ori,"results/BaiduII_OV60/8. Mutation/fvsm_compare_ov60_sig_ori.csv",col.names = F,row.names = F,quote = F)
## check if pathway genes exist in FvsM sig
pathway_gene_list_fvsm_ov60 <- pathway_gene_list
for (n in 1:10){
  pathway_gene_list_fvsm_ov60[[n]] <- intersect(pathway_gene_list[[n]],fvsm_compare_ov60_sig )
}
pathway_gene_list_fvsm_ov60

fvsm_compare_ov60_sig_ori[fvsm_compare_ov60_sig_ori$Hugo_Symbol %in% unlist(pathway_gene_list_fvsm_ov60),]

## plot sig pathway (PI3K MYC)
pdf("results/ov60_fvsm_PI3K_pathway_mutation.pdf")
PlotOncogenicPathways(maf = mut_female_combine_ov60, pathways = "PI3K")
PlotOncogenicPathways(maf = mut_male_combine_ov60, pathways = "PI3K")
dev.off()
pdf("results/ov60_fvsm_MYC_pathway_mutation.pdf")
PlotOncogenicPathways(maf = mut_female_combine_ov60, pathways = "MYC")
PlotOncogenicPathways(maf = mut_male_combine_ov60, pathways = "MYC")
dev.off()

## forest plot
# pathway
pdf("results/ov60_fvsm_pathway_forest_plot.pdf")
forestPlot(mafCompareRes = fvsm_compare_ov60_pathway, pVal = 1,
                              color = c('royalblue', 'maroon'), geneFontSize = 0.8)
dev.off()

# seperate pathways
for (i in 1:10){
  fvsm_compare_ov60_pathway <- fvsm_compare_ov60
  fvsm_compare_ov60_pathway$results <- fvsm_compare_ov60_pathway$results[fvsm_compare_ov60_pathway$results$Hugo_Symbol %in% pathway_gene_list[[i]],]
  pdf(paste0("results/pathway/ov60/ov60_fvsm",names(pathway_gene_list)[i],"_forest_plot.pdf"),height = 10)
  forestPlot(mafCompareRes = fvsm_compare_ov60_pathway, pVal = 1,
             color = c('royalblue', 'maroon'), geneFontSize = 0.8)
  dev.off()
  
}


#### 3. Import BL60 data ######
pre_mut_female_combine_bl60 <- read.table("rawdata/pre_mut_female_combine_bl60.xls",header = T)
pre_mut_male_combine_bl60 <- read.table("rawdata/pre_mut_male_combine_bl60.xls",header = T)

new_variant_types <- unique(pre_mut_female_combine_bl60$Variant_Classification)
new_variant_types <- new_variant_types[!new_variant_types %in% c("IGR","Intron")]

mut_female_combine_bl60 <- read.maf(maf=pre_mut_female_combine_bl60,vc_nonSyn = new_variant_types)
mut_male_combine_bl60 <- read.maf(maf=pre_mut_male_combine_bl60,vc_nonSyn = new_variant_types)
### get gene mutation rate
female_bl60 <- get_gene_mutation_rate(mut_female_combine_bl60)
male_bl60 <- get_gene_mutation_rate(mut_male_combine_bl60)

dim(female_bl60)
dim(male_bl60)

male_bl60$load <- -male_bl60$load
merged_bl60 <- merge(female_bl60,male_bl60,by = "Hugo_Symbol",all = T)
merged_bl60 <- merged_bl60[!duplicated(merged_bl60),]
merged_bl60 <- na.replace(merged_bl60,0)
colnames(merged_bl60) <- c("Gene","Female","Male")
dim(merged_bl60)

#write.table(merged_bl60,"results/new_variant_type_mutation_rate_fandm_bl60.txt",sep = "\t",quote = F,row.names = F)
## merge ov60 and bl60
merged_both <- merge(merged_ov60,merged_bl60, by = "Gene")
colnames(merged_both) <- c("Gene","OV60_female","OV60_male","BL60_female","BL60_male")

#write.table(merged_both,"results/new_variant_type_mutation_rate_fandm_both.txt",sep = "\t",quote = F,row.names = F)

### compare FvsM
fvsm_compare_bl60 <- mafCompare(m1 = mut_female_combine_bl60, m2 = mut_male_combine_bl60, m1Name = 'Female', m2Name = 'Male', minMut = 1,useCNV = F)
fvsm_compare_bl60_pathway <- mafCompare(m1 = mut_female_combine_bl60, m2 = mut_male_combine_bl60,
                                        m1Name = 'Female', m2Name = 'Male', minMut = 1,
                                        useCNV = F,pathways = T)

fvsm_compare_bl60_sig_ori <- fvsm_compare_bl60$results[fvsm_compare_bl60$results$pval <= 0.05,]
fvsm_compare_bl60_sig  <- fvsm_compare_bl60_sig_ori$Hugo_Symbol

fvsm_compare_bl60_pathway
## check if pathway genes exist in FvsM sig
pathway_gene_list_fvsm_bl60 <- pathway_gene_list
for (n in 1:10){
  pathway_gene_list_fvsm_bl60[[n]] <- intersect(pathway_gene_list[[n]],fvsm_compare_bl60_sig )
}
pathway_gene_list_fvsm_bl60

fvsm_compare_bl60_sig_ori[fvsm_compare_bl60_sig_ori$Hugo_Symbol %in% unlist(pathway_gene_list_fvsm_bl60),]

## forest plot
# pathway
pdf("results/pathway/bl60/bl60_fvsm_pathway_forest_plot.pdf")
forestPlot(mafCompareRes = fvsm_compare_bl60_pathway, pVal = 1,
           color = c('royalblue', 'maroon'), geneFontSize = 0.8)
dev.off()

# seperate pathways
for (i in 1:10){
  fvsm_compare_bl60_pathway <- fvsm_compare_bl60
  fvsm_compare_bl60_pathway$results <- fvsm_compare_bl60_pathway$results[fvsm_compare_bl60_pathway$results$Hugo_Symbol %in% pathway_gene_list[[i]],]
  pdf(paste0("results/pathway/bl60/bl60_fvsm",names(pathway_gene_list)[i],"_forest_plot.pdf"),height = 10)
  forestPlot(mafCompareRes = fvsm_compare_bl60_pathway, pVal = 1,
             color = c('royalblue', 'maroon'), geneFontSize = 0.8)
  dev.off()
  
}
#### 4. Pathway significant genes coOncoplot #####
mutation_cell_cycle <- c("TP53","CDKN2A","MDM4","FBXW7","CCND1","MDM2","ATM","CCNE1","CDKN1A")
mutation_WNT <- c("GSK3B","APC","CTNNB1","AMER1","TLE1","TLE4","TCF7L2","MYC")
mutation_PI3K <- c("ERBB4","PIK3CB","PIK3CA","KIT","IRS2","MAPK1","IGF1R","AKT3","ERBB3","PTEN","ARAF",
                   "RET")
mutation_all <- c(mutation_cell_cycle, mutation_WNT,mutation_PI3K)
write.table(mutation_all,"results/mutation_pathway_all.txt",row.names = F, quote = F)

pdf("results/cellcyle_mutation_coOncoplot_pathway_draft2.pdf",width = 5)
coOncoplot(m1 = mut_female_combine_ov60, m2 = mut_male_combine_ov60, m1Name = 'Female', m2Name = 'Male',
           genes = mutation_cell_cycle,colors = color_order,keepGeneOrder = F,geneNamefont = 0.6)
coOncoplot(m1 = mut_female_combine_bl60, m2 = mut_male_combine_bl60, m1Name = 'Female', m2Name = 'Male',
           genes = mutation_cell_cycle,colors = color_order,keepGeneOrder = F,geneNamefont = 0.6)

dev.off()

pdf("results/WNT_mutation_coOncoplot_pathway.pdf",width = 5)
coOncoplot(m1 = mut_female_combine_ov60, m2 = mut_male_combine_ov60, m1Name = 'Female', m2Name = 'Male',
           genes = mutation_WNT,colors = color_order,keepGeneOrder = F,geneNamefont = 0.6)
coOncoplot(m1 = mut_female_combine_bl60, m2 = mut_male_combine_bl60, m1Name = 'Female', m2Name = 'Male',
           genes = mutation_WNT,colors = color_order,keepGeneOrder = F,geneNamefont = 0.6)

dev.off()

pdf("results/PI3K_mutation_coOncoplot_pathway.pdf",width = 5)
coOncoplot(m1 = mut_female_combine_ov60, m2 = mut_male_combine_ov60, m1Name = 'Female', m2Name = 'Male',
           genes = mutation_PI3K,colors = color_order,keepGeneOrder = F,geneNamefont = 0.6)
coOncoplot(m1 = mut_female_combine_bl60, m2 = mut_male_combine_bl60, m1Name = 'Female', m2Name = 'Male',
           genes = mutation_PI3K,colors = color_order,keepGeneOrder = F,geneNamefont = 0.6)

dev.off()





