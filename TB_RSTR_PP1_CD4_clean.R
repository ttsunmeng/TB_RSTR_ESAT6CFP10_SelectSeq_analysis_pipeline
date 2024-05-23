# # install.packages('remotes')
# remotes::install_version("Seurat", version = "3.2.3")
# library(Seurat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(xlsx)
library(purrr)
library("readxl")
library(igraph)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggrepel)
library(forcats)
library(scales)
library(ggrepel)
logFC_cutoff <- 0.5
p_cutoff <- 0.05
library("org.Hs.eg.db")
hs <- org.Hs.eg.db
library(reshape2)

setwd("/Volumes/GoogleDrive/My Drive/Stanford/RNA-seq/data_analysis/TB_RSTR/")
############# data inputting #############################################
TB_shrink.seurat <- readRDS('TB_RSTR_PP1_CD4_lognorm.rds')
TB_shrink.seurat <- NormalizeData(TB_shrink.seurat, normalization.method = "LogNormalize", scale.factor = 3e6)
TB_shrink.seurat <- ScaleData(TB_shrink.seurat)
TB_shrink.seurat <- RunPCA(TB_shrink.seurat, features = seurat_gene_list, npcs = 50)
n_pca_selected <- 30
TB_shrink.seurat <- RunUMAP(TB_shrink.seurat, reduction = "pca", dims = 1:n_pca_selected)

cols_clusters <- c('#984EA3','#4DAF4A')
names(cols_clusters) <- c('1','0')

DimPlot(TB_shrink.seurat, group.by = 'batch', label = F, reduction = "umap")
dev.print(pdf, 'TB_RSTR_PP1_CD4_lognormscale_cleaned_umap_batch.pdf',width = 4.3, height = 4)
DimPlot(TB_shrink.seurat, group.by = 'Group', label = F, reduction = "umap",cols = cols_group)
dev.print(pdf, 'TB_RSTR_PP1_CD4_lognormscale_cleaned_umap_Group.pdf',width = 5, height = 4)
DimPlot(TB_shrink.seurat, group.by = 'Donor', label = F, reduction = "umap")
dev.print(pdf, 'TB_RSTR_PP1_CD4_lognormscale_cleaned_umap_Donor.pdf',width = 4.3, height = 4)
DimPlot(TB_shrink.seurat, group.by = 'Group_DonorID', label = F, reduction = "umap")
dev.print(pdf, 'TB_RSTR_PP1_CD4_lognormscale_cleaned_umap_Group_DonorID.pdf',width = 5.5, height = 4)
DimPlot(TB_shrink.seurat, group.by = 'Sex', label = F, reduction = "umap")
dev.print(pdf, 'TB_RSTR_PP1_CD4_lognormscale_cleaned_umap_Sex.pdf',width = 4.7, height = 4)

TB_shrink.seurat <- FindNeighbors(TB_shrink.seurat, reduction = "pca", dims = 1:n_pca_selected)
TB_shrink.seurat <- FindClusters(TB_shrink.seurat, resolution = 0.4)

DimPlot(TB_shrink.seurat, label = F, reduction = "umap", group.by = 'seurat_clusters',cols = cols_clusters)
dev.print(pdf, 'TB_RSTR_PP1_CD4_lognorm_cleaned_umap_clusters.pdf',width = 4.3, height = 4)

##### DEG calculation #####################################
library(ggrepel)
logFC_cutoff <- 0.5
p_cutoff <- 0.05
library("org.Hs.eg.db")
hs <- org.Hs.eg.db

DefaultAssay(TB_shrink.seurat) <- "RNA"
CD4_RSTR_LTBI.markers <- FindMarkers(TB_shrink.seurat, group.by = "Group",ident.1 = 'RSTR',ident.2 = 'LTBI',logfc.threshold = 0,min.pct = 0)
CD4_RSTR_LTBI.markers$gene <- rownames(CD4_RSTR_LTBI.markers)
CD4_RSTR_LTBI.markers$significant <- 'no'
CD4_RSTR_LTBI.markers$significant[((CD4_RSTR_LTBI.markers$avg_log2FC <= -logFC_cutoff) | (CD4_RSTR_LTBI.markers$avg_log2FC >= logFC_cutoff)) & (CD4_RSTR_LTBI.markers$p_val_adj <= p_cutoff)] <- 'yes'
geneID <- AnnotationDbi::select(hs, 
                 keys = CD4_RSTR_LTBI.markers$gene,
                 columns = c("ENTREZID", "SYMBOL"),
                 keytype = "SYMBOL")
# There is a repeat of gene Entrez, pick the later one
# if using the ENTREZID does not work because there is NA value in this column.
p <- match(CD4_RSTR_LTBI.markers$gene,geneID$SYMBOL)
sum((p[2:length(p)] - p[1:length(p) - 1]) != 1)
which((p[2:length(p)] - p[1:length(p) - 1]) != 1)
View(geneID)
# # TEC & MEMO1
geneID <- geneID[((rownames(geneID) != 25385) & (rownames(geneID) != 22838)),]
CD4_RSTR_LTBI.markers$geneID <- geneID$ENTREZID

CD4_RSTR_LTBI.markers$labels <- rownames(CD4_RSTR_LTBI.markers)
CD4_RSTR_LTBI.markers$labels[CD4_RSTR_LTBI.markers$significant == 'no'] <- ''
CD4_RSTR_LTBI.markers$pct_diff <- CD4_RSTR_LTBI.markers$pct.1 - CD4_RSTR_LTBI.markers$pct.2
CD4_RSTR_LTBI.markers <- CD4_RSTR_LTBI.markers %>% arrange(desc(significant),desc(avg_log2FC))
write.xlsx(CD4_RSTR_LTBI.markers,'TB_RSTR_PP1_CD4_lognorm_cleaned_DE_RSTR_vs_LTBI_all.xlsx')

graphics.off()
ggplot(CD4_RSTR_LTBI.markers, aes(x=avg_log2FC, y=-log10(p_val_adj),color = significant,legend = significant)) +
  geom_point(size = 1,alpha = 0.5) +
  scale_colour_manual(values = c("yes" = "red", "no" = "darkgrey")) +
  ggtitle('DE CD4 Tcell RSTR vs LTBI') + #PP1') +
  xlab("log2 Fold Change") + ylab("-log10 adjusted p-value") +
  theme(text = element_text(size = 20)) +
  theme_bw() +
  geom_text_repel(aes(x=avg_log2FC, y=-log10(p_val_adj),label = labels), color = "black", size = 3) +
  geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff), linetype="dotted") +
  geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")
dev.print(pdf, 'TB_RSTR_PP1_CD4_lognorm_cleaned_DE_RSTR_vs_LTBI.pdf',width = 5.5, height = 5)

##### DE GO DAVID visualization #######
# run DEG GO in DAVID website
# https://david.ncifcrf.gov/tools.jsp
table1 <- 'TB_resistance_ESAT6_CD4_lognorm_RSTR_BP_DIRECT.xlsx'
table1_DE_data <- read_excel(table1, sheet = "Sheet1")
RSTR_signficant_genes <- RSTR_signficant_gene_table[!is.na(RSTR_signficant_gene_table$geneID),][,c('gene','geneID')]
RSTR_signficant_genes$geneID <- as.numeric(RSTR_signficant_genes$geneID)
rownames(RSTR_signficant_genes) <- RSTR_signficant_genes$geneID
RSTR_signficant_genes <- RSTR_signficant_genes %>% arrange(desc(geneID))

table1_DE_data$Term <- gsub('.*\\~','',table1_DE_data$Term)
table1_DE_data <- table1_DE_data %>% dplyr::filter(FDR <= 0.05)
BP_direct_select <- c('translation','mitochondrial electron transport, NADH to ubiquinone',
                        'cell-cell adhesion','movement of cell or subcellular component','Wnt signaling pathway, planar cell polarity pathway')

# BP_direct_select <- c('NIK/NF-kappaB signaling','T cell receptor signaling pathway','MAPK cascade',
# 'tumor necrosis factor-mediated signaling pathway')
table1_DE_data <- table1_DE_data[table1_DE_data$Term %in% BP_direct_select,]
table1_DE_data %>%
  mutate(Term = fct_reorder(Term, Count)) %>%
  ggplot( aes(x=Term, y=Count, fill = FDR)) +
  geom_bar(stat="identity", colour="black") +
  coord_flip() +
  # ggtitle('metabolic and basic activities') +
  ggtitle('T cell activation') +
  xlab(expression(paste('GO BP terms of RSTR',symbol('\255')))) + ylab('#genes') + 
  theme_bw() +
  scale_fill_gradientn(colors = c('khaki','#984EA3'),limits = c(0,0.0063)) +
  scale_x_discrete(labels = wrap_format(35)) +
  theme(panel.grid = element_blank(), 
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 10),
              plot.title = element_text(size = 12, face = "bold"))
dev.print(pdf, paste('TB_RSTR_PP1_CD4_lognorm_cleaned_DE_RSTR_DAVID_GO_BF_DIRECT_differentiation.pdf',sep = ''),width = 4.5, height = 2.5)

# LTBI GO
table1 <- 'TB_RSTR_PP1_CD4_stringent_lognorm_RSTR_vs_LTBI_downregulated_sig_DAVID_BP_DIRECT.xlsx'
table1_DE_data <- read_excel(table1, sheet = "Sheet1")
RSTR_signficant_genes <- RSTR_signficant_gene_table[!is.na(RSTR_signficant_gene_table$geneID),][,c('gene','geneID')]
RSTR_signficant_genes$geneID <- as.numeric(RSTR_signficant_genes$geneID)
rownames(RSTR_signficant_genes) <- RSTR_signficant_genes$geneID
RSTR_signficant_genes <- RSTR_signficant_genes %>% arrange(desc(geneID))

table1_DE_data$Term <- gsub('.*\\~','',table1_DE_data$Term)
table1_DE_data <- table1_DE_data %>% dplyr::filter(FDR <= 0.05)

BP_direct_select <- c('apoptotic process','negative regulation of gene expression',
                      'immune response','positive regulation of interferon-gamma production',
'inflammatory response','positive regulation of I-kappaB kinase/NF-kappaB signaling')
table1_DE_data <- table1_DE_data[table1_DE_data$Term %in% BP_direct_select,]
table1_DE_data %>%
  mutate(Term = fct_reorder(Term, Count)) %>%
  ggplot( aes(x=Term, y=Count, fill = FDR)) +
  geom_bar(stat="identity", colour="black") +
  coord_flip() +
  # ggtitle('metabolic and basic activities') +
  xlab(expression(paste('GO BP terms of LTBI',symbol('\255')))) + ylab('#genes') + 
  theme_bw() +
  scale_fill_gradientn(colors = c('khaki','#984EA3'),limits = c(0,0.05)) +
  scale_x_discrete(labels = wrap_format(35)) +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 10),
        plot.title = element_text(size = 12, face = "bold"))
# theme(axis.text = element_text(size = 20),plot.title = element_text(size = 25, face = "bold"))
dev.print(pdf, paste('TB_RSTR_PP1_CD4_lognorm_cleaned_DE_LTBI_DAVID_GO_BF_DIRECT_selected.pdf',sep = ''),width = 4.5, height = 2.5)


#### GSEA ########################
library(msigdbr)
# msigdbr_species()
# H: hallmark gene sets
# C1: positional gene sets
# C2: curated gene sets
# C3: motif gene sets
# C4: computational gene sets
# C5: GO gene sets
# C6: oncogenic signatures
# C7: immunologic signatures
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)

CD4_RSTR_LTBI.markers <- CD4_RSTR_LTBI.markers %>% arrange(desc(avg_log2FC))

CD4_RSTR_LTBI.markers_filtered <- CD4_RSTR_LTBI.markers[!is.na(CD4_RSTR_LTBI.markers$geneID),]
geneList <- CD4_RSTR_LTBI.markers_filtered$avg_log2FC
names(geneList) <- CD4_RSTR_LTBI.markers_filtered$geneID
em2 <- GSEA(geneList, TERM2GENE = m_t2g)
temp_output <- em2@result
rownames(temp_output) <- NULL
temp_output$Description <- NULL
all_gene_table <- CD4_RSTR_LTBI.markers
all_gene_table <- all_gene_table[!is.na(all_gene_table$geneID),][,c('gene','geneID')]
all_gene_table$geneID <- as.numeric(all_gene_table$geneID)
rownames(all_gene_table) <- all_gene_table$geneID
all_gene_table <- all_gene_table %>% arrange(desc(geneID))
for (geneID_index in rownames(all_gene_table)) {
  temp_output$core_enrichment <- gsub(paste('\\/',geneID_index,'\\/',sep = ''),paste('\\/',all_gene_table[geneID_index,][,c('gene')],'\\/',sep = ''),temp_output$core_enrichment)
  temp_output$core_enrichment <- gsub(paste('^',geneID_index,'\\/',sep = ''),paste(all_gene_table[geneID_index,][,c('gene')],'\\/',sep = ''),temp_output$core_enrichment)
  temp_output$core_enrichment <- gsub(paste('\\/',geneID_index,'$',sep = ''),paste('\\/',all_gene_table[geneID_index,][,c('gene')],sep = ''),temp_output$core_enrichment)
  }
temp_output$core_enrichment <- gsub('\\/',', ',temp_output$core_enrichment)
temp_output <- temp_output %>% arrange(desc(NES))
write.xlsx(temp_output,'TB_RSTR_PP1_CD4_lognorm_cleaned_DE_RSTR_vs_LTBI_GSEA_msigdbr_C5.xlsx')


########### gene visualization violin ###############################################
RNAseq_dataframe_log1p <- data.frame(t(as.matrix(TB_shrink.seurat@assays[["RNA"]]@data)),check.names = F)
RNAseq_dataframe_log1p$Group <- TB_shrink.seurat$Group

library(reshape2)
# generating memory figure
temp <- RNAseq_dataframe_log1p[,c('CCR7','TCF7','FOXP1','Group')]
temp2 <- melt(temp, id.vars = c("Group"), variable.name = "genes", value.name = "expression")
temp2$Group <- factor(temp2$Group,levels = c('RSTR','LTBI'))
temp2$genes <- factor(temp2$genes, levels = c('CCR7','TCF7','FOXP1'))
temp1 <- data.frame(genes = c('CCR7','TCF7','FOXP1'),label = c('***','**','*')) 
temp1$genes <- factor(temp1$genes,levels = c('CCR7','TCF7','FOXP1'))
ggplot(temp2, aes(y = expression, x = Group)) +
  geom_violin(aes(fill = Group),position = position_dodge(width = 1),scale = 'width') +
  geom_point(aes(fill = Group),position = position_jitterdodge(dodge.width = 0.9),size = 0.01) +
  scale_fill_manual(values = cols_group) + theme_bw() + theme(panel.grid = element_blank(), 
                                                              strip.text = element_text(face = "bold"), strip.background = element_rect(fill = NA, color = NA), 
                                                              # axis.text = element_text(color = "black"), 
                                                              # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
                                                              # legend.key.height = unit(0.8, "lines"),
                                                              text = element_text(size = 13),plot.title = element_text(size = 13, face = "bold")) +
  xlab('') + ylab('log(counts)') + RotatedAxis() + ylim(0,10.5) + 
  geom_text(data = temp1, aes(x = 1.5,y = 10, label = label),size = 4) +
  facet_wrap(~genes, scales = 'fixed')
dev.print(pdf, paste('TB_RSTR_PP1_CD4_lognorm_cleaned_violin_Group_differentiation_panel.pdf',sep = ''),width = 4, height = 2.5)

# generating activation figure
gene_list <- c('FAS','IL2RB','TNFRSF4','ICOS','CTLA4','IKZF1')
temp <- RNAseq_dataframe_log1p[,c(gene_list,'Group')]
temp2 <- melt(temp, id.vars = c("Group"), variable.name = "genes", value.name = "expression")
temp2$Group <- factor(temp2$Group,levels = c('RSTR','LTBI'))
temp2$genes <- factor(temp2$genes, levels = gene_list)
temp1 <- data.frame(genes = gene_list,label = c('***','***','***','***','**','***')) 
temp1$genes <- factor(temp1$genes,levels = gene_list)

ggplot(temp2, aes(y = expression, x = Group)) +
  geom_violin(aes(fill = Group),position = position_dodge(width = 1),scale = 'width') +
  # geom_jitter(position=position_jitter(0.01)) +
  geom_point(aes(fill = Group),position = position_jitterdodge(dodge.width = 0.9),size = 0.01) +
  scale_fill_manual(values = cols_group) + theme_bw() + theme(panel.grid = element_blank(), 
                                                              strip.text = element_text(face = "bold"), strip.background = element_rect(fill = NA, color = NA), 
                                                              # axis.text = element_text(color = "black"), 
                                                              # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
                                                              # legend.key.height = unit(0.8, "lines"),
                                                              text = element_text(size = 13),plot.title = element_text(size = 13, face = "bold")) +
  xlab('') + ylab('log(counts)') + RotatedAxis() + ylim(0,10.5) + 
  geom_text(data = temp1, aes(x = 1.5,y = 10, label = label),size = 4) +
  facet_wrap(~genes, scales = 'fixed',ncol = 2) #+ scale_x_continuous(limits=c(10,35))# + scale_y_continuous(limits=c(0,400))
dev.print(pdf, paste('TB_RSTR_PP1_CD4_lognorm_cleaned_scale_violin_Group_actviation_panel.pdf',sep = ''),width = 4, height = 5)

# Th17 figure
gene_list <- c('BATF','RORA','TBX21','RORC')
temp <- RNAseq_dataframe_log1p[,c(gene_list,'Group')]
temp2 <- melt(temp, id.vars = c("Group"), variable.name = "genes", value.name = "expression")
temp2$Group <- factor(temp2$Group,levels = c('RSTR','LTBI'))
temp2$genes <- factor(temp2$genes, levels = gene_list)
temp1 <- data.frame(genes = gene_list,label = c('***','***','***','***')) 
temp1$genes <- factor(temp1$genes,levels = gene_list)

ggplot(temp2, aes(y = expression, x = Group)) +
  geom_violin(aes(fill = Group),position = position_dodge(width = 1),scale = 'width') +
  # geom_jitter(position=position_jitter(0.01)) +
  geom_point(aes(fill = Group),position = position_jitterdodge(dodge.width = 0.9),size = 0.01) +
  scale_fill_manual(values = cols_group) + theme_bw() + theme(panel.grid = element_blank(), 
                                                              strip.text = element_text(face = "bold"), strip.background = element_rect(fill = NA, color = NA), 
                                                              # axis.text = element_text(color = "black"), 
                                                              # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
                                                              # legend.key.height = unit(0.8, "lines"),
                                                              text = element_text(size = 13),plot.title = element_text(size = 13, face = "bold")) +
  xlab('') + ylab('log(counts)') + ylim(0,10.5) + 
  geom_text(data = temp1, aes(x = 1.5,y = 10, label = label),size = 4) +
  facet_wrap(~genes, scales = 'fixed',ncol = 2) #+ scale_x_continuous(limits=c(10,35))# + scale_y_continuous(limits=c(0,400))
dev.print(pdf, paste('TB_RSTR_PP1_CD4_lognorm_cleaned_scale_violin_Group_Th17_panel.pdf',sep = ''),width = 4, height = 4)

# Flynn_NHM_T1T17pop1
gene_list <- c('TXNIP','CCR6','PDE4D','FYB')
temp <- RNAseq_dataframe_log1p[,c(gene_list,'Group')]
temp2 <- melt(temp, id.vars = c("Group"), variable.name = "genes", value.name = "expression")
temp2$Group <- factor(temp2$Group,levels = c('RSTR','LTBI'))
temp2$genes <- factor(temp2$genes, levels = gene_list)
temp1 <- data.frame(genes = gene_list,label = c('***','***','***','**')) 
temp1$genes <- factor(temp1$genes,levels = gene_list)

ggplot(temp2, aes(y = expression, x = Group)) +
  geom_violin(aes(fill = Group),position = position_dodge(width = 1),scale = 'width') +
  # geom_jitter(position=position_jitter(0.01)) +
  geom_point(aes(fill = Group),position = position_jitterdodge(dodge.width = 0.9),size = 0.01) +
  scale_fill_manual(values = cols_group) + theme_bw() + theme(panel.grid = element_blank(), 
                                                              strip.text = element_text(face = "bold"), strip.background = element_rect(fill = NA, color = NA), 
                                                              # axis.text = element_text(color = "black"), 
                                                              # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
                                                              # legend.key.height = unit(0.8, "lines"),
                                                              text = element_text(size = 13),plot.title = element_text(size = 13, face = "bold")) +
  xlab('') + ylab('log(counts)') + ylim(0,10.5) + 
  # annotate("text", x = 1.5, y = 10,label = c('***','***','**'),size = 5) +
  geom_text(data = temp1, aes(x = 1.5,y = 10, label = label),size = 5) +
  facet_wrap(~genes, scales = 'fixed',ncol = 4) #+ scale_x_continuous(limits=c(10,35))# + scale_y_continuous(limits=c(0,400))
dev.print(pdf, paste('TB_RSTR_PP1_CD4_violin_Flynn_NHM_genes.pdf',sep = ''),width = 6, height = 2.5)

######## gene visualization heatmap #####################################
RNAseq_dataframe_log1p_scale <- data.frame(t(TB_shrink.seurat@assays[["RNA"]]@scale.data),check.names = F)
RNAseq_dataframe_log1p_scale$DonorID <- TB_shrink.seurat$DonorID
RNAseq_dataframe_log1p_scale$Donor <- TB_shrink.seurat$Donor
RNAseq_dataframe_log1p_scale$Group <- TB_shrink.seurat$Group
RNAseq_dataframe_log1p_scale$Group_DonorID <- TB_shrink.seurat$Group_DonorID

cytokine_global_table <- read_excel('../Global landscape of cytokines Supplementary Table S1.xlsx', sheet = "Cytokines")
cytokine_global_list <- unique(cytokine_global_table$`HGNC symbol`)
cytokine_global_list <- c(cytokine_global_list,'GZMB','PRF1','MIF')
seurat_cytokine_list <- cytokine_global_list[cytokine_global_list %in% seurat_gene_list]
seurat_costimulation_list <- c('IL2RA', 'IL7R', 'IL2RG','CD28','CD274','CD27', 'PDCD1', 'LAG3', 'KLRG1', 'B3GAT1', 'KLRB1', 'CD38', 'HLA-DRA','MKI67','IL2RB', 'FAS', 'ICOS','TNFRSF4','TNFRSF9','CD69','CD40LG')
seurat_receptor_list <- c('CD74','ITGAL','CD58','ITGA1','ITGAE','ICAM1','TCF7','LEF1','SELL','CD27','CD28','IL7R','CCR7','CCR4','CXCR3','CXCR4','CCR5','CXCR5','CCR3','PTGDR2','CCR6','CCR10','SELPLG')
seurat_exhaustion_list <- c('TCF7','PDCD1', 'CTLA4', 'TIGIT', 'HAVCR2','LAG3','NFATC1', 'IRF4', 'BATF', 'NR4A1','TOX','BTLA')
Treg_gene_list <- c('FOXP3','IL2RA','IL2RB','TNFRSF4','IKZF2','CTLA4','GPR83','CAPG','IZUMO1R','CHCHD10','TNFRSF18','DAPL1','IGFBP4')
GO_Th17_gene_list <- c('BATF','MALT1','SLAMF6','RC3H1','IL23R','TBX21','IL2','IL4','IL6','IL6R','IL12B','IL12RB1','IRF4','JUNB','LY9','MIR21','SMAD7','ZBTB7B','IL23A','RC3H2','ENTPD7','RORA','RORC','NFKBIZ','STAT3','ZC3H12A','LOXL3','NFKBID','TNFSF18','CCR4','CCR6','CXCR3')
# GO Th17 gene list remove FOXP3 but add CCR4, CCR6 and CXCR3.
GO_Th17_gene_list <- GO_Th17_gene_list[GO_Th17_gene_list %in% seurat_gene_list]

Flynn_NHM_Th1_Th17_subcluster1 <- c('IL7R','PDE4D','ZFP36L2','ITGB1','TXNIP','FGL2','FYB','MBNL1','LMNA','DAPP1','FAM46C',
                    'CCR6','ST6GAL1','S1PR1','RPS18','RGL2','ITGA6','C4H6orf62')
Flynn_NHM_stem_like <- c('PLK2','TCF7','PECAM1','FOS','TMEM123','JUNB','CCR7','TXNIP','S1PR1','IL7R','ITGB1','LEF1','SELL',
                         'IL6R','FGL2','ZFP36L2','ZFP36','CR1','CYSLTR1','ITGA6')
Flynn_NHM_metallothionine <- c('RPS18','LMNA','SPINK2','NACA2','TOB1','MAFA-F','IL7R','BCL2L11','CD52','UBE2J1','SKIL',
                               'FAM174B','GPCPD1','RPL36AL','BCAT1','RPL3')
Flynn_NHM_proliferating <- c('MKI67','STMN1','TOP2A','RRM2','HMGB2','TUBB','KIAA0101','TUBA1B','ASPM','CCNA2','KIF11','SMC2',
                             'TMPO','CENPE','BIRC5','SMC4','TPX2','NUSAP1','HMGB1')
seurat_migration_list <- c('ITGAL','ITGA1','ITGAE','ICAM1')
seurat_naive_list <- c('CCR7','TCF7','FOXP1','LEF1','SELL','CD27','IL7R')
seurat_receptor_list <- c('CD74','CD58','CD28','CCR4','CXCR4','CCR5','CXCR5','SELPLG')
temp <- CD4_RSTR_LTBI.markers[Flynn_NHM_stem_like,]
temp$category <- 'receptor'
temp$category[temp$gene %in% seurat_cytokine_list] <- 'cytokine'
temp$category[temp$gene %in% seurat_migration_list] <- 'migration'
temp$category[temp$gene %in% seurat_naive_list] <- 'differentiation'

temp <- temp[(temp$pct.1 + temp$pct.2) > 0.1,]
temp <- temp %>% arrange(desc(avg_log2FC))

list_interest <- temp$gene
list_interest <- list_interest[!is.na(list_interest)]
TB_shrink.seurat$Group <- factor(TB_shrink.seurat$Group, levels = c('LTBI','RSTR'))

cytokine_lognormmean <- data.frame(t(RNAseq_dataframe_log1p_scale %>% group_by(Donor,Group_DonorID) %>% summarize_at(vars(list_interest),funs(mean))),check.names = F)
temp <- data.frame(t(cytokine_lognormmean),check.names = F)

rownames(temp) <- temp$Donor
temp <- temp[sort(rownames(temp)),]
temp$Donor <- NULL
temp$Group_DonorID <- NULL
cytokine_lognormmean <- data.frame(t(temp),check.names = F)
cytokine_lognormmean <- mutate_all(cytokine_lognormmean, function(x) as.numeric(x))
cytokine_lognormmean$wilcoxon_p <- CD4_RSTR_LTBI.markers[list_interest,]$p_val_adj
cytokine_lognormmean$RSTR_pct <- CD4_RSTR_LTBI.markers[list_interest,]$pct.1
cytokine_lognormmean$LTBI_pct <- CD4_RSTR_LTBI.markers[list_interest,]$pct.2
cytokine_lognormmean$pct_diff <- cytokine_lognormmean$RSTR_pct - cytokine_lognormmean$LTBI_pct
cytokine_lognormmean$ttest_p <- 0
cytokine_lognormmean$mean_diff <- 0
for (gene_name in list_interest){
    ttest_result <- t.test(cytokine_lognormmean[,c('1','2','3')][rownames(cytokine_lognormmean) == gene_name,],
                           cytokine_lognormmean[,c('4','5','6','7')][rownames(cytokine_lognormmean) == gene_name,])
    cytokine_lognormmean$ttest_p[rownames(cytokine_lognormmean) == gene_name] <- ttest_result$p.value
    cytokine_lognormmean$mean_diff[rownames(cytokine_lognormmean) == gene_name] <- ttest_result[["estimate"]][["mean of x"]] - ttest_result[["estimate"]][["mean of y"]]
}
cytokine_lognormmean$ttest_p_adj <- p.adjust(cytokine_lognormmean$ttest_p,method = 'fdr')

# cytokine_lognormmean_selected <- cytokine_lognormmean[(cytokine_lognormmean$RSTR_pct >= 0.1) | (cytokine_lognormmean$RSTR_pct - cytokine_lognormmean$LTBI_pct > 0.05) | (cytokine_lognormmean$wilcoxon_p <= 0.05),]
# cytokine_lognormmean_selected <- cytokine_lognormmean[rowMaxs(as.matrix(cytokine_lognormmean[,cytokine_lognormmean_colnames])) >= 1,]
cytokine_lognormmean_selected <- cytokine_lognormmean#[(cytokine_lognormmean$LTBI_pct > 0.05) | (cytokine_lognormmean$RSTR_pct > 0.05),]

# cytokine_lognormmean_notselected <- cytokine_lognormmean[!rowMaxs(as.matrix(cytokine_lognormmean[,cytokine_lognormmean_colnames])) >= 0.5,]
# cytokine_select_rownames <- cytokine_lognormmean_selected$gene
cytokine_lognormmean_matrix <- cytokine_lognormmean_selected[,c('1','2','3','4','5','6','7')]
# rownames(cytokine_lognormmean_matrix) <- cytokine_select_rownames
cytokine_lognormmean_matrix <- as.matrix(t(cytokine_lognormmean_matrix))

# the vertical layout
library(circlize)
# col_fun = colorRamp2(c(8, 4, 0), c("red", "white", "blue"))
col_fun = colorRamp2(c(1, 0, -1), c("red", "white", "blue"))
# col_fun = colorRamp2(c(1, 0, -1), c("yellow", "black", "purple"))

ha_row = HeatmapAnnotation(group = c('RSTR','RSTR','RSTR','LTBI','LTBI','LTBI','LTBI'),
                           col = list(group = cols_group), which = 'column')
column_split <- data.frame(c('LTBI','LTBI','LTBI','RSTR','RSTR','RSTR','RSTR'))
# row_split <- data.frame(c(rep('cytokine',22),rep('differentiation',7),rep('migration',4),rep('receptor',6)))

Heatmap(t(cytokine_lognormmean_matrix), cluster_columns = F,col = col_fun, cluster_rows = F, 
        top_annotation = ha_row,name = "mean\nnormalized\nexpression", show_column_names = F,
        # row_split = row_split, 
        show_row_names = T,
        row_title = "genes",
        row_gap = unit(3, "mm"),
        column_split = column_split,
        column_title = "RSTR    LTBI   ", column_title_side = "top",row_title_side = "left",row_names_side = 'left')
dev.print(pdf, 'TB_RSTR_PP1_CD4_lognormmean_donor_heatmap_Flynn_stem_like.pdf',width = 4, height = 5)

########## module calculation  at single cell level ####################
Flynn_NHM_Th1_Th17_subcluster1 <- c('IL7R','PDE4D','ZFP36L2','ITGB1','TXNIP','FGL2','FYB','MBNL1','LMNA','DAPP1','FAM46C',
                                    'CCR6','ST6GAL1','S1PR1','RPS18','RGL2','ITGA6','C4H6orf62')
Flynn_NHM_stem_like <- c('PLK2','TCF7','PECAM1','FOS','TMEM123','JUNB','CCR7','TXNIP','S1PR1','IL7R','ITGB1','LEF1','SELL',
                         'IL6R','FGL2','ZFP36L2','ZFP36','CR1','CYSLTR1','ITGA6')
temp <- CD4_RSTR_LTBI.markers[Flynn_NHM_stem_like,]
temp <- temp[(temp$pct.1 + temp$pct.2) > 0.1,]
list_interest <- temp$gene
list_interest <- list_interest[!is.na(list_interest)]

TB_shrink.seurat <- AddModuleScore(TB_shrink.seurat,features = list(list_interest),name = 'Flynn_NHM_T1_T17_pop1')
wilcox_test <- wilcox.test(TB_shrink.seurat$Flynn_NHM_T1_T17_pop11[TB_shrink.seurat$Group == 'RSTR'],
                           TB_shrink.seurat$Flynn_NHM_T1_T17_pop11[TB_shrink.seurat$Group == 'LTBI'])

temp2 <- data.frame(stem_score = TB_shrink.seurat$Flynn_NHM_T1_T17_pop11,Group = TB_shrink.seurat$Group)
temp2$Group <- factor(temp2$Group,levels = c('RSTR','LTBI'))

ggplot(temp2, aes(y = stem_score, x = Group)) +
  geom_violin(aes(fill = Group),position = position_dodge(width = 1),scale = 'width') +
  geom_point(aes(fill = Group),position = position_jitterdodge(dodge.width = 0.9),size = 0.01) +
  scale_fill_manual(values = cols_group) + theme_bw() + theme(panel.grid = element_blank(), 
                                                              strip.text = element_text(face = "bold"), strip.background = element_rect(fill = NA, color = NA), 
                                                              # axis.text = element_text(color = "black"), 
                                                              # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
                                                              # legend.key.height = unit(0.8, "lines"),
                                                              text = element_text(size = 13),plot.title = element_text(size = 13, face = "bold")) +
  xlab('') + ylab('expression level') + ylim(-2.5,4.5) + ggtitle('T1-T17 pop1 cell score') +
  annotate("text", x = 1.5,y=4,label = 'p = 3.1E-5',size = 5)
dev.print(pdf, paste('TB_RSTR_PP1_CD4_violin_Flynn_NHM_T1T17pop1_score.pdf',sep = ''),width = 4, height = 4)

#### stringdb network ##########################
# BiocManager::install("STRINGdb")browseVignettes("STRINGdb")
# browseVignettes("STRINGdb")

# install.packages('STRINGdb_2.6.1.tar.gz', repos = NULL, type="source")

library(STRINGdb)
immune_reference_table <- read_excel('GO0002376_immune_system_process.xlsx', sheet = "Sheet1",col_names = F)
GO0002376_immune_system_process <- toupper(unique(immune_reference_table$...1))
CD4_RSTR_LTBI.markers_immune <- CD4_RSTR_LTBI.markers[toupper(CD4_RSTR_LTBI.markers$gene) %in% GO0002376_immune_system_process,]

string_db <- STRINGdb$new(version="11.5", species=9606,score_threshold=200, input_directory="", protocol="http") # start instantiating the STRINGdb reference class. 
# In the constructor of the class you can also define the STRING version to be used and a threshold for the combined scores of the interactions, 
# such that any interaction below that threshold is not loaded in the object (by default the score threshold is set to 400)
stringent_p_cutoff <- .05#1e-7#
CD4_RSTR_LTBI.markers_immune_sig <- CD4_RSTR_LTBI.markers_immune %>% dplyr::filter(p_val_adj <= stringent_p_cutoff)

CD4_RSTR_LTBI.markers_immune_sig <- CD4_RSTR_LTBI.markers_immune_sig %>% dplyr::filter(avg_log2FC >= logFC_cutoff)  %>% dplyr::select(gene,geneID,avg_log2FC,p_val_adj)# %>% arrange(p_val_adj) #
# CD4_gene_lognormfraction_gene_list <- CD4_gene_lognormfraction$gene[(CD4_gene_lognormfraction$pos_diff > 0) & (CD4_gene_lognormfraction$ttest_p <= 0.05)]
# CD4_gene_lognormmean_gene_list <- CD4_gene_lognormmean$gene[(CD4_gene_lognormmean$mean_diff > 0) & (CD4_gene_lognormmean$ttest_p <= 0.05)]
# CD4_RSTR_LTBI.markers_immune_sig <- CD4_RSTR_LTBI.markers_immune_sig[CD4_RSTR_LTBI.markers_immune_sig$gene %in% intersect(CD4_gene_lognormmean_gene_list,CD4_gene_lognormfraction_gene_list),]
CD4_RSTR_LTBI.markers_immune_sig <- data.frame(CD4_RSTR_LTBI.markers_immune_sig)

stringdb_significant_RSTR_mapped <- string_db$map(CD4_RSTR_LTBI.markers_immune_sig, "gene", removeUnmappedRows = TRUE)
# stringdb_significant_RSTR_mapped <- stringdb_significant_RSTR_mapped[stringdb_significant_RSTR_mapped$STRING_id != '9606.ENSP00000383715',]
# stringdb_significant_RSTR_mapped <- stringdb_significant_RSTR_mapped[stringdb_significant_RSTR_mapped$STRING_id != '9606.ENSP00000349967',]

hits <- stringdb_significant_RSTR_mapped$STRING_id#[1:200]
string_db$plot_network(hits)
# filter by p-value and add a color column
# (i.e. green down-regulated gened and red for up-regulated genes)
selected_gene_set1_org_mapped_pval05 <- string_db$add_diff_exp_color(stringdb_significant_RSTR_mapped,logFcColStr="avg_log2FC" )
# post payload information to the STRING server
payload_id <- string_db$post_payload( selected_gene_set1_org_mapped_pval05$STRING_id,colors=selected_gene_set1_org_mapped_pval05$color )
# display a STRING network pdf with the "halo"
string_db$plot_network(hits, payload_id=payload_id )
dev.print(pdf, paste('TB_RSTR_PP1_CD4_lognorm_cleaned_PP1_DE_RSTR_vs_LTBI_upregulated_stringdb.pdf',sep = ''),width = 10, height = 9)

# string_db$plot_network(stringdb_significant_RSTR_mapped$STRING_id)
clustersList_significant_RSTR <- string_db$get_clusters(stringdb_significant_RSTR_mapped$STRING_id)
for (cluster_name in c(1:3)){
  graphics.off()
  mapped_table <- stringdb_significant_RSTR_mapped[stringdb_significant_RSTR_mapped$STRING_id %in% clustersList_significant_RSTR[[cluster_name]],]
  input_graph <- string_db$get_subnetwork(mapped_table$STRING_id)
  Degree <- igraph::degree(input_graph)
  Eig <- evcent(input_graph)$vector
  Hub <- hub.score(input_graph)$vector
  Authority <- authority.score(input_graph)$vector
  Closeness <- closeness(input_graph)
  Reach_2 <- (ego_size(input_graph, 2)-1)/(vcount(input_graph)-1)
  ## Reach at k=3
  Reach_3 <- (ego_size(input_graph, 3)-1)/(vcount(input_graph)-1)
  Betweenness <- betweenness.estimate(input_graph,cutoff=3)
  centralities_dataframe <- data.frame(cbind(Degree, Eig, Hub, Authority, Closeness, Reach_2, Reach_3, Betweenness))
  centralities_dataframe$string_id <- rownames(centralities_dataframe)
  p <- match(centralities_dataframe$string_id,mapped_table$STRING_id)
  centralities_dataframe$gene <- mapped_table$gene[p]
  centralities_dataframe$avg_log2FC <- mapped_table$avg_log2FC[p]
  centralities_dataframe$p_val_adj <- mapped_table$p_val_adj[p]
  write.xlsx(centralities_dataframe %>% arrange(desc(Hub)),paste('TB_RSTR_PP1_CD4_lognorm_cleaned_RSTR_immune_only_stringdb_',stringent_p_cutoff,'_cluster',cluster_name,'.xlsx',sep = ''))
  
  lay=layout_with_fr(input_graph)
  stringdb_id_list <- V(input_graph)$name
  p <- match(stringdb_id_list,mapped_table$STRING_id)
  gene_list <- mapped_table$gene[p]
  input_graph_genename <- set.vertex.attribute(input_graph, "name", value=gene_list)

  p <- match(V(input_graph_genename)$name,mapped_table$gene)
  mapped_table <- mapped_table[p,]
  V(input_graph_genename)$color <- mapped_table$avg_log2FC
  input_graph_genename$palette=rev(colorRampPalette(c("red","white"))(max(mapped_table$avg_log2FC)))
  # plot <- plot(input_graph_genename,vertex.size=Hub/sum(Hub)*500,layout=lay,rescale=F,main = paste('significantly upregulated immune genes in RSTR cluster',cluster_name,'\nnode size: Hub\nnode color: log2FC',sep = ''))

  plot <- plot(input_graph_genename,vertex.size=Hub/sum(Hub)*500,layout=lay,xlim=range(lay[,1]),ylim=range(lay[,2]),rescale=F,main = paste('significantly upregulated immune genes in RSTR cluster',cluster_name,'\nnode size: Hub\nnode color: log2FC',sep = ''))
  print(plot)
  dev.print(pdf, paste('TB_RSTR_PP1_CD4_lognorm_cleaned_PP1_RSTR_immune_only_stringdb_cluster',cluster_name,'_Hub_log2FC.pdf',sep = ''),width = 11, height = 10)
}


###### SCENIC gene regulatory network analysis ######################
vignette("SCENIC_Running")
exprMat <- as.matrix(TB_shrink.seurat@assays[["RNA"]]@data)
library(SCENIC)

dbs <- defaultDbNames[["hgnc"]]
dbs['500bp'] <- 'hg19-500bp-upstream-10species.mc9nr.feather'
dbs['10kb'] <- 'hg19-tss-centered-10kb-10species.mc9nr.feather'
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir="TB_RSTR_PP1_only_CD4_stringent_cutoff_lognorm_Seurat/SCENIC/database", 
                                  dbs = dbs,
                                  nCores=8) #library(parallel) #detectCores()
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]

# add into loom 
library(SCopeLoomR)
library(SCENIC)
library(loomR)

# https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scenic.html
add_cell_annotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}


loom <- build_loom("TB_RSTR_CD4_lognorm_FilterLQCells.loom", dgem=exprMat_filtered)
cellInfo <- TB_shrink.seurat@meta.data
cellInfo$CellID <- NULL

loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

#move to the terminal, open the virtual environment

#### cell cycle ###########
s.genes <- cc.genes$s.genes
s.genes <- s.genes[s.genes %in% seurat_gene_list]
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- g2m.genes[g2m.genes %in% seurat_gene_list]
TB_shrink.seurat <- CellCycleScoring(TB_shrink.seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
TB_shrink.seurat$Phase
DimPlot(TB_shrink.seurat, group.by = "Phase") 
dev.print(pdf, paste('TB_RSTR_PP1_CD4_lognorm_Cycle_Phase.pdf',sep = ''),width = 4.3, height = 4.5)

temp1 <- dplyr::count(TB_shrink.seurat@meta.data,Group,Phase)#,DonorID)
temp1 <- temp1 %>% group_by(Group) %>% dplyr::mutate(total = sum(n))
temp1$percentage <- temp1$n/temp1$total*100
# Filling zero values for the condition that has no counts!!
for (donor_name in unique(temp1$Group)){
  temp_Tcount <- temp1[(temp1$Group == donor_name),][1,]
  for (temp2 in c('S','G1','G2M')){
    if_row <- (temp1$Group == donor_name)
    if (sum(if_row) == 0){
      temp_Tcount$Phase <- temp_Tcount
      temp_Tcount$n <- 0
      temp_Tcount$percentage <- 0
      temp1 <- rbind(temp1,temp_Tcount)
    }
  }  
}

cols_TCR_specificity <- c("black","bisque4","grey")
names(cols_TCR_specificity) <- c('G2M','S','G1')
temp1$Group <- factor(temp1$Group, levels = c('RSTR','LTBI'))
ggplot(temp1,aes(x=Group, y=percentage,fill = Phase)) +  ggtitle('Cell Cycle') +
  geom_bar(position="stack", stat="identity") +
  theme(text = element_text(size = 15),plot.title = element_text(size = 15, face = "bold")) + RotatedAxis() +
  ylab('activated CD4 cells (%)') + xlab('Group') +
  scale_fill_manual(values = cols_TCR_specificity)
dev.print(pdf, paste('TB_RSTR_PP1_CD4_Cycle_Phase_fraction_Group.pdf',sep = ''),width = 3.5, height = 3.5)

temp1 <- dplyr::count(TB_shrink.seurat@meta.data,Group,Phase,DonorID)
temp1 <- temp1 %>% group_by(Group,DonorID) %>% dplyr::mutate(total = sum(n))
temp1$percentage <- temp1$n/temp1$total*100

phase_name <- 'G2M'
wilcox_result <- wilcox.test(temp1$percentage[(temp1$Phase == phase_name) & (temp1$Group == 'RSTR')], temp1$percentage[(temp1$Phase == phase_name) & (temp1$Group == 'LTBI')])
temp1$Group <- factor(temp1$Group, levels = c('RSTR','LTBI'))
ggplot(temp1[temp1$Phase == phase_name,],aes(x=Group, y=percentage)) +
  geom_point(shape = 2) +
  theme(text = element_text(size = 15),plot.title = element_text(size = 15, face = "bold")) + RotatedAxis() +
  ylab('activated CD4 cells (%)') + xlab('Group') + 
  ggtitle(paste('Cell Cycle',phase_name,'phase\nwilcoxon p:',format(wilcox_result$p.value,digits = 3)))
dev.print(pdf, paste('TB_RSTR_PP1_CD4_Cycle_Phase_fraction_Group_',phase_name,'.pdf',sep = ''),width = 3, height = 3.5)

