library(dplyr)
library(ggplot2)
library(csv)
library(stringr)
library("readxl")
library(flowCore)
library(CATALYST)
library(cowplot)
library(FlowSOM)
library(xlsx)
library(Seurat)
library(ggridges)

ESAT6Table_cleaned <- read.csv('TB_RSTR_PP1_index_sort_PCR_TCR_table.csv',check.names = F)
flow_channels_full <- colnames(ESAT6Table_cleaned)[grepl('idx',colnames(ESAT6Table_cleaned))]
ESAT6Table_cleaned <- ESAT6Table_cleaned[!is.na(ESAT6Table_cleaned[,flow_channels_full[1]]),]
ESAT6Table_cleaned$donorGrp[ESAT6Table_cleaned$donorGrp == 'LTB'] <- 'LTBI'
ESAT6Table_cleaned$`age at study baseline` <- str_pad(ESAT6Table_cleaned$`age at study baseline`, 2, pad = "0")
ESAT6Table_cleaned$donorGrp_AgeBase_Sex_donorID <- paste(ESAT6Table_cleaned$donorGrp,' ',ESAT6Table_cleaned$`age at study baseline`,'yrs ',ESAT6Table_cleaned$SEX,' ',ESAT6Table_cleaned$donorId,sep = '')

### finding TCR clones #####################
library(tidyverse)
library(pheatmap)
library(plyr)
library(tidyverse)
library(dplyr)

CDR3a_unique_table_CrossDonor <- data.frame(table(ESAT6Table_cleaned$CDR3a)) %>% arrange(desc(Freq))
CDR3a_unique_table_CrossDonor$clone_group <- 1:dim(CDR3a_unique_table_CrossDonor)[1]
CDR3b_unique_table_CrossDonor <- data.frame(table(ESAT6Table_cleaned$CDR3b)) %>% arrange(Freq)
CDR3b_unique_table_CrossDonor$clone_group <- 1:dim(CDR3b_unique_table_CrossDonor)[1]

ESAT6Table_cleaned$CDR3a_clone_number_WithinDonor <- NaN # -1 is na, 0 is not calculated
ESAT6Table_cleaned$CDR3a_clone_number_CrossDonor <- NaN # -1 is na, 0 is not calculated
ESAT6Table_cleaned$CDR3a_clone_group <- NaN
ESAT6Table_cleaned$CDR3b_clone_number_WithinDonor <- NaN # -1 is na, 0 is not calculated
ESAT6Table_cleaned$CDR3b_clone_number_CrossDonor <- NaN # -1 is na, 0 is not calculated
ESAT6Table_cleaned$CDR3b_clone_group <- NaN

CDR3a_unique_table_WithinDonor_All <- data.frame(Group_DonorID = numeric(0), stim = numeric(0), CDR3a = numeric(0), CDR3a_clone_group  = numeric(0), CDR3a_clone_number_WithinDonor  = numeric(0), CDR3a_clone_number_CrossDonor  = numeric(0))
CDR3b_unique_table_WithinDonor_All <- data.frame(Group_DonorID = numeric(0), stim = numeric(0), CDR3b = numeric(0), CDR3b_clone_group  = numeric(0), CDR3b_clone_number_WithinDonor  = numeric(0), CDR3b_clone_number_CrossDonor  = numeric(0))

CDR3a_clone_index <- 1
CrossDonor_CDR3b_clone_number <- 1
DonorID_list <- unique(ESAT6Table_cleaned$Group_DonorID)

for (donorID in DonorID_list){
  print(donorID)
  PerDonorTable <- ESAT6Table_cleaned[ESAT6Table_cleaned$Group_DonorID == donorID,]
  # CDR3a
  CDR3a_unique_table_WithinDonor <- data.frame(table(PerDonorTable$CDR3a)) %>% arrange(desc(Freq))
  p <- match(CDR3a_unique_table_WithinDonor$Var1,CDR3a_unique_table_CrossDonor$Var1)
  temp <- CDR3a_unique_table_CrossDonor[p,]
  CDR3a_unique_table_WithinDonor$Freq_CrossDonor <- temp$Freq
  CDR3a_unique_table_WithinDonor$clone_group <- temp$clone_group
  CDR3a_unique_table_WithinDonor$Group_DonorID <- donorID
  CDR3a_unique_table_WithinDonor$stim <- PerDonorTable$stim[1]
  colnames(CDR3a_unique_table_WithinDonor) <- c('CDR3a','CDR3a_clone_number_WithinDonor','CDR3a_clone_number_CrossDonor','CDR3a_clone_group','Group_DonorID','stim')
  CDR3a_unique_table_WithinDonor <- CDR3a_unique_table_WithinDonor[,c('CDR3a','CDR3a_clone_number_CrossDonor','CDR3a_clone_group','CDR3a_clone_number_WithinDonor','Group_DonorID','stim')]
  CDR3a_unique_table_WithinDonor_All <- rbind(CDR3a_unique_table_WithinDonor_All,CDR3a_unique_table_WithinDonor)
  
  p <- match(PerDonorTable$CDR3a,CDR3a_unique_table_WithinDonor$CDR3a)
  temp <- CDR3a_unique_table_WithinDonor[p,]
  PerDonorTable$CDR3a_clone_number_WithinDonor <- temp$CDR3a_clone_number_WithinDonor
  PerDonorTable$CDR3a_clone_number_CrossDonor <- temp$CDR3a_clone_number_CrossDonor
  PerDonorTable$CDR3a_clone_group <- temp$CDR3a_clone_group
  
  
  # CDR3b
  CDR3b_unique_table_WithinDonor <- data.frame(table(PerDonorTable$CDR3b)) %>% arrange(desc(Freq))
  p <- match(CDR3b_unique_table_WithinDonor$Var1,CDR3b_unique_table_CrossDonor$Var1)
  temp <- CDR3b_unique_table_CrossDonor[p,]
  CDR3b_unique_table_WithinDonor$Freq_CrossDonor <- temp$Freq
  CDR3b_unique_table_WithinDonor$clone_group <- temp$clone_group
  CDR3b_unique_table_WithinDonor$Group_DonorID <- donorID
  CDR3b_unique_table_WithinDonor$stim <- PerDonorTable$stim[1]
  colnames(CDR3b_unique_table_WithinDonor) <- c('CDR3b','CDR3b_clone_number_WithinDonor','CDR3b_clone_number_CrossDonor','CDR3b_clone_group','Group_DonorID','stim')
  CDR3b_unique_table_WithinDonor <- CDR3b_unique_table_WithinDonor[,c('CDR3b','CDR3b_clone_number_CrossDonor','CDR3b_clone_group','CDR3b_clone_number_WithinDonor','Group_DonorID','stim')]
  CDR3b_unique_table_WithinDonor_All <- rbind(CDR3b_unique_table_WithinDonor_All,CDR3b_unique_table_WithinDonor)
  
  p <- match(PerDonorTable$CDR3b,CDR3b_unique_table_WithinDonor$CDR3b)
  temp <- CDR3b_unique_table_WithinDonor[p,]
  PerDonorTable$CDR3b_clone_number_WithinDonor <- temp$CDR3b_clone_number_WithinDonor
  PerDonorTable$CDR3b_clone_number_CrossDonor <- temp$CDR3b_clone_number_CrossDonor
  PerDonorTable$CDR3b_clone_group <- temp$CDR3b_clone_group
  
  ESAT6Table_cleaned[ESAT6Table_cleaned$Group_DonorID == donorID,] <- PerDonorTable
  
}
rm(CDR3a_unique_table_WithinDonor)
rm(CDR3b_unique_table_WithinDonor)
rm(PerDonorTable)

ESAT6Table_cleaned$CDR3a_clone_number_WithinDonor[is.na(ESAT6Table_cleaned$CDR3a_clone_number_WithinDonor)] <- 0
ESAT6Table_cleaned$CDR3a_clone_number_CrossDonor[is.na(ESAT6Table_cleaned$CDR3a_clone_number_CrossDonor)] <- 0
ESAT6Table_cleaned$CDR3a_clone_group[is.na(ESAT6Table_cleaned$CDR3a_clone_group)] <- -1
ESAT6Table_cleaned$CDR3b_clone_number_WithinDonor[is.na(ESAT6Table_cleaned$CDR3b_clone_number_WithinDonor)] <- 0
ESAT6Table_cleaned$CDR3b_clone_number_CrossDonor[is.na(ESAT6Table_cleaned$CDR3b_clone_number_CrossDonor)] <- 0
ESAT6Table_cleaned$CDR3b_clone_group[is.na(ESAT6Table_cleaned$CDR3b_clone_group)] <- -1

ESAT6Table_cleaned$if_cloned <- 0
ESAT6Table_cleaned$if_cloned[(ESAT6Table_cleaned$CDR3a_clone_number_WithinDonor > 1) | (ESAT6Table_cleaned$CDR3b_clone_number_WithinDonor > 1)] <- 1

############ quantify cloanl expansion proportion #######################
ESAT6Table_cleaned$if_CD4 <- ESAT6Table_cleaned$`idx_CD8_UV515-A`<1000
ESAT6Table_cleaned <- ESAT6Table_cleaned[,c('donorId','if_CD4','if_cloned','stim','donorGrp','CDR3a','CDR3a_clone_group','CDR3a_clone_number_WithinDonor','CDR3a_clone_number_CrossDonor','CDR3b','CDR3b_clone_group','CDR3b_clone_number_WithinDonor','CDR3b_clone_number_CrossDonor')]#'Va','Ja','Vb','Jb',

sum((ESAT6Table_cleaned$CDR3b_clone_number_WithinDonor > 1) | (ESAT6Table_cleaned$CDR3a_clone_number_WithinDonor > 1))
sum((ESAT6Table_cleaned$CDR3b_clone_number_WithinDonor > 1))
sum((ESAT6Table_cleaned$CDR3b_clone_number_WithinDonor == 0) & (ESAT6Table_cleaned$CDR3a_clone_number_WithinDonor == 0))

ESAT6Table_count_TCR <- unique(ESAT6Table_cleaned)
ESAT6Table_count_TCRb <- unique(ESAT6Table_count_TCR[,c('donorId','stim','donorGrp','CDR3b','CDR3b_clone_group','CDR3b_clone_number_WithinDonor','CDR3b_clone_number_CrossDonor')])
ESAT6Table_count_TCRb <- ESAT6Table_count_TCRb[ESAT6Table_count_TCRb$CDR3b_clone_number_WithinDonor > 0,]

ESAT6Table_count_TCR <- dplyr::count(ESAT6Table_cleaned,donorId,stim,donorGrp,CDR3a,CDR3b,CDR3a_clone_group,CDR3b_clone_group)
ESAT6Table_count_TCR <- ESAT6Table_count_TCR[(ESAT6Table_count_TCR$CDR3a_clone_group != -1) | (ESAT6Table_count_TCR$CDR3b_clone_group != -1),]
ESAT6Table_count_TCR_noNA <- ESAT6Table_count_TCR[(ESAT6Table_count_TCR$CDR3a_clone_group != -1) & (ESAT6Table_count_TCR$CDR3b_clone_group != -1),]

RSTR_table <- ESAT6Table_count_TCRb[ESAT6Table_count_TCRb$stim == 'ESAT-6/CFP-10 6 hrs',]

cols_group <- c('#984EA3','#4DAF4A')
names(cols_group) <- c('RSTR','LTBI')

RSTR_table$Donor <- ''
RSTR_table$Donor[RSTR_table$donorId == '92527-1-02'] <- '1'
RSTR_table$Donor[RSTR_table$donorId == '93506-1-02'] <- '2'
RSTR_table$Donor[RSTR_table$donorId == '93774-1-05'] <- '3'
RSTR_table$Donor[RSTR_table$donorId == '84165-1-06'] <- '4'
RSTR_table$Donor[RSTR_table$donorId == '91512-1-05'] <- '5'
RSTR_table$Donor[RSTR_table$donorId == '92422-1-03'] <- '6'
RSTR_table$Donor[RSTR_table$donorId == '93334-1-02'] <- '7'

RSTR_table$donorGrp <- factor(RSTR_table$donorGrp, levels = c('RSTR','LTBI'))
ggplot(RSTR_table, aes(Donor, log2(RSTR_table$CDR3b_clone_number_WithinDonor), color=donorGrp)) +
  geom_jitter(size=log2(RSTR_table$CDR3b_clone_number_WithinDonor), position=position_jitter(width=.2)) +
  scale_y_continuous(name="log2(TCRb clone size)") +#, limits=c(-0.2, 3)
  facet_grid(~ donorGrp, scales = "free") +
  theme_bw(base_size=16) + 
  scale_color_manual(values = cols_group) + 
  theme(text = element_text(size = 20),plot.title = element_text(size = 30, face = "bold"), legend.position = "none")
dev.print(pdf, paste('TB_RSTR_PP1_Tcell_TCRb_clone_size.pdf',sep = ''),width = 5, height = 4)

ESAT6Table_cleaned$TCRb_if_cloned <- 0
ESAT6Table_cleaned$TCRb_if_cloned[ESAT6Table_cleaned$CDR3b_clone_number_WithinDonor > 1] <- 1
ESAT6Table_cleaned$TCRb_if_cloned[ESAT6Table_cleaned$CDR3b_clone_number_WithinDonor == 0] <- -1

RSTR_CD4_CD8_count <- dplyr::count(ESAT6Table_cleaned,stim,donorGrp,donorId,if_CD4,TCRb_if_cloned)
RSTR_CD4_CD8_count <- RSTR_CD4_CD8_count %>% dplyr::group_by(stim, donorGrp, donorId,if_CD4) %>% dplyr::mutate(total = sum(n))
RSTR_CD4_CD8_count$percentage <- RSTR_CD4_CD8_count$n/RSTR_CD4_CD8_count$total*100
# Filling zero values for the condition that has no counts!!
for (donor_name in unique(RSTR_CD4_CD8_count$donorId)){
  # print(donor_name)
  temp_Tcount <- RSTR_CD4_CD8_count[(RSTR_CD4_CD8_count$donorId == donor_name),][1,]
  for (cluster_name in c(1,0)){
    for (temp_clone in c(-1,0,1)){
      if_row <- ((RSTR_CD4_CD8_count$donorId == donor_name) & (RSTR_CD4_CD8_count$if_CD4 == cluster_name) & (RSTR_CD4_CD8_count$TCRb_if_cloned == temp_clone))
      if (sum(if_row) == 0){
        temp_Tcount$if_CD4 <- cluster_name
        temp_Tcount$TCRb_if_cloned <- temp_clone
        temp_Tcount$n <- 0
        temp_Tcount$percentage <- 0
        RSTR_CD4_CD8_count <- rbind(RSTR_CD4_CD8_count,temp_Tcount)
      }
    }
  }
}

cols_group <- c('#984EA3','#4DAF4A')
names(cols_group) <- c('RSTR','LTBI')

print(cluster_name)
temp_cell_cluster <- RSTR_CD4_CD8_count[(RSTR_CD4_CD8_count$if_CD4 == 1) & (RSTR_CD4_CD8_count$TCRb_if_cloned == 1) & (RSTR_CD4_CD8_count$stim == 'ESAT-6/CFP-10 6 hrs'),]
ttest_result1 <- wilcox.test(temp_cell_cluster$percentage[(temp_cell_cluster$donorGrp =='LTBI')],
                             temp_cell_cluster$percentage[(temp_cell_cluster$donorGrp =='RSTR')])
temp_cell_cluster$donorGrp <- factor(temp_cell_cluster$donorGrp,levels = c('RSTR','LTBI'))
ggplot(temp_cell_cluster,aes(x=donorGrp, y=percentage,color = donorGrp,group = donorGrp)) +  ggtitle(cluster_name) +
  geom_boxplot() +
  geom_jitter(size=1,shape = 9) +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) + RotatedAxis() +
  ylab(paste('activated CD4 cells %')) +
  scale_colour_manual(values = cols_group) +
  ggtitle('TCRb clonal expansion') + 
  theme_bw() +
  ylim(0,60)

dev.print(pdf, paste('TB_RSTR_PP1_CD4_TCRb_clone_ratio.pdf',sep = ''),width = 3, height = 3)

######## csv outputing #####################
ESAT6Table_cleaned <- read.csv('TB_RSTR_PP1_index_sort_PCR_TCR_table.csv',check.names = F)
flow_channels_full <- colnames(ESAT6Table_cleaned)[grepl('idx',colnames(ESAT6Table_cleaned))]
ESAT6Table_cleaned <- ESAT6Table_cleaned[!is.na(ESAT6Table_cleaned[,flow_channels_full[1]]),]
ESAT6Table_cleaned$donorGrp[ESAT6Table_cleaned$donorGrp == 'LTB'] <- 'LTBI'
ESAT6Table_cleaned$`age at study baseline` <- str_pad(ESAT6Table_cleaned$`age at study baseline`, 2, pad = "0")
ESAT6Table_cleaned$donorGrp_AgeBase_Sex_donorID <- paste(ESAT6Table_cleaned$donorGrp,' ',ESAT6Table_cleaned$`age at study baseline`,'yrs ',ESAT6Table_cleaned$SEX,' ',ESAT6Table_cleaned$donorId,sep = '')

ESAT6_flow_dataframe <- ESAT6Table_cleaned[,flow_channels_full[1:(length(flow_channels_full) - 1)]]
flow_panel_name_full <- colnames(ESAT6_flow_dataframe)
flow_panel_name_short <- gsub('_.*','',substring(flow_panel_name_full,first = 5,last = 30))
colnames(ESAT6_flow_dataframe) <- paste(flow_panel_name_short,':',flow_panel_name_short,sep = '')

ESAT6_PCR_dataframe <- ESAT6Table_cleaned[,c(43:59)]
PCR_panel_name_short <- colnames(ESAT6_PCR_dataframe)
colnames(ESAT6_PCR_dataframe) <- paste(PCR_panel_name_short,':',PCR_panel_name_short,sep = '')

ESAT6_flow_PCR_dataframe <- cbind(ESAT6_flow_dataframe,ESAT6_PCR_dataframe)
ESAT6_flow_PCR_dataframe$'if_cloned:if_cloned' <- ESAT6Table_cleaned$if_cloned
ESAT6_flow_PCR_dataframe$'cell_index:cell_index' <- 1:dim(ESAT6_flow_PCR_dataframe)[1]

# everything is shifted up to avoid negative values
min_flow_raw <- - min(ESAT6_flow_PCR_dataframe)
ESAT6_flow_PCR_dataframe_nonzero <- ESAT6_flow_PCR_dataframe + min_flow_raw

for (donorID in unique(ESAT6Table_cleaned$donorGrp_AgeBase_Sex_donorID)){
  ESAT6_flow_PCR_dataframe_nonzero_per_sample <- ESAT6_flow_PCR_dataframe_nonzero[ESAT6Table_cleaned$donorGrp_AgeBase_Sex_donorID == donorID,]
  write.csv(ESAT6_flow_PCR_dataframe_nonzero_per_sample,gsub(' ','_',paste('flow_PCR_raw_nonzero_',donorID,'.csv',sep = '')),row.names = F)
}

##### from .fcs to flowSet ###############
# converting csv2fcs using .sh and .jar tools and inputing .fcs 
# csv2fcs does not take negative values and will overspill!!!
# sce_CyTOF <- readRDS('tonsil_organoid_PMA_sce.rds')
# fcs_CyTOF <- readRDS('tonsil_organoid_PMA_fcs.rds')
fcs_list_flow_PCR_raw_nonzero <- dir('flow_PCR_raw_nonzero_042822', full.names=T, pattern="*.fcs$")
fcs_flow_PCR_raw_nonzero <- read.flowSet(fcs_list_flow_PCR_raw_nonzero, transformation = FALSE, truncate_max_range = FALSE)
# summary(fcs_flow_PCR_raw_nonzero)
panel_fcs <- pData(parameters(fcs_flow_PCR_raw_nonzero[[1]]))
all_marker_list <- as.character(panel_fcs$desc)
# shift back to the orginial data values!
fcs_flow_PCR_raw <- fsApply(fcs_flow_PCR_raw_nonzero, function(x){
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- expr[, c(panel_fcs$desc)] - min_flow_raw 
  exprs(x) <- expr
  x
})

##### creating sce object in CATALYST ###############
# meta data, must be these column names and cannot recognize others.
meta_form_flow_PCR_raw <- data.frame(file_name = paste('flow_PCR_raw_nonzero_',gsub(' ','_',unique(ESAT6Table_cleaned$donorGrp_AgeBase_Sex_donorID)),'.fcs',sep = ''))
meta_form_flow_PCR_raw$patient_id <- gsub('_',' ',gsub('flow_PCR_raw_nonzero_','',gsub('.fcs','',meta_form_flow_PCR_raw$file_name)))
meta_form_flow_PCR_raw$patient_id <- gsub('non-RSTR', 'LTBI',meta_form_flow_PCR_raw$patient_id)
meta_form_flow_PCR_raw$condition <- gsub(' .*','',meta_form_flow_PCR_raw$patient_id)
meta_form_flow_PCR_raw$sample_id <- paste(meta_form_flow_PCR_raw$condition,substring(gsub('.*yrs ','',meta_form_flow_PCR_raw$patient_id),first = 3,last = 30))
meta_form_flow_PCR_raw <- meta_form_flow_PCR_raw %>% arrange(sample_id)
# panel names
panel_flow_PCR_raw <- pData(parameters(fcs_flow_PCR_raw[[1]]))
panel_flow_PCR_raw$fcs_colname <- panel_flow_PCR_raw$desc
panel_flow_PCR_raw$antigen <- panel_flow_PCR_raw$name # was marker originally
panel_flow_PCR_raw <- panel_flow_PCR_raw[!is.na(panel_flow_PCR_raw$antigen),]
panel_flow_PCR_raw$marker_class <- 'type'
panel_flow_PCR_raw$marker_class[panel_flow_PCR_raw$name %in% c('AQUA-CD14-CD19','if_cloned','cell_index','FSC-A','SSC-A','CD3',"TCRab",'IL5')] <- 'state'
panel_input_flow_PCR_raw <- panel_flow_PCR_raw[,c('fcs_colname','antigen','marker_class')]
marker_list <- panel_flow_PCR_raw$antigen[panel_flow_PCR_raw$marker_class == 'type']

# converting to sce object
# here it can only import patient_id, sample_id and condition
# I manually add batch and Group
sce_flow_PCR_raw <- prepData(fcs_flow_PCR_raw, panel_input_flow_PCR_raw, meta_form_flow_PCR_raw, features = panel_flow_PCR_raw$fcs_colname,FACS = T,transform = T)
sce_flow_PCR_raw@colData@listData[["batch"]] <- '040919'
sce_flow_PCR_raw@colData@listData[["batch"]][sce_flow_PCR_raw@colData@listData[["patient_id"]] == 'LTBI 12yrs M 84165-1-06'] <- '040919'
sce_flow_PCR_raw@colData@listData[["batch"]][sce_flow_PCR_raw@colData@listData[["patient_id"]] %in% c('LTBI 08yrs F 91512-1-05','RSTR 12yrs M 93774-1-05')] <- '041019'
sce_flow_PCR_raw@colData@listData[["batch"]][sce_flow_PCR_raw@colData@listData[["patient_id"]] %in% c('LTBI 10yrs F 92422-1-03','RSTR 25yrs F 93506-1-02')] <- '041119'
sce_flow_PCR_raw@colData@listData[["batch"]][sce_flow_PCR_raw@colData@listData[["patient_id"]] %in% c('LTBI 29yrs F 93334-1-02','RSTR 62yrs M 92527-1-02')] <- '041219'
sce_flow_PCR_raw@colData@listData[["Group"]] <- sce_flow_PCR_raw@colData@listData[["condition"]]
sce_flow_PCR_raw@colData@listData[["cell_index"]] <- sce_flow_PCR_raw@assays@data@listData[["counts"]]['cell_index',]
sce_flow_PCR_raw@colData@listData[["clonal_expanded"]] <- sce_flow_PCR_raw@assays@data@listData[["counts"]]['if_cloned',]
sce_flow_PCR_raw@colData@listData[["clonal_expanded"]] <- gsub('1','yes',gsub('0','no',sce_flow_PCR_raw@colData@listData[["clonal_expanded"]]))

flow_PCR_raw_table_sce <- sce_flow_PCR_raw@assays@data@listData[["counts"]]
sce_flow_PCR_raw@assays@data@listData[["exprs"]] <- flow_PCR_raw_table_sce

# Here transformed data is flow arcsinh5 and count log1p
flow_plus_PCR_transformed_table_sce <- flow_PCR_raw_table_sce
flow_plus_PCR_transformed_table_sce[flow_panel_name_short,] <- asinh((flow_plus_PCR_transformed_table_sce[flow_panel_name_short,] + min_flow_raw)/5)
flow_plus_PCR_transformed_table_sce[PCR_panel_name_short,] <- log1p(flow_plus_PCR_transformed_table_sce[PCR_panel_name_short,])
sce_flow_plus_PCR_transformed <- sce_flow_PCR_raw
sce_flow_plus_PCR_transformed@assays@data@listData[["exprs"]] <- flow_plus_PCR_transformed_table_sce

flow_plus_PCR_transformed_table_sce_scale <- t(scale(t(flow_plus_PCR_transformed_table_sce)))

sce_flow_plus_PCR_transformed_scale <- sce_flow_PCR_raw
sce_flow_plus_PCR_transformed_scale@assays@data@listData[["exprs"]] <- flow_plus_PCR_transformed_table_sce_scale

flow_plus_transformed_PCR_binary_table_sce <- flow_plus_PCR_transformed_table_sce
flow_plus_transformed_PCR_binary_table_sce[PCR_panel_name_short,][flow_PCR_raw_table_sce[PCR_panel_name_short,] <= 5] <- 0
flow_plus_transformed_PCR_binary_table_sce[PCR_panel_name_short,][flow_PCR_raw_table_sce[PCR_panel_name_short,] > 5] <- 1
flow_plus_transformed_PCR_binary_table_sce_scale <- t(scale(t(flow_plus_transformed_PCR_binary_table_sce)))
flow_plus_transformed_PCR_binary_scale <- sce_flow_PCR_raw
flow_plus_transformed_PCR_binary_scale@assays@data@listData[["exprs"]] <- flow_plus_transformed_PCR_binary_table_sce_scale

########### clustering and UMAP ####################
# cluster tSNE and UMAP
flow_plus_transformed_PCR_binary_scale <- runDR(flow_plus_transformed_PCR_binary_scale, "TSNE", cells = 280, features = "type")

# To make the tSNE consistent, we read in the pre-saved rds
flow_plus_transformed_PCR_binary_scale <- readRDS('TB_RSTR_PP1_flow_PCR_SingleCellExperiment_object.rds')
library(gridExtra)
new_marker_list <- c("CD4","CD8","CD69","CD154","CD137",
                     "CD25","CD127","CD38","HLA-DR","CD45RA",
                     "CXCR3","CCR6","FOXP3","TGFB","RUNX1",
                     "TBET","IFNG","IL2","TNF","RUNX3",
                     "RORC","IL17A","GZMB","PERF","IL21",
                     "GATA3","IL13","IL4")
plots <- lapply(new_marker_list, function(.x) plotDR(flow_plus_transformed_PCR_binary_scale, "TSNE", color_by = .x) + 
         theme(text = element_text(size = 20)))

do.call(grid.arrange,plots)
dev.print(pdf, paste('TB_RSTR_PP1_only_flow_plus_transformed_PCR_binary_scale_tSNE_allmarkers.pdf',sep = ''),width = 25, height = 25)

# FlowSom clustering
n_clusters <- 40
grid_size <- 15
flow_plus_transformed_PCR_binary_scale <- cluster(flow_plus_transformed_PCR_binary_scale, xdim = grid_size, ydim = grid_size, features = "type", maxK = n_clusters,seed = 1)

plotExprHeatmap(flow_plus_transformed_PCR_binary_scale, features = "type", by = "cluster_id", k = paste("meta",n_clusters,sep = ''), bars = TRUE, perc = TRUE)
dev.print(pdf, paste('TB_RSTR_PP1_only_flow_plus_transformed_PCR_binary_scale_size',grid_size,'_cluster',n_clusters,'_heatmap.pdf',sep = ''),width = 8, height = 5)

plotDR(flow_plus_transformed_PCR_binary_scale, "TSNE", color_by = paste("meta",n_clusters,sep = '')) 
dev.print(pdf, paste('TB_RSTR_PP1_only_flow_plus_transformed_PCR_binary_scale_size',grid_size,'_cluster',n_clusters,'_tSNE.pdf',sep = ''),width = 5, height = 4)

##### cluster annotation #####################################
# To make the tSNE consistent, we read in the pre-saved rds
flow_plus_transformed_PCR_binary_scale <- readRDS('TB_RSTR_PP1_flow_PCR_SingleCellExperiment_object.rds')

merging_table2 <- data.frame(original_cluster = c(1:n_clusters), new_cluster = '')
merging_table2$new_cluster[merging_table2$original_cluster %in% c(36,37,39)] <- '1' 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(28,30)] <- '2' 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(16,20,23)] <- '3' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(35,38)] <- '4' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(40)] <- '5' 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(31)] <- '6' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(34)] <- '7' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(11,32)] <- '8' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(9)] <- '9' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(6)] <- '10' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(3,7,8,12,16)] <- '11' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(5)] <- '12' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(4)] <- '13' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(22)] <- '14'
merging_table2$new_cluster[merging_table2$original_cluster %in% c(14,15,26)] <- '15' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(33)] <- '16' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(17,25)] <- '17' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(10,13)] <- '18' # 
merging_table2$new_cluster[merging_table2$original_cluster %in% c(1,2)] <- '19'

merging_table2$new_cluster <- factor(merging_table2$new_cluster, levels = as.character(c(1:19)))

flow_plus_transformed_PCR_binary_scale <- mergeClusters(flow_plus_transformed_PCR_binary_scale, k = paste("meta",n_clusters,sep = ''), table = merging_table2, id = "T_cell_subsets",overwrite = T)

plotDR(flow_plus_transformed_PCR_binary_scale, "TSNE", color_by = 'T_cell_subsets',facet_by = 'Group', ncol = 2) + xlab('TSNE 1') + ylab('TSNE 2') +
  guides(fill=guide_legend(title="T cell\nsubsets"))
dev.print(pdf, paste('TB_RSTR_PP1_only_flow_plus_transformed_PCR_binary_scale_size',grid_size,'_cluster',n_clusters,'_annotated_subsets_TSNE_Group_vertical.pdf',sep = ''),width = 7, height = 4.5)

plotExprHeatmap(flow_plus_transformed_PCR_binary_scale, features = "type", row_clust = F, by = "cluster_id", k = "T_cell_subsets",bars = T, perc = F)
dev.print(pdf, paste('TB_RSTR_PP1_only_flow_plus_transformed_PCR_binary_scale_size',grid_size,'_cluster',n_clusters,'_heatmap_annotated_subsets.pdf',sep = ''),width = 8, height = 5)

merging_table <- data.frame(original_cluster = c(1:n_clusters), new_cluster = '')
merging_table$new_cluster[merging_table$original_cluster %in% c(36,37,39)] <- 'Treg' # the difference is CD25 and CD137
merging_table$new_cluster[merging_table$original_cluster %in% c(30,28)] <- 'RORC+IL13+IL4+CD4' #
merging_table$new_cluster[merging_table$original_cluster %in% c(16,20,23)] <- 'RORC+IL13+CD4' #
merging_table$new_cluster[merging_table$original_cluster %in% c(35,38)] <- 'RORC+GATA3+CD4' #
merging_table$new_cluster[merging_table$original_cluster %in% c(40)] <- 'RORC+GATA3+IL17A+CD4' #
merging_table$new_cluster[merging_table$original_cluster %in% c(31)] <- 'RORC+IL17A+CD4' #
merging_table$new_cluster[merging_table$original_cluster %in% c(34)] <- 'RORC+CD4' #
merging_table$new_cluster[merging_table$original_cluster %in% c(32,11)] <- 'RORC+TNF+CD4' #
merging_table$new_cluster[merging_table$original_cluster %in% c(9)] <- 'RORC+IFNG+IL2+CD4' #
merging_table$new_cluster[merging_table$original_cluster %in% c(6)] <- 'RORC+TBET+CD4' #
merging_table$new_cluster[merging_table$original_cluster %in% c(3,7,8,12,16)] <- 'RORC+IFNG+cytolytic CD4' #Here 7 is different
merging_table$new_cluster[merging_table$original_cluster %in% c(5)] <- 'TBET+CD4' #
merging_table$new_cluster[merging_table$original_cluster %in% c(4)] <- 'CXCR3+TNF+CD4' #IFNG+
merging_table$new_cluster[merging_table$original_cluster %in% c(22)] <- 'IL21+CD4'
merging_table$new_cluster[merging_table$original_cluster %in% c(14,15,26)] <- 'TNF+CD4' #
merging_table$new_cluster[merging_table$original_cluster %in% c(33)] <- 'CD137+CD4' # partially CXCR3+
merging_table$new_cluster[merging_table$original_cluster %in% c(25,17)] <- 'CD45RA-CD4' # partially CXCR3+
merging_table$new_cluster[merging_table$original_cluster %in% c(10,13)] <- 'CD45RA+CD4'
merging_table$new_cluster[merging_table$original_cluster %in% c(1,2)] <- 'CD8' # the resolution is so low for CD8 subset to be seperable!

merging_table$new_cluster <- factor(merging_table$new_cluster,
                                     levels = c("Treg",'RORC+IL13+IL4+CD4',"RORC+IL13+CD4","RORC+GATA3+CD4",'RORC+GATA3+IL17A+CD4',"RORC+IL17A+CD4","RORC+CD4",'RORC+TNF+CD4',"RORC+IFNG+IL2+CD4","RORC+TBET+CD4","RORC+IFNG+cytolytic CD4",
                                                'TBET+CD4',"CXCR3+TNF+CD4", "IL21+CD4", 'TNF+CD4','CD137+CD4',"CD45RA-CD4","CD45RA+CD4", "CD8"))

flow_plus_transformed_PCR_binary_scale <- mergeClusters(flow_plus_transformed_PCR_binary_scale, k = paste("meta",n_clusters,sep = ''), table = merging_table, id = "merging1",overwrite = T)

##### filtering only CD4 and subset visualization ####################
CD4_flow_plus_transformed_PCR_binary_scale <- filterSCE(flow_plus_transformed_PCR_binary_scale, k = "merging1",cluster_id != 'CD8')

cluster_dataframe <- data.frame(cluster = cluster_ids(CD4_flow_plus_transformed_PCR_binary_scale,'T_cell_subsets'),
                                Group = CD4_flow_plus_transformed_PCR_binary_scale@colData@listData[["Group"]],
                                sample_id = CD4_flow_plus_transformed_PCR_binary_scale@colData@listData[["sample_id"]])
cluster_dataframe <- mutate_all(cluster_dataframe, function(x) as.character(x))
cluster_group_sample_count <- dplyr::count(cluster_dataframe, cluster, Group,sample_id)
cluster_group_sample_count <- cluster_group_sample_count %>% dplyr::group_by(sample_id) %>% dplyr::mutate(total = sum(n))
cluster_group_sample_count$percentage <- cluster_group_sample_count$n/cluster_group_sample_count$total*100
# Filling zero values for the Group that has no counts!!
for (sample_name in unique(cluster_dataframe$sample_id)){
  temp_Bcount <- cluster_group_sample_count[(cluster_group_sample_count$sample_id == sample_name),][1,]
    for (cluster_name in unique(cluster_dataframe$cluster)){
      if_row <- ((cluster_group_sample_count$sample_id == sample_name) & (cluster_group_sample_count$cluster == cluster_name) )
      if (sum(if_row) == 0){
        temp_Bcount$cluster <- cluster_name
        temp_Bcount$n <- 0
        temp_Bcount$percentage <- 0
        cluster_group_sample_count[nrow(cluster_group_sample_count) + 1,] <- temp_Bcount
      }
    }
}

cols_group <- c('#984EA3','#4DAF4A')
names(cols_group) <- c('RSTR','LTBI')

cluster_group_sample_count <- cluster_group_sample_count %>% dplyr::group_by(cluster,Group) %>% dplyr::mutate(median = median(percentage))
temp <- unique(cluster_group_sample_count[,c('cluster','Group','median')])
temp <- temp %>% pivot_wider(names_from = Group, values_from = median)
temp$diff <- temp$RSTR - temp$LTBI
temp$cluster <- as.character(temp$cluster)
temp <- temp %>% arrange(desc(diff))
cluster_group_sample_count$cluster <- factor(cluster_group_sample_count$cluster,levels = temp$cluster)
ggplot(cluster_group_sample_count,aes(x=cluster, y=percentage))  +
  geom_boxplot(position=position_dodge(1),aes(fill = Group)) +
  theme(text = element_text(size = 15),plot.title = element_text(size = 20, face = "bold")) +
  ylab('activated CD4 %') +
  geom_jitter(position=position_dodge(1),size = 1,aes(fill = Group),shape = 2) +
  ggtitle('antigen-specific CD4 T cell subsets') +
  scale_fill_manual(values = cols_group) +
  theme_bw()
dev.print(pdf, paste('TB_RSTR_PP1_only_CD4_flow_plus_transformed_PCR_binary_scale_Group_CD4_subsets.pdf',sep = ''),width = 10, height = 2.5)

###### gating ###############################
# trace("plotAbundances", edit=TRUE)
# # take the last p of shape out of the aes_string. shape = 2

es <- assay(CD4_flow_plus_transformed_PCR_binary_scale, "exprs")
es <- es[type_markers(CD4_flow_plus_transformed_PCR_binary_scale), ]

cs_FOXP3 <- split(seq_len(ncol(CD4_flow_plus_transformed_PCR_binary_scale)), cut(es["FOXP3", ], nk <- 2))
cs_CD25 <- split(seq_len(ncol(CD4_flow_plus_transformed_PCR_binary_scale)), cut(es["CD25", ], nk <- 2))
cs <- list()
cs[['FOXP3+CD25+CD4']] <- intersect(cs_FOXP3$`(1.34,3.02]`,cs_CD25$`(0.61,4.39]`)
cs[['FOXP3-CD25+CD4']] <- intersect(cs_FOXP3$`(-0.335,1.34]`,cs_CD25$`(0.61,4.39]`)
cs[['FOXP3+CD25-CD4']] <- intersect(cs_FOXP3$`(1.34,3.02]`,cs_CD25$`(-3.17,0.61]`)
cs[['FOXP3-CD25-CD4']] <- intersect(cs_FOXP3$`(-0.335,1.34]`,cs_CD25$`(-3.17,0.61]`)
kids <- lapply(seq_len(4), function(i) {
  rep(i, length(cs[[i]]))
})
kids <- factor(unlist(kids))

# store cluster IDs in cell metadata & codes in metadata
CD4_flow_plus_transformed_PCR_binary_scale$cluster_id <- as.character(CD4_flow_plus_transformed_PCR_binary_scale$cluster_id)
CD4_flow_plus_transformed_PCR_binary_scale$cluster_id[unlist(cs)] <- unlist(kids)
CD4_flow_plus_transformed_PCR_binary_scale$cluster_id <- as.factor(CD4_flow_plus_transformed_PCR_binary_scale$cluster_id)
metadata(CD4_flow_plus_transformed_PCR_binary_scale)$cluster_codes <- data.frame(
  custom = factor(levels(kids), levels = levels(kids)))

merging_table <- data.frame(original_cluster = c(1:4), new_cluster = '')
merging_table$new_cluster[merging_table$original_cluster %in% c(1)] <- 'FOXP3+CD25+CD4' # the difference is CD25 and CD137
merging_table$new_cluster[merging_table$original_cluster %in% c(2)] <- 'FOXP3-CD25+CD4' # the difference is CD25 and CD137
merging_table$new_cluster[merging_table$original_cluster %in% c(3)] <- 'FOXP3+CD25-CD4' # the difference is CD25 and CD137
merging_table$new_cluster[merging_table$original_cluster %in% c(4)] <- 'FOXP3-CD25-CD4' # the difference is CD25 and CD137

CD4_flow_plus_transformed_PCR_binary_scale <- mergeClusters(CD4_flow_plus_transformed_PCR_binary_scale, k =  'custom', table = merging_table, id = "FOXP3_CD25_gate",overwrite = T)

plotAbundances(CD4_flow_plus_transformed_PCR_binary_scale, k = "FOXP3_CD25_gate", by = "cluster_id") + 
  scale_color_manual(values = cols_group) + 
  scale_fill_manual(values = cols_group) + 
  ylab('activated CD4 %') + 
  theme(text = element_text(size = 13),plot.title = element_text(size = 13, face = "bold")) +
  ylim(0,10) + 
  facet_wrap('cluster_id',ncol = 2)
dev.print(pdf, paste('TB_RSTR_PP1_only_CD4_flow_plus_transformed_PCR_binary_scale_CD25_FOXP3_gate_Group.pdf',sep = ''),width = 5.2, height = 5)

# RORC and TBET
cs_RORC <- split(seq_len(ncol(CD4_flow_plus_transformed_PCR_binary_scale)), cut(es["RORC", ], nk <- 2))
cs_TBET <- split(seq_len(ncol(CD4_flow_plus_transformed_PCR_binary_scale)), cut(es["TBET", ], nk <- 2))
cs <- list()
cs[['RORC+TBET+CD4']] <- intersect(cs_RORC$`(0.0297,1.03]`,cs_TBET$`(0.845,2.16]`)
cs[['RORC-TBET+CD4']] <- intersect(cs_RORC$`(-0.973,0.0297]`,cs_TBET$`(0.845,2.16]`)
cs[['RORC+TBET-CD4']] <- intersect(cs_RORC$`(0.0297,1.03]`,cs_TBET$`(-0.467,0.845]`)
cs[['RORC-TBET-CD4']] <- intersect(cs_RORC$`(-0.973,0.0297]`,cs_TBET$`(-0.467,0.845]`)
kids <- lapply(seq_len(4), function(i) {
  rep(i, length(cs[[i]]))
})
kids <- factor(unlist(kids))

# store cluster IDs in cell metadata & codes in metadata
CD4_flow_plus_transformed_PCR_binary_scale$cluster_id <- as.character(CD4_flow_plus_transformed_PCR_binary_scale$cluster_id)
CD4_flow_plus_transformed_PCR_binary_scale$cluster_id[unlist(cs)] <- unlist(kids)
CD4_flow_plus_transformed_PCR_binary_scale$cluster_id <- as.factor(CD4_flow_plus_transformed_PCR_binary_scale$cluster_id)
metadata(CD4_flow_plus_transformed_PCR_binary_scale)$cluster_codes <- data.frame(
  custom = factor(levels(kids), levels = levels(kids)))

merging_table <- data.frame(original_cluster = c(1:4), new_cluster = '')
merging_table$new_cluster[merging_table$original_cluster %in% c(1)] <- 'RORC+TBET+CD4' # the difference is CD25 and CD137
merging_table$new_cluster[merging_table$original_cluster %in% c(2)] <- 'RORC-TBET+CD4' # the difference is CD25 and CD137
merging_table$new_cluster[merging_table$original_cluster %in% c(3)] <- 'RORC+TBET-CD4' # the difference is CD25 and CD137
merging_table$new_cluster[merging_table$original_cluster %in% c(4)] <- 'RORC-TBET-CD4' # the difference is CD25 and CD137

CD4_flow_plus_transformed_PCR_binary_scale <- mergeClusters(CD4_flow_plus_transformed_PCR_binary_scale, k =  'custom', table = merging_table, id = "RORC_TBET_gate",overwrite = T)

plotAbundances(CD4_flow_plus_transformed_PCR_binary_scale, k = "RORC_TBET_gate", by = "cluster_id") + 
  scale_color_manual(values = cols_group) + 
  scale_fill_manual(values = cols_group) + 
  ylab('activated CD4 %') + 
  theme(text = element_text(size = 13),plot.title = element_text(size = 13, face = "bold")) +
  ylim(0,75) + 
  facet_wrap('cluster_id',ncol = 2)
dev.print(pdf, paste('TB_RSTR_PP1_only_CD4_flow_plus_transformed_PCR_binary_scale_TBET_RORC_gate_Group.pdf',sep = ''),width = 5.2, height = 5)

##### clonal expansion per Group visualization ################################
cluster_clonal_expanded_dataframe <- data.frame(cluster = cluster_ids(flow_plus_transformed_PCR_binary_scale,'T_cell_subsets'),clonally_expanded = sce_flow_PCR_raw@colData@listData[["clonal_expanded"]],Group = flow_plus_transformed_PCR_binary_scale@colData@listData[["Group"]])
cluster_clonal_expanded_group_sample_count <- dplyr::count(cluster_clonal_expanded_dataframe, cluster, clonally_expanded,Group)
cluster_clonal_expanded_group_sample_count <- cluster_clonal_expanded_group_sample_count %>% dplyr::group_by(Group,cluster) %>% dplyr::mutate(total = sum(n))
cluster_clonal_expanded_group_sample_count$percentage <- cluster_clonal_expanded_group_sample_count$n/cluster_clonal_expanded_group_sample_count$total*100
# Filling zero values for the Group that has no counts!!
for (cluster_name in unique(cluster_clonal_expanded_dataframe$cluster)){
  temp_Bcount <- cluster_clonal_expanded_group_sample_count[cluster_clonal_expanded_group_sample_count$cluster == cluster_name,][1,]
  for (temp in unique(cluster_clonal_expanded_dataframe$Group)) {
    for (clonal_expanded_type in c('yes','no')){
      if_row <- ((cluster_clonal_expanded_group_sample_count$cluster == cluster_name) & 
                   (cluster_clonal_expanded_group_sample_count$clonally_expanded == clonal_expanded_type) & 
                   (cluster_clonal_expanded_group_sample_count$Group == temp))
      if (sum(if_row) == 0){
        temp_Bcount$cluster <- cluster_name
        temp_Bcount$clonally_expanded <- clonal_expanded_type
        temp_Bcount$Group <- temp
        temp_Bcount$n <- 0
        temp_Bcount$percentage <- 0
        cluster_clonal_expanded_group_sample_count[nrow(cluster_clonal_expanded_group_sample_count) + 1,] <- temp_Bcount
      }
    }
  }
}  
ggplot(cluster_clonal_expanded_group_sample_count,aes(x=cluster, y=percentage,fill = clonally_expanded))  +
  geom_bar(stat="identity") +
  facet_grid(Group ~ .) +
  theme(text = element_text(size = 20,family = 'Arial'),plot.title = element_text(size = 20, face = "bold")) + RotatedAxis()+
  ylab('% T cell clusters') +
  theme_bw() + RotatedAxis() +
  scale_fill_manual(values = c("yes" = "red", "no" = "darkgrey")) +
  guides(fill=guide_legend(title="clonal\nexpansion"))
dev.print(pdf, paste('TB_RSTR_PP1_only_flow_plus_transformed_PCR_binary_scale_Group_allcluster_clonal_expanded.pdf',sep = ''),width = 5, height = 4)
