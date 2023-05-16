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
library(dplyr)

#https://www.bioconductor.org/packages/release/bioc/html/CATALYST.html
setwd("/Volumes/GoogleDrive/My\ Drive/Stanford/RNA-seq/data_analysis/TB_resistance/")

# The dataset was 2023 Musvosvi et al T cell receptor repertoires associated with control and disease progression following Mycobacterium tuberculosis infection
# supplementary table 2.
progressor_table <- read.csv('./data_input/2023_supplementary_table2_T cell receptor repertoires associated with control and disease progression following Mycobacterium tuberculosis infection.csv',check.names = F)
progressor_table$days <- as.numeric(gsub('D','',progressor_table$visit))
progressor_table$days[progressor_table$donorGrp == 'Progressor'] <- as.numeric(progressor_table$Days.To.TB[progressor_table$donorGrp == 'Progressor'])

progressor_table <- progressor_table %>% group_by(donorId) %>% mutate(max_day = max(days))
progressor_table$days[progressor_table$donorGrp == 'Non Progressor'] <- 
  progressor_table$days[progressor_table$donorGrp == 'Non Progressor'] - 
  progressor_table$max_day[progressor_table$donorGrp == 'Non Progressor']
CD4_progressor_table <- progressor_table[(progressor_table$Cell.Type == 'CD4') & (!(is.na(progressor_table$Cell.Type))),]
# colnames(CD4_progressor_table)

CD4_progressor_table$MAIT.Classification[is.na(CD4_progressor_table$MAIT.Classification)] <- 'NA'
CD4_progressor_table <- CD4_progressor_table[CD4_progressor_table$MAIT.Classification == 'Non.MAIT',]
CD4_progressor_table$binary <- (CD4_progressor_table$FOXP3 > 5)
# (CD4_progressor_table$TBET <= 5) & (CD4_progressor_table$RORC <= 5) 

temp <- data.frame(cbind(CD4_progressor_table$donorId,CD4_progressor_table$days,CD4_progressor_table$visit,CD4_progressor_table$Days.To.TB,CD4_progressor_table$donorGrp,as.numeric(CD4_progressor_table$binary)),check.names = F)
colnames(temp) <- c('donor_ID','days','visit','days_TB','Group','binary')
CD4_progressor_table$days <- as.numeric(CD4_progressor_table$days)
sapply(temp,class)
CD4_donor_cell_count <- dplyr::count(temp, donor_ID,days, visit,days_TB, Group,binary)
CD4_donor_cell_count <- CD4_donor_cell_count %>% group_by(donor_ID, days,visit, days_TB,Group) %>% mutate(total = sum(n))
CD4_donor_cell_count$percentage <- CD4_donor_cell_count$n/CD4_donor_cell_count$total*100
for (donor_name in unique(CD4_donor_cell_count$donor_ID)){
  for (dayname in unique(CD4_donor_cell_count$days[CD4_donor_cell_count$donor_ID == donor_name])){
    temp <- CD4_donor_cell_count[(CD4_donor_cell_count$donor_ID == donor_name) & 
                                   (CD4_donor_cell_count$days == dayname),]
    for (cluster_name in unique(CD4_donor_cell_count$binary)){
      if_row <- ((CD4_donor_cell_count$donor_ID == donor_name) & 
                   (CD4_donor_cell_count$binary == cluster_name) & 
                   (CD4_donor_cell_count$days == dayname))
      if (sum(if_row) == 0){
        temp$binary <- cluster_name
        temp$n <- 0
        temp$percentage <- 0
        CD4_donor_cell_count <- rbind(CD4_donor_cell_count,temp)
      }
    }
  }
  
}
CD4_donor_cell_count <- CD4_donor_cell_count[CD4_donor_cell_count$binary == '1',]
CD4_donor_cell_count$days <- as.numeric(CD4_donor_cell_count$days)

wilcoxon_test_result <- wilcox.test(CD4_donor_cell_count$percentage[(CD4_donor_cell_count$Group == 'Progressor')],
            CD4_donor_cell_count$percentage[(CD4_donor_cell_count$Group != 'Progressor')])

cols_group <- c('red','darkgrey')
names(cols_group) <- c('Progressor','Non Progressor')
ggplot(CD4_donor_cell_count,aes(x=Group, y=percentage,fill = Group))  +
  geom_boxplot(position=position_dodge(1), outlier.size = 1) +
  theme(text = element_text(size = 14),axis.text = element_text(size = 16),plot.title = element_text(size = 16, face = "bold")) + RotatedAxis() +
  ylab('activated CD4 %') +
  ggtitle('FOXP3+CD4') +
  geom_jitter(size = 1,aes(fill = Group)) +
  geom_text(aes(x = 1.5,y = 20, label = paste('p =',sprintf(wilcoxon_test_result$p.value, fmt = '%#.2f'))),size = 3) +
  scale_fill_manual(values = cols_group) + 
  theme_bw()
dev.print(pdf, paste('TB_2023_ACS_progressor_Treg_nonMAIT_CD4_alldays.pdf',sep = ''),width = 4, height = 4)
