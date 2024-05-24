ReadMe File

1. Index sort and targeted PCR FlowSom analysis
TB_RSTR_PP1_flow_PCR_TCR_analysis_clean.R

Line 14 Input: TB_RSTR_PP1_index_sort_PCR_TCR_table.csv to calculate clonality and prepare the data for .fcs in FlowSom

Line 224 Input: flow_PCR_raw_nonzero_042822/*.fcs files into FlowSom

Lines 297 and 324: need to input TB_RSTR_PP1_flow_PCR_SingleCellExperiment_object.rds to make tSNE and clustering consistent


2. Select-Seq analysis
TB_RSTR_PP1_CD4_clean.R

Line 23 Input TB_RSTR_PP1_CD4_raw.tsv
Line 27 Input 20200212_SI34_b1_to_b7_Subset_T_cloneVSsingleton_wCP.xlsx

Line 67 Need to Input TB_RSTR_PP1_CD4_Seurat_object.rds to make UMAP consistent

Line 145 input GO BP get from DAVID GO website

Line 345 Get the cytokine list

Line 472 Get the immune related genes for stringdb analysis

Line 542 Need to download the database from SCENIC website and save in the proper folder. 

Line 579 Generate .loom file and run it in pyscenic_pipeline.sh. Need to download the proper database from SCENIC



3. ACS validation
ACS_validation_scRNAseq_clean.R

Unzip 2023_supplementary_table2_T cell receptor repertoires associated with control and disease progression following Mycobacterium tuberculosis infection.csv.zip 

Input 2023_supplementary_table2_T cell receptor repertoires associated with control and disease progression following Mycobacterium tuberculosis infection.csv
