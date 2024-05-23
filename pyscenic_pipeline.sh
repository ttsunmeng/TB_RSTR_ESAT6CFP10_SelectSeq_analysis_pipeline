# following 2020
# gets to the virtual environment
conda activate scenic_protocol


pyscenic grn TB_RSTR_CD4_lognorm_FilterLQCells.loom hs_hgnc_tfs.txt --num_workers 8 --output adj.csv --method grnboost2

pyscenic ctx adj.csv hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname TB_RSTR_CD4_lognorm_FilterLQCells.loom --output reg.csv --num_workers 8 --mask_dropouts

pyscenic aucell TB_RSTR_CD4_lognorm_FilterLQCells.loom reg.csv --output TB_RSTR_CD4_lognorm_FilterLQCells_SCENIC.loom --num_workers 8

#visualize the TF RSS in juypter notebook
Jupyter notebook
