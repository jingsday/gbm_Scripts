import pandas as pd

import cstarpy
import os
import numpy as np
from cstarpy.separation import CellStateTransition


df=pd.read_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_annotated/gbm_OUTPUT_sct_full_annotated.pkl')
lm_df=pd.read_csv('/home/jing/Phd_project/project_UCD_blca/blca_Scripts/perturbations/01_outputs_2020/L1000_Data_norm_data.csv',index_col=0)

primary_df = df[df['Stage']=='Primary']
common_genes=primary_df.columns.intersection(lm_df.index)

primary_tumor=primary_df[primary_df['Tumor_Normal_annotation']=='Tumor']
primary_tumor=primary_tumor[common_genes]

primary_normal=primary_df[primary_df['Tumor_Normal_annotation']=='Normal']
primary_normal=primary_normal[common_genes]


print('Documents loaded')


cst = CellStateTransition('p_norm_tum', primary_tumor, primary_normal)
print('Cells separated')
dpd_scores = cst.get_dpd()

norm_s_df = pd.DataFrame(np.stack([cst.n, cst.s], axis=1), index=cst.svm_input.data.columns, columns=["n", "s"])
norm_s_df['n']

pd.to_pickle(dpd_scores,"/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_svm/gbm_OUTPUT_svm_sctranform/gbm_OUTPUT_sct_P_dpd_lm_norm_tum.pkl")
pd.to_pickle(norm_s_df,"/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_svm/gbm_OUTPUT_svm_sctranform/gbm_OUTPUT_sct_P_stv_lm_norm_tum.pkl")