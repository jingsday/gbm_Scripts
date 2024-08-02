import pandas as pd

import cstarpy
import os
import numpy as np
from cstarpy.separation import CellStateTransition


df=pd.read_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_sctransform/gbm_OUTPUT_sct_annotated/gbm_OUTPUT_sct_full_annotated.pkl')

primary_df = df[df['Stage']=='Primary']

primary_tumor=primary_df[primary_df['Tumor_Normal_annotation']=='Tumor']

primary_normal=primary_df[primary_df['Tumor_Normal_annotation']=='Normal']


primary_tumor.drop(columns=['Patients_id',	'Barcode',	'Tumor_Normal_annotation',	'Stage'],inplace=True)
primary_normal.drop(columns=['Patients_id',	'Barcode',	'Tumor_Normal_annotation',	'Stage'],inplace=True)

print('Documents loaded')


cst = CellStateTransition('p_norm_tum', primary_tumor, primary_normal)
print('Cells separated')
dpd_scores = cst.get_dpd()

norm_s_df = pd.DataFrame(np.stack([cst.n, cst.s], axis=1), index=cst.svm_input.data.columns, columns=["n", "s"])


pd.to_pickle(dpd_scores,"/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_svm/gbm_OUTPUT_svm_sctranform/gbm_OUTPUT_sct_P_dpd_norm_tum.pkl")
pd.to_pickle(norm_s_df,"/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_svm/gbm_OUTPUT_svm_sctranform/gbm_OUTPUT_sct_P_stv_norm_tum.pkl")