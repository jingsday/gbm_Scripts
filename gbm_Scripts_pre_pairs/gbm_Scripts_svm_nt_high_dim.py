import anndata
import pandas as pd

import cstarpy
import os
import numpy as np
from cstarpy.separation import CellStateTransition


gbm_recurrent_tumor_df=pd.read_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_high_dim/gbm_recurrent_tumor_high_dim.pkl')
gbm_recurrent_normal_df=pd.read_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_high_dim/gbm_recurrent_normal_high_dim.pkl')
gbm_primary_tumor_df=pd.read_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_high_dim/gbm_primary_tumor_high_dim.pkl')
gbm_primary_normal_df=pd.read_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_high_dim/gbm_primary_normal_high_dim.pkl')


tumor=pd.concat([gbm_recurrent_tumor_df,gbm_primary_tumor_df],axis=0)
normal=pd.concat([gbm_primary_normal_df,gbm_recurrent_normal_df],axis=0)
tumor.set_index('index_clean',inplace=True)
normal.set_index('index_clean',inplace=True)
tumor.drop(columns=['Tumor_Normal_annotation'],inplace=True)
normal.drop(columns=['Tumor_Normal_annotation'],inplace=True)

print('Documents loaded')


cst = CellStateTransition('tmr_nml', tumor, normal)
print('Cells separated')
dpd_scores = cst.get_dpd()

norm_s_df = pd.DataFrame(np.stack([cst.n, cst.s], axis=1), index=cst.svm_input.data.columns, columns=["n", "s"])


pd.to_pickle(dpd_scores,"/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_svm_dpd_stvs/tumor_normal/dpd_tmr_nml.pkl")
pd.to_pickle(norm_s_df,"/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_svm_dpd_stvs/tumor_normal/stv_tmr_nml.pkl")




