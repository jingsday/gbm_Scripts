import anndata
import pandas as pd

import cstarpy
import os
import numpy as np
from cstarpy.separation import CellStateTransition
import seaborn as sns

gbm_recurrent_tumor_df=pd.read_pickle('~/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scvi_high_dim_df/gbm_recurrent_tumor_high_dim.pkl')
gbm_primary_tumor_df=pd.read_pickle('~/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scvi_high_dim_df/gbm_primary_tumor_high_dim.pkl')


gbm_primary_tumor_df.drop(columns='Tumor_Normal_annotation',inplace=True)
gbm_primary_tumor_df.set_index('index_clean',inplace=True)

gbm_recurrent_tumor_df.drop(columns='Tumor_Normal_annotation',inplace=True)
gbm_recurrent_tumor_df.set_index('index_clean',inplace=True)


cst = CellStateTransition('p_r_tmr', gbm_primary_tumor_df, gbm_recurrent_tumor_df)
print('Cells separated')
dpd_scores = cst.get_dpd()

norm_s_df = pd.DataFrame(np.stack([cst.n, cst.s], axis=1), index=cst.svm_input.data.columns, columns=["n", "s"])


pd.to_pickle(dpd_scores,"/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scvi_svm_dpd_stvs/primary_recurrent/PR_dpd_tmr_nml.pkl")
pd.to_pickle(norm_s_df,"/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scvi_svm_dpd_stvs/primary_recurrent/PR_stv_tmr_nml.pkl")




