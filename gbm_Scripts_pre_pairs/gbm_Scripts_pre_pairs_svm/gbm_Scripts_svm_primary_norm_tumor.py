import pandas as pd

import cstarpy
import os
import numpy as np
from cstarpy.separation import CellStateTransition

gbm_primary_tumor=pd.read_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scvi_high_dim_df/gbm_primary_tumor_high_dim.pkl')


gbm_primary_normal=pd.read_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scvi_high_dim_df/gbm_primary_normal_high_dim.pkl')

gbm_primary_tumor.set_index('index_clean',inplace=True)
gbm_primary_tumor.drop(columns='Tumor_Normal_annotation',inplace=True)

gbm_primary_normal.set_index('index_clean',inplace=True)
gbm_primary_normal.drop(columns='Tumor_Normal_annotation',inplace=True)

cst = CellStateTransition('pnml_ptmr', gbm_primary_tumor, gbm_primary_normal)

dpd_scores = cst.get_dpd()

norm_s_df = pd.DataFrame(np.stack([cst.n, cst.s], axis=1), index=cst.svm_input.data.columns, columns=["n", "s"])

print(cst.h)

dpd_scores.to_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scvi_svm_dpd_stvs/primary_recurrent/dpd_primay_high_normal_tumor.pkl')

norm_s_df.to_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scvi_svm_dpd_stvs/primary_recurrent/stv_primay_high_normal_tumor.pkl')