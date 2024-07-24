import anndata
import pandas as pd

import cstarpy
import os
import numpy as np
from cstarpy.separation import CellStateTransition

gbm_primary_df=pd.read_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_primary_df.pkl')
gbm_recurrent_df=pd.read_pickle('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_recurrent_df.pkl')
print('Documents loaded')

cst = CellStateTransition('pmr_rct', gbm_primary_df, gbm_recurrent_df)
print('Cells separated')
dpd_scores = cst.get_dpd()

norm_s_df = pd.DataFrame(np.stack([cst.n, cst.s], axis=1), index=cst.svm_input.data.columns, columns=["n", "s"])


pd.to_pickle(dpd_scores,"/home/jing/Phd_project/project_GBM/gbm_OUTPUT/dpd_pmr_rct.pkl")
pd.to_pickle(norm_s_df,"/home/jing/Phd_project/project_GBM/gbm_OUTPUT/stv_pmr_rct.pkl")
