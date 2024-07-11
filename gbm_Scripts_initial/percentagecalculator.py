import os
import numpy 
import pandas 
import scanpy 
import loompy 
from MulticoreTSNE import MulticoreTSNE 
import seaborn 
import matplotlib.pyplot 
import scipy
import csv
import gzip
import anndata 
from pathlib import Path
import glob

def calculator(sample):
    percentage = 3

    # STEP 1 -> set corresponding folders which the results of each specimen are written to
    output_folder = "/Users/lidiayung/github/notebooks/GBM/specimens"
    resource_folder = "/Users/lidiayung/project/resource/GSE174554_RAW/"
 
    matrix_path = glob.glob(f"{resource_folder}/GSM*_{sample}_matrix.mtx.gz")[0]
    features_path = glob.glob(f"{resource_folder}/GSM*_{sample}_features.tsv.gz")[0]
    barcodes_path = glob.glob(f"{resource_folder}/GSM*_{sample}_barcodes.tsv.gz")[0]


    output_path = os.path.join(output_folder, sample)
    # End of Step 1

    # STEP 2 -> creat output documents
    # path to unfiltered loom file (this will be created in the optional steps below)
    f_loom_path_unfilt = "unfiltered.loom" # test dataset, n=500 cells
    f_loom_path_scenic = "filtered_scenic.loom"
    f_anndata_path = "anndata.h5ad"
    f_pyscenic_output = "output.loom"
    f_final_loom = 'scenic_integrated-output.loom'

    sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_versions()
    sc.set_figure_params(dpi=150, fontsize=10, dpi_save=600)
    # End of Step 2

    # Step 3-> read data into anndata
    mat = scipy.io.mmread(matrix_path)
    feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
    gene_names = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
    feature_types = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
    barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]

    # convert the index and columns to DataFrame objects
    obs_df = matrix.index.to_frame(index=False)
    var_df = matrix.columns.to_frame(index=False)

    adata = ad.AnnData(X=matrix.values, obs=obs_df, var=var_df)


    row_attrs = {  "Gene": np.array(var_df[0]) ,}
    col_attrs = { "CellID":  np.array(matrix.index) , 
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
    }
    lp.create( f_loom_path_unfilt, adata.X.transpose(), row_attrs, col_attrs )   

    # End of Step 3

    #Step 4-> Filtering cells that have >2.5% mitochondrial read counts and <200 expressed genes

    adata = sc.read_loom( f_loom_path_unfilt )
    nCountsPerGene = np.sum(adata.X, axis=0)
    nCellsPerGene = np.sum(adata.X>0, axis=0)
        
    nCells=adata.X.shape[0]

    mito_genes = adata.var_names.str.startswith('MT-')
    # for each cell compute fraction of counts in mito genes vs. all genes
    adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / numpy.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    
    sc.pp.filter_cells(adata, min_genes=200 )
    adata = adata[adata.obs['percent_mito'] <=0.025, :]
    #End of Step 4

    #Step 5 -> Read filtered data
    adata.write( f_anndata_path )
    # create basic row and column attributes for the loom file:
    row_attrs = {    "Gene": np.array(adata.var_names) ,}
    col_attrs = {
     "CellID": np.array(adata.obs_names) ,
     "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
     "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
     }
    lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)
    info = adata.X.transpose()
    print(info.shape)
    #End of Step 5

    #Step 6 -> integrate metadata
    file = "/Users/lidiayung/project/resource/GSM5319572_SF4324/GSE174554_Tumor_normal_metadata.txt"
    metadata= pd.read_csv(file,sep=' ')
    metadata.head()
    #End of Step 6

    #Step 7 -> Calculate percentage and return the value from both pre and post data
    new_df = metadata[metadata["Sample#"] == sample].copy()
    new_df['Barcode'] = new_df['Barcode'].astype(str) + '-1'
    intersection_barcodes = set(new_df['Barcode']).intersection(adata.obs.index)
    tumor = new_df[new_df['Barcode'].isin(intersection_barcodes) & (new_df['Tumor_Normal_annotation'] == 'Tumor')]
    percentage = len(tumor)/len(adata.obs.nUMI)
    print("Percentage: {:.2%}".format(len(tumor)/len(adata.obs.nUMI)))


    #End of Step 7


    print(sample)
    print("Post-filtering tumor/normal: " + f"{len(tumor)}/{len(adata.obs.nUMI) - len(tumor)}")
    print("Percentage:" +percentage)

    return percentage