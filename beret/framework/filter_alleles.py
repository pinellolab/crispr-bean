from scipy.stats import fisher_exact
import pandas as pd
from tqdm.notebook import tqdm
import multiprocessing
import mlcrate as mlc
from beret.framework.Edit import Allele


def get_edit_significance_to_ctrl(sample_adata, ctrl_adata):
    sample_columns = sample_adata.condit['index'].tolist()
    edit_counts_ctrl = ctrl_adata.uns['edit_counts']
    edit_counts_sample = sample_adata.uns['edit_counts'][['guide', 'edit'] + sample_columns]
    edit_counts_merged = pd.merge(edit_counts_sample, edit_counts_ctrl, on = ['guide', 'edit'], 
                                  how = 'outer').fillna(0)
    guide_count_ctrl = ctrl_adata._get_allele_norm(edit_counts_merged, thres = 0)[:,0]
    guide_count_sample = sample_adata._get_allele_norm(edit_counts_merged, thres = 0)
    edit_counts_merged = edit_counts_merged.set_index(['guide', 'edit'])

    sample_edit = edit_counts_merged[sample_columns]
    ctrl_edit = edit_counts_merged['LDLRCDS_plasmid']
    odds_ratios_tbl = pd.DataFrame().reindex_like(sample_edit)
    p_value_tbl = pd.DataFrame().reindex_like(sample_edit)
    
    def fisher_test_single_sample(j):
        odds_ratio_series = pd.Series(dtype='float64').reindex_like(ctrl_edit)
        p_value_series = pd.Series(dtype='float64').reindex_like(ctrl_edit)

        for i in tqdm(range(len(edit_counts_merged))):
            fisher_tbl = [[sample_edit.iloc[i, j], guide_count_sample[i, j] - sample_edit.iloc[i, j]],
                            [ctrl_edit.iloc[i], guide_count_ctrl[i] - ctrl_edit.iloc[i]]]
            odds_ratio_series[i], p_value_series[i] = fisher_exact(fisher_tbl, alternative = "greater")
        return(odds_ratio_series, p_value_series)
    
    pool = multiprocessing.Pool(len(sample_columns))
    odds_ratios, p_values = zip(*pool.map(fisher_test_single_sample, range(len(sample_columns))))
    
    p_value_tbl = pd.concat(p_values, axis=1)
    p_value_tbl.columns = sample_columns
    odds_ratio_tbl = pd.concat(odds_ratios, axis=1)
    odds_ratio_tbl.columns = sample_columns
    
    edit_counts_sample = edit_counts_sample.set_index(['guide', 'edit'])
    
    q_bonf_tbl = pd.DataFrame().reindex_like(p_value_tbl)
    for j in range(len(sample_columns)):
        q_bonf_tbl.iloc[:, j] = p_value_tbl.iloc[:,j] * \
            (np.sum(edit_counts_sample.iloc[:,j] > 0) - np.sum(np.isnan(p_value_tbl.iloc[:, j])))                    
    q_bonf_tbl[q_bonf_tbl > 1] = 1
    return(odds_ratio_tbl, q_bonf_tbl) 

#https://stackoverflow.com/questions/8804830/python-multiprocessing-picklingerror-cant-pickle-type-function
def _filter_single_allele(allele, sig_edits):
    for edit in allele.edits.copy():
        if not edit in sig_edits:
            allele.edits.remove(edit)
    return(allele)

def _row_filter_allele(row, sample_guide_to_sig_edit_dict):
    try:
        sig_edits = sample_guide_to_sig_edit_dict[row.guide]
    except KeyError:
        return(Allele())
    filtered_allele = _filter_single_allele(row.allele, sig_edits)
    return(filtered_allele)

def _filter_allele_sample(sample_allele_tbl, sample_guide_to_sig_edit_dict):
    sample_allele_tbl = sample_allele_tbl.reset_index()
    sample_allele_tbl['filtered_allele_str'] = sample_allele_tbl.progress_apply(
        lambda row: _row_filter_allele(row, sample_guide_to_sig_edit_dict), axis=1).map(str)
    sample_filtered_allele_tbl = sample_allele_tbl.groupby(['guide', 'filtered_allele_str']).sum().reset_index()
    sample_filtered_allele_tbl = sample_filtered_allele_tbl.loc[
        sample_filtered_allele_tbl.filtered_allele_str != '', :]
    sample_filtered_allele_tbl = sample_filtered_allele_tbl.set_index(['guide', 'filtered_allele_str'], drop=True)
    return(sample_filtered_allele_tbl)

def _filter_alleles(allele_df, edit_significance_tbl):
    allele_df = allele_df.copy().set_index(['guide', 'allele'])
    def _filter_allele_sample_multiproc(i):
        sample_allele_df = allele_df.iloc[:,i]
        sample_sig_edit = edit_significance_tbl.iloc[:, i]
        sample_sig_edit = sample_sig_edit.loc[sample_sig_edit]
        sample_guide_to_sig_edit_dict = sample_sig_edit.index.to_frame().reset_index(drop=True).groupby('guide')['edit'].apply(list).to_dict()
        return(_filter_allele_sample(sample_allele_df, sample_guide_to_sig_edit_dict))
    filtered_allele_dfs = []
#     for j in tqdm(range(len(allele_df.columns))):
#         filtered_allele_dfs.append(_filter_allele_sample_multiproc(j))
    pool = mlc.SuperPool(len(allele_df.columns))
    filtered_allele_dfs = zip(*pool.map(_filter_allele_sample_multiproc, 
                                        range(len(allele_df.columns))))
    filtered_alleles = pd.DataFrame().join(filtered_allele_dfs, how="outer").reset_index()
    filtered_alleles.insert(1, 'allele', 
                                  filtered_alleles.filtered_allele_str.map(
                                      lambda s: Allele.from_str(s)))
    filtered_alleles = filtered_alleles.drop('filtered_allele_str', axis=1)
    return(filtered_alleles)

def filter_alleles(sample_adata, ctrl_adata, q_thres = 0.05):
    odds_ratio_tbl, q_bonf_tbl = get_edit_significance_to_ctrl(sample_adata, ctrl_adata)
    print("Done calculating significance")
    edit_significance_tbl = q_bonf_tbl < q_thres
    print("Filtering alleles for those containing significant edits (q < {})...".format(q_thres))
    filtered_alleles = _filter_alleles(sample_adata.uns['allele_counts'], edit_significance_tbl)
    return(q_bonf_tbl, filtered_alleles)