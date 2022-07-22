from copy import deepcopy
from functools import reduce
import multiprocessing
import numpy as np
from scipy.stats import fisher_exact
from tqdm.auto import tqdm
import pandas as pd
import mlcrate as mlc
from beret.framework.Edit import Allele

def sum_column_groups(mat, column_index_list):
    cols = []
    for colidx in column_index_list:
        cols.append(mat[:,np.array(colidx)].sum(axis=1, keepdims = True))
        
    group_sums = np.concatenate(cols, axis=1)
    return(group_sums)

def get_edit_significance_to_ctrl(sample_adata, ctrl_adata, aggregate_cond = None):
    if 'index' in sample_adata.condit.columns:
        sample_columns = sample_adata.condit['index'].tolist()
    else:
        sample_columns = sample_adata.condit.index.tolist()
    n_samples = len(sample_columns)
    if not "edit_counts" in sample_adata.uns.keys():
        sample_adata.uns['edit_counts'] = sample_adata.get_edit_from_allele()
    if not "edit_counts" in ctrl_adata.uns.keys():
        ctrl_adata.uns['edit_counts'] = ctrl_adata.get_edit_from_allele()

    edit_counts_ctrl = ctrl_adata.uns['edit_counts']

    if not aggregate_cond is None:
        conds = sample_adata.condit.reset_index().groupby(aggregate_cond).groups
        edit_counts_raw = sample_adata.uns['edit_counts'][
            ['guide', 'edit'] + sample_columns].set_index(['guide', 'edit'])
        edit_counts_sample = sample_adata.uns['edit_counts'][['guide', 'edit']].set_index(['guide', 'edit'])
        
        for k, col_idx in conds.items():
            edit_counts_sample[k] = edit_counts_raw.iloc[:,np.array(col_idx.tolist())].sum(axis=1)
            
            
        n_samples = len(conds.keys())
        sample_columns = conds.keys()
    else:
        edit_counts_sample = sample_adata.uns['edit_counts'][['guide', 'edit'] + sample_columns]

    edit_counts_merged = pd.merge(edit_counts_sample, edit_counts_ctrl, on = ['guide', 'edit'], 
                                  how = 'left').fillna(0)
    guide_count_ctrl = ctrl_adata._get_allele_norm(edit_counts_merged, thres = 0)[:,0]
    guide_count_sample = sample_adata._get_allele_norm(edit_counts_merged, thres = 0)
    edit_counts_merged = edit_counts_merged.set_index(['guide', 'edit'])
    if not aggregate_cond is None:
        guide_count_sample = sum_column_groups(guide_count_sample, conds.values())  
    sample_edit = edit_counts_merged[sample_columns]
    ctrl_edit = edit_counts_merged.iloc[:,-1]
    
    def fisher_test_single_sample(j):
        odds_ratio_series = pd.Series(dtype='float64').reindex_like(ctrl_edit)
        p_value_series = pd.Series(dtype='float64').reindex_like(ctrl_edit)
        for i in tqdm(range(len(edit_counts_merged))):
            if sample_edit.iloc[i, j] == 0: 
                odds_ratio_series[i], p_value_series[i] = 0, 1
                continue
            fisher_tbl = [[sample_edit.iloc[i, j], guide_count_sample[i, j] - sample_edit.iloc[i, j]],
                            [ctrl_edit.iloc[i], guide_count_ctrl[i] - ctrl_edit.iloc[i]]]
            odds_ratio_series[i], p_value_series[i] = fisher_exact(fisher_tbl, alternative = "greater")
        return(odds_ratio_series, p_value_series)
    
    print("Running Fisher's exact test to get significant edits compared to control...")
    with multiprocessing.Pool(min(n_samples, 30)) as pool:
        odds_ratios, p_values = zip(*pool.map(fisher_test_single_sample, range(n_samples)))
    
    p_value_tbl = pd.concat(p_values, axis=1)
    p_value_tbl.columns = sample_columns
    odds_ratio_tbl = pd.concat(odds_ratios, axis=1)
    odds_ratio_tbl.columns = sample_columns
    
    q_bonf_tbl = pd.DataFrame().reindex_like(p_value_tbl)
    for j in range(len(sample_columns)):
        q_bonf_tbl.iloc[:, j] = p_value_tbl.iloc[:,j] * \
            (np.sum(edit_counts_sample[sample_columns].iloc[:,j] > 0) - np.sum(np.isnan(p_value_tbl.iloc[:, j])))                    
    q_bonf_tbl[q_bonf_tbl > 1] = 1
    return(odds_ratio_tbl, q_bonf_tbl) 

#https://stackoverflow.com/questions/8804830/python-multiprocessing-picklingerror-cant-pickle-type-function
def _filter_single_allele(allele: Allele, sig_edits):
    filtered_allele = deepcopy(allele)
    for edit in allele.edits:
        if not edit in sig_edits:
            filtered_allele.edits.remove(edit)
    return(filtered_allele)

def _row_filter_allele(row, sample_guide_to_sig_edit_dict):
    try:
        sig_edits = sample_guide_to_sig_edit_dict[row.guide]
    except KeyError:
        return(Allele())
    filtered_allele = _filter_single_allele(row.allele, sig_edits)
    return(filtered_allele)

def _filter_allele_sample(sample_allele_tbl, sample_guide_to_sig_edit_dict):
    sample_allele_tbl = sample_allele_tbl.reset_index()
    sample_allele_tbl['filtered_allele_str'] = sample_allele_tbl.apply(
        lambda row: _row_filter_allele(row, sample_guide_to_sig_edit_dict), axis=1).map(str)
    sample_filtered_allele_tbl = sample_allele_tbl.groupby(['guide', 'filtered_allele_str']).sum().reset_index()
    sample_filtered_allele_tbl = sample_filtered_allele_tbl.loc[
        sample_filtered_allele_tbl.filtered_allele_str != '', :]
    sample_filtered_allele_tbl = sample_filtered_allele_tbl.set_index(['guide', 'filtered_allele_str'], drop=True)
    return(sample_filtered_allele_tbl)

def _filter_allele_sample_multiproc(i):
    sample_allele_df = allele_df.iloc[:,i]
    sample_allele_df = sample_allele_df[sample_allele_df > 0]
    if filter_each_sample:
        sample_sig_edit = edit_significance_tbl.iloc[:, i] < q_thres
    else:
        sample_sig_edit = (edit_significance_tbl < q_thres).any(axis=1)
    sample_sig_edit = sample_sig_edit.loc[sample_sig_edit]
    sample_guide_to_sig_edit_dict = sample_sig_edit.index.to_frame().reset_index(drop=True).groupby('guide')['edit'].apply(list).to_dict()
    return(_filter_allele_sample(sample_allele_df, sample_guide_to_sig_edit_dict))

def _filter_allele_sample_loop(i, allele_df, 
                          edit_significance_tbl, q_thres, filter_each_sample = False):
    sample_allele_df = allele_df.iloc[:,i]
    sample_allele_df = sample_allele_df[sample_allele_df > 0]
    if filter_each_sample:
        sample_sig_edit = edit_significance_tbl.iloc[:, i] < q_thres
    else:
        sample_sig_edit = (edit_significance_tbl < q_thres).any(axis=1)
    sample_sig_edit = sample_sig_edit.loc[sample_sig_edit]
    sample_guide_to_sig_edit_dict = sample_sig_edit.index.to_frame().reset_index(drop=True).groupby('guide')['edit'].apply(list).to_dict()
    return(_filter_allele_sample(sample_allele_df, sample_guide_to_sig_edit_dict))

def _filter_alleles(allele_df, edit_significance_tbl, q_thres, 
                    n_threads = None, 
                    filter_each_sample = False, multiprocess = False):
    '''
    - args
        filter_each_sample (bool) : filter out edits that are insignificant in each of the sample. 
        If False, preserve the edit if called significant in any of the sample.
    '''
    if n_threads is None:
        n_threads = len(allele_df.columns)
    allele_df = allele_df.set_index(['guide', 'allele'])
    
    def child_initialize(_allele_df, _edit_significance_tbl, _filter_each_sample, _q_thres):
        #https://stackoverflow.com/questions/25825995/python-multiprocessing-only-one-process-is-running
        global allele_df, edit_significance_tbl, filter_each_sample, q_thres
        allele_df = _allele_df
        edit_significance_tbl = _edit_significance_tbl
        filter_each_sample = _filter_each_sample
        q_thres = _q_thres

    if multiprocess:
        print("Running {} parallel processes to filter alleles...".format(n_threads))
        with multiprocessing.Pool(n_threads, 

            initializer = child_initialize, 
            initargs=(allele_df, edit_significance_tbl, filter_each_sample, q_thres)) as pool:
            filtered_allele_dfs = pool.map(_filter_allele_sample_multiproc, list(range(len(allele_df.columns))))
    else:
        print("Filtering alleles...".format(n_threads))
        filtered_allele_dfs = []
        for i in tqdm(range(len(allele_df.columns))):
            filtered_allele_dfs.append(_filter_allele_sample_loop(
                i, allele_df, 
                edit_significance_tbl, 
                q_thres, 
                filter_each_sample))
        
    print("Done filtering alleles, merging the result...")

    try:
        filtered_alleles = reduce(lambda left, right: left.join(right, how="outer"), 
                                  filtered_allele_dfs).reset_index()
        filtered_alleles.insert(1, 'allele', 
                                    filtered_alleles.filtered_allele_str.map(
                                        lambda s: Allele.from_str(s)))
        filtered_alleles = filtered_alleles.drop('filtered_allele_str', axis=1)
        return(filtered_alleles)
    except Exception as e:
        print(e)
        return(filtered_allele_dfs)

def filter_alleles(sample_adata, ctrl_adata, q_thres = 0.05, aggregate_cond = None, 
                   filter_each_sample = False, edit_sig_tbl = None, n_threads = 30,
                  parallel_filter = False):
    if not aggregate_cond is None:
        sample_tested = sample_adata.condit.groupby(aggregate_cond).ngroups
    elif not filter_each_sample:
        sample_tested = len(sample_adata.condit)
    else:
        sample_tested = 1
        
    q_thres /= sample_tested
    
    if edit_sig_tbl is None:
        odds_ratio_tbl, q_bonf_tbl = get_edit_significance_to_ctrl(
            sample_adata, ctrl_adata, aggregate_cond)
        print("Done calculating significance.\n\n")
    else:
        q_bonf_tbl = edit_sig_tbl
        print("Using provided edit significance table.\n")    
    print("Filtering alleles for those containing significant edits (q < {})...".format(q_thres))
    filtered_alleles = _filter_alleles(sample_adata.uns['allele_counts'], q_bonf_tbl, q_thres,              
        filter_each_sample=filter_each_sample, n_threads = n_threads, multiprocess = parallel_filter)
    print("Done!")
    return(q_bonf_tbl, filtered_alleles)

