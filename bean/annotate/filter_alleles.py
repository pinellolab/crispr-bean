from copy import deepcopy
from functools import reduce
import multiprocessing
import numpy as np
from scipy.stats import fisher_exact
from tqdm.auto import tqdm
import pandas as pd
from scipy.special import softmax
from ..framework.Edit import Allele
from ..framework.AminoAcidEdit import CodingNoncodingAllele
from ._supporting_fn import map_alleles_to_filtered


def sum_column_groups(mat, column_index_list):
    cols = [
        mat[:, np.array(colidx)].sum(axis=1, keepdims=True)
        for colidx in column_index_list
    ]
    return np.concatenate(cols, axis=1)


def fisher_test_single_sample_multiproc(j):
    odds_ratio_series = pd.Series(dtype="float64").reindex_like(ctrl_edit)
    p_value_series = pd.Series(dtype="float64").reindex_like(ctrl_edit)
    for i in tqdm(range(len(sample_edit)), leave=False):
        if sample_edit.iloc[i, j] == 0:
            odds_ratio_series[i], p_value_series[i] = 0, 1
            continue
        fisher_tbl = [
            [sample_edit.iloc[i, j], guide_count_sample[i, j] - sample_edit.iloc[i, j]],
            [ctrl_edit.iloc[i], guide_count_ctrl[i] - ctrl_edit.iloc[i]],
        ]
        odds_ratio_series[i], p_value_series[i] = fisher_exact(
            fisher_tbl, alternative="greater"
        )
    return (odds_ratio_series, p_value_series)


def fisher_test_single_sample(
    j, sample_edit, ctrl_edit, guide_count_sample, guide_count_ctrl
):
    odds_ratio_series = pd.Series(dtype="float64").reindex_like(ctrl_edit)
    p_value_series = pd.Series(dtype="float64").reindex_like(ctrl_edit)
    for i in tqdm(range(len(sample_edit)), leave=False):
        if sample_edit.iloc[i, j] == 0:
            odds_ratio_series[i], p_value_series[i] = 0, 1
            continue
        fisher_tbl = [
            [sample_edit.iloc[i, j], guide_count_sample[i, j] - sample_edit.iloc[i, j]],
            [ctrl_edit.iloc[i], guide_count_ctrl[i] - ctrl_edit.iloc[i]],
        ]
        assert (np.array(fisher_tbl) >= 0).all(), fisher_tbl
        odds_ratio_series[i], p_value_series[i] = fisher_exact(
            fisher_tbl, alternative="greater"
        )
    return (odds_ratio_series, p_value_series)


def get_edit_significance_to_ctrl(
    sample_adata,
    ctrl_adata,
    aggregate_cond=None,
    run_parallel=False,
    allele_counts_key="allele_counts",
):
    """
    Calculate edit counts for all edits in sample_adata, regardless of the edit presence in the sample.
    """
    if "index" in sample_adata.condit.columns:
        sample_columns = sample_adata.condit["index"].tolist()
    else:
        sample_columns = sample_adata.condit.index.tolist()
    n_samples = len(sample_columns)
    sample_adata.uns["edit_counts"] = sample_adata.get_edit_from_allele(
        allele_count_key=allele_counts_key, return_result=True
    )
    ctrl_adata.uns["edit_counts"] = ctrl_adata.get_edit_from_allele(
        allele_count_key=allele_counts_key, return_result=True
    )
    edit_counts_ctrl = ctrl_adata.uns["edit_counts"]

    if aggregate_cond is not None:
        conds = sample_adata.condit.reset_index().groupby(aggregate_cond).groups
        edit_counts_raw = sample_adata.uns["edit_counts"][
            ["guide", "edit"] + sample_columns
        ].set_index(["guide", "edit"])
        edit_counts_sample = sample_adata.uns["edit_counts"][
            ["guide", "edit"]
        ].set_index(["guide", "edit"])
        for k, col_idx in conds.items():
            edit_counts_sample[k] = edit_counts_raw.iloc[
                :, np.array(col_idx.tolist())
            ].sum(axis=1)
        n_samples = len(conds.keys())
        sample_columns = conds.keys()
    else:
        edit_counts_sample = sample_adata.uns["edit_counts"][
            ["guide", "edit"] + sample_columns
        ]
    edit_counts_merged = pd.merge(
        edit_counts_sample, edit_counts_ctrl, on=["guide", "edit"], how="left"
    ).fillna(0)
    guide_count_ctrl = ctrl_adata._get_allele_norm(edit_counts_merged, thres=0)[:, 0]
    guide_count_sample = sample_adata._get_allele_norm(edit_counts_merged, thres=0)
    edit_counts_merged = edit_counts_merged.set_index(["guide", "edit"])
    if aggregate_cond is not None:
        guide_count_sample = sum_column_groups(guide_count_sample, conds.values())
    sample_edit = edit_counts_merged[sample_columns]
    ctrl_edit = edit_counts_merged.iloc[:, -1]

    print("Running Fisher's exact test to get significant edits compared to control...")

    def child_initialize2(
        _sample_edit, _ctrl_edit, _guide_count_sample, _guide_count_ctrl
    ):
        global sample_edit, ctrl_edit, guide_count_sample, guide_count_ctrl
        sample_edit = _sample_edit
        ctrl_edit = _ctrl_edit
        guide_count_sample = _guide_count_sample
        guide_count_ctrl = _guide_count_ctrl

    if run_parallel:
        with multiprocessing.Pool(
            min(n_samples, 30),
            initializer=child_initialize2,
            initargs=(sample_edit, ctrl_edit, guide_count_sample, guide_count_ctrl),
        ) as pool:
            odds_ratios, p_values = zip(
                *pool.map(fisher_test_single_sample_multiproc, range(n_samples))
            )
    else:
        odds_ratios, p_values = zip(
            *map(
                lambda i: fisher_test_single_sample(
                    i, sample_edit, ctrl_edit, guide_count_sample, guide_count_ctrl
                ),
                range(n_samples),
            )
        )
    p_value_tbl = pd.concat(p_values, axis=1)
    p_value_tbl.columns = sample_columns
    odds_ratio_tbl = pd.concat(odds_ratios, axis=1)
    odds_ratio_tbl.columns = sample_columns

    q_bonf_tbl = pd.DataFrame().reindex_like(p_value_tbl)
    for j in range(len(sample_columns)):
        q_bonf_tbl.iloc[:, j] = p_value_tbl.iloc[:, j] * (
            np.sum(edit_counts_sample[sample_columns].iloc[:, j] > 0)
            - np.sum(np.isnan(p_value_tbl.iloc[:, j]))
        )
    q_bonf_tbl[q_bonf_tbl > 1] = 1
    return (odds_ratio_tbl, q_bonf_tbl)


# https://stackoverflow.com/questions/8804830/python-multiprocessing-picklingerror-cant-pickle-type-function
def _filter_single_allele(allele: Allele, sig_edits):
    filtered_allele = deepcopy(allele)
    for edit in allele.edits:
        if edit not in sig_edits:
            filtered_allele.edits.remove(edit)
    return filtered_allele


def _row_filter_allele(row, sample_guide_to_sig_edit_dict):
    try:
        sig_edits = sample_guide_to_sig_edit_dict[row.guide]
    except KeyError:
        return Allele()
    return _filter_single_allele(row.allele, sig_edits)


def _filter_allele_sample(sample_allele_tbl, sample_guide_to_sig_edit_dict):
    sample_allele_tbl = sample_allele_tbl.reset_index()
    sample_allele_tbl["filtered_allele_str"] = sample_allele_tbl.apply(
        lambda row: _row_filter_allele(row, sample_guide_to_sig_edit_dict), axis=1
    ).map(str)
    sample_filtered_allele_tbl = (
        sample_allele_tbl.groupby(["guide", "filtered_allele_str"])
        .sum(numeric_only=True)
        .reset_index()
    )
    sample_filtered_allele_tbl = sample_filtered_allele_tbl.loc[
        sample_filtered_allele_tbl.filtered_allele_str != "", :
    ]
    return sample_filtered_allele_tbl.set_index(
        ["guide", "filtered_allele_str"], drop=True
    )


def _filter_allele_sample_multiproc(i):
    sample_allele_df = allele_df.iloc[:, i]
    sample_allele_df = sample_allele_df[sample_allele_df > 0]
    if filter_each_sample:
        sample_sig_edit = edit_significance_tbl.iloc[:, i] < q_thres
    else:
        sample_sig_edit = (edit_significance_tbl < q_thres).any(axis=1)
    sample_sig_edit = sample_sig_edit.loc[sample_sig_edit]
    sample_guide_to_sig_edit_dict = (
        sample_sig_edit.index.to_frame()
        .reset_index(drop=True)
        .groupby("guide")["edit"]
        .apply(list)
        .to_dict()
    )
    return _filter_allele_sample(sample_allele_df, sample_guide_to_sig_edit_dict)


def _filter_allele_sample_loop(
    i, allele_df, edit_significance_tbl, q_thres, filter_each_sample=False
):
    sample_allele_df = allele_df.iloc[:, i]
    sample_allele_df = sample_allele_df[sample_allele_df > 0]
    if filter_each_sample:
        sample_sig_edit = edit_significance_tbl.iloc[:, i] < q_thres
    else:
        sample_sig_edit = (edit_significance_tbl < q_thres).any(axis=1)
    sample_sig_edit = sample_sig_edit.loc[sample_sig_edit]
    sample_guide_to_sig_edit_dict = (
        sample_sig_edit.index.to_frame()
        .reset_index(drop=True)
        .groupby("guide")["edit"]
        .apply(list)
        .to_dict()
    )
    return _filter_allele_sample(sample_allele_df, sample_guide_to_sig_edit_dict)


def _filter_alleles(
    allele_df,
    edit_significance_tbl,
    q_thres,
    n_threads=None,
    filter_each_sample=False,
    run_parallel=False,
):
    """
    - args
        filter_each_sample (bool) : filter out edits that are insignificant in each of the sample.
        If False, preserve the edit if called significant in any of the sample.
    """
    if n_threads is None:
        n_threads = len(allele_df.columns)
    allele_df = allele_df.set_index(["guide", "allele"]).copy().fillna(0)

    def child_initialize(
        _allele_df, _edit_significance_tbl, _filter_each_sample, _q_thres
    ):
        # https://stackoverflow.com/questions/25825995/python-multiprocessing-only-one-process-is-running
        global allele_df, edit_significance_tbl, filter_each_sample, q_thres
        allele_df = _allele_df
        edit_significance_tbl = _edit_significance_tbl
        filter_each_sample = _filter_each_sample
        q_thres = _q_thres

    if run_parallel:
        print("Running {} parallel processes to filter alleles...".format(n_threads))
        with multiprocessing.Pool(
            n_threads,
            initializer=child_initialize,
            initargs=(allele_df, edit_significance_tbl, filter_each_sample, q_thres),
        ) as pool:
            filtered_allele_dfs = pool.map(
                _filter_allele_sample_multiproc, list(range(len(allele_df.columns)))
            )
    else:
        print("Filtering alleles...")
        filtered_allele_dfs = []
        for i in tqdm(range(len(allele_df.columns))):
            filtered_allele_dfs.append(
                _filter_allele_sample_loop(
                    i, allele_df, edit_significance_tbl, q_thres, filter_each_sample
                )
            )

    print("Done filtering alleles, merging the result...")

    try:
        filtered_alleles = reduce(
            lambda left, right: left.join(right, how="outer"), filtered_allele_dfs
        ).reset_index()
        filtered_alleles.insert(
            1,
            "allele",
            filtered_alleles.filtered_allele_str.map(lambda s: Allele.from_str(s)),
        )
        filtered_alleles = filtered_alleles.drop("filtered_allele_str", axis=1)
        return filtered_alleles
    except Exception as e:
        print(e)
        return filtered_allele_dfs


def filter_alleles(
    sample_adata,
    ctrl_adata,
    allele_counts_key="allele_counts",
    q_thres=0.05,
    OR_thres=2,
    aggregate_cond=None,
    filter_each_sample=False,
    edit_sig_tbl=None,
    n_threads=30,
    run_parallel=False,
    map_to_filtered=True,
):
    if aggregate_cond is not None:
        sample_tested = sample_adata.condit.groupby(aggregate_cond).ngroups
    elif not filter_each_sample:
        sample_tested = len(sample_adata.condit)
    else:
        sample_tested = 1

    q_thres /= sample_tested

    if edit_sig_tbl is None:
        odds_ratio_tbl, q_bonf_tbl = get_edit_significance_to_ctrl(
            sample_adata,
            ctrl_adata,
            aggregate_cond,
            run_parallel=run_parallel,
            allele_counts_key=allele_counts_key,
        )
        print("Done calculating significance.\n\n")
    else:
        q_bonf_tbl = edit_sig_tbl
        print("Using provided edit significance table.\n")
    print(
        f"Filtering alleles for those containing significant edits (q < {q_thres})..."
    )
    filtered_alleles = _filter_alleles(
        sample_adata.uns[allele_counts_key],
        q_bonf_tbl,
        q_thres,
        filter_each_sample=filter_each_sample,
        n_threads=n_threads,
        run_parallel=run_parallel,
    )
    if map_to_filtered:
        filtered_alleles = map_alleles_to_filtered(
            sample_adata.uns[allele_counts_key], filtered_alleles
        )
    print("Done!")
    return (q_bonf_tbl, filtered_alleles)


def filter_allele_prop(
    adata,
    allele_uns_key: str,
    allele_prop_thres: float = 0.05,
    sample_prop_thres: float = 0.1,
    map_to_filtered=True,
    retain_max=True,
    allele_col: str = "allele",
    jaccard_threshold=0.5,
    aa_jaccard_threshold=0.5,
    nt_jaccard_threshold=0.5,
    distribute=False,
):
    """
    Filter allele based on the proportion of that allele among the guides.
    Arguments
    -- allele_prop_thres: Proportion of allele among the barcode matched guide counts to filter for
    -- sample_prop_thres: Proportion of samples where allele proportion exceeds the allele_prop_thres.
    -- map_to_filtered: If True, map the allele counts that are filtered out to the closest and most abundant allele that is retained after filtering. If False, discard the read counts that are filtered out.
    -- retain_max: If True, in case there is no allele that is retained during the filtering step, retain one allele with maximum median allele proportion across samples.
    -- alelle_col: String key for the allele column in adata.uns[allele_uns_key]
    """
    alleles = adata.uns[allele_uns_key].copy()
    aa_prop = adata.get_normalized_allele_counts(alleles).set_index(
        ["guide", allele_col]
    )
    aa_prop_filtered = aa_prop.loc[
        (np.nan_to_num(aa_prop) > allele_prop_thres).mean(axis=1) >= sample_prop_thres,
        :,
    ]
    # if no allele is left for the guide, retain the maximum frequency allele
    if retain_max:
        rows_to_add = []
        for guide in (
            aa_prop_filtered.reset_index().groupby("guide")[allele_col].count() == 0
        ).index:
            guide_df = aa_prop.loc[
                aa_prop.index.get_level_values("guide") == guide,
            ]
            max_frequency_idx = np.argmax(guide_df.median(axis=1))
            rows_to_add.append(guide_df.iloc[max_frequency_idx, :])
        aa_prop_filtered = pd.concat(
            [aa_prop_filtered] + rows_to_add, ignore_index=True
        ).reset_index(drop=True)
    alleles = alleles.set_index(["guide", allele_col])
    aa_filtered = alleles.loc[alleles.index.isin(aa_prop_filtered.index)]
    alleles = alleles.reset_index()
    aa_filtered = aa_filtered.reset_index()
    if map_to_filtered:
        print("Mapping filtered alleles...")
        if distribute:
            aa_filtered = _distribute_alleles_to_filtered(
                alleles,
                aa_filtered,
                jaccard_threshold=jaccard_threshold,
                allele_col=allele_col,
            )
        else:
            aa_filtered = _map_alleles_to_filtered(
                alleles,
                aa_filtered,
                jaccard_threshold=jaccard_threshold,
                aa_jaccard_threshold=aa_jaccard_threshold,
                nt_jaccard_threshold=nt_jaccard_threshold,
                allele_col=allele_col,
            )
    return aa_filtered


def _map_alleles_to_filtered(
    raw_allele_counts: pd.DataFrame,
    filtered_allele_counts: pd.DataFrame,
    jaccard_threshold=0.5,
    aa_jaccard_threshold=0.5,
    nt_jaccard_threshold=0.5,
    allele_col="allele",
):
    """
    Map pre-filtering alleles to the closest post-filtering alleles and aggregate the read count.
    """
    mapped_allele_counts = []
    is_cn_allele = isinstance(
        raw_allele_counts.reset_index()[allele_col][0], CodingNoncodingAllele
    )
    for guide, guide_raw_counts in tqdm(
        raw_allele_counts.groupby("guide"),
        desc=f"Mapping alleles to closest filtered alleles with min Jaccard index (nt: {nt_jaccard_threshold}, aa: {aa_jaccard_threshold})",
    ):

        guide_filtered_allele_counts = filtered_allele_counts.loc[
            filtered_allele_counts.guide == guide, :
        ].set_index(allele_col)
        guide_filtered_alleles = guide_filtered_allele_counts.index.tolist()
        if len(guide_filtered_alleles) != 0:
            # prioritize allele with most mean counts
            merge_priority = guide_filtered_allele_counts.mean(
                axis=1, numeric_only=True
            )
            if is_cn_allele:
                guide_raw_counts["allele_mapped"] = guide_raw_counts[allele_col].map(
                    lambda allele: allele.map_to_closest(
                        guide_filtered_alleles,
                        aa_jaccard_threshold=aa_jaccard_threshold,
                        nt_jaccard_threshold=nt_jaccard_threshold,
                        merge_priority=merge_priority,
                    )
                )
            else:
                guide_raw_counts["allele_mapped"] = guide_raw_counts[allele_col].map(
                    lambda allele: allele.map_to_closest(
                        guide_filtered_alleles,
                        jaccard_threshold=jaccard_threshold,
                        merge_priority=merge_priority,
                    )
                )
            guide_raw_counts = (
                guide_raw_counts.drop(allele_col, axis=1)
                .rename(columns={"allele_mapped": allele_col})
                .groupby(["guide", allele_col])
                .sum()
            )

            mapped_allele_counts.append(guide_raw_counts)
    res = pd.concat(mapped_allele_counts).reset_index()
    res = res.loc[res[allele_col].map(bool), :]
    return res


def _distribute_alleles_to_filtered(
    raw_allele_counts: pd.DataFrame,
    filtered_allele_counts: pd.DataFrame,
    jaccard_threshold=0.1,
    allele_col="allele",
):
    """
    Distrubute pre-filtering allele counts to the post-filtering alleles based on
    similarity.
    """
    mapped_allele_counts = []
    is_cn_allele = isinstance(
        raw_allele_counts.reset_index()[allele_col][0], CodingNoncodingAllele
    )
    for guide, guide_filtered_counts in tqdm(
        filtered_allele_counts.groupby("guide"),
        desc="Mapping alleles to closest filtered alleles",
    ):
        guide_filtered_counts = guide_filtered_counts.set_index(["guide", allele_col])
        guide_raw_counts = raw_allele_counts.loc[
            raw_allele_counts.guide == guide, :
        ].set_index(["guide", allele_col])
        guide_filtered_alleles = guide_filtered_counts.index.get_level_values(
            allele_col
        ).tolist()
        res = guide_filtered_counts.values
        for i, (g, allele) in enumerate(guide_raw_counts.index):
            # for allelel in raw_counts, distribute based on similarity if not present in filtered_counts
            if allele in guide_filtered_alleles:
                continue
            if is_cn_allele:
                aa_jaccards, nt_jaccards = allele.get_jaccards(guide_filtered_alleles)
                jaccards = (aa_jaccards * 3 + nt_jaccards) / 4
            else:
                jaccards = allele.get_jaccards(guide_filtered_alleles)

            weights = np.zeros(len(jaccards))
            if not (jaccards > jaccard_threshold).any():
                # print(f"map {allele} to {guide_filtered_alleles} failed.")
                continue
            weights[jaccards > jaccard_threshold] = softmax(
                jaccards[jaccards > jaccard_threshold]
            )
            weights = np.minimum(weights, jaccards)
            dist_values = np.rint(
                weights[:, None] * guide_raw_counts.values[i, :][None, :]
            ).astype(int)
            if not (dist_values.sum(axis=0) <= guide_raw_counts.values[i, :]).all():
                pass  # rounding process can inflate total number of allele counts
            res += dist_values
        added_counts = pd.DataFrame(
            res,
            index=guide_filtered_counts.index,
            columns=guide_filtered_counts.columns,
        )
        mapped_allele_counts.append(added_counts)
    res = pd.concat(mapped_allele_counts).reset_index()
    res = res.loc[res[allele_col].map(bool), :]
    return res
