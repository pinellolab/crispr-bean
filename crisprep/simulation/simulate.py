from typing import List
import numpy as np
import pandas as pd
import scipy.stats

"""
class SimulatedScreen:
    """
    #Simulate a screen.

    """
    def __init__(self, 
        total_num_cells
        edit_count_stats_file = None,
        n_guides,
        n_targets,
        n_guides_per_target = 5,
        total_cells,
        vars_per_mu = 10,
        mu_steps = 10,
        ):
        self.count_edit_stats = pd.read_csv(edit_count_stats_file)
        self.n_guides = n_guides
        self.n_targets = n_targets
        self.n_guides_per_target = n_guides_per_target
        self.total_cells = total_cells
        self.vars_per_mu = vars_per_mu
        self.mu_steps = mu_steps
"""

COUNT_EDIT_STATS = pd.read_csv("../data/0721_ANBE_data/guide_counts_lenient/LDLvar_bulk_count_stats.csv")

def _sample_all_from_data(n_guides, total_num_cells):
    row = COUNT_EDIT_STATS.sample(n_guides, replace = True).reset_index()
    return((np.floor(row.guide_count / COUNT_EDIT_STATS.guide_count.sum() * total_num_cells), 
            row.edit_rate))

def _get_effect_sizes(sdev = 1, 
                      n_targets = 700, 
                      n_guides_per_target = 5, 
                      total_cells = 8*10**6, 
                      vars_per_mu = 10, 
                      mu_steps = 10):
    no_effect_targets = n_targets - mu_steps*vars_per_mu
    mu_delta = np.linspace(sdev/mu_steps, sdev, num = mu_steps)
    coverage = np.floor(total_cells/n_targets/n_guides_per_target)
    effect_sizes = np.concatenate((np.zeros(no_effect_targets), np.repeat(mu_delta, vars_per_mu)))
    
    guide_info = pd.DataFrame({"target_id" : np.repeat(list(range(n_targets)), n_guides_per_target), 
                          "guide_num" : np.tile(list(range(n_guides_per_target)), n_targets), 
                          "effect_size" : np.repeat(effect_sizes, n_guides_per_target)})

    return((effect_sizes, guide_info))


def simulate_phenotypes(sdev = 1,
                    n_targets = 700,
                    n_guides_per_target = 5,
                    total_cells = 10**6*8,
                    edit_efficiency = 0.1,
                    sort_population = 0.3,
                    vars_per_mu = 10,
                    mu_steps = 10):
    
    effect_sizes, guide_info = _get_effect_sizes(sdev, n_targets, n_guides_per_target, total_cells, vars_per_mu, mu_steps)
    
    covs, edit_rates = _sample_all_from_data(n_targets * n_guides_per_target, total_cells)
    guide_info["coverage"] = covs.astype(int)
    guide_info["edit_rate"] = edit_rates
    
    phenotypes = guide_info.reindex(guide_info.index.repeat(guide_info.coverage)).reset_index()
    phenotypes['edited'] = scipy.stats.bernoulli.rvs(phenotypes.edit_rate)
    phenotypes['reporter_edited'] = scipy.stats.bernoulli.rvs(phenotypes.edit_rate)
    phenotypes['phenotype_mean'] = phenotypes['edited']*phenotypes['effect_size']
    phenotypes['phenotype'] = scipy.stats.norm.rvs(loc = phenotypes.phenotype_mean, scale = 1)
    return(phenotypes)

def _aggregate_counts(df, prefix = None, group_by = ["target_id", "guide_num"]):
    counts = df.groupby(group_by).agg(
        {"edited":["count", "sum"], "reporter_edited":"sum"})
    
    counts.columns = ["guide", "target_edit", "reporter_edit"]
    if not prefix is None:
        counts.columns = counts.columns.map(lambda s: "{}_".format(prefix) + s)
    return(counts)

def get_sorted_cell_counts(phenotype_df, 
                           n_bulk_cells, 
                           top_quantile = 0.7, 
                           bot_quantile = 0.3, 
                          ):
    """Sort cells (phenotypes) into designated sorting scheme: bulk, top, bottom"""
    phenotype_df["sorted_top"] = phenotype_df.phenotype > phenotype_df.phenotype.quantile(top_quantile)
    phenotype_df["sorted_bot"] = phenotype_df.phenotype < phenotype_df.phenotype.quantile(bot_quantile)
    bulk = phenotype_df.sample(n = n_bulk_cells, replace = True)
    guide_info = phenotype_df[["target_id", "guide_num"]].drop_duplicates()

    phenotype_top = phenotype_df.loc[phenotype_df.sorted_top, :]
    top_counts = _aggregate_counts(phenotype_top, "top")

    phenotype_bot = phenotype_df.loc[phenotype_df.sorted_bot, :]
    bot_counts = _aggregate_counts(phenotype_bot, "bot")
    
    bulk_counts = _aggregate_counts(bulk, "bulk")
    
    sorted_counts = guide_info.merge(
        top_counts, on = ["target_id", "guide_num"]).merge(
        bot_counts, on = ["target_id", "guide_num"]).merge(
    bulk_counts, on = ["target_id", "guide_num"])
    sorted_counts = sorted_counts.fillna(0)
    return(sorted_counts)

def get_sorted_cell_counts2(phenotype_df, 
                            n_bulk_cells,
                            quantile_width = 0.2,
                            bin_labels = None,
                          ):
    """Sort cells (phenotypes) into designated sorting scheme: bulk, 1/2/3/4"""
    q = np.arange(0, 1.001, quantile_width)
    if bin_labels is None:
        bin_labels = np.arange(len(q) - 1)
    else:
        if len(bin_labels) != len(q) - 1: 
            raise ValueError("provided length of bin_labels doesn't match the number of bins made by quantile")      

    phenotype_df["sorted_bin"] = pd.qcut(phenotype_df.phenotype,
                                        q = q,
                                        labels = bin_labels)

    guide_info = phenotype_df[["target_id", "guide_num"]].drop_duplicates()

    counts = _aggregate_counts(phenotype_df, 
                               group_by = ["target_id", "guide_num", "sorted_bin"]
                              ).reset_index()
    
    bulk = phenotype_df.sample(n = n_bulk_cells, replace = True)
    bulk_counts = _aggregate_counts(bulk, "bulk")

    sorted_counts = guide_info.merge(counts, on = ["target_id", "guide_num"])
    sorted_counts = sorted_counts.pivot(index = ["target_id", "guide_num"], 
                       columns = "sorted_bin",
                       values = ["guide", "target_edit", "reporter_edit"])
    sorted_counts.columns = ['_'.join(col[::-1]).strip() for col in sorted_counts.columns.values]
    sorted_counts = sorted_counts.reset_index().merge(bulk_counts, 
                                                      on = ["target_id", "guide_num"])
    sorted_counts = sorted_counts.fillna(0)
    return(sorted_counts)

def get_sorted_read_counts(cell_counts : pd.DataFrame, 
                           n_read_sample = 2*10**6, 
                          sorting_bins = ["top", "bot", "bulk"]):
    measures = ["guide", "target_edit", "reporter_edit"]
    for sorting_bin in sorting_bins:
        for measure in measures:
            cell_counts["{}_p".format(sorting_bin)] = cell_counts[
                "{}_{}".format(sorting_bin, measure)
            ] / cell_counts["{}_{}".format(sorting_bin, measure)].sum()
            cell_counts.loc[:,"{}_{}".format(sorting_bin, measure)] = scipy.stats.binom.rvs(
                n_read_sample, cell_counts["{}_p".format(sorting_bin)])        
    return(cell_counts[
        ["target_id", "guide_num"] + 
        ["{}_{}".format(sb, m) for sb in sorting_bins for m in measures]]
          )
    
def simulate_reps(total_cells_per_rep = 10**5, 
                  nreps = 4, 
                  nreads_per_sample = 10**6, 
                  sorting_mode = "topbot"
                 ) -> List[pd.DataFrame]:
    res = []
    for rep in range(nreps):
        phenotype = simulate_phenotypes(total_cells = total_cells_per_rep)
        if sorting_mode == "topbot":
            cell_counts = get_sorted_cell_counts(
                phenotype, 
                n_bulk_cells = int(total_cells_per_rep / 3))
            read_counts = get_sorted_read_counts(cell_counts, nreads_per_sample)
        elif sorting_mode == "bins":
            sorting_bins = ["top", "high", "mid", "low", "bot"]
            cell_counts = get_sorted_cell_counts2(
                phenotype, 
                n_bulk_cells = int(total_cells_per_rep / 3),
                bin_labels = sorting_bins)
            read_counts = get_sorted_read_counts(
                cell_counts, 
                nreads_per_sample, 
                sorting_bins = sorting_bins
            )
        else:
            raise ValueError("Unknown sorting scheme provided.")
        
        #log2fc = get_log2fc(read_counts)
        res.append(read_counts)
    return(res)  