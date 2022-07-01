import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np


def collect_guide_edit_counts(reps, conditions, guide_filename_pattern, reporter_filename_pattern, guide_info_filename, 
                             ref_base = "A", sample_base = "G"):
    counts_df = pd.read_csv(guide_info_filename)
    counts_df = counts_df.set_index("name")
    
    for rep in reps:
        for cond in conditions:
            res = pd.read_csv(guide_filename_pattern.format(r = rep, c = cond), delimiter="\t")
            res = res.set_index("guide_name")
            mapped_reads = res["read_counts"].sum()
            res["guide_RPM"] = res["read_counts"] / mapped_reads * 10**6
            res.columns = res.columns.map(lambda x: "{r}_{c}_".format(r = rep, c = cond) + x)
            assert all([x in counts_df.index.tolist() for x in res.index.tolist()])
            cols_to_join = res.columns[[0, 2, 3]]
            counts_df = counts_df.join(res[cols_to_join], how = "left")

            res = pd.read_csv(reporter_filename_pattern.format(r = rep, c = cond), delimiter="\t")
            res = res.set_index("name")
            all_edits = res.groupby('name')["count"].sum()
            all_edits.name = "{r}_{c}_all_edits".format(r = rep, c = cond)
            abe_edits = res.loc[(res.ref_base == ref_base) & (res.sample_base == sample_base),:].groupby('name')["count"].sum()
            abe_edits.name = "{r}_{c}_AG_edits".format(r = rep, c = cond)
            
            # count exact match in position
            all_reporter_edits_df = counts_df.join(res, on = "name", how = "left")
            all_reporter_edits_df = all_reporter_edits_df.reset_index()
            exact_edits = all_reporter_edits_df.loc[(all_reporter_edits_df["Target base position in gRNA"] == (all_reporter_edits_df.pos + 1)) & (all_reporter_edits_df.ref_base == ref_base) & (all_reporter_edits_df.sample_base == sample_base)][["name", "count"]].set_index("name", drop = True)["count"]
            exact_edits.name = "{r}_{c}_exact_edits".format(r = rep, c = cond)

            counts_df = counts_df.join(all_edits, how = "left")
            counts_df = counts_df.join(abe_edits, how = "left")
            counts_df = counts_df.join(exact_edits, how = "left")
            
            counts_df["{r}_{c}_edit_RPM".format(r = rep, c = cond)] = counts_df["{r}_{c}_exact_edits".format(r = rep, c = cond)] / mapped_reads * 10**6
            counts_df["{r}_{c}_edit_rate".format(r = rep, c = cond)] = counts_df["{r}_{c}_exact_edits".format(r = rep, c = cond)] / counts_df["{r}_{c}_read_counts".format(r = rep, c = cond)]
        
    counts_df = counts_df.fillna(0)
    return(counts_df)


def write_count_stats(count_file_name, reps, outfile_name, condition = "bulk", stat_fn = np.median):
    counts_df = pd.read_csv(count_file_name)
    edit_rates = counts_df.set_index("name")[["{r}_bulk_edit_rate".format(r = r) for r in reps]]
    edit_rates["edit_rate"] = edit_rates.apply(stat_fn, axis = 1)
    guide_counts = counts_df.set_index("name")[["{r}_bulk_read_counts".format(r = r) for r in reps]]
    guide_counts["guide_count"] = guide_counts.apply(stat_fn, axis = 1)
    edit_rates[["edit_rate"]].join(guide_counts[["guide_count"]]).to_csv(outfile_name)
    

def _log2FC(x, y):
    with np.errstate(divide='ignore', invalid='ignore'):
        return(np.log2((x+0.5)/(y+0.5)))

    
def _get_log2FC(df, prefix, cond1, cond2, measure, count_thres = 5, count_column = "read_counts"):
    res = _log2FC(df["{p}_{c}_{m}".format(p = prefix, c = cond1, m = measure)], 
                df["{p}_{c}_{m}".format(p = prefix, c = cond2, m = measure)])
    ignore_idx = np.where(np.logical_and(df["{p}_{c}_{cn}".format(p = prefix, c = cond1, cn = count_column)] < count_thres,
                         df["{p}_{c}_{cn}".format(p = prefix, c = cond2, cn = count_column)] < count_thres))
    res.iloc[ignore_idx] = np.nan
    return(res)
    
    
def get_log2FC(count_file_name, reps, cond1 = "top", cond2 = "bot", guide_count_measure = "RPM", 
              edit_count_measure = "edit_RPM", guide_count_filter_measure = "read_counts", 
              edit_count_filter_measure = "exact_edits"):
    
    counts_df = pd.read_csv(count_file_name)
    guide_enrichment = counts_df[["name"]]
    
    for rep in reps:
        top_vs_bulk = _get_log2FC(counts_df, rep, "top", "bulk", guide_count_measure, count_column = "read_counts")
        bot_vs_bulk = _get_log2FC(counts_df, rep, "bot", "bulk", guide_count_measure, count_column = "read_counts")
        bot_vs_top = _get_log2FC(counts_df, rep, "bot", "top", guide_count_measure, count_column = "read_counts")
        top_vs_bulk.name = rep + "_guide_top_vs_bulk"
        bot_vs_bulk.name = rep + "_guide_bot_vs_bulk"
        bot_vs_top.name = rep + "_guide_bot_vs_top"

        guide_enrichment = pd.concat((guide_enrichment, top_vs_bulk, bot_vs_bulk, bot_vs_top), axis = 1)

        # Reporter enrichment
        top_vs_bulk = _get_log2FC(counts_df, rep, "top", "bulk", "edit_RPM", count_column = edit_count_filter_measure)
        bot_vs_bulk = _get_log2FC(counts_df, rep, "bot", "bulk", "edit_RPM", count_column = edit_count_filter_measure)
        bot_vs_top = _get_log2FC(counts_df, rep, "bot", "top", "edit_RPM", count_column = edit_count_filter_measure)
        top_vs_bulk.name = rep + "_edit_top_vs_bulk"
        bot_vs_bulk.name = rep + "_edit_bot_vs_bulk"
        bot_vs_top.name = rep + "_edit_bot_vs_top"

        guide_enrichment = pd.concat((guide_enrichment, top_vs_bulk, bot_vs_bulk, bot_vs_top), axis = 1)
    return(guide_enrichment)
        
    

def plot_count_dist_ranked(count_file_name, reps, conditions, sample = "",
    figsize = None):
    if figsize is None:
        figsize = (3*len(conditions), 9)
    plt.rcParams.update({'font.size': 15})

    counts_df = pd.read_csv(count_file_name)

    figs, axs = plt.subplots(3, len(conditions), figsize=figsize)
    for j, c in enumerate(conditions):
        counts_df["RPM_sum"] = 0

        for r in reps:
            prefix = "{r}_{c}_".format(r = r, c = c)
            counts_df.RPM_sum += counts_df[prefix + "RPM"]
        counts_df = counts_df.sort_values("RPM_sum", ascending = False)
        counts_df = counts_df.reset_index(drop = True)

        # Panel 1: guide counts
        for r in reps:
            prefix = "{r}_{c}_".format(r = r, c = c)
            p = axs[0, j].plot(counts_df.index, counts_df[prefix + "read_counts"], ".", label = r, markersize = 0.5, alpha = 0.5)
        #lgnd = axs[0, j].legend()

        axs[0, j].set_title(c)
        axs[0, j].set_ylabel("guide counts")
        axs[0, j].set_xlabel("rank by guide count")

            # Panel 2: edit rates sorted by guide counts
        for r in reps:
            prefix = "{r}_{c}_".format(r = r, c = c)
            axs[1, j].plot(counts_df.index, counts_df[prefix + "edit_rate"], ".", markersize = 0.5, alpha = 0.5, label =r)
        #axs[1, j].set_title("LDLRCDS edit efficiency, sorted by total guide RPM ({})".format(c))
        axs[1, j].set_ylabel("edit rate")
        axs[1, j].set_xlabel("guide rank")
        axs[1, j].set_xlabel("rank by guide count")

        # Panel 3: Edit rates
        counts_df["mean_edit_rate"] = 0
        for r in reps:
            prefix = "{r}_{c}_".format(r = r, c = c)
            counts_df.mean_edit_rate += counts_df[prefix + "exact_edits"]
        counts_df.mean_edit_rate /= 4
        counts_df = counts_df.sort_values("mean_edit_rate", ascending = False)
        counts_df = counts_df.reset_index(drop = True)
        for r in reps:
            prefix = "{r}_{c}_".format(r = r, c = c)
            axs[2, j].plot(counts_df.index, counts_df[prefix + "exact_edits"], ".", markersize = 0.5, alpha = 0.5, label = r)
        #axs[2, j].set_title("LDLRCDS edit efficiency, sorted by mean edit rate ({})".format(c))
        axs[2, j].set_ylabel("exact edit counts")
        axs[2, j].set_xlabel("rank by edit count")

    if len(reps) > 1:
        lgnd = axs[1, len(conditions)-1].legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
        for i in range(len(reps)):
            lgnd.legendHandles[i]._legmarker.set_markersize(6)
    plt.suptitle("{} guide / edit counts".format(sample))
    plt.tight_layout()

    
def plot_count_dist(count_file_name, reps, conditions, sample = "",
    figsize = None):
    if figsize is None:
        figsize = (3*len(conditions), 9)
    plt.rcParams.update({'font.size': 15})
    counts_df = pd.read_csv(count_file_name)

    figs, axs = plt.subplots(3, len(conditions), figsize=figsize)
    for j, c in enumerate(conditions):
        counts_df["RPM_sum"] = 0

        for r in reps:
            prefix = "{r}_{c}_".format(r = r, c = c)
            counts_df.RPM_sum += counts_df[prefix + "RPM"]
        
        for r in reps:
            prefix = "{r}_{c}_".format(r = r, c = c)
            sns.kdeplot(ax = axs[0, j], x = counts_df[prefix + "read_counts"], label = r, alpha = 0.5).set(
                xlim=(0, max(counts_df[prefix + "read_counts"])))
        #lgnd = axs[0, j].legend()

        axs[0, j].set_title(c)
        axs[0, j].set_xlabel("guide counts")
        axs[0, 0].set_ylabel("Density")

            # Panel 2: edit rates sorted by guide counts
        for r in reps:
            prefix = "{r}_{c}_".format(r = r, c = c)
            sns.kdeplot(ax = axs[1, j], x = counts_df[prefix + "edit_rate"], alpha = 0.5, label =r).set(xlim=(0, 1))
        #axs[1, j].set_title("LDLRCDS edit efficiency, sorted by total guide RPM ({})".format(c))
        axs[1, j].set_xlabel("edit rate")
        axs[1, 0].set_ylabel("Density")

        # Panel 3: Edit rates
        counts_df["mean_edit_rate"] = 0
        for r in reps:
            prefix = "{r}_{c}_".format(r = r, c = c)
            counts_df.mean_edit_rate += counts_df[prefix + "exact_edits"]
        counts_df.mean_edit_rate /= 4
        counts_df = counts_df.sort_values("mean_edit_rate", ascending = False)
        counts_df = counts_df.reset_index(drop = True)
        for r in reps:
            prefix = "{r}_{c}_".format(r = r, c = c)
            sns.kdeplot(ax = axs[2, j], x = counts_df[prefix + "exact_edits"], alpha = 0.5, label = r).set(
                xlim=(0, max(counts_df[prefix + "exact_edits"])))
        axs[2, j].set_xlabel("exact edit counts")
        axs[2, 0].set_ylabel("Density")

    if len(reps) > 1:
        lgnd = axs[1, len(conditions)-1].legend(bbox_to_anchor=(1.04,0.5), 
            loc="center left", borderaxespad=0)

    plt.suptitle("{} count / edit rate distribution".format(sample))
    plt.tight_layout()

    
def plot_counts_vs_edit_rate(count_file_name, reps, conditions, 
    title = "Guide counts vs edit rate", outfile = None, figsize = None):
    if figsize is None:
        figsize = (3*len(conditions), 3)
    counts_df = pd.read_csv(count_file_name)
    fig, axs = plt.subplots(1, len(conditions), figsize = figsize)
    for j, c in enumerate(conditions):
        read_counts = pd.Series(dtype = 'float')
        edit_rates = pd.Series(dtype = 'float')
        for r in reps:
            prefix = "{r}_{c}".format(r = r, c = c)
            axs[j].plot(counts_df[prefix + "_read_counts"], counts_df[prefix + "_edit_rate"], ".", 
                label = prefix, markersize = 1, alpha = 0.1)
            read_counts = pd.concat((read_counts, counts_df[prefix + "_read_counts"]))
            edit_rates = pd.concat((edit_rates, counts_df[prefix + "_edit_rate"]))
        cor = round(read_counts.corr(edit_rates, "spearman"), 2)
        axs[j].text(max(read_counts)*0.4, max(edit_rates)*0.7, "cor_Sp={}".format(cor), fontsize = 10)
        axs[j].set_title(c)
        axs[j].set_ylabel("edit rate")
        axs[j].set_xlabel("guide counts")
    
    if len(reps) > 1:
        lgnd = axs[len(conditions)-1].legend(
            bbox_to_anchor=(1.04,0.45), 
            loc="center left", 
            borderaxespad=0)
        for i in range(len(lgnd.legendHandles)):
            lgnd.legendHandles[i]._legmarker.set_markersize(6)
            lgnd.legendHandles[i]._legmarker.set_alpha(1)
        
    plt.suptitle(title)
    plt.tight_layout()
    plt.show()
    if outfile: plt.savefig(outfile)
        

def _plot_replicate_pairplot(df, reps, measure = "read_counts"):
    count_df = df.iloc[:, np.append(np.where(df.columns == "name")[0], np.where(df.columns.map(lambda x: measure in x))[0])] # select columns with measure in their name
    count_df = pd.wide_to_long(count_df, reps, "name", "condition", suffix = r"(.*_{})".format(measure), sep = "_")
    count_df = count_df.reset_index()
    return(sns.pairplot(count_df, kind = "scatter", hue = "condition", plot_kws={"s": 0.5}))


def _plot_replicate_pairplot_foldchange(df, reps, measure = "guide_bot_vs_top"):
    count_df = df.iloc[:, np.append(np.where(df.columns == "name")[0], np.where(df.columns.map(lambda x: measure in x))[0])]
    count_df.columns = count_df.columns.map(lambda s: s.split("_" + measure)[0])
    count_df = count_df.set_index("name", drop = True)
    assert len(count_df.columns) == len(reps), "{}".format(count_df.columns)
    return(sns.pairplot(count_df, kind = "scatter"))


def _annotate_pearsonr(x, y, **kws):
    r, _ = stats.pearsonr(x, y)
    ax = plt.gca()
    # count how many annotations are already present
    n = len([c for c in ax.get_children() if isinstance(c, matplotlib.text.Annotation)])
    pos = (.1, .9 - .1*n)
    ax.annotate("{}:r={:.2f}".format(kws['label'].split("_")[0],r),
                xy=pos, xycoords=ax.transAxes)

    
def _annotate_pearsonr_foldchange(x, y, **kws):
    ignore_idx = np.isnan(x) | np.isnan(y)
    x = x[~ignore_idx]
    y = y[~ignore_idx]
    r, _ = stats.pearsonr(x, y)
    ax = plt.gca()
    
    # count how many annotations are already present
    n = len([c for c in ax.get_children() if isinstance(c, matplotlib.text.Annotation)])
    pos = (.1, .9 - .1*n)
    
    ax.annotate("r={:.2f}".format(r),
                xy=pos, xycoords=ax.transAxes)
    
    
def plot_replicate_read_count_correlation(reporter_file_name: str, reps, outfile_name = None, title = ""):
    read_count_df = pd.read_csv(reporter_file_name)
    pairplot_obj = _plot_replicate_pairplot(read_count_df, measure = "read_counts", reps = reps)
    pairplot_obj.map_lower(_annotate_pearsonr)
    #pairplot_obj.fig.suptitle(title)
    #pairplot_obj.fig.tight_layout()
    if outfile_name:
        plt.savefig(outfile_name)
        

def plot_replicate_enrichment_correlation(enrichment_df, count, reps, outfile_name = None, title = ""):
    assert count in ["edit", "guide"], "argument 'count' must be one of 'guide' or 'edit'."
    pairplot_obj = _plot_replicate_pairplot_foldchange(enrichment_df, measure = "{}_bot_vs_top".format(count), reps = reps)
    pairplot_obj.map_lower(_annotate_pearsonr_foldchange)
    pairplot_obj.fig.suptitle(title)
    pairplot_obj.fig.tight_layout()
    