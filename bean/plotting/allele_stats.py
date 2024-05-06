"""Plotting functions to describe allele/guide/edit stats for allele count information stored in ReporterScreen.uns[allele_df_key]."""
import numpy as np
import matplotlib.pyplot as plt


def plot_n_alleles_per_guide(
    bdata, allele_df_key: str, allele_col: str = "allele", ax=None
):
    guide_to_allele = dict(
        list(
            bdata.uns[allele_df_key][["guide", allele_col]].groupby("guide")[allele_col]
        )
    )
    lens = [
        len(guide_to_allele[g]) if g in guide_to_allele else 0
        for g in bdata.guides.index
    ]
    ax.hist(lens, bins=np.arange(min(lens) - 0.5, max(lens) + 0.5))
    ax.set_xlabel("# alleles per guide")
    ax.set_title(f"n_alleles={len(bdata.uns[allele_df_key])}")
    ax.set_ylabel("# guides")
    return ax


def plot_n_guides_per_edit(
    bdata, allele_df_key: str, allele_col: str = "aa_allele", ax=None
):
    a = bdata.uns[allele_df_key][["guide", allele_col]].copy()
    if allele_col == "aa_allele":
        a["edits"] = a[allele_col].map(
            lambda a: list(a.aa_allele.edits)
            + list(map(lambda e: e.get_abs_edit(), a.nt_allele.edits))
        )
    elif allele_col == "allele":
        a["edits"] = a[allele_col].map(lambda a: list(a.edits))
    edits_df = a.explode("edits")[["guide", "edits"]]
    edits_df["edits"] = edits_df.edits.map(str)
    edits_df = edits_df.loc[~edits_df.edits.str.startswith("-"), :].drop_duplicates()
    n_guides = edits_df.groupby("edits")["guide"].count()
    if ax is None:
        fig, ax = plt.subplots()
    if len(n_guides) > 0:
        ax.hist(n_guides, bins=np.arange(min(n_guides) - 0.5, max(n_guides) + 0.5))
    ax.set_xlabel("# guides per variant")
    ax.set_title(f"n_variants={len(edits_df)}")
    ax.set_ylabel("# edits")
    return ax


def plot_allele_stats(bdata, allele_df_keys, plot_save_path):
    fig = plt.figure(constrained_layout=True, figsize=(6, 3 * len(allele_df_keys)))
    fig.suptitle("Allele stats")

    subfigs = fig.subfigures(nrows=len(allele_df_keys), ncols=1)
    for row, subfig in enumerate(subfigs):
        key = allele_df_keys[row]
        subfig.suptitle(f"{key}")

        # create 1x3 subplots per subfig
        ax = subfig.subplots(nrows=1, ncols=2)
        plot_n_alleles_per_guide(bdata, key, bdata.uns[key].columns[1], ax[0])
        plot_n_guides_per_edit(bdata, key, bdata.uns[key].columns[1], ax[1])

    #plt.tight_layout()
    fig.savefig(plot_save_path, bbox_inches="tight")
