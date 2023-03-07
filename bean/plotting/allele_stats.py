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
    ax.set_xlabel("# alleles")
    ax.set_title(f"# Alleles per guide, raw (n={len(bdata.uns[allele_df_key])})")
    ax.set_ylabel("# guides")
    return ax


def plot_n_guides_per_edit(
    bdata, allele_df_key: str, allele_col: str = "aa_allele", ax=None
):
    a = bdata.uns[allele_df_key][["guide", allele_col]].copy()
    a["edits"] = a[allele_col].map(
        lambda a: list(a.aa_allele.edits)
        + list(map(lambda e: e.get_abs_edit(), a.nt_allele.edits))
    )
    edits_df = a.explode("edits")[["guide", "edits"]]
    edits_df["edits"] = edits_df.edits.map(str)
    edits_df = edits_df.loc[~edits_df.edits.str.startswith("-"), :].drop_duplicates()
    n_guides = edits_df.groupby("edits")["guide"].count()
    if ax is None:
        fig, ax = plt.subplots()
    ax.hist(n_guides, bins=np.arange(min(n_guides) - 0.5, max(n_guides) + 0.5))
    ax.set_xlabel("# guides")
    ax.set_title(f"# guides per edit (n={len(bdata.uns[allele_df_key])})")
    ax.set_ylabel("# edits")
    return ax
