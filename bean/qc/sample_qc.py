"""Calculate sample quality"""

from typing import Literal
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

linestyles = ["solid", "dotted", "dashed", "dashdot"]


def plot_guide_edit_rates(
    bdata, ax=None, figsize=(5, 3), title="", n_bins=30, plot_normed: bool = True
):
    """Plot guide edit rates
    Args:
        plot_normed: Plot normalized edit rates. If False, plot total edit rates in editing window for tiling ReporterScreen with `.tiling == True`.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    if "edit_rate_norm" in bdata.guides.columns and plot_normed:
        sns.histplot(bdata.guides.edit_rate_norm, bins=n_bins)
    else:
        sns.histplot(bdata.guides.edit_rate, bins=n_bins)
    ax.set_title(title)
    ax.set_xlabel("Editing rate")
    return ax


def plot_sample_edit_rates(bdata, ax=None, figsize=(5, 3), title="", agg_method="mean"):
    if agg_method == "median":
        agg = np.nanmedian
    elif agg_method == "mean":
        agg = np.nanmean
    else:
        raise ValueError(
            "Invalid aggregation method provided. Please pass one of 'median' or 'mean'."
        )
    set_sample_edit_rates(bdata, agg_method)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    for i, sample in enumerate(bdata.samples.index):
        sns.kdeplot(
            bdata.layers["edit_rate"][:, i],
            label=f"{sample}({agg_method} {agg(bdata.layers['edit_rate'][:, i]):.2f})",
            linestyle=linestyles[(i // 10) % 4],
            clip=(0, 1),
        )
    sns.kdeplot(
        bdata.guides.edit_rate,
        linewidth=2,
        c="black",
        label=f"Bulk {agg_method}",
        clip=(0, 1),
    )
    ax.set_xlabel("Editing rate")
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1, 1))
    return ax


def set_sample_edit_rates(
    bdata, agg_method: Literal["mean", "median"] = "median"
) -> None:
    if "edit_rate" not in bdata.layers or bdata.layers["edit_rate"].max() == 0:
        raise ValueError(
            "ReporterScreen object doesn't have .edit_rate. Calculate target site edits first using .get_edit_rate(...)"
        )
    if agg_method == "median":
        bdata.samples[f"{agg_method}_editing_rate"] = np.nanmedian(
            bdata.layers["edit_rate"].copy(), axis=0
        )
    if agg_method == "mean":
        bdata.samples[f"{agg_method}_editing_rate"] = np.nanmean(
            bdata.layers["edit_rate"].copy(), axis=0
        )
