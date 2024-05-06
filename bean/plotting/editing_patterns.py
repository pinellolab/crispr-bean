from copy import deepcopy
from typing import Optional, Sequence, Dict, List, Union, Tuple
import re
from bean import ReporterScreen
from bean.framework.Edit import Edit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logomaker
from tqdm.auto import tqdm

BASES = ["A", "T", "C", "G"]
REVCOMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
other_from = {"A": "C", "T": "G", "C": "A", "G": "T"}


def _add_absent_edits(
    bdata,
    guide: str,
    edit_tbl: pd.DataFrame,
    target_base_edits: Dict[str, str],
):
    """If edit is not observed in editable position, add into the edit rate table."""

    if guide not in edit_tbl.guide.tolist():
        observed_rel_pos = []
    else:
        edited_db = edit_tbl.loc[edit_tbl.guide == guide, :]
        observed_rel_pos = edited_db.rel_pos.tolist()
    edits = []
    positions = []
    for edited_base, target_alt in target_base_edits.items():
        editable_positions = np.where(
            (np.array(list(bdata.guides.loc[guide, "reporter"])) == edited_base)
        )[0]
        for editable_pos in editable_positions:
            if editable_pos not in observed_rel_pos:
                edits.append(
                    Edit(
                        editable_pos,
                        edited_base,
                        target_alt,
                    )
                )
                positions.append(editable_pos)
    edit_tbl = pd.concat(
        [
            edit_tbl,
            pd.DataFrame(
                {
                    "guide": guide,
                    "edit": edits,
                    "rel_pos": positions,
                    "rep_median": 0.0,
                    "rep_mean": 0.0,
                },
            ),
        ],
        ignore_index=True,
    )
    return edit_tbl


def get_edit_rates(
    bdata,
    edit_count_key="edit_counts",
    adjust_spacer_pos: bool = True,
    reporter_column: str = "reporter",
):
    """
    Obtain position- and context-wise editing rate (context: base preceding the target base position).
    Args
    bdata (bean.reporterScreen): input ReporterScreen object to be used for analyzing posiiton-wide editing rate.
    edit_count_key: key in bdata.uns with edit count DataFrame to be used for analyzing posiiton-wide editing rate.
    """
    edit_rates = bdata.get_normalized_allele_counts(bdata.uns[edit_count_key])
    edit_rates_agg = edit_rates[["guide", "edit"]]
    edit_rates_agg = edit_rates_agg.loc[
        edit_rates_agg.guide.map(lambda s: "CONTROL" not in s)
    ].reset_index(drop=True)
    edit_rates_agg["rep_median"] = edit_rates.iloc[:, 2:].median(axis=1)
    edit_rates_agg["rep_mean"] = edit_rates.iloc[:, 2:].mean(axis=1)
    edit_rates_agg["rel_pos"] = edit_rates_agg.edit.map(lambda e: e.rel_pos).astype(int)

    for guide in tqdm(
        edit_rates_agg.guide.unique(),
        desc="Calibrating edits in editable positions...",
    ):
        edit_rates_agg = _add_absent_edits(
            bdata,
            guide,
            edit_rates_agg,
            target_base_edits=bdata.target_base_changes,
        )

    if adjust_spacer_pos:
        if "guide_len" not in bdata.guides.columns:
            bdata.guides["guide_len"] = bdata.guides.sequence.map(len)
        start_pad = (
            32 - 6 - bdata.guides.guide_len[edit_rates_agg.guide].reset_index(drop=True)
        )
        edit_rates_agg["spacer_pos"] = (edit_rates_agg.rel_pos - start_pad).astype(int)
        edit_rates_agg = edit_rates_agg[
            (edit_rates_agg.spacer_pos >= 0) & (edit_rates_agg.spacer_pos < 20)
        ].reset_index(drop=True)
        edit_rates_agg.loc[:, "spacer_pos"] = edit_rates_agg["spacer_pos"] + 1
    else:
        edit_rates_agg["spacer_pos"] = edit_rates_agg.rel_pos + 1
    edit_rates_agg["base_change"] = edit_rates_agg.edit.map(
        lambda e: e.get_base_change()
    )
    edit_rates_agg.rel_pos = edit_rates_agg.rel_pos.astype(int)
    edit_rates_agg["context"] = edit_rates_agg.apply(
        lambda row: bdata.guides.loc[row.guide, reporter_column][
            row.rel_pos - 1 : row.rel_pos + 1
        ],
        axis=1,
    )
    return edit_rates_agg


def plot_by_pos_context(
    edit_rates_df: pd.DataFrame,
    target_base="A",
):
    target_base_alt = {"A": "G", "C": "T"}[target_base]
    context_to_offset_map = {
        f"A{target_base}": -0.4,
        f"G{target_base}": -0.2,
        f"C{target_base}": 0,
        f"T{target_base}": 0.2,
    }
    edit_rates_df["spacer_pos_ctxt"] = (
        edit_rates_df.spacer_pos + edit_rates_df.context.map(context_to_offset_map)
    )
    fig, ax = plt.subplots(figsize=(6, 3))
    sns.scatterplot(
        edit_rates_df.loc[
            (edit_rates_df.base_change == f"{target_base}>{target_base_alt}")
        ],
        x="spacer_pos_ctxt",
        y="rep_mean",
        hue="context",
        alpha=0.3,
        ax=ax,
        s=5,
        rasterized=True,
    )
    ax.legend(bbox_to_anchor=(1.02, 0.5), loc="center left", title="Context")
    ax.set_xlabel("Protospacer position")
    ax.set_ylabel(f"{target_base}>{target_base_alt} editing rate")
    ax.set_ylim((0, 1))
    return ax


def _get_complementary_base_change(base_change_str: str):
    if not re.fullmatch(r"[ATCG]>[ATCG]", base_change_str):
        raise ValueError(
            f"Input argument {base_change_str} doesn't conform to the valid format ex. 'A>G'"
        )
    base_from = base_change_str[0]
    base_to = base_change_str[-1]
    return f"{REVCOMP[base_from]}>{REVCOMP[base_to]}"


def _combine_complementary_base_changes(edit_rate_df: pd.DataFrame, show_list=None):
    """Combine complementary list"""

    edit_rate_df_reduced = edit_rate_df[show_list].copy()
    for base_change in show_list:
        rev_base_change = _get_complementary_base_change(base_change)
        if rev_base_change in edit_rate_df.columns:
            edit_rate_df_reduced.loc[:, base_change] += (
                edit_rate_df[base_change] + edit_rate_df[rev_base_change]
            )
        # edit_rate_df_reduced.loc[:,base_change] = edit_rate_df_reduced.loc[:,base_change]/2
    return edit_rate_df_reduced


def _get_norm_rates_df(
    bdata,
    edit_rates_df=None,
    edit_count_key="edit_counts",
    base_changes: Optional[Sequence[str]] = None,
):
    change_by_pos = pd.pivot(
        edit_rates_df[["base_change", "spacer_pos", "rep_mean"]]
        .groupby(["base_change", "spacer_pos"])
        .sum()["rep_mean"]
        .reset_index(),
        index="spacer_pos",
        columns="base_change",
        values="rep_mean",
    ).fillna(0)

    norm_matrix = pd.DataFrame(index=change_by_pos.index, columns=BASES)
    for pos in norm_matrix.index:
        pos_base = bdata.guides.reporter.map(
            lambda s: s[pos] if pos < len(s) else " "
        ).values
        for b in BASES:
            norm_matrix.loc[pos, b] = (pos_base == b).sum()
    ref_bases = change_by_pos.columns.map(lambda s: s.split(">")[0])
    change_by_pos = change_by_pos.loc[:, ref_bases.isin(BASES)]
    norm_rate = (
        change_by_pos
        / norm_matrix.loc[
            :, change_by_pos.columns.map(lambda s: s.split(">")[0])
        ].values
    )
    # norm_rate_reduced = _combine_complementary_base_changes(norm_rate).astype(float)
    return norm_rate.astype(float)[base_changes]  # _reduced


def _get_possible_changes_from_target_base(target_basechange: str) -> List[str]:
    """Return base changes strings (ex. A>C) for given the same reference base to edit from to the input target_basechange. ex) returns ['A>C', 'A>T'] given input 'A>G'."""
    if not re.fullmatch(r"[ATCG]>[ATCG]", target_basechange):
        raise ValueError(
            f"Input argument {target_basechange} doesn't conform to the valid format ex. 'A>G'"
        )
    base_from = target_basechange[0]
    base_to = target_basechange[-1]
    return [f"{base_from}>{b}" for b in BASES if b not in [base_from, base_to]]


def plot_by_pos_behive(
    norm_rates_df: Optional[pd.DataFrame] = None,
    bdata=None,
    edit_count_key: str = "edit_counts",
    target_basechanges: Optional[Dict[str, str]] = None,
    normalize=False,
):
    """Plot position-wise editing pattern as in BE-Hive.
    Args
    bdata (bean.ReporterScreen): input ReporterScreen object to be used for analyzing posiiton-wide editing rate.
    edit_count_key: key in bdata.uns with edit count DataFrame to be used for analyzing posiiton-wide editing rate.
    df_to_draw: input DataFrame with position-wise editing rates (n_position by n_basechange)
    norm_rates_df: pd.DataFrame with normalized editing rate with columns ['guide', 'edit', 'median/mean editing rate', 'spacer_pos']
    normalize: Normalize the editing rate relative to the max editing rate by position (as 100).

    """
    if target_basechanges is None:
        if bdata is None:
            raise ValueError("Target base change not provided.")
        target_basechange = bdata.target_base_changes
    fig, axes = plt.subplots(
        1,
        len(target_basechanges.keys()),
        figsize=(3 * len(target_basechanges.keys()), 7),
    )
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    dfs = []
    for i, (edited_base, alt_base) in enumerate(target_basechanges.items()):  # type: ignore
        target_basechange = f"{edited_base}>{alt_base}"
        print(target_basechange)
        if target_basechange == "A>G":
            nonref_base_changes = ["C>T", "C>G"]
        elif target_basechange == "C>T":
            nonref_base_changes = ["G>A", "G>C"]
        else:
            print("No non-ref base changes specified. not drawing them")
            nonref_base_changes = []

        ref_other_changes = _get_possible_changes_from_target_base(target_basechange)

        df_to_draw = _get_norm_rates_df(
            bdata,
            norm_rates_df,
            edit_count_key,
            base_changes=ref_other_changes
            + [
                target_basechange,
            ]
            + nonref_base_changes,
        )

        vmax = df_to_draw.max().max()
        if normalize:
            df_to_draw = df_to_draw / vmax * 100
            vmax = 100

        target_data = df_to_draw.copy()
        target_data.loc[:, target_data.columns != target_basechange] = np.nan
        sns.heatmap(
            target_data,
            ax=axes[i],
            annot=True,
            cmap="Reds",
            vmax=vmax,
            cbar=False,
            vmin=-0.03,
            fmt=".0f" if normalize else ".2g",
            annot_kws={"fontsize": 8},
        )

        ref_data = df_to_draw.copy()
        ref_data.loc[
            :,
            ~ref_data.columns.isin(ref_other_changes),
        ] = np.nan
        sns.heatmap(
            ref_data,
            ax=axes[i],
            annot=True,
            cmap="Blues",
            vmax=vmax,
            cbar=False,
            fmt=".1g",
            vmin=-0.03,
            annot_kws={"fontsize": 8},
        )

        nonref_data = df_to_draw.copy()
        nonref_data.loc[:, ~nonref_data.columns.isin(nonref_base_changes)] = np.nan
        sns.heatmap(
            nonref_data,
            ax=axes[i],
            annot=True,
            cmap="Greys",
            vmax=vmax,
            cbar=False,
            fmt=".1g",
            vmin=-0.03,
            annot_kws={"fontsize": 8},
        )
        axes[i].set_ylabel("Protospacer position")
        dfs.append(df_to_draw)
    df = pd.concat(dfs, axis=1)
    df = df.loc[:, ~df.columns.duplicated()].copy()
    return axes, df


def get_position_by_pam_rates(
    bdata, edit_rates_df: pd.DataFrame, target_base_change: str, pam_col="5-nt PAM"
):
    edit_rates_df["pam"] = bdata.guides.loc[edit_rates_df.guide, pam_col].reset_index(
        drop=True
    )
    edit_rates_df["pam23"] = edit_rates_df.pam.map(lambda s: s[1:3])
    return pd.pivot(
        edit_rates_df.loc[
            (edit_rates_df.base_change == target_base_change),
            ["rep_mean", "pam23", "spacer_pos"],
        ]
        .groupby(["pam23", "spacer_pos"])
        .mean()
        .reset_index(),
        index="pam23",
        columns="spacer_pos",
        values="rep_mean",
    )


def plot_by_pos_pam(
    bdata: ReporterScreen,
    edit_rates_df: pd.DataFrame,
    target_base_changes: List[str],
    pam_col: str = "5-nt PAM",
    axes=None,
    figsize: Tuple[float, float] = (6, 4),
    save_fig: bool = False,
    save_path: Optional[str] = None,
) -> pd.DataFrame:
    """Plot editing efficiency for each position x PAM combination."""
    if axes is None:
        fig, axes = plt.subplots(
            1,
            len(target_base_changes),
            figsize=(figsize[0] * len(target_base_changes), figsize[1]),
        )
    if len(target_base_changes) == 1:
        axes = np.array([axes])
    pos_by_pam_basechanges = []
    for i, target_base_change in enumerate(target_base_changes):
        edit_rates_df_base = edit_rates_df.loc[
            edit_rates_df.base_change == target_base_change, :
        ].reset_index(drop=True)
        pos_by_pam = get_position_by_pam_rates(
            bdata, edit_rates_df_base, target_base_change, pam_col
        )
        sns.heatmap(pos_by_pam, ax=axes[i], cmap="Blues")
        axes[i].set_yticklabels(axes[i].get_yticklabels(), rotation=0)
        axes[i].set_title(target_base_change)
        pos_by_pam["base_change"] = target_base_change
        pos_by_pam_basechanges.append(pos_by_pam)
    if save_fig:
        if save_path is None:
            save_path = "bean_profile_pos_by_pam.pdf"
        plt.savefig(save_path, bbox_inches="tight")
    return pd.concat(pos_by_pam_basechanges, axis=0)


def get_pam_preference(
    edit_rates_df: pd.DataFrame,
    edit_start_pos: int = 3,
    edit_end_pos: int = 8,
    norm=True,
):
    window_edits = deepcopy(
        edit_rates_df.loc[
            (edit_rates_df.spacer_pos >= edit_start_pos)
            & (edit_rates_df.spacer_pos < edit_end_pos)
            & (edit_rates_df.base_change == "A>G")
        ]
    )
    if norm:
        ctxt_mean = window_edits.groupby("context")["rep_mean"].mean()
        ctxt_mean = ctxt_mean / ctxt_mean.max()
        pos_mean = window_edits.groupby("spacer_pos")["rep_mean"].mean()
        pos_mean = pos_mean / pos_mean.max()
        print("Normalizing per context mean")
        window_edits["norm_rate"] = (
            window_edits["rep_mean"]
            / ctxt_mean.loc[window_edits.context].values
            / pos_mean.loc[window_edits.spacer_pos].values
        )
    else:
        window_edits["norm_rate"] = window_edits.rep_mean
    pam_pref = pd.pivot(
        window_edits.groupby(["pam2", "pam3"]).mean().reset_index(),
        columns="pam3",
        index="pam2",
        values="norm_rate",
    )
    return pam_pref


def plot_pam_preference(
    edit_rates_df,
    bdata=None,
    pam_col="5-nt PAM",
    edit_start_pos: int = 3,
    edit_end_pos: int = 8,
    ax=None,
    norm=True,
):
    """ """
    if "pam2" not in edit_rates_df.columns or "pam3" not in edit_rates_df.columns:
        edit_rates_df["pam"] = bdata.guides.loc[
            edit_rates_df.guide, pam_col
        ].reset_index(drop=True)
        edit_rates_df["pam2"] = edit_rates_df.pam.map(lambda s: s[1])
        edit_rates_df["pam3"] = edit_rates_df.pam.map(lambda s: s[2])
    pam_pref = get_pam_preference(
        edit_rates_df, edit_start_pos, edit_end_pos, norm=norm
    )
    fig, ax = plt.subplots(figsize=(3, 3))
    sns.heatmap(data=pam_pref, annot=True, cmap="Blues")
    ax.set_box_aspect(1)
    ax.set_xticklabels([f"NN{t.get_text()}" for t in ax.get_xticklabels()])
    ax.set_yticklabels([f"N{t.get_text()}N" for t in ax.get_yticklabels()])
    ax.set(xlabel=None, ylabel=None)
    ax.set_title("PAM preference")
    return ax


def plot_context_specificity(
    bdata: ReporterScreen,
    edit_rates_df: pd.DataFrame,
    target_base_changes: List[str],
    window: Tuple[int, int],
    axes=None,
    figsize: Tuple[float, float] = (3, 5),
    save_fig: bool = False,
    save_path: Optional[str] = None,
) -> pd.DataFrame:
    window_start, window_end = window
    if axes is None:
        fig, axes = plt.subplots(
            1,
            len(target_base_changes),
            figsize=(figsize[0] * len(target_base_changes), figsize[1]),
        )
    if len(target_base_changes) == 1:
        axes = np.array([axes])
    save_tbls = []
    for j, target_base_change in enumerate(target_base_changes):
        edit_rates_df_base = edit_rates_df.loc[
            edit_rates_df.base_change == target_base_change, :
        ].reset_index(drop=True)
        edit_rates_df_base_window = edit_rates_df_base.loc[
            (edit_rates_df_base.spacer_pos >= window_start)
            & (edit_rates_df_base.spacer_pos <= window_end)
        ].copy()
        edit_rates_df_base_window["context"] = edit_rates_df_base_window.apply(
            lambda row: bdata.guides.loc[row.guide, "reporter"][
                row.rel_pos - 1 : row.rel_pos + 2
            ],
            axis=1,
        )
        context_mean_edit_rate = {}
        bases = ["A", "C", "G", "T"]

        for i in range(3):
            edit_rates_df_base_window[f"context_{i}"] = (
                edit_rates_df_base_window.context.map(lambda s: s[i])
            )
            context_mean_edit_rate[i] = edit_rates_df_base_window.groupby(
                f"context_{i}"
            )["rep_mean"].mean()
            context_mean_edit_rate[i] = (
                context_mean_edit_rate[i].reindex(bases).fillna(0)
            )

        target_df = pd.DataFrame(
            [int(b == target_base_change[0]) for b in bases], index=bases
        )
        ic_tbl = pd.concat(
            [
                context_mean_edit_rate[0] / context_mean_edit_rate[0].sum(),
                target_df,
                context_mean_edit_rate[2] / context_mean_edit_rate[2].sum(),
            ],
            axis=1,
        ).T
        save_tbl = ic_tbl.transpose()
        save_tbl["base_change"] = target_base_change
        save_tbls.append(save_tbl)
        ic_tbl.index = [-1, 0, 1]
        logomaker.Logo(ic_tbl.astype(float), ax=axes[j])
        axes[j].set_ylabel("Relative frequency")
        axes[j].set_title(target_base_change)
    if save_fig:
        if save_path is None:
            save_path = (
                f"bean_profile_context_preference_{window_start}_{window_end}.pdf"
            )
        fig.savefig(save_path, bbox_inches="tight")
    return pd.concat(save_tbls, axis=0)
