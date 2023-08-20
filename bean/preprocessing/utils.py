from typing import Dict
import torch
import numpy as np
import pyBigWig
import pandas as pd
import anndata as ad
import bean as be


class Alias:
    def __init__(self, source_name):
        self.source_name = source_name

    def __get__(self, obj, objtype=None):
        if obj is None:
            # Class lookup, return descriptor
            return self
        return getattr(obj, self.source_name)

    def __set__(self, obj, value):
        setattr(obj, self.source_name, value)


def _get_accessibility_single(
    pos: int,
    track,
    chrom: str = "chr19",
    guide_start_pos=0,
    half_window_size: int = 100,
):
    """Obtain mean log accessibility signal for a given genomic position, padded by `half_window_size`.

    Args
    --
    pos: genomic position
    track: pyBigWig.bigWigFile object with accessibility signal information
    chrom: chromosome name in the bigWig file.
    guide_start_pos: offset to the genomic position to be added
    half_window_size: half-size of the window centered at pos + guide_start_pos of which the signal would be taken
    """
    if half_window_size < 0:
        raise ValueError("Window size must be non-negative.")
    if pos == "control" or np.isnan(pos):
        return np.nan
    try:
        return np.exp(
            np.nanmean(
                np.log(
                    np.array(
                        track.values(
                            chrom,
                            int(guide_start_pos + pos - half_window_size),
                            int(guide_start_pos + pos + half_window_size),
                        )
                    )
                    + 1.0
                )
            )
        )
    except Exception as e:
        print(e)
        return np.nan


def get_accessibility_guides(
    accessibility_bw_path: str, guide_info: pd.DataFrame, half_window_size: int = 100
) -> torch.Tensor:
    """Gets guide accessibility given a dataframe with `genomic_pos`, `chr` column. Returns the mean logged signal of the guide position padded by `half_window_size`.

    Args:
    accessibility_bw_path: bigWig file path with ATAC-seq signal.
    guide_info: screen.guides dataframe with `genomic_pos`, `chr` column.
    half_window_size: For each guide, accessibility signal its position padded by this half_window_size upstream and downstream
    """
    acc = pyBigWig.open(accessibility_bw_path)
    guide_accessibility = torch.as_tensor(
        guide_info.apply(
            lambda row: _get_accessibility_single(
                row.genomic_pos,
                acc,
                chrom=row.chr,
                guide_start_pos=0,
                half_window_size=half_window_size,
            ),
            axis=1,
        )
    )
    guide_accessibility[torch.isnan(guide_accessibility)] = torch.nanmedian(
        guide_accessibility
    )
    return guide_accessibility


def get_edit_to_index_dict(cnalleles: pd.Series) -> Dict[str, int]:
    """Returns mapping (dictionary) of edit to unique index
    Returns Dict[edit(str) -> index(int)] for edits
    ---
    Arguments
    cnalleles -- pd.Series object of CodingNoncodingAllele objects.
    """
    edit_lists = cnalleles.map(
        lambda a: list(a.aa_allele.edits) + list(a.nt_allele.edits)
    )
    edits = pd.Series(
        pd.Series(
            [e.get_abs_edit() for l in edit_lists.tolist() for e in l], dtype="object"
        ).unique()
    )
    edits_to_index = pd.Series(edits)
    edits_to_index = (
        edits_to_index.reset_index()
        .rename(columns={0: "edit"})
        .set_index("edit")["index"]
        .to_dict()
    )
    return edits_to_index


def _insert_row_(row_number: int, df: pd.DataFrame, row):
    # https://www.geeksforgeeks.org/insert-row-at-given-position-in-pandas-dataframe/
    # Function to insert row in the dataframe
    df1 = df.iloc[0:row_number]
    df2 = df.iloc[row_number:]
    return pd.concat([df1, row, df2], axis=0)


def _insert_row_to_obs(adata: be.ReporterScreen, row_number: int, row_value, index):
    obs = _insert_row_(row_number, adata.guides, row_value)
    X = np.insert(adata.X, row_number, 0, axis=0)
    layers = {}
    for layer in adata.layers.keys():
        layers[layer] = np.insert(adata.layers[layer], row_number, 0, axis=0)
    if "repguide_mask" in adata.uns:
        adata.uns["repguide_mask"] = pd.concat(
            [
                adata.uns["repguide_mask"].iloc[:row_number, :],
                pd.DataFrame(
                    columns=adata.uns["repguide_mask"].columns,
                    index=[index],
                ).fillna(0),
                adata.uns["repguide_mask"].iloc[row_number:, :],
            ],
            axis=0,
        )
    adata = be.ReporterScreen(
        X=X,
        guides=obs,
        samples=adata.samples,
        uns=adata.uns,
        obsm=adata.obsm,
        varm=adata.varm,
        layers=layers,
    )
    return adata


def check_consecutive_targets(target_list, guide_per_target_counts=5):
    item_count = 0
    prev_item = ""
    for item in target_list:
        if item_count == 0:
            item_count += 1
            prev_item = item
        else:
            if item == prev_item:
                item_count += 1
            else:
                assert (
                    item_count == guide_per_target_counts
                ), f"not all targets are consecutive in len {guide_per_target_counts}"
                item_count = 1
                prev_item = item


def _assign_rep_ids_and_sort(
    screen: be.ReporterScreen, rep_col: str, condition_column: str = None
) -> be.ReporterScreen:
    """Assign replicate IDs to samples and sort them accordingly."""
    for i, rep in enumerate(sorted(screen.samples[rep_col].unique())):
        screen.samples.loc[screen.samples[rep_col] == rep, f"{rep_col}_id"] = i
        if condition_column is None:
            sort_key = f"{rep_col}_id"
        else:
            sort_key = [f"{rep_col}_id", f"{condition_column}_id"]
        screen = screen[
            :,
            screen.samples.sort_values(sort_key).index,
        ]
    return screen


def _obtain_effective_edit_rate(ndata, count_thres=10):
    """Calculates effective editing rate of each variant for tiling screen.

    Allele count
    Args
    --
    count_thres: Per-replicate allele editing rate is ignored if X_bcmatch is lower than this threshold.
    """
    allele_rates = ndata.allele_counts_control / ndata.X_bcmatch_control[:, :, :, None]
    allele_rates[
        (ndata.X_bcmatch_control < count_thres)[:, :, :, None].expand(
            allele_rates.shape
        )
    ] = torch.nan
    mean_allele_rates = allele_rates.nanmean(axis=(0, 1))[:, 1:]
    assert mean_allele_rates.shape == (ndata.n_guides, ndata.n_max_alleles - 1)
    allele_to_edit_normed = (
        ndata.allele_to_edit / ndata.allele_to_edit.nansum(axis=-1)[:, :, None]
    )

    allele_to_edit_normed.shape == (
        ndata.n_guides,
        ndata.n_max_alleles - 1,
        ndata.n_edits,
    )
    return (mean_allele_rates[:, :, None] * allele_to_edit_normed).nansum(axis=(0, 1))


def _obtain_n_guides_alleles_per_variant(ndata):
    return (ndata.allele_to_edit.sum(axis=1) > 0).sum(axis=0)
