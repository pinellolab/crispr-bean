from copy import deepcopy
from typing import Collection, Iterable, List, Optional, Union, Sequence, Literal
import re
import anndata as ad
import numpy as np
import pandas as pd
from perturb_tools import Screen

from ..annotate._supporting_fn import (
    filter_allele_by_base,
    filter_allele_by_pos,
    get_aa_alleles,
)
from .AminoAcidEdit import AminoAcidEdit, CodingNoncodingAllele
from .Edit import Allele, Edit
from ..annotate.filter_alleles import (
    _distribute_alleles_to_filtered,
    _map_alleles_to_filtered,
)


def _get_counts(
    filename_pattern: str,
    guides: pd.DataFrame,
    reps: Iterable[str],
    samples: Iterable[str],
) -> pd.DataFrame:
    """Read multiple .csv file with "guide_name,read_counts" columns across replicates and samples.
    File names are dictated by filename_pattern and reps, condiitons argument provided.

    Args
    --
    filename_pattern: format string with 'r' as replicate and 'c' as condition. Will be formatted as filename_pattern.format(r=replicate, c=condition).
    guides: Guide information that will be used as reference
    reps: List of replicate string to be used for formatting
    samples: List of condition strings to be used for formatting
    """
    for rep in reps:
        for spl in samples:
            res = pd.read_csv(filename_pattern.format(r=rep, c=spl), delimiter="\t")
            res = res.set_index("guide_name")[["read_counts"]]
            res.columns = res.columns.map(lambda x: "{rep}_{spl}")
            guides = guides.join(res, how="left")
        guides = guides.fillna(0)
        return guides


def _get_edits(
    filename_pattern: str,
    guide_info: pd.DataFrame,
    reps: Iterable[str],
    samples: Iterable[str],
    count_exact: bool = True,
) -> pd.DataFrame:
    """ """
    edits = pd.DataFrame(index=(guide_info.index))
    for rep in reps:
        for spl in samples:
            res = pd.read_csv(filename_pattern.format(r=rep, c=spl), delimiter="\t")
            res = res.set_index("name")
            if count_exact:
                res_info = guide_info.join(res, how="left").reset_index()
                this_edit = res_info.loc[
                    (
                        (res_info["Target base position in gRNA"] == res_info.pos + 1)
                        & (res_info.ref_base == "A")
                        & (res_info.sample_base == "G")
                    )
                ][["name", "count"]].set_index("name", drop=True)["count"]
            else:
                this_edit = (
                    res.loc[(res.ref_base == "A") & (res.sample_base == "G"), :]
                    .groupby("name")["count"]
                    .sum()
                )
            this_edit.name = "{r}_{c}".format(r=rep, c=spl)
            edits = edits.join(this_edit, how="left")
        edits = edits.fillna(0)
        return edits


class ReporterScreen(Screen):
    def __init__(
        self,
        X=None,
        X_edit=None,
        X_bcmatch=None,
        target_base_change: Optional[str] = None,
        replicate_label: str = "rep",
        condition_label: str = "bin",
        tiling: Optional[bool] = None,
        *args,
        **kwargs,
    ):
        (super().__init__)(X=X, *args, **kwargs)
        if X_edit is not None:
            self.layers["edits"] = X_edit
        if X_bcmatch is not None:
            self.layers["X_bcmatch"] = X_bcmatch
        for k, df in self.uns.items():
            if not isinstance(df, pd.DataFrame):
                if k == "sample_covariates" and not isinstance(df, list):
                    self.uns[k] = df.tolist()
                continue
            if "guide" in df.columns and len(df) > 0:
                if (
                    "allele" in df.columns
                    and (not isinstance(df["allele"].iloc[0], Allele))
                    and any(df["allele"].map(Allele.match_str))
                ):
                    self.uns[k].loc[:, "allele"] = self.uns[k].allele.map(
                        lambda s: Allele.from_str(s)
                    )
                if (
                    "edit" in df.columns
                    and (not isinstance(df["edit"].iloc[0], Edit))
                    and any(df["edit"].map(Edit.match_str))
                ):
                    self.uns[k].loc[:, "edit"] = self.uns[k].edit.map(
                        lambda s: Edit.from_str(s)
                    )
                if "reporter_allele" in df.columns and "guide_allele" in df.columns:
                    self.uns[k].loc[:, "reporter_allele"] = self.uns[
                        k
                    ].reporter_allele.map(lambda s: Allele.from_str(s))
                    self.uns[k].loc[:, "guide_allele"] = self.uns[k].guide_allele.map(
                        lambda s: Allele.from_str(s)
                    )
                if "aa_allele" in df.columns and any(
                    df["aa_allele"].map(CodingNoncodingAllele.match_str)
                ):
                    self.uns[k].loc[:, "aa_allele"] = self.uns[k].aa_allele.map(
                        lambda s: CodingNoncodingAllele.from_str(s)
                    )
        if target_base_change is not None:
            if not re.fullmatch(r"[ACTG]>[ACTG]", target_base_change):
                raise ValueError(
                    f"target_base_change {target_base_change} doesn't match the allowed base change. Feed in valid base change string ex) 'A>G', 'C>T'"
                )
            self.uns["target_base_change"] = target_base_change
        if tiling is not None:
            self.uns["tiling"] = tiling

    @property
    def X_edits(self):
        return self.layers["edits"]

    @property
    def X_bcmatch(self):
        return self.layers["X_bcmatch"]

    @property
    def edit_tables(self):
        return {k: self.uns[k] for k in self.uns.keys() if "edit" in k}

    @property
    def allele_tables(self):
        return {k: self.uns[k] for k in self.uns.keys() if "allele" in k}

    @property
    def base_edited_from(self):
        return self.uns["target_base_change"][0]

    @property
    def base_edited_to(self):
        return self.uns["target_base_change"][-1]

    @property
    def target_base_change(self):
        return self.uns["target_base_change"]

    @property
    def tiling(self):
        return self.uns["tiling"]

    @classmethod
    def from_file_paths(
        cls,
        reps: List[str] = None,
        samples: List[str] = None,
        guide_info_file_name: str = None,
        guide_count_filenames: str = None,
        guide_bcmatched_count_filenames: str = None,
        edit_count_filenames: str = None,
    ):
        guide_info = pd.read_csv(guide_info_file_name).set_index("name")
        guides = pd.DataFrame(index=(pd.read_csv(guide_info_file_name)["name"]))
        guides_lenient = _get_counts(guide_count_filenames, guides, reps, samples)
        edits_ag = _get_edits(
            edit_count_filenames, guide_info, reps, samples, count_exact=False
        )
        edits_exact = _get_edits(edit_count_filenames, guide_info, reps, samples)
        if guide_bcmatched_count_filenames is not None:
            guides_bcmatch = _get_counts(
                guide_bcmatched_count_filenames, guides, reps, samples
            )
            repscreen = cls(
                guides_lenient, edits_exact, X_bcmatch=guides_bcmatch, guides=guide_info
            )
        else:
            repscreen = cls(guides_lenient, edits_exact, guides=guide_info)
        repscreen.layers["edits_ag"] = edits_ag
        repscreen.samples["replicate"] = np.repeat(reps, len(samples))
        repscreen.samples["sort"] = np.tile(samples, len(reps))
        repscreen.uns["replicates"] = reps
        repscreen.uns["samples"] = samples
        return repscreen

    @classmethod
    def from_adata(cls, adata):
        X_bcmatch = adata.layers["X_bcmatch"] if "X_bcmatch" in adata.layers else None
        edits = adata.layers["edits"] if "edits" in adata.layers else None
        return cls(
            (adata.X),
            X_edit=edits,
            X_bcmatch=X_bcmatch,
            guides=(adata.obs),
            samples=(adata.var),
            uns=(adata.uns),
            layers=(adata.layers),
        )

    def copy(self):
        adata = super().copy()
        return type(self).from_adata(adata)

    def rename(self, new_index: Sequence[str], axis: Literal[0, 1], keep_old=False):
        """Rename index"""
        if axis == 0:
            if len(new_index) != len(self.guides):
                raise ValueError(
                    f"New index of length {len(new_index)} doesn't match original index length of {len(self.guides)}."
                )
            index_subs_map = {
                self.guides.index[i]: new_index[i] for i in range(len(new_index))
            }
            if keep_old:
                self.guides["old_index"] = self.guides.index
            self.guides.index = new_index
            for k in self.uns.keys():
                if "guide" in self.uns[k].columns:
                    self.uns[k].guide = self.uns[k].guide.map(index_subs_map)
        else:
            if len(new_index) != len(self.samples):
                raise ValueError(
                    f"New index of length {len(new_index)} doesn't match original index length of {len(self.samples)}."
                )
            index_subs_map = {
                self.samples.index[i]: new_index[i] for i in range(len(new_index))
            }
            if keep_old:
                self.samples["old_index"] = self.samples.index
            self.samples.index = new_index
            for k in self.uns.keys():
                if "guide" in self.uns[k].columns and (
                    "edit" in self.uns[k].columns[1]
                    or "allele" in self.uns[k].columns[1]
                ):
                    self.uns[k].columns = (
                        self.uns[k].columns[:2].tolist()
                        + self.uns[k].columns[2:].map(index_subs_map).tolist()
                    )

    def __add__(self, other):
        if all(self.guides.index == other.guides.index) and all(
            self.samples.index == other.samples.index
        ):
            added_uns = {}
            for k in self.uns.keys():
                if k not in other.uns.keys():
                    continue
                if k == "allele_counts":
                    index_pair = ["guide", "allele"]
                elif k in ["guide_edit_counts", "edit_counts"]:
                    index_pair = ["guide", "edit"]
                elif k == "guide_reporter_allele_counts":
                    index_pair = ["guide", "reporter_allele", "guide_allele"]
                else:
                    continue
                self_df = self.uns[k].set_index(index_pair, drop=True)
                add_df = other.uns[k].set_index(index_pair, drop=True)
                add_df = add_df.loc[:, self_df.columns]
                added_uns[k] = (
                    self_df.add(add_df, fill_value=0).astype(int).reset_index()
                )

            return (
                ReporterScreen(
                    (self.X + other.X),
                    (self.layers["edits"] + other.layers["edits"]),
                    (self.layers["X_bcmatch"] + other.layers["X_bcmatch"]),
                    guides=(self.guides),
                    samples=(self.samples),
                    uns=added_uns,
                )
                if "X_bcmatch" in self.layers and "X_bcmatch" in other.layers
                else ReporterScreen(
                    (self.X + other.X),
                    (self.layers["edits"] + self.layers["edits"]),
                    guides=(self.guides),
                    samples=(self.samples),
                    uns=added_uns,
                )
            )
        raise ValueError("Guides/sample description mismatch")

    def __getitem__(self, index):
        """TODO: currently the samples names are in ['index'] column. Making it to be the index will
        allow the subsetting by condition names.
        """
        oidx, vidx = self._normalize_indices(index)
        if isinstance(oidx, (int, np.integer)):
            if not (-self.n_obs <= oidx < self.n_obs):
                raise IndexError(f"Observation index `{oidx}` is out of range.")
            oidx += self.n_obs * (oidx < 0)
            oidx = slice(oidx, oidx + 1, 1)
        if isinstance(vidx, (int, np.integer)):
            if not (-self.n_vars <= vidx < self.n_vars):
                raise IndexError(f"Variable index `{vidx}` is out of range.")
            vidx += self.n_vars * (vidx < 0)
            vidx = slice(vidx, vidx + 1, 1)
        guides_include = self.guides.iloc[oidx].index.tolist()
        condit_include = self.samples.iloc[vidx].index.tolist()
        adata = super().__getitem__(index)
        new_uns = deepcopy(self.uns)
        for k, df in adata.uns.items():
            if k.startswith("repguide_mask"):
                if "sample_covariates" in adata.uns:
                    adata.var["_rc"] = adata.var[
                        ["replicate"] + list(adata.uns["sample_covariates"])
                    ].values.tolist()
                    adata.var["_rc"] = adata.var["_rc"].map(
                        lambda slist: ".".join(slist)
                    )
                    new_uns[k] = df.loc[guides_include, adata.var._rc.unique()]
                    # adata.var.pop("_rc")
                else:
                    try:
                        new_uns[k] = df.loc[
                            guides_include, adata.var["replicate"].unique()
                        ]
                    except KeyError as e:
                        raise ValueError(
                            f"Replicate column should be `replicate`. Modify your ReporterScreen.samples. Current columns: {adata.var.columns.tolist()}"
                        ) from e
            if not isinstance(df, pd.DataFrame):
                if k == "sample_covariates":
                    new_uns[k] = df
                continue
            if "guide" in df.columns:
                if "allele" in df.columns:
                    key_col = ["guide", "allele"]
                elif "edit" in df.columns:
                    key_col = ["guide", "edit"]
                elif "aa_allele" in df.columns:
                    key_col = ["guide", "aa_allele"]
                else:
                    continue
                df_new = df.loc[df.guide.isin(guides_include), key_col + condit_include]
                # df_new = df_new.loc[df_new.loc[:, condit_include].sum(axis=1) > 0, :]
                new_uns[k] = df_new.reset_index(drop=True)
        adata.uns = new_uns
        return type(self).from_adata(adata)

    def remove_zero_allele_counts(
        self, key: Optional[Union[str, Iterable[str]]] = None
    ):
        condit_include = self.samples.index.tolist()

        def _remove_zero_count(key):
            df = self.uns[key]
            df = df.loc[df.loc[:, condit_include].sum(axis=1) > 0, :]
            self.uns[key] = df

        if key is None:
            for k in self.uns.keys():
                if "count" in k:
                    _remove_zero_count(k)
        elif isinstance(key, str):
            _remove_zero_count(key)
        elif isinstance(key, Iterable):
            for k in key:
                _remove_zero_count(k)

    def get_edit_mat_from_uns(
        self,
        ref_base: Optional[str] = None,
        alt_base: Optional[str] = None,
        match_target_position: Optional[bool] = None,
        rel_pos_start=0,
        rel_pos_end=np.Inf,
        rel_pos_is_reporter=False,
        target_pos_col="target_pos",
        edit_count_key="edit_counts",
    ):
        """
        Get the edit matrix from `.uns[edit_count_key]` to store the result in `.layers["edits"]`
        and returns the old `.layers["edits"]`.

        Args
        --
        ref_base: reference base of editing event to count. If not provided, default value in `.base_edited_from` is used.
        alt_base: alternate base of editing event to count. If not provided, default value in `.base_edited_to` is used.
        match_target_position: If `True`, edits with `.rel_pos` that matches `.guides[target_pos_col]` are counted.
        If `False`, edits with `.rel_pos` in range [rel_pos_start, rel_pos_end) is counted.
        rel_pos_start: Position from which edits will be counted when `match_target_position` is `False`.
        rel_pos_end: Position until where edits will be counted when `match_target_position` is `False`.
        rel_pos_is_reporter: `rel_pos_start` and `rel_pos_end` is relative to the reporter position. If `False`, those are treated to be relative of spacer position.
        target_pos_col: Column name in `.guides` DataFrame that has target position information when `match_target_position` is `True`.
        edit_count_key: Key of the edit counts DataFrame to be used to count the edits (`.uns[edit_count_key]`).
        """
        if ref_base is None:
            ref_base = self.base_edited_from
        if alt_base is None:
            alt_base = self.base_edited_to
        if match_target_position is None:
            match_target_position = not self.tiling
        if edit_count_key not in self.uns or len(self.uns[edit_count_key]) == 0:
            raise ValueError(
                "Edit count isn't calculated. "
                + "Call .get_edit_from_allele(allele_count_key, allele_key)"
            )
        edits = self.uns[edit_count_key].copy()
        if "edits" in self.layers:
            old_edits = self.layers["edits"].copy()
        else:
            old_edits = None
        self.layers["edits"] = np.zeros_like(
            self.X,
        ).astype(float)
        edits["ref_base"] = edits.edit.map(lambda e: e.ref_base)
        edits["alt_base"] = edits.edit.map(lambda e: e.alt_base)
        edits = edits.loc[
            (edits.ref_base == ref_base) & (edits.alt_base == alt_base), :
        ].reset_index()
        guide_len = self.guides.sequence.map(len)
        guide_name_to_idx = self.guides.reset_index().reset_index().set_index("name")
        edits["guide_idx"] = guide_name_to_idx.loc[edits.guide, "index"].reset_index(
            drop=True
        )
        edits["guide_start_pos"] = (
            32 - 6 - guide_len[edits.guide_idx].reset_index(drop=True)
        )
        if not match_target_position:
            edits["rel_pos"] = edits.edit.map(lambda e: e.rel_pos)
            edits["in_rel_pos_range"] = (
                edits.rel_pos
                >= rel_pos_start + (0 if rel_pos_is_reporter else edits.guide_start_pos)
            ) & (
                edits.rel_pos
                < rel_pos_end + (0 if rel_pos_is_reporter else edits.guide_start_pos)
            )
            good_edits = edits.loc[
                edits.in_rel_pos_range, ["guide", "edit"] + self.samples.index.tolist()
            ]
        else:
            edits["rel_pos"] = edits.edit.map(lambda e: e.rel_pos)
            edits["target_pos_guides"] = self.guides.iloc[edits.guide_idx, :][
                target_pos_col
            ].reset_index(drop=True)
            edits["target_pos_matches"] = edits.rel_pos == edits.target_pos_guides
            good_edits = edits.loc[
                edits.target_pos_matches,
                ["guide", "edit"] + self.samples.index.tolist(),
            ]

        good_guide_idx = guide_name_to_idx.loc[good_edits.guide, "index"].astype(int)
        for gidx, eidx in zip(good_guide_idx, good_edits.index):
            self.layers["edits"][gidx, :] += good_edits.loc[
                eidx, self.samples.index.tolist()
            ].astype(int)
        print("New edit matrix saved in .layers['edits']. Returning old edits.")
        return old_edits

    def get_guide_edit_rate(
        self,
        normalize_by_editable_base: Optional[bool] = None,
        edited_base: Optional[str] = None,
        editable_base_start=3,
        editable_base_end=8,
        bcmatch_thres=1,
        prior_weight: float = None,
        return_result=False,
        count_layer="X_bcmatch",
        edit_layer="edits",
        unsorted_condition_label=None,
    ):
        """
        prior_weight:
        Considering the edit rate to have prior of beta distribution with mean 0.5,
        prior weight to use when calculating posterior edit rate.
        unsorted_condition_label: Editing rate is calculated only for the samples that have this string in the sample index.
        """
        if edited_base is None:
            edited_base = self.base_edited_from
        if normalize_by_editable_base is None:
            normalize_by_editable_base = self.tiling
        if self.layers[count_layer] is None or self.layers[edit_layer] is None:
            raise ValueError("edits or barcode matched guide counts not available.")
        num_targetable_sites = 1.0
        if normalize_by_editable_base:
            if edited_base not in ["A", "C", "T", "G"]:
                raise ValueError("Specify the correct edited_base")
            num_targetable_sites = self.guides.sequence.map(
                lambda s: s[editable_base_start:editable_base_end].count(edited_base)
            )
        if unsorted_condition_label is not None:
            bulk_idx = np.where(
                self.samples.index.astype(str).map(
                    lambda s: unsorted_condition_label in s
                )
            )[0]
        else:
            bulk_idx = np.arange(0, len(self.samples)).astype(int)

        if prior_weight is None:
            prior_weight = 1
        n_edits = self.layers[edit_layer].copy()[:, bulk_idx].sum(axis=1)
        n_counts = self.layers[count_layer].copy()[:, bulk_idx].sum(axis=1)
        edit_rate = (n_edits + prior_weight / 2) / (
            (n_counts * num_targetable_sites) + prior_weight / 2
        )
        edit_rate[n_counts < bcmatch_thres] = np.nan
        if normalize_by_editable_base:
            print("normalize by editable counts")
            edit_rate[num_targetable_sites == 0] = np.nan
        if return_result:
            return edit_rate
        else:
            self.guides["edit_rate"] = edit_rate
            print(self.guides.edit_rate)

    def get_edit_rate(
        self,
        normalize_by_editable_base=None,
        edited_base=None,
        editable_base_start=3,
        editable_base_end=8,
        bcmatch_thres=1,
        prior_weight: float = None,
        return_result=False,
        count_layer="X_bcmatch",
        edit_layer="edits",
    ):
        """
        prior_weight:
        Considering the edit rate to have prior of beta distribution with mean 0.5,
        prior weight to use when calculating posterior edit rate.
        """
        if normalize_by_editable_base is None:
            normalize_by_editable_base = self.tiling
        if edited_base is None:
            edited_base = self.base_edited_from
        if self.layers[count_layer] is None or self.layers[edit_layer] is None:
            raise ValueError("edits or barcode matched guide counts not available.")
        num_targetable_sites = 1.0
        if normalize_by_editable_base:
            if edited_base not in ["A", "C", "T", "G"]:
                raise ValueError("Specify the correct edited_base")
            num_targetable_sites = self.guides.sequence.map(
                lambda s: s[editable_base_start:editable_base_end].count(edited_base)
            )
        if prior_weight is None:
            prior_weight = 1
        n_edits = self.layers[edit_layer]
        n_counts = self.layers[count_layer].astype(float)
        n_counts[n_counts < bcmatch_thres] = np.nan
        edit_rate = n_edits / n_counts
        if normalize_by_editable_base:
            edit_rate[num_targetable_sites == 0, :] = np.nan
        if return_result:
            return edit_rate
        else:
            self.layers["edit_rate"] = edit_rate

    def get_edit_from_allele(
        self, allele_count_key="allele_counts", allele_key="allele", return_result=False
    ):
        if allele_count_key not in self.uns or len(self.uns[allele_count_key]) == 0:
            raise ValueError(f"No allele information stored: {self.uns}")

        df = self.uns[allele_count_key].copy()
        df = df.loc[df[allele_key].map(str) != "", :]
        df["edits"] = df[allele_key].map(lambda a: str(a).split(","))
        df = (
            df[["guide", "edits"] + self.samples.index.tolist()]
            .explode("edits")
            .groupby(["guide", "edits"])
            .sum()
        )
        df = df.reset_index().rename(columns={"edits": "edit"})

        def try_objectify(s):
            try:
                return Edit.from_str(s)
            except ValueError:
                return AminoAcidEdit.from_str(s)

        df["edit"] = df.edit.map(try_objectify)
        if return_result:
            return df
        self.uns["edit_counts"] = df

    def _get_allele_norm(self, allele_count_df=None, thres=10):
        """
        Get bcmatched count to normalize allele counts
        thres: normalizing count below this number would be disregarded.
        """
        if allele_count_df is None:
            allele_count_df = self.uns["allele_counts"]
        guide_to_idx = self.guides.reset_index().reset_index().set_index("name")
        guide_idx = guide_to_idx.loc[allele_count_df.guide, "index"]
        norm_counts = self.layers["X_bcmatch"][guide_idx.values.astype(int), :].copy()
        norm_counts[norm_counts < thres] = np.nan
        return norm_counts

    def get_normalized_allele_counts(self, allele_count_df=None, norm_thres=10):
        if allele_count_df is None:
            allele_count_df = self.uns["allele_counts"]
        alleles_norm = allele_count_df.copy()
        count_columns = self.samples.index.tolist()
        norm_counts = self._get_allele_norm(
            allele_count_df=allele_count_df, thres=norm_thres
        )
        alleles_norm.loc[:, count_columns] = (
            alleles_norm.loc[:, count_columns] / norm_counts
        )
        return alleles_norm

    def filter_allele_counts_by_pos(
        self,
        rel_pos_start: int = 0,
        rel_pos_end: int = 32,
        rel_pos_is_reporter: bool = True,
        allele_uns_key: str = "allele_counts",
        filter_rel_pos: bool = True,
        map_to_filtered: bool = True,
        distribute: bool = False,
        jaccard_threshold: float = 0.1,
    ):
        """
        Filter alleles based on barcode matched counts, allele counts,
        or proportion

        Args
        --
        map_to_filtered: Map allele to the closest filtered allele to preserve total allele count. Ignores the case where there is no alleles filtered.
        rel_pos_start: rel_pos to start including (inclusive)
        rel_pos_end: rel_pos to end including (exclusive)
        rel_pos_is_reporter: rel_pos_start and rel_pos_end is 0-based relative to reporter sequence start.
        jaccard_threshold (float) --
        """
        allele_count_df = self.uns[allele_uns_key].copy()
        if rel_pos_is_reporter:
            filtered_allele, filtered_edits = zip(
                *allele_count_df.allele.map(
                    lambda a: filter_allele_by_pos(
                        a, rel_pos_start, rel_pos_end, filter_rel_pos
                    )
                )
            )
        else:
            if "guide_len" not in self.guides.columns.tolist():
                self.guides["guide_len"] = self.guides.sequence.map(len)
            guide_start_pos = (
                32 - 6 - self.guides.loc[allele_count_df.guide, "guide_len"].values
            )
            filtered_allele, filtered_edits = zip(
                *[
                    filter_allele_by_pos(
                        allele_count_df.allele[i],
                        guide_start_pos[i] + rel_pos_start,
                        guide_start_pos[i] + rel_pos_end,
                        filter_rel_pos,
                    )
                    for i in range(len(allele_count_df))
                ]
            )
        allele_count_df.loc[:, "allele"] = filtered_allele
        # Hashing on Allele object messes up the order. Converting it to str and back to allele for groupby.
        allele_count_df["str_allele"] = allele_count_df.allele.map(str)
        allele_count_df = (
            allele_count_df.groupby(["guide", "str_allele"])
            .sum(numeric_only=True)
            .reset_index()
        )
        allele_count_df.insert(
            1, "allele", allele_count_df.str_allele.map(lambda s: Allele.from_str(s))
        )
        allele_count_df = allele_count_df.drop("str_allele", axis=1)
        print(
            f"{sum(filtered_edits)} edits filtered from {len(filtered_edits)} alleles."
        )
        if map_to_filtered:
            print("mapping filtered alleles ...")
            if distribute:
                allele_count_df = _distribute_alleles_to_filtered(
                    self.uns[allele_uns_key],
                    allele_count_df,
                    jaccard_threshold=jaccard_threshold,
                )
            else:
                allele_count_df = _map_alleles_to_filtered(
                    self.uns[allele_uns_key],
                    allele_count_df,
                    jaccard_threshold=jaccard_threshold,
                )
        return allele_count_df

    def filter_allele_counts_by_base(
        self,
        ref_base: Union[List, str] = "A",
        alt_base: Union[List, str] = "G",
        allele_uns_key="allele_counts",
        map_to_filtered=True,
        jaccard_threshold: float = 0.5,
    ):
        """
        Filter alleles based on base change.

        Keyword arguments:
        map_to_filtered -- Map allele to the closest filtered allele to preserve total allele count. Ignores the case where there is no alleles filtered.
        """
        allele_count_df = self.uns[allele_uns_key].copy()
        filtered_allele, filtered_edits = zip(
            *allele_count_df.allele.map(
                lambda a: filter_allele_by_base(
                    a, allowed_ref_base=ref_base, allowed_alt_base=alt_base
                )
            )
        )
        allele_count_df.loc[:, "allele"] = filtered_allele
        # Hashing on Allele object messes up the order. Converting it to str and back to allele for groupby.
        allele_count_df["str_allele"] = allele_count_df.allele.map(str)
        allele_count_df = (
            allele_count_df.groupby(["guide", "str_allele"])
            .sum(numeric_only=True)
            .reset_index()
        )
        allele_count_df.insert(
            1, "allele", allele_count_df.str_allele.map(lambda s: Allele.from_str(s))
        )
        allele_count_df = allele_count_df.drop("str_allele", axis=1)
        print(
            f"{sum(filtered_edits)} edits filtered from {len(filtered_edits)} alleles."
        )
        if map_to_filtered:
            allele_count_df = _map_alleles_to_filtered(
                self.uns[allele_uns_key],
                allele_count_df,
                jaccard_threshold=jaccard_threshold,
            )
        return allele_count_df

    def collapse_allele_by_target(
        self, ref_base, alt_base, target_pos_column="target_pos"
    ):
        if target_pos_column not in self.guides.columns:
            raise ValueError(
                f"The .guides have to have '{target_pos_column}' specifying the relative position of target edit."
            )
        df = self.uns["allele_counts"].copy().reset_index(drop=True)
        df["target_pos"] = self.guides.loc[df.guide, target_pos_column].reset_index(
            drop=True
        )
        df["has_target"] = df.apply(
            lambda row: row.allele.has_edit(ref_base, alt_base, rel_pos=row.target_pos),
            axis=1,
        )
        df["has_nontarget"] = df.apply(
            lambda row: row.allele.has_other_edit(
                ref_base, alt_base, rel_pos=row.target_pos
            ),
            axis=1,
        )
        df = df.drop("target_pos", axis=1)
        return df.groupby(["guide", "has_target", "has_nontarget"]).sum()

    def collapse_allele_by_nedit(self, ref_base, alt_base):
        df = self.uns["allele_counts"].copy().reset_index(drop=True)
        df["n_edits"] = df.allele.map(
            lambda a: sum(
                (
                    e.ref_base == ref_base
                    and e.alt_base == alt_base
                    and e.rel_pos >= 6
                    and e.rel_pos < 6 + 20
                )
                for e in a.edits
            )
        )
        df.loc[df.n_edits > 4, "n_edits"] = 4
        return df.groupby(["guide", "n_edits"]).sum()

    def translate_alleles(self):
        if self.uns["allele_counts"] is None:
            print("No allele information. Run crisrpep-count with -a option.")
            return
        self.uns["allele_counts"]["aa_allele"] = self.uns["allele_counts"].allele.map(
            lambda s: get_aa_alleles(s)
        )

    def log_norms(self, **kwargs):
        """Calculate log normalized guide counts and log normalized edit counts"""
        super().log_norm(**kwargs)
        super().log_norm(output_layer="lognorm_edits", read_count_layer="edits")

    def log_fold_changes(self, cond1, cond2, return_result=False):
        """Calculate log fold changes in normalized guide and edit counts."""
        guide_fc = super().log_fold_change(cond1, cond2, return_result=return_result)
        edit_fc = super().log_fold_change(
            cond1,
            cond2,
            lognorm_counts_key="lognorm_edits",
            out_guides_suffix="edit_lfc",
            return_result=return_result,
        )
        if return_result:
            return (guide_fc, edit_fc)

    def log_fold_change_aggregates(
        self,
        cond1,
        cond2,
        aggregate_condit="replicate",
        compare_condit="sort",
        aggregate_fn="median",
        return_result=False,
        keep_per_replicate=False,
    ):
        guide_fc_agg = super().log_fold_change_aggregate(
            cond1,
            cond2,
            lognorm_counts_key="lognorm_counts",
            aggregate_condit=aggregate_condit,
            compare_condit=compare_condit,
            out_guides_suffix="lfc",
            aggregate_fn=aggregate_fn,
            return_result=return_result,
            keep_per_replicate=keep_per_replicate,
        )
        edit_fc_agg = super().log_fold_change_aggregate(
            cond1,
            cond2,
            lognorm_counts_key="lognorm_edits",
            aggregate_condit=aggregate_condit,
            compare_condit=compare_condit,
            out_guides_suffix="edit_lfc",
            aggregate_fn=aggregate_fn,
            return_result=return_result,
            keep_per_replicate=keep_per_replicate,
        )
        if return_result:
            return (guide_fc_agg, edit_fc_agg)

    def allele_log_fold_change_aggregate(
        self,
        cond1,
        cond2,
        aggregate_condit="replicate",
        compare_condit="sort",
        aggregate_fn="median",
        return_result=False,
        keep_per_replicate=False,
    ):
        pass

    def write(self, out_path):
        """
        Write .h5ad
        """
        adata = self.copy()
        for k, v in adata.uns.items():
            if isinstance(v, pd.DataFrame) and len(v) > 0:
                if "edit" in v.columns and isinstance(
                    adata.uns[k]["edit"].iloc[0], Edit
                ):
                    adata.uns[k].edit = adata.uns[k].edit.map(str)
                try:
                    for c in [
                        colname for colname in v.columns if "allele" in str(colname)
                    ]:
                        if isinstance(v[c].iloc[0], (Allele, CodingNoncodingAllele)):
                            adata.uns[k].loc[:, c] = adata.uns[k][c].map(str)
                except TypeError as e:
                    raise TypeError(f"error with {e}: {k, v} cannot be written")
        super(ReporterScreen, adata).write(out_path)


def _convert_obj_column_to_str(df, obj_column):
    if obj_column not in df.columns:
        return df
    df = df.rename(columns={obj_column: "obj_col"})
    df[obj_column] = df["obj_col"].map(str)
    df = df.drop("obj_col", axis=1)
    return df


def concat(screens: Collection[ReporterScreen], *args, axis=1, **kwargs):
    """Concatenate multiple Screen objects."""
    # TODO: var/obs info not concated if doesn't overlap
    # TODO: concat 2 times: allele_counts merging doesn't work?
    if axis == 1 and not all(
        screen.guides.index.equals(screens[0].guides.index) for screen in screens
    ):
        raise ValueError("Guide index doesn't match.")

    adata = ad.concat(screens, *args, axis=axis, **kwargs)
    if axis == 1:
        adata.obs = screens[0].guides
    elif axis == 0:
        adata.var = screens[0].samples
    else:
        raise ValueError(f"axis must be 0 or 1. Provided of invalid axis {axis}.")
    keys = set(screens[0].uns.keys())
    for screen in screens:
        keys.intersection(screen.uns.keys())

    if axis == 0:
        for k in keys:
            if k in ["target_base_change", "tiling", "sample_covariates"]:
                adata.uns[k] = screens[0].uns[k]
                continue
            elif "edit" not in k and "allele" not in k:
                continue
            adata.uns[k] = pd.concat([screen.uns[k] for screen in screens])

    if axis == 1:
        # If combining multiple samples, edit/allele tables should be merged.
        for k in keys:
            if k in ["target_base_change", "tiling", "sample_covariates"]:
                adata.uns[k] = screens[0].uns[k]
                continue
            elif "edit" not in k and "allele" not in k:
                continue
            screen_uns_dfs = [screen.uns[k] for screen in screens]

            if max(len(screen_uns_df) for screen_uns_df in screen_uns_dfs) == 0:
                continue
            if "guide_reporter_allele" in k:
                merge_on = screen[0].uns[k].columns[:3].tolist()
            else:
                merge_on = screen[0].uns[k].columns[:2].tolist()
            try:
                edit_cls = type(screens[0].uns[k][merge_on[1]][0])
            except IndexError:
                print(merge_on)
                print(screen[0].uns[k])
                exit(1)
            screen_uns_dfs = [
                _convert_obj_column_to_str(df, merge_on[1]) for df in screen_uns_dfs
            ]

            for i, screen in enumerate(screens):
                if screen.uns[k][merge_on].duplicated().any():
                    raise ValueError(
                        f"{i}-th ReporterScreen .uns['{k}'] has duplicated row. Aborting.: {screen.uns[k][merge_on].loc[screen.uns[k][merge_on].duplicated()]}"
                    )
                if i == 0:
                    adata.uns[k] = screen.uns[k]
                else:
                    adata.uns[k] = pd.merge(
                        adata.uns[k], screen.uns[k], on=merge_on, how="outer"
                    )
            adata.uns[k] = adata.uns[k].fillna(0)

            df = adata.uns[k].rename(columns={merge_on[1]: "_str"})
            df.insert(1, merge_on[1], df["_str"].map(lambda s: edit_cls.from_str(s)))
            adata.uns[k] = df.drop("_str", axis=1)

            float_col = adata.uns[k].select_dtypes(include=["float64"])
            for col in float_col.columns.values:
                adata.uns[k][col] = adata.uns[k][col].astype("int64")

    return ReporterScreen.from_adata(adata)


def read_h5ad(filename):
    adata = ad.read_h5ad(filename)
    return ReporterScreen.from_adata(adata)
