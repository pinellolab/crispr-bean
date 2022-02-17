# uncompyle6 version 3.8.1.dev0
# Python bytecode 3.8.0 (3413)
# Decompiled from: Python 3.8.11 (default, Aug  3 2021, 15:09:35)
# [GCC 7.5.0]
# Embedded file name: /data/pinello/PROJECTS/2021_08_ANBE/software/crisprep/crisprep/ReporterScreen.py
# Compiled at: 2021-12-30 21:34:00
# Size of source mod 2**32: 8519 bytes
from functools import reduce
from perturb_tools import Screen
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Union, Collection
from anndata import AnnData
import anndata as ad
from .Edit import Edit
from ..annotate.translate_allele import CDS


def edit_alleles(allele_str):
    ldlr_cds = CDS()
    ldlr_cds.edit_allele(allele_str)

    return ldlr_cds.get_aa_change(True)


def _get_counts(filename_pattern, guides, reps, conditions):
    for rep in reps:
        for cond in conditions:
            res = pd.read_csv(filename_pattern.format(r=rep, c=cond), delimiter="\t")
            res = res.set_index("guide_name")[["read_counts"]]
            res.columns = res.columns.map(lambda x: "{r}_{c}".format(r=rep, c=cond))
            guides = guides.join(res, how="left")
        else:
            guides = guides.fillna(0)
            return guides


def _get_edits(filename_pattern, guide_info, reps, conditions, count_exact=True):
    edits = pd.DataFrame(index=(guide_info.index))
    for rep in reps:
        for cond in conditions:
            res = pd.read_csv(filename_pattern.format(r=rep, c=cond), delimiter="\t")
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
            this_edit.name = "{r}_{c}".format(r=rep, c=cond)
            edits = edits.join(this_edit, how="left")
        else:
            edits = edits.fillna(0)
            return edits


class ReporterScreen(Screen):
    def __init__(self, X=None, X_edit=None, X_bcmatch=None, *args, **kwargs):
        (super().__init__)(X, *args, **kwargs)
        self.layers["edits"] = X_edit
        self.layers["X_bcmatch"] = X_bcmatch

    def get_edit_rate(
        self, 
        normalize_by_editable_base = False, 
        edited_base = None,
        editable_base_start = 3,
        editable_base_end = 8,
        bcmatch_thres = 1,
        return_result = False
        ):
        if self.layers["X_bcmatch"] is not None and self.layers["edits"] is not None:
            num_targetable_sites = 1.0
            if normalize_by_editable_base:
                if not edited_base in ["A", "C", "T", "G"]: 
                    raise ValueError("Specify the correct edited_base")
                num_targetable_sites = self.guides.sequence.map(
                    lambda s: s[editable_base_start:editable_base_end].count(edited_base))
            bulk_idx = np.where(
                self.condit.reset_index()["index"].map(lambda s: "bulk" in s)
            )[0]
            self.layers["_edit_rate"] = (self.layers["edits"] + 0.5) / (
                self.layers["X_bcmatch"] + 0.5
            )
            self.layers["_edit_rate"][self.layers["X_bcmatch"] < bcmatch_thres] = np.nan
            self.guides["edit_rate"] = self.layers["_edit_rate"][:, bulk_idx].mean(
                axis=1
            ) / num_targetable_sites
            
        else:
            raise ValueError("edits or barcode matched guide counts not available.")

        
    @classmethod
    def from_file_paths(
        cls,
        reps: List[str] = None,
        conditions: List[str] = None,
        guide_info_file_name: str = None,
        guide_count_filenames: str = None,
        guide_bcmatched_count_filenames: str = None,
        edit_count_filenames: str = None,
    ):
        guide_info = pd.read_csv(guide_info_file_name).set_index("name")
        guides = pd.DataFrame(index=(pd.read_csv(guide_info_file_name)["name"]))
        guides_lenient = _get_counts(guide_count_filenames, guides, reps, conditions)
        edits_ag = _get_edits(
            edit_count_filenames, guide_info, reps, conditions, count_exact=False
        )
        edits_exact = _get_edits(edit_count_filenames, guide_info, reps, conditions)
        if guide_bcmatched_count_filenames is not None:
            guides_bcmatch = _get_counts(
                guide_bcmatched_count_filenames, guides, reps, conditions
            )
            repscreen = cls(
                guides_lenient, edits_exact, X_bcmatch=guides_bcmatch, guides=guide_info
            )
        else:
            repscreen = cls(guides_lenient, edits_exact, guides=guide_info)
        repscreen.layers["edits_ag"] = edits_ag
        repscreen.condit["replicate"] = np.repeat(reps, len(conditions))
        repscreen.condit["sort"] = np.tile(conditions, len(reps))
        repscreen.uns["replicates"] = reps
        repscreen.uns["conditions"] = conditions
        return repscreen

    @classmethod
    def from_adata(cls, adata):
        if "X_bcmatch" in adata.layers:
            X_bcmatch = adata.layers["X_bcmatch"]
        else:
            X_bcmatch = None
        repscreen = cls(
            (adata.X),
            (adata.layers["edits"]),
            X_bcmatch=X_bcmatch,
            guides=(adata.obs),
            condit=(adata.var),
            uns=(adata.uns),
        )
        return repscreen

    def __add__(self, other):
        if all(self.guides.index == other.guides.index):
            if all(self.condit.index == other.condit.index):
                added_uns = dict()
                for k in self.uns.keys():
                    if k not in other.uns.keys(): continue
                    if k == "guide_edit_counts":
                        index_pair = ["guide", "edit"]
                    elif k == "edit_counts":
                        index_pair = ["guide", "edit"]
                    elif k == "allele_counts":
                        index_pair = ["guide", "allele"]
                    elif k == "guide_reporter_allele_counts":
                        index_pair = ["guide", "reporter_allele", "guide_allele"]
                    else:
                        continue
                    self_df = self.uns[k].set_index(index_pair, drop = True)
                    add_df = other.uns[k].set_index(index_pair, drop = True)
                    add_df = add_df.loc[:,self_df.columns]
                    added_uns[k] = self_df.add(add_df, fill_value = 0).astype(int).reset_index()

                if "X_bcmatch" in self.layers and "X_bcmatch" in other.layers:
                
                    added = ReporterScreen(
                        (self.X + other.X),
                        (self.layers["edits"] + other.layers["edits"]),
                        (self.layers["X_bcmatch"] + other.layers["X_bcmatch"]),
                        guides=(self.guides),
                        condit=(self.condit),
                        uns = added_uns
                    )
                else:
                    added = ReporterScreen(
                        (self.X + other.X),
                        (self.layers["edits"] + self.layers["edits"]),
                        guides=(self.guides),
                        condit=(self.condit),
                        uns = added_uns
                    )
                return added
        raise ValueError("Guides/sample description mismatch")

    def get_edit_mat_from_uns(
        self, 
        ref_base, 
        alt_base, 
        match_target_position=True,
        rel_pos_start = 0,
        rel_pos_end = np.Inf
        ):
        edits = self.uns["edit_counts"]
        if not self.layers["edits"] is None:
            old_edits = self.layers["edits"].copy()
            self.layers["edits"] = np.zeros_like(
                self.X,
            )
        else:
            old_edits = None
        for i in edits.index:
            edit = Edit.from_str(edits.edit[i])
            if edit.rel_pos < rel_pos_start or edit.rel_pos >= rel_pos_end: continue
            guide_idx = np.where(edits.guide[i] == self.guides.reset_index().name)[
                0
            ].item()
            if (
                match_target_position
                and edit.rel_pos
                == self.guides.reset_index().loc[guide_idx, "target_pos"]
            ) or not match_target_position:
                if edit.ref_base == ref_base and edit.alt_base == alt_base:
                    self.layers["edits"][guide_idx, :] += edits.iloc[i, 2:].astype(int)
        return old_edits

    def allele_df_to_edits(self):
        pass

    def translate_allele(self):
        if self.uns["allele_counts"] is None:
            print("No allele information. Run crisrpep-count with -a option.")
            return
        self.uns["allele_counts"]["aa_allele"] = self.uns["allele_counts"].allele.map(
            lambda s: edit_alleles(s)
        )

    def log_norm(self):
        super().log_norm()
        super().log_norm(output_layer = "lognorm_edits", read_count_layer="edits")

    def log_fold_changes(self, cond1, cond2, return_result=False):
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

    def log_fold_change_aggregate(
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


def concat(screens: Collection[ReporterScreen], *args, **kwargs):
    adata = ad.concat(screens, *args, **kwargs)
    keys = set()

    # Merging uns
    for screen in screens:
        keys.update(screen.uns.keys())

    for k in keys:
        if k == "guide_edit_counts":
            merge_on = ["guide", "edit"]
        elif k == "edit_counts":
            merge_on = ["guide", "edit"]
        elif k == "allele_counts":
            merge_on = ["guide", "allele"]
        elif k == "guide_reporter_allele_counts":
            merge_on = ["guide", "reporter_allele", "guide_allele"]
        else:
            continue
        for i, screen in enumerate(screens):
            if i == 0:
                adata.uns[k] = screen.uns[k]
            else:
                try:
                    adata.uns[k] = pd.merge(
                        adata.uns[k], screen.uns[k], on=merge_on, how="outer"
                    )
                except:
                    print(k)
                    print(screen.uns[k])
        adata.uns[k] = adata.uns[k].fillna(0)
        float_col = adata.uns[k].select_dtypes(include=["float64"])
        for col in float_col.columns.values:
            adata.uns[k][col] = adata.uns[k][col].astype("int64")

    return ReporterScreen.from_adata(adata)


# okay decompiling ReporterScreen.cpython-38.pyc


def read_h5ad(filename):
    adata = ad.read_h5ad(filename)
    return ReporterScreen.from_adata(adata)
