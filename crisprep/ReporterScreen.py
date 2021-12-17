from perturb_tools import Screen
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Union
from anndata import AnnData
import anndata as ad

def _get_counts(filename_pattern, guides, reps, conditions):
    for rep in reps:
        for cond in conditions:
            res = pd.read_csv(filename_pattern.format(r = rep, c = cond), delimiter="\t")
            res = res.set_index("guide_name")[["read_counts"]]
            res.columns = res.columns.map(lambda x: "{r}_{c}".format(r = rep, c = cond))
            guides = guides.join(res, how = "left")
    guides = guides.fillna(0)
    return(guides)

def _get_edits(filename_pattern, guide_info, reps, conditions, count_exact = True):
    edits = pd.DataFrame(index = guide_info.index)
    for rep in reps:
        for cond in conditions:
            res = pd.read_csv(filename_pattern.format(r = rep, c = cond), delimiter="\t")
            res = res.set_index("name")
            if count_exact:
                res_info = guide_info.join(res, how = "left").reset_index()
                this_edit = res_info.loc[(res_info["Target base position in gRNA"] == (res_info.pos + 1)) & 
                                            (res_info.ref_base == "A") & 
                                            (res_info.sample_base == "G")][["name", "count"]].set_index("name", drop = True)["count"]
            else:
                this_edit = res.loc[(res.ref_base == "A") & (res.sample_base == "G"),:].groupby('name')["count"].sum()
            this_edit.name = "{r}_{c}".format(r = rep, c = cond)
            edits = edits.join(this_edit, how = "left")
    edits = edits.fillna(0)
    return(edits)


class ReporterScreen(Screen):
    def __init__(self, X = None, X_edit = None, X_bcmatch = None, *args, **kwargs):
        super().__init__(X, *args, **kwargs)

        self.layers["edits"] = X_edit
        self.layers["X_bcmatch"] = X_bcmatch

    def get_edit_rate(self):
        if (not X_bcmatch is None) and (not X_edit is None):
            self.layers["edit_rate"] = self.layers["edits"] / self.layers["X_bcmatch"]
        else:
            raise ValueError("edits or barcode matched guide counts not available.")


    @classmethod
    def from_file_paths(cls, 
         reps: List[str] = None, 
         conditions: List[str] = None, 
         guide_info_file_name: str = None,
         guide_count_filenames: str = None,
         guide_bcmatched_count_filenames: str = None,
         edit_count_filenames: str = None
        ):
            
        guide_info = pd.read_csv(guide_info_file_name).set_index('name')
        guides = pd.DataFrame(index = pd.read_csv(guide_info_file_name)['name'])
        
        guides_lenient = _get_counts(guide_count_filenames, guides, reps, conditions)

        edits_ag = _get_edits(edit_count_filenames, guide_info, reps, conditions, count_exact = False)
        edits_exact = _get_edits(edit_count_filenames, guide_info, reps, conditions)
        

        if not guide_bcmatched_count_filenames is None:
            guides_bcmatch = _get_counts(guide_bcmatched_count_filenames, guides, reps, conditions)
            repscreen = cls(guides_lenient, edits_exact, X_bcmatch = guides_bcmatch, guides = guide_info)
        else:
            repscreen = cls(guides_lenient, edits_exact, guides = guide_info)
            
        repscreen.layers["edits_ag"] = edits_ag
        repscreen.condit["replicate"] = np.repeat(reps, len(conditions))
        repscreen.condit["sort"] = np.tile(conditions, len(reps))
        repscreen.uns["replicates"] = reps
        repscreen.uns["conditions"] = conditions

        return(repscreen)

    @classmethod
    def from_adata(cls, adata):
        if "X_bcmatch" in adata.layers:
            X_bcmatch = adata.layers["X_bcmatch"]
        else:
            X_bcmatch = None
        repscreen = cls(adata.X, adata.layers["edits"], X_bcmatch = X_bcmatch, guides = adata.obs, condit = adata.var, uns = adata.uns)
        return(repscreen)


    def __add__(self, other):
        if all(self.guides.index == other.guides.index) and all(self.condit.index == other.condit.index):
            if "X_bcmatch" in self.layers and "X_bcmatch" in other.layers:
                added = ReporterScreen(self.X + other.X, self.layers["edits"] + self.layers["edits"],
                                        self.layers["X_bcmatch"] + other.layers["X_bcmatch"], 
                                        guides = self.guides, condit = self.condit)
            else:
                added = ReporterScreen(self.X + other.X, self.layers["edits"] + self.layers["edits"],
                                        guides = self.guides, condit = self.condit)
            return(added)
        else:
            raise ValueError("Guides/sample description mismatch")



    def log_norm(self):
        super().log_norm()
        super().log_norm(read_count_layer = "edits")


    def log_fold_changes(self, cond1, cond2, return_result = False):
        guide_fc = super().log_fold_change(cond1, cond2, return_result = return_result)
        edit_fc = super().log_fold_change(cond1, cond2, lognorm_counts_key = "lognorm_edits", out_guides_suffix = "edit_lfc", return_result = return_result)
        if return_result :
            return((guide_fc, edit_fc))

    def log_fold_change_aggregates(
            self, 
            cond1, 
            cond2, 
            aggregate_condit = "replicate", 
            compare_condit = "sort", 
            aggregate_fn = "median", 
            return_result = False, 
            keep_per_replicate = False):
        guide_fc_agg = super().log_fold_change_aggregate(
                cond1, 
                cond2, 
                lognorm_counts_key = "lognorm_counts", 
                aggregate_condit = aggregate_condit,
                compare_condit = compare_condit,
                out_guides_suffix = "lfc", 
                aggregate_fn = aggregate_fn,
                return_result = return_result,
                keep_per_replicate = keep_per_replicate
                )

        edit_fc_agg = super().log_fold_change_aggregate(
                cond1, 
                cond2, 
                lognorm_counts_key = "lognorm_edits", 
                aggregate_condit = aggregate_condit,
                compare_condit = compare_condit,
                out_guides_suffix = "edit_lfc", 
                aggregate_fn = aggregate_fn,
                return_result = return_result,
                keep_per_replicate = keep_per_replicate
                )

        if return_result:
            return((guide_fc_agg, edit_fc_agg))


def concat(screens, *args, **kwargs):
    adata = ad.concat(screens, *args, **kwargs)
    keys = set()
    for i in range(len(screens)):
        keys.add(adata.uns.keys())
    for k in keys:
        adata.uns[k] = pd.concat([screen.uns[k] for screen in screens], axis = 1)

    return(ReporterScreen.from_adata(adata))

