from perturb_tools import Screen
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import List
from anndata import AnnData

class ReporterScreen(Screen):
    def __init__(self, 
                 reps:List[str], 
                 conditions:List[str], 
                 guide_info_file_name:str,
                 guide_count_filenames:str = None,
                 guide_bcmatched_count_filenames: str = None,
                 edit_count_filenames:str = None
                ):
        guide_info = pd.read_csv(guide_info_file_name).set_index('name')
        guides = pd.DataFrame(index = pd.read_csv(guide_info_file_name)['name'])
        
        guides_lenient = self.get_counts(guide_count_filenames, guides, reps, conditions)
        guides_bcmatch = self.get_counts(guide_bcmatched_count_filenames, guides, reps, conditions)

        edits_ag = self.get_edits(edit_count_filenames, guide_info, reps, conditions, count_exact = False)
        edits_exact = self.get_edits(edit_count_filenames, guide_info, reps, conditions)
        
        super().__init__(guides_lenient, guide_info)
        self.condit["replicate"] = np.repeat(reps, len(conditions))
        self.condit["sort"] = np.tile(conditions, len(reps))
        self.condit["mapped_reads"] = self.X.sum(axis = 0)
        
        self.layers["guide_RPM"] = np.divide(self.X, self.condit.mapped_reads[np.newaxis, :])*10**6
        self.layers["X_bcmatch"] = guides_bcmatch
        self.layers["edits_ag"] = edits_ag
        self.layers["edits"] = edits_exact
        self.layers["edit_rate"] = self.layers["edits"] / self.layers["X_bcmatch"]

        self.uns["replicates"] = reps
        self.uns["conditions"] = conditions


    def get_counts(self, filename_pattern, guides, reps, conditions):
        for rep in reps:
            for cond in conditions:
                res = pd.read_csv(filename_pattern.format(r = rep, c = cond), delimiter="\t")
                res = res.set_index("guide_name")[["read_counts"]]
                res.columns = res.columns.map(lambda x: "{r}_{c}".format(r = rep, c = cond))
                guides = guides.join(res, how = "left")
        guides = guides.fillna(0)
        return(guides)

    def get_edits(self, filename_pattern, guide_info, reps, conditions, count_exact = True):
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
