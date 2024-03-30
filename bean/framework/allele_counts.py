# functionalities for manipulating allele counts
from typing import Dict
from .ReporterScreen import ReporterScreen
from .Edit import Edit
import pandas as pd

def get_edit_from_allele(allele_count_df: pd.DataFrame, allele_col: str = 'allele'):
    """
    Get edit counts from allele counts pd.DataFrame
    -----
    Arguments
    allele_count_df -- DataFrame with guide, allele, and sample columns
    allele_col -- column name of allele_count_df that specifies allele.
    """
    allele_count_df = self.uns[allele_count_key].copy()
    allele_count_df["edits"] = allele_count_df[allele_col].map(lambda a: str(a).split(","))
    allele_count_df = allele_count_df.explode("edits").groupby(["guide", "edits"]).sum()
    allele_count_df = allele_count_df.reset_index().rename(columns={"edits":"edit"})
    def str_to_edit(s):
        try:
            return(Edit.from_str(s))
        except:
            return(AminoAcidEdit.from_str(s))
        # except Exception as e:
        #     print(f"Can't make edit {s} as Edit or AminoAcidEdit")
        #     print(e)
        #     return(s)
    allele_count_df["edit"] = allele_count_df.edit.map(str_to_edit(s))
    return(df)

def simulate_allele_counts(bdata: ReporterScreen, allele_counts_df: pd.DataFrame, allele_col: str = 'allele'):
    edit_counts = get_edit_from_allele(allele_counts_df, allele_col)
    nt_edit_rates = bdata._get_allele_norm(edit_counts_df)
    
    # For each guide, calculate allele proportion. 



# def get_nt_edit_rates(bdata, edit_counts_df: pd.DataFrame, norm_layer_key = "X_bcmatch") -> Dict[str, Edit]:
#     """
#     Get nucleotide level editing rates per guide
#     """
#     nt_edit_rates = bdata._get_allele_norm(edit_counts_df)
#     nt_edit_rates.
    