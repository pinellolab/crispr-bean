from copy import deepcopy
from typing import List, Tuple
from .Edit import Edit, Allele
from ..annotate.translate_allele import CDS, RefBaseMismatchException


def filter_allele_by_pos(
    allele: Allele,
    pos_start: int = None, 
    pos_end: int = None,
    filter_rel_pos = True
    ):
    '''
    Filter alleles based on position and return the filtered allele and
    number of filtered edits.
    '''
    filtered_edits = 0
    allele_filtered = deepcopy(allele)
    if not (pos_start is None and pos_end is None):
        if filter_rel_pos:
            for edit in allele.edits:
                if not (edit.rel_pos >= pos_start and edit.rel_pos < pos_end): 
                    filtered_edits += 1
                    allele_filtered.edits.remove(edit)
        else:
            for edit in allele.edits:
                if not (edit.pos >= pos_start and edit.pos < pos_end): 
                    filtered_edits += 1
                    allele_filtered.edits.remove(edit)
    else:
        print("No threshold specified") # TODO: warn
    return(allele_filtered, filtered_edits)

def filter_allele_by_base(
    allele: Allele,
    allowed_base_changes: List[Tuple] = None,
    allowed_ref_base: str = None, allowed_alt_base: str = None
):
    '''
    Filter alleles based on position and return the filtered allele and
    number of filtered edits.
    '''
    filtered_edits = 0
    if not (allowed_ref_base is None and allowed_alt_base is None) + (allowed_base_changes is None) == 1:
        print("No filters specified or misspecified filters.")
    elif not allowed_base_changes is None:
        for edit in allele.edits.copy():
            if not (edit.ref_base, edit.alt_base) in allowed_base_changes: 
                filtered_edits += 1
                allele.edits.remove(edit)
    elif not allowed_ref_base is None:
        for edit in allele.edits.copy():
            if edit.ref_base != allowed_ref_base: 
                filtered_edits += 1
                allele.edits.remove(edit)
            elif not allowed_alt_base is None and edit.alt_base != allowed_alt_base:
                filtered_edits += 1
                allele.edits.remove(edit)
    else:
        for edit in allele.edits.copy():
            if edit.alt_base != allowed_alt_base: 
                filtered_edits += 1
                allele.edits.remove(edit)
    return(allele, filtered_edits)

def get_aa_alleles(allele_str, include_synonymous = True):
    ldlr_cds = CDS()
    try:
        ldlr_cds.edit_allele(allele_str)
        ldlr_cds.get_aa_change(include_synonymous)
    except RefBaseMismatchException as e:
        print(e)
        return("ref mismatch")
    return ldlr_cds.get_aa_change(True)