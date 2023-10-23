from __future__ import annotations
from typing import Iterable, Optional
from enum import IntEnum
import warnings
import numpy as np
import re
from ..framework.Edit import Allele, Edit
from ..utils.arithmetric import jaccard

AA_SET = {
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "*",
    "/",
}


class MutationType(IntEnum):
    NO_CHANGE = -1
    SYNONYMOUS = 0
    MISSENSE = 1
    NONSENSE = 2


class AminoAcidEdit(Edit):
    def __init__(self, pos: int, ref: str, alt: str, gene: str = None):
        self.gene = gene
        self.pos = pos
        assert ref in AA_SET, f"Invalid ref aa: {ref}"
        assert alt in AA_SET, f"Invalid alt aa: {alt}"
        self.ref = ref
        self.alt = alt

    @classmethod
    def from_str(cls, edit_str):
        if len(edit_str.split(":")) == 2:
            pos, aa_change = edit_str.split(":")
            gene = None
        else:
            gene, pos, aa_change = edit_str.split(":")
        ref, alt = aa_change.split(">")
        return cls(pos, ref, alt, gene=gene)

    @classmethod
    def match_str(cls, edit_str):
        if isinstance(edit_str, AminoAcidEdit):
            return True
        if edit_str == "":
            return True
        pattern = r"(\w+:)?-?\d+:[+-]:[A-Z*-]>[A-Z*-]"
        return re.fullmatch(pattern, edit_str)

    def get_abs_edit(self):
        return f"{f'{self.gene}:' if self.gene else ''}A{int(self.pos)}:{self.ref}>{self.alt}"

    def __repr__(self):
        return f"{f'{self.gene}:' if self.gene else ''}{int(self.pos)}:{self.ref}>{self.alt}"

    def __hash__(self):
        return hash(self.__repr__())

    def _severity(self):
        if self.ref == self.alt:
            return MutationType.SYNONYMOUS
        if self.alt == "*":
            return MutationType.NONSENSE
        if self.alt in AA_SET:
            return MutationType.MISSENSE
        else:
            raise ValueError(f"Alt base invalid:{self.alt}")

    def __eq__(self, other):
        return (
            self.gene == other.gene
            and self.pos == other.pos
            and self.ref == other.ref
            and self.alt == other.alt
        )

    def __lt__(self, other):  # Implemented for pandas compatibility
        if not self.gene and other.gene:
            return True
        if self.gene and not other.gene:
            return False
        if self.gene == other.gene:
            return self.pos < other.pos
        if self.gene < other.gene:
            return True
        if self.gene > other.gene:
            return False

    def __gt__(self, other):  # Implemented for pandas compatibility
        if self.gene and not other.gene:
            return True
        if not self.gene and other.gene:
            return False
        if self.gene == other.gene:
            return self.pos > other.pos
        if self.gene > other.gene:
            return True
        if self.gene < other.gene:
            return False


class AminoAcidAllele(Allele):
    def __init__(self, edits: Iterable[AminoAcidEdit] = None, gene=None):
        self.gene = None
        self.edits = set() if edits is None else set(edits)

    @classmethod
    def match_str(cls, allele_str):
        if isinstance(allele_str, AminoAcidAllele):
            return True
        # try:
        return all(map(AminoAcidEdit.match_str, allele_str.split(",")))

    @classmethod
    def from_str(cls, allele_str):  # pos:strand:start>end
        if type(allele_str) is AminoAcidAllele:
            return allele_str
        edits = set()
        try:
            for edit_str in allele_str.split(","):
                edit = AminoAcidEdit.from_str(edit_str)
                edits.add(edit)
        except ValueError:
            if allele_str.strip() == "":
                return cls(None)
        return cls(edits)

    def get_most_severe(self):
        if len(self.edits) == 0:
            return MutationType.NO_CHANGE
        severities = [edit._severity() for edit in self.edits]
        return max(severities)

    def get_most_severe_edit(self):
        if len(self.edits) == 0:
            return None
        severities = []
        edits = []
        for edit in self.edits:
            severities.append(edit._severity())
            edits.append(edit)
        return edits[np.argmax(severities)]


class CodingNoncodingAllele(Allele):
    def __init__(
        self,
        aa_edits: Iterable[AminoAcidEdit] = None,
        base_edits: Iterable[Edit] = None,
        unique_identifier=None,
    ):
        self.aa_allele = AminoAcidAllele(aa_edits)
        self.nt_allele = Allele(base_edits)
        self.uid = unique_identifier
        if self.uid is not None:
            self.set_uid(self.uid)

    @classmethod
    def from_str(cls, allele_str):
        if isinstance(allele_str, cls):
            return allele_str
        try:
            aa_allele_str, nt_allele_str = allele_str.split("|")
            aa_allele = AminoAcidAllele.from_str(aa_allele_str)
            nt_allele = Allele.from_str(nt_allele_str)
        except ValueError as e:
            print(e)
            if allele_str.strip() == "":
                return cls(None)
            print(allele_str)
            exit(1)
        return cls(aa_allele.edits, nt_allele.edits, nt_allele.get_uid())

    @classmethod
    def from_alleles(
        cls,
        aa_allele: Optional[AminoAcidAllele] = None,
        nt_allele: Optional[Allele] = None,
    ):
        if not aa_allele:
            aa_allele = AminoAcidAllele()
        if not nt_allele:
            nt_allele = Allele()
        return cls(aa_allele.edits, nt_allele.edits, nt_allele.get_uid())

    @classmethod
    def match_str(cls, allele_str):
        if isinstance(allele_str, CodingNoncodingAllele):
            return True
        if allele_str.count("|") != 1:
            return False
        aa_allele, nt_allele = allele_str.split("|")
        aa_match = AminoAcidAllele.match_str(aa_allele)
        nt_match = Allele.match_str(nt_allele)
        return aa_match and nt_match

    def get_most_severe(self):
        aa_sev = self.aa_allele.get_most_severe()
        return max(aa_sev, 0.1) if len(self.nt_allele.edits) > 0 else aa_sev

    def get_most_severe_edit(self):
        aa_sev = self.aa_allele.get_most_severe()
        sev_aa = self.aa_allele.get_most_severe_edit()
        if len(self.nt_allele.edits) > 0 and aa_sev > 0.1:
            return next(iter(self.nt_allele.edits))
        return sev_aa

    def set_uid(self, uid):
        self.uid = uid
        self.nt_allele.edits = {e.set_uid(uid) for e in self.nt_allele.edits}
        # self.aa_allele.edits = set([e.set_uid(uid) for e in self.aa_allele.edits])

    def get_jaccard(self, other):
        aa_jaccard = jaccard(self.aa_allele.edits, other.aa_allele.edits)
        nt_jaccard = jaccard(self.nt_allele.edits, other.nt_allele.edits)
        return (aa_jaccard, nt_jaccard)

    def get_jaccards(self, allele_list: Iterable[CodingNoncodingAllele]):
        aa_jaccards, nt_jaccards = zip(
            *list(map(lambda o: self.get_jaccard(o), allele_list))
        )
        aa_jaccards = np.array(aa_jaccards)
        nt_jaccards = np.array(nt_jaccards)
        return (aa_jaccards, nt_jaccards)

    def map_to_closest(
        self,
        allele_list,
        aa_jaccard_threshold=0.5,
        nt_jaccard_threshold=0.5,
        merge_priority: np.ndarray = None,
    ):
        """
        Arguments
        merge_priority -- Priority on which allele to merge if the jaccard index is the same.
        """
        if len(allele_list) == 0:
            return CodingNoncodingAllele()
        aa_jaccards, nt_jaccards = self.get_jaccards(allele_list)
        with warnings.catch_warnings():  # catch all NaN warning
            warnings.simplefilter("ignore")
            has_aa_jac = bool(self.aa_allele)
            has_nt_jac = bool(self.nt_allele)
            if has_aa_jac:
                aa_max_idx = np.where(aa_jaccards == np.nanmax(aa_jaccards))[0]
            else:
                aa_max_idx = np.where(
                    map(lambda cn_allele: not cn_allele.aa_allele, allele_list)
                )[0]
            if has_nt_jac:
                nt_max_idx = np.where(nt_jaccards == np.nanmax(nt_jaccards))[0]
            else:
                nt_max_idx = np.where(
                    map(lambda cn_allele: not cn_allele.nt_allele, allele_list)
                )[0]
            both_max_idx = np.intersect1d(aa_max_idx, nt_max_idx)

            def get_max_idx(max_idces):
                if len(max_idces) > 1 and merge_priority is not None:
                    if len(merge_priority) != len(allele_list):
                        raise ValueError(
                            "merge_priority length {} is not the same as allele_list length {}".format(
                                len(merge_priority), len(allele_list)
                            )
                        )
                    return max_idces[np.argmax(merge_priority[max_idces])]
                elif len(max_idces) > 0:
                    return max_idces[0]
                else:
                    return -1

            both_max_idx = get_max_idx(both_max_idx)
            if both_max_idx >= 0:
                if (
                    has_aa_jac and aa_jaccards[both_max_idx] >= aa_jaccard_threshold
                ) or (has_nt_jac and nt_jaccards[both_max_idx] >= nt_jaccard_threshold):
                    return allele_list[both_max_idx.item()]
            elif has_aa_jac:
                aa_max_idx = get_max_idx(aa_max_idx)
                if aa_max_idx >= 0:
                    return allele_list[aa_max_idx.item()]
            elif has_nt_jac:
                nt_max_idx = get_max_idx(nt_max_idx)
                if nt_max_idx >= 0:
                    return allele_list[nt_max_idx.item()]
            return CodingNoncodingAllele()

    def __bool__(self):
        return bool(self.aa_allele) or bool(self.nt_allele)

    def __len__(self):
        return len(self.aa_allele) + len(self.nt_allele)

    def __repr__(self):
        return f"{str(self.aa_allele)}|{str(self.nt_allele)}"

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        return (
            self.uid == other.uid
            and self.aa_allele == other.aa_allele
            and self.nt_allele == other.nt_allele
        )

    def __lt__(self, other):
        if self.aa_allele < other.aa_allele:
            return True
        if self.aa_allele == other.aa_allele and self.nt_allele < other.nt_allele:
            return True
        return False

    def __gt__(self, other):
        if self.aa_allele > other.aa_allele:
            return True
        if self.aa_allele == other.aa_allele and self.nt_allele > other.nt_allele:
            return True
        return False

    def __hash__(self):
        return hash(self.__repr__())
