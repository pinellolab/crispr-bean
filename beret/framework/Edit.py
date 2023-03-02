from __future__ import annotations
from enum import unique
from typing import Iterable
import numpy as np
import re
from ..utils.arithmetric import jaccard


class Edit:
    reverse_map = {"A": "T", "C": "G", "T": "A", "G": "C", "-": "-"}
    strand_map = {"+": 1, "-": -1}

    def __init__(
        self,
        rel_pos: int,
        ref_base: chr,
        alt_base: chr,
        offset: int = None,
        strand: int = 1,
        unique_identifier=None,
    ):
        assert strand in [+1, -1]
        strand_to_symbol = {1: "+", -1: "-"}
        self.rel_pos = rel_pos
        self.ref_base = ref_base  # TODO make it ref / alt instead of ref_base and alt_base for AAEdit comp. or make abstract class
        self.alt_base = alt_base
        self.uid = unique_identifier
        if type(strand) == int:
            self.strand = strand_to_symbol[strand]
        else:
            assert strand in ["+", "-"]
            self.strand = strand
        if not offset is None:
            self.pos = offset + self.rel_pos * strand
        else:
            self.pos = self.rel_pos

    @classmethod
    def from_str(cls, edit_str):  # pos:strand:start>end
        if type(edit_str) is Edit:
            return edit_str
        if not cls.match_str(edit_str):
            raise ValueError(f"{edit_str} doesn't match with Edit string format.")
        uid = None
        if "!" in edit_str:
            uid, edit_str = edit_str.split("!")

        pos, rel_pos, strand, base_change = edit_str.split(":")
        pos = int(pos)
        rel_pos = int(rel_pos)
        assert strand in ["+", "-"]
        strand = cls.strand_map[strand]
        offset = pos - rel_pos * strand
        ref_base, alt_base = base_change.split(">")
        return cls(
            rel_pos,
            ref_base,
            alt_base,
            offset=offset,
            strand=strand,
            unique_identifier=uid,
        )

    @classmethod
    def match_str(cls, edit_str):
        if isinstance(edit_str, Edit):
            return True
        pattern = r"-?\d+:-?\d+:[+-]:[A-Z*-]>[A-Z*-]"
        pattern2 = r"[\w*]!-?\d+:-?\d+:[+-]:[A-Z*-]>[A-Z*-]"
        return re.fullmatch(pattern, edit_str) or re.fullmatch(pattern2, edit_str)

    def get_abs_edit(self):
        """
        Returns absolute edit representation regardless of the relative edit position by guide.
        """
        if self.strand == "-":
            ref_base = type(self).reverse_map[self.ref_base]
            alt_base = type(self).reverse_map[self.alt_base]
        else:
            ref_base = self.ref_base
            alt_base = self.alt_base
        if not self.uid is None:
            return "{}!{}:{}>{}".format(self.uid, int(self.rel_pos), ref_base, alt_base)
        return "{}:{}>{}".format(int(self.pos), ref_base, alt_base)

    def set_uid(self, uid):
        if "!" in uid:
            raise ValueError("Cannot use special character `!` in uid.")
        self.uid = uid
        return self

    def get_abs_base_change(self):
        if self.strand == "-":
            ref_base = type(self).reverse_map[self.ref_base]
            alt_base = type(self).reverse_map[self.alt_base]
        else:
            ref_base = self.ref_base
            alt_base = self.alt_base
        return f"{ref_base}>{alt_base}"

    def get_base_change(self):
        ref_base = self.ref_base
        alt_base = self.alt_base
        return f"{ref_base}>{alt_base}"

    def __eq__(self, other):
        if self.__repr__() == other.__repr__():
            return True
        return False

    def __lt__(self, other):
        if isinstance(other, Edit):
            if self.pos != other.pos:
                return self.pos < other.pos
        return self.__repr__() < str(other)

    def __gt__(self, other):
        if isinstance(other, Edit):
            if self.pos != other.pos:
                return self.pos > other.pos
        return self.__repr__() > str(other)

    def __hash__(self):
        return hash(self.__repr__())

    def __repr__(self):
        if self.uid is None:
            return "{}:{}:{}:{}>{}".format(
                int(self.pos),
                int(self.rel_pos),
                self.strand,
                self.ref_base,
                self.alt_base,
            )

        return "{}!{}:{}:{}:{}>{}".format(
            self.uid,
            int(self.pos),
            int(self.rel_pos),
            self.strand,
            self.ref_base,
            self.alt_base,
        )


class Allele:
    # pos, ref, alt
    def __init__(self, edits: Iterable[Edit] = None):
        if edits is None:
            self.edits = set()
        else:
            self.edits = set(edits)

    @classmethod
    def from_str(cls, allele_str):  # pos:strand:start>end
        if type(allele_str) is Allele:
            return allele_str
        edits = set()
        try:
            for edit_str in allele_str.split(","):
                edit = Edit.from_str(edit_str)
                edits.add(edit)
        except ValueError:
            if allele_str.strip() == "":
                return cls(None)
        return cls(edits)

    @classmethod
    def match_str(cls, allele_str):
        if isinstance(allele_str, Allele):
            return True
        if allele_str == "":
            return True
        return all(map(Edit.match_str, allele_str.split(",")))

    def has_edit(self, ref_base, alt_base, pos=None, rel_pos=None):
        if not (pos is None) + (rel_pos is None):
            raise ValueError("Either pos or rel_pos should be specified")

        for e in self.edits:
            if e.ref_base == ref_base and e.alt_base == alt_base:
                if not pos is None:
                    if e.pos == pos:
                        return True
                elif e.rel_pos == rel_pos:
                    return True
        return False

    def has_other_edit(self, ref_base, alt_base, pos=None, rel_pos=None):
        """
        Has other edit than specified.
        """
        if len(self.edits) == 0:
            return False
        if not (pos is None) + (rel_pos is None):
            raise ValueError("Either pos or rel_pos should be specified")
        for e in self.edits:
            if e.ref_base == ref_base and e.alt_base == alt_base:
                if not pos is None:
                    if e.pos != pos:
                        return True
                elif e.rel_pos != rel_pos:
                    return True
            else:
                return True
        return False

    def get_jaccard(self, other):
        return jaccard(set(map(str, self.edits)), set(map(str, other.edits)))

    def get_jaccards(self, allele_list: Iterable[Allele]):
        return np.array(list(map(lambda o: self.get_jaccard(o), allele_list)))

    def map_to_closest(
        self, allele_list, jaccard_threshold=0.5, merge_priority: np.ndarray = None
    ):
        """
        Arguments
        merge_priority -- Priority on which allele to merge if the jaccard index is the same.
        """
        if len(allele_list) == 0:
            return Allele()
        nt_jaccards = np.array(list(map(lambda o: self.get_jaccard(o), allele_list)))
        if not np.isnan(np.nanmax(nt_jaccards)):
            nt_max_idx = np.where(nt_jaccards == np.nanmax(nt_jaccards))[0]
            if len(nt_max_idx) > 0:
                if len(nt_max_idx) > 1 and not merge_priority is None:
                    if not len(merge_priority) == len(allele_list):
                        raise ValueError(
                            "merge_priority length {} is not the same as allele_list length {}".format(
                                len(merge_priority), len(allele_list)
                            )
                        )
                    nt_max_idx = nt_max_idx[np.argmax(merge_priority[nt_max_idx])]
                else:
                    nt_max_idx = nt_max_idx[0]
                if nt_jaccards[nt_max_idx] > jaccard_threshold:
                    return allele_list[nt_max_idx.item()]
        return Allele()

    def __bool__(self):
        return len(self.edits) > 0

    def __len__(self):
        return len(self.edits)

    def __eq__(self, other):
        if self.__repr__() == other.__repr__():
            return True
        return False

    def __lt__(self, other):  # Implemented for pandas compatibility
        return len(self.edits) < len(other.edits)

    def __hash__(self):
        return hash(self.__repr__())

    def add(self, edit: Edit):
        self.edits.add(edit)  # TBD: adding to set?

    def update(self, edits: Iterable[Edit]):
        self.edits.update(edits)

    def __repr__(self):
        if len(self.edits) == 0:
            return ""
        list_edits = sorted(list(self.edits.copy()))
        list_edits = list(map(lambda s: str(s), list_edits))
        return ",".join(list_edits)
