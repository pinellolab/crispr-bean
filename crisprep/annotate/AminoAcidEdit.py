from typing import Iterable
from enum import IntEnum
from ..framework.Edit import Allele

AA_SET = {'A', 'C', 'D', 'E', 'F', 
'G', 'H', 'I', 'K', 'L', 
'M', 'N', 'P', 'Q', 'R',
 'S', 'T', 'V', 'W', 'Y', '*'}

class MutationType(IntEnum):
    NO_CHANGE = -1
    SYNONYMOUS = 0
    MISSENSE = 1
    NONSENSE = 2

class AminoAcidEdit:
    def __init__(self, pos: int, ref: str, alt: str):
        self.pos = pos
        assert ref in AA_SET, "Invalid ref aa: {}".format(ref)
        assert alt in AA_SET, "Invalid alt aa: {}".format(alt)
        self.ref = ref
        self.alt = alt

    @classmethod
    def from_str(cls, edit_str):
        pos, aa_change = edit_str.split(":")
        ref, alt = aa_change.split(">")
        return(cls(pos, ref, alt))

    def __repr__(self):
        return("{}:{}>{}".format(int(self.pos),
        self.ref, self.alt))

    def __hash__(self):
        return(hash((self.pos, self.ref, self.alt)))

    def _severity(self):
        if self.ref == self.alt: return(MutationType.SYNONYMOUS)
        if self.alt == "*": return(MutationType.NONSENSE)
        if self.alt in AA_SET: return(MutationType.MISSENSE)
        else:
            raise ValueError("Alt base invalid:{}".format(self.alt))

class AminoAcidAllele(Allele):
    def __init__(self, edits: Iterable[AminoAcidEdit] = None):
        if edits is None:
            self.edits = set()
        else:
            self.edits = set(edits)

    def get_most_severe(self):
        if len(self.edits) == 0:
            return(MutationType.NO_CHANGE)
        severities = []
        for edit in self.edits:
            severities.append(edit._severity())
        return(max(severities))


