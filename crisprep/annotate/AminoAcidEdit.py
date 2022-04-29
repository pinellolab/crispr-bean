from typing import Iterable
from enum import IntEnum
from ..framework.Edit import Allele, Edit

AA_SET = {'A', 'C', 'D', 'E', 'F', 
'G', 'H', 'I', 'K', 'L', 
'M', 'N', 'P', 'Q', 'R',
 'S', 'T', 'V', 'W', 'Y', '*', '/'}

class MutationType(IntEnum):
    NO_CHANGE = -1
    SYNONYMOUS = 0
    MISSENSE = 1
    NONSENSE = 2

class AminoAcidEdit(Edit):
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
        return(hash(self.__repr__()))

    def _severity(self):
        if self.ref == self.alt: return(MutationType.SYNONYMOUS)
        if self.alt == "*": return(MutationType.NONSENSE)
        if self.alt in AA_SET: return(MutationType.MISSENSE)
        else:
            raise ValueError("Alt base invalid:{}".format(self.alt))

    def __eq__(self, other):
        if (
            self.pos == other.pos
            and self.ref == other.ref
            and self.alt == other.alt
        ):
            return True
        return False

    def __lt__(self, other): # Implemented for pandas compatibility
        return self.pos < other.pos

    def __gt__(self, other): # Implemented for pandas compatibility
        return self.pos > other.pos

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


class CodingNoncodingAllele():
    def __init__(self, aa_edits: Iterable[AminoAcidEdit] = None, base_edits: Iterable[Edit] = None):
        self.aa_allele = AminoAcidAllele(aa_edits)
        self.nt_allele = Allele(base_edits)

    def get_most_severe(self):
        aa_sev = self.aa_allele.get_most_severe()
        if len(self.nt_allele.edits) > 0: return max(aa_sev, 0.1)
        return aa_sev

    def __repr__(self):
        return str(self.aa_allele) + "|" + str(self.nt_allele)

    def __eq__(self, other):
        if not isinstance(other, type(self)): return False
        return(self.aa_allele == other.aa_allele and self.nt_allele == other.nt_allele)
    
    def __lt__(self, other):
        if self.aa_allele < other.aa_allele: return True
        if self.aa_allele == other.aa_allele:
            if self.nt_allele < other.nt_allele: return True
        return False

    def __gt__(self, other):
        if self.aa_allele > other.aa_allele: return True
        if self.aa_allele == other.aa_allele:
            if self.nt_allele > other.nt_allele: return True
        return False

    def __hash__(self):
        return(hash(self.__repr__()))