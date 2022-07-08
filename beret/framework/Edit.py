from enum import unique
from typing import Iterable

class Edit:
    reverse_map = {"A":"T", "C":"G", "T":"A", "G":"C", "-":"-"}
    strand_map = {"+":1, "-":-1}
    def __init__(self, rel_pos: int, ref_base: chr, alt_base: chr, offset: int = None, strand: int = 1, 
    unique_identifier = None):
        assert strand in [+1, -1]
        strand_to_symbol = {1:'+', -1:'-'}
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
    def from_str(cls, edit_str): #pos:strand:start>end
        if type(edit_str) is Edit: return(edit_str)
        pos, rel_pos, strand, base_change = edit_str.split(":")
        pos = int(pos)
        rel_pos = int(rel_pos)
        assert strand in ["+", "-"]
        strand = cls.strand_map[strand]
        offset = pos - rel_pos*strand
        ref_base, alt_base = base_change.split(">")
        
        # if strand == "-":
        #     start = cls.reverse_map[start]
        #     end = cls.reverse_map[end]
        return(cls(rel_pos, ref_base, alt_base, offset = offset, strand = strand))

    def get_abs_edit(self):
        '''
        Returns absolute edit representation regardless of the relative edit position by guide.
        '''
        if self.strand == '-':
            ref_base = type(self).reverse_map[self.ref_base]
            alt_base = type(self).reverse_map[self.alt_base]
        else:
            ref_base = self.ref_base
            alt_base = self.alt_base
        if not self.uid is None:
            return("{}!{}:{}>{}".format(self.uid, int(self.pos), ref_base, alt_base))
        return("{}:{}>{}".format(int(self.pos), ref_base, alt_base))

    def set_uid(self, uid):
        self.uid = uid
        return(self)

    def __eq__(self, other):
        if (
            self.pos == other.pos
            and self.ref_base == other.ref_base
            and self.alt_base == other.alt_base
        ):
            return True
        return False
    
    def __lt__(self, other): # Implemented for pandas compatibility
        return self.pos < other.pos

    def __gt__(self, other): # Implemented for pandas compatibility
        if isinstance(other, Edit):
            return self.pos > other.pos
        else:
            return self.__repr__() > str(other)

    def __hash__(self):
        # Note that this doesn't include relative bases. 
        # This wouldn't matter if we assign edit to each guide.
        return(hash((self.pos, self.ref_base, self.alt_base)))

    def __repr__(self):
        return(("{}:{}:{}:{}>{}".format(
            int(self.pos), 
            int(self.rel_pos),
            self.strand, 
            self.ref_base, 
            self.alt_base)))


class Allele:
    # pos, ref, alt
    def __init__(self, edits: Iterable[Edit] = None):
        if edits is None:
            self.edits = set()
        else:
            self.edits = set(edits)

    @classmethod
    def from_str(cls, allele_str): #pos:strand:start>end
        if type(allele_str) is Allele: return(allele_str)
        edits = set()
        try:
            for edit_str in allele_str.split(","):
                edit = Edit.from_str(edit_str)
                edits.add(edit)
        except ValueError:
            if allele_str.strip() == '':
                return(cls(None))
        return(cls(edits))

    def has_edit(self, ref_base, alt_base, pos = None, rel_pos = None):
        if not (pos is None) + (rel_pos is None): 
            raise ValueError("Either pos or rel_pos should be specified")
        
        for e in self.edits:
            if e.ref_base == ref_base and e.alt_base == alt_base:
                if not pos is None:
                    if e.pos == pos: return True
                elif e.rel_pos == rel_pos: return True
        return False

    def has_other_edit(self, ref_base, alt_base, pos = None, rel_pos = None):
        '''
        Has other edit than specified.
        '''
        if len(self.edits) == 0: return False
        if not (pos is None) + (rel_pos is None): 
            raise ValueError("Either pos or rel_pos should be specified")
        for e in self.edits:
            if e.ref_base == ref_base and e.alt_base == alt_base:
                if not pos is None:
                    if e.pos != pos: return True
                elif e.rel_pos != rel_pos: return True
            else: return True
        return False

    def __eq__(self, other):
        if self.__repr__() == other.__repr__():
            return True
        return False

    def __lt__(self, other): # Implemented for pandas compatibility
        return len(self.edits) < len(other.edits)

    def __hash__(self):
        return(hash(self.__repr__()))

    def add(self, edit: Edit):
        self.edits.add(edit)  # TBD: adding to set?

    def update(self, edits: Iterable[Edit]):
        self.edits.update(edits)

    def __repr__(self):
        if len(self.edits) == 0: return ""
        list_edits = sorted(list(self.edits.copy()))
        list_edits = list(map(lambda s: str(s), list_edits))
        return(",".join(list_edits))
