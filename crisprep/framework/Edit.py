from typing import Iterable, Literal

class Edit:
    reverse_map = {"A":"T", "C":"G", "T":"A", "G":"C", "-":"-"}
    strand_map = {"+":1, "-":-1}
    def __init__(self, rel_pos: int, ref_base: chr, alt_base: chr, offset: int = None, strand: Literal[1, -1] = 1):
        strand_to_symbol = {1:'+', -1:'-'}
        self.rel_pos = rel_pos
        self.ref_base = ref_base
        self.alt_base = alt_base
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

    def __eq__(self, other):
        if (
            self.pos == other.pos
            and self.ref_base == other.ref_base
            and self.alt_base == other.alt_base
        ):
            return True
        return False

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
        if self.edits == other.edits:
            return True
        return False

    def __lt__(self, other): # Implemented for pandas compatibility
        return self.edits < other.edits

    def __hash__(self):
        return(hash(frozenset(self.edits)))

    def add(self, edit: Edit):
        self.edits.add(edit)  # TBD: adding to set?

    def __repr__(self):
        list_edits = list(map(lambda s: str(s), self.edits))
        list_edits.sort()
        return(",".join(list_edits))
