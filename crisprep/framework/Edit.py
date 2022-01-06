from typing import Iterable, Literal

class Edit:
    reverse_map = {"A":"T", "C":"G", "T":"A", "G":"C"}
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
        pos, strand, base_change = edit_str.split(":")
        pos = int(pos)
        start, end = base_change.split(">")
        assert strand in ["+", "-"]
        # if strand == "-":
        #     start = cls.reverse_map[start]
        #     end = cls.reverse_map[end]
        return(cls(pos, start, end, strand = strand))

    def __eq__(self, other):
        if (
            self.pos == other.pos
            and self.ref_base == other.ref_base
            and self.alt_base == other.alt_base
        ):
            return True
        return False

    def __hash__(self):
        return(hash((self.pos, self.ref_base, self.alt_base)))

    def __repr__(self):
        return(("{}:{}:{}>{}".format(
            int(self.pos), 
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
        edits = set()
        for edit_str in allele_str.split(","):
            pos, strand, base_change = edit_str.split(":")
            pos = int(pos)
            start, end = base_change.split(">")
            assert strand in ["+", "-"]
            if strand == "-":
                start = cls.reverse_map[start]
                end = cls.reverse_map[end]
            edits.update(cls(pos, start, end))
        return(cls(edits))


    def __eq__(self, other):
        if self.edits == other.edits:
            return True
        return False

    def __hash__(self):
        return(hash(frozenset(self.edits)))

    def add(self, edit: Edit):
        self.edits.add(edit)  # TBD: adding to set?

    def __repr__(self):
        list_edits = list(map(lambda s: str(s), self.edits))
        list_edits.sort()
        return(",".join(list_edits))

