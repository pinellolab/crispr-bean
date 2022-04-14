from typing import List
import numpy as np
from Bio import SeqIO
from crisprep.framework.Edit import Edit, Allele
from crisprep.annotate.AminoAcidEdit import AminoAcidEdit, AminoAcidAllele
import logging
import sys

logging.basicConfig(level=logging.INFO,
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info

BASE_SET = {"A", "C", "T", "G"}
reverse_map = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}

codon_map = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I',
         'ATA':'I', 'ATG':'M', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
         'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A',
         'GCA':'A', 'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
         'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'TGT':'C', 'TGC':'C',
         'TGA':'*', 'TGG':'W', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
         'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

# This function is adopted from https://github.com/gpp-rnd/be-validation-pipeline/blob/main/notebooks/01_BEV_allele_frequencies.ipynb
def _translate(seq, codon_map):
#     if not seq: # if remove_introns returned False -> possible splice site mutation
    if seq == 'UTR':
        return 'Possible UTR mutation'
    if seq == 'intron':
        return 'Possible intronic mutation'
    aa = ''
    i = 0
    frame = 1
    while i < len(seq):
        substring = ''
        while frame <= 3:
            if i<len(seq):
                if seq[i] == '-': #print('deletion')
                    i += 1
                else:
                    substring += seq[i]  
                    i += 1
                    frame+=1
            else: # if reached end of the sequence and frame still <=3, complete codon sequence based on last_codon
                substring += last_codon[frame-1]
                i += 1
                frame+=1
        if len(substring) == 3:
            frame = 1 # reset frame 
            if ('N' in substring):
                aa = aa + '-'
            else:
                aa = aa + codon_map[substring] # translate codon
        else:
            frame = 1
    return aa

def _parse_range(chr_range):
    chrom, gen_range = chr_range.split(":")
    start, end = gen_range.split('-')
    start = int(start)
    end = int(end)
    return(chrom, start, end)
    
def _parse_description(desc_str):
    sstr = desc_str.split(" ")
    for s in sstr:
        if s.startswith("range="):
            return(_parse_range(s[6:]))

def _get_seq_pos_from_fasta(fasta_file_name: str):
    exons = list(SeqIO.parse(fasta_file_name, "fasta"))
    translated_seq = []
    genomic_pos = []
    for exon in exons:
        _, exon_start, exon_end = _parse_description(exon.description)
        for i, nt in enumerate(str(exon.seq)):
            if nt.islower(): continue
            translated_seq.append(nt)
            genomic_pos.append(exon_start + i)
    return(translated_seq, genomic_pos)

def _translate_single_codon(nt_seq_string: str, aa_pos: int) -> str:
    codon = ''.join(nt_seq_string[aa_pos*3:(aa_pos*3+3)])
    try:
        aa = codon_map[codon]
        return(aa)
    except KeyError:
        if codon[-1] == "N" and codon[0] in BASE_SET and codon[1] in BASE_SET:
            aa_set = set()
            for N in BASE_SET:
                aa_set.add(codon_map[codon[:2] + N])
            if len(aa_set) == 1:
                return(next(iter(aa_set)))
            else:
                print("warning: no matching aa with codon {}".format(codon))
                return("_")
        else:
            #raise ValueError("Cannot translate codon due to ambiguity: {}".format(codon))
            print("Cannot translate codon due to ambiguity: {}".format(codon))
            return("_")
    

class CDS():  
    fasta_file_name = "/data/pinello/PROJECTS/2021_08_ANBE/data/LDLR/exons.fa"
    translated_seq, genomic_pos = _get_seq_pos_from_fasta(fasta_file_name)
    nt = translated_seq
    pos = np.array(genomic_pos)
    
    def __init__(self):
        self.edited_nt = type(self).nt.copy()
        self.edited_aa_pos = set()
    
    @classmethod
    def set_exon_fasta_name(cls, fasta_file_name: str):
        cls.fasta_file_name = fasta_file_name
        cls.translated_seq, cls.genomic_pos = _get_seq_pos_from_fasta(fasta_file_name)
        cls.nt = cls.translated_seq
        cls.pos = np.array(cls.genomic_pos)
    
    def translate(self):
        self.aa = _translate(self.edited_nt, codon_map)
    
    def _get_relative_nt_pos(self, absolute_pos):
        nt_relative_pos = np.where(type(self).pos == absolute_pos)[0]
        assert len(nt_relative_pos) <= 1, nt_relative_pos
        if not nt_relative_pos: return(-1)
        return(nt_relative_pos.astype(int).item())
    
    def _edit_pos_to_aa_pos(self, edit_pos):
        nt_relative_pos = self._get_relative_nt_pos(edit_pos)
        if nt_relative_pos != -1 : 
            self.edited_aa_pos.add(nt_relative_pos // 3)
        return(nt_relative_pos)
    
    def edit_single(self, edit_str):
        edit = Edit.from_str(edit_str)
        rel_pos= self._edit_pos_to_aa_pos(edit.pos)
        if rel_pos == -1 : return # do nothing
        if edit.strand == '-': 
            ref_base = reverse_map[edit.ref_base]
            alt_base = reverse_map[edit.alt_base]
        else:
            ref_base = edit.ref_base
            alt_base = edit.alt_base
        if type(self).nt[rel_pos] != ref_base:
            raise ValueError("Ref base mismatch: {},{},{}".format(type(self).nt[rel_pos], edit, rel_pos))
        self.edited_nt[rel_pos] = alt_base
        
    def edit_allele(self, allele_str):
        if isinstance(allele_str, Allele):
            edit_strs = allele_str.edits
        else: edit_strs = allele_str.split(",")
        for edit_str in edit_strs:
            self.edit_single(edit_str)
    
    def get_aa_change(self, include_synonymous = True) -> List[str]:
        aa_mutations = AminoAcidAllele()

        for edited_aa_pos in self.edited_aa_pos:
            ref_aa = _translate_single_codon(type(self).nt, edited_aa_pos)
            mt_aa = _translate_single_codon(self.edited_nt, edited_aa_pos)
            if mt_aa == "_": 
                return("translation error")
            if not include_synonymous and ref_aa == mt_aa:
                continue
            aa_mutations.add(AminoAcidEdit(
                edited_aa_pos + 1, ref_aa, mt_aa))
        return(aa_mutations)
   