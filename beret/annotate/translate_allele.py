from contextlib import suppress
import os
from typing import List, Iterable
from copy import deepcopy
import numpy as np
import pandas as pd
from Bio import SeqIO
import beret as be
from beret.framework.Edit import Edit, Allele
from beret.framework.AminoAcidEdit import AminoAcidEdit, AminoAcidAllele, CodingNoncodingAllele
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
reverse_map = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N", "-":"-"}

codon_map = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I',
         'ATA':'I', 'ATG':'M', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
         'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A',
         'GCA':'A', 'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
         'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'TGT':'C', 'TGC':'C',
         'TGA':'*', 'TGG':'W', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
         'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

class RefBaseMismatchException(Exception):
    pass

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
    if len(codon) != 3:
        print("reached the end of CDS, frameshift.")
        return("/")
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
            print("Cannot translate codon due to ambiguity: {}".format(codon))
            return("_")
    

class CDS():  
    def __init__(self, fasta_file_name=None, suppressMessage = True):
        if fasta_file_name is None: 
            if not suppressMessage:
                print("No fasta file provided as reference: using LDLR")
            fasta_file_name=os.path.dirname(be.__file__) + "/annotate/ldlr_exons.fa"
        try:
            type  (self).set_exon_fasta_name(fasta_file_name)
        except FileNotFoundError as e:
            print(os.getcwd())
            raise e
        self.edited_nt = type(self).nt.copy()
        self.edited_aa_pos = set()
        self.edits_noncoding = set()
        self.set_exon_fasta_name(fasta_file_name)
        
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
        if edit.strand == '-': 
            ref_base = reverse_map[edit.ref_base]
            alt_base = reverse_map[edit.alt_base]
        else:
            ref_base = edit.ref_base
            alt_base = edit.alt_base
        if rel_pos == -1 : # position outside CDS
            self.edits_noncoding.add(edit)
        elif type(self).nt[rel_pos] != ref_base:
            if ref_base != "-":
                raise RefBaseMismatchException("ref:{} at pos {}, got edit {}".format(type(self).nt[rel_pos], rel_pos, edit))
        else:
            self.edited_nt[rel_pos] = alt_base
        if alt_base == '-': #frameshift
            #self.edited_nt.pop(rel_pos)
            self.edited_aa_pos.update(list(range(rel_pos, len(CDS.nt) // 3)))
        
    def edit_allele(self, allele_str):
        if isinstance(allele_str, Allele):
            edit_strs = allele_str.edits
        else: edit_strs = allele_str.split(",")
        for edit_str in edit_strs:
            self.edit_single(edit_str)
        if '-' in self.edited_nt:
            self.edited_nt.remove('-')
    
    def get_aa_change(self, include_synonymous = True) -> List[str]:
        mutations = CodingNoncodingAllele()
        mutations.nt_allele.update(self.edits_noncoding)
        for edited_aa_pos in self.edited_aa_pos:
            ref_aa = _translate_single_codon(type(self).nt, edited_aa_pos)
            mt_aa = _translate_single_codon(self.edited_nt, edited_aa_pos)
            if mt_aa == "_": 
                return("translation error")
            if not include_synonymous and ref_aa == mt_aa:
                continue
            mutations.aa_allele.add(AminoAcidEdit(
                edited_aa_pos + 1, ref_aa, mt_aa))
        
        return(mutations)

def get_allele_aa_change(allele_str, include_synonymous = True, fasta_file=None):
    """
    Obtain amino acid changes
    """
    cds = CDS(fasta_file_name=fasta_file)
    cds.edit_allele(allele_str)
    return(cds.get_aa_change(include_synonymous))

def translate_allele(allele: Allele, include_synonymouns = True, allow_ref_mismatch = True):
    try:
        return(get_allele_aa_change(allele, include_synonymous=include_synonymouns))
    except RefBaseMismatchException as e:
        if allow_ref_mismatch:
            print(e)
            return("ref mismatch")
        else:
            raise e

def translate_allele_df(allele_df, include_synonymouns = True, allow_ref_mismatch = True):
    allele_df = allele_df.copy()
    translated_alleles = allele_df.allele.map(lambda a: be.translate_allele(a, include_synonymouns=include_synonymouns, allow_ref_mismatch=allow_ref_mismatch))
    allele_df.insert(2, "cn_allele", translated_alleles)
    allele_df['cn_allele'] = allele_df.cn_allele.map(str)
    allele_df= allele_df.loc[~allele_df.cn_allele.map(str).isin(["ref mismatch", "|", "", "translation error"]), :]
    allele_df = allele_df.drop('allele', axis=1).groupby(['guide', 'cn_allele']).sum()
    allele_df = allele_df.reset_index().rename(columns={"cn_allele":"aa_allele"})
    allele_df.aa_allele = allele_df.aa_allele.map(lambda s:be.CodingNoncodingAllele.from_str(s))
    return(allele_df)


def filter_nt_allele(cn_allele: CodingNoncodingAllele, pos_include: Iterable[int]):
    """
    For CodingNoncodingAllele object, retain all amino acid mutation while filtering the nt_allele based on edit position.
    """
    cn_allele = deepcopy(cn_allele)
    edit_list = []
    for e in cn_allele.nt_allele.edits:
        if e.pos in pos_include: 
            edit_list.append(e)
    cn_allele.nt_allele = be.Allele(edit_list)
    return cn_allele

def filter_nt_alleles(cn_allele_df: pd.DataFrame, pos_include: Iterable[int]):
    """
    For CodingNoncodingAllele object, retain all amino acid mutation while filtering the nt_allele based on edit position.
    Arguments
    -- cn_allele_df (pd.DataFrame): Allele dataframe that has 'guide', 'aa_allele' as columns.
    """
    splice_only = cn_allele_df.aa_allele.map(lambda a: be.an.translate_allele.filter_nt_allele(a, pos_include))
    alleles = cn_allele_df.copy()
    alleles = alleles.drop('aa_allele', axis=1, inplace = False)
    alleles.insert(1, "aa_allele", splice_only)
    alleles = alleles.groupby(['guide', 'aa_allele']).sum().reset_index()
    alleles = alleles.loc[alleles.aa_allele.map(bool),:]
    return(alleles)