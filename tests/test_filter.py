import pytest
import subprocess
from bean.annotate.translate_allele import CDS
from bean.annotate.translate_allele import (
    get_cds_seq_pos_from_gene_name,
    get_cds_seq_pos_from_fasta,
)


@pytest.mark.order(311)
def test_filter_varscreen():
    cmd = "bean filter tests/data/var_mini_screen_masked.h5ad -o tests/data/var_mini_screen_annotated -s 0 -e 19 -w -b -t -ap 0.1 -sp 0.3"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(312)
def test_filter_tiling_screen():
    cmd = "bean filter tests/data/tiling_mini_screen_masked.h5ad -o tests/data/tiling_mini_screen_annotated -s 0 -e 19 -w -b -t -ap 0.1 -sp 0.3"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(313)
def test_filter_tiling_screen_translate_genename():
    cmd = "bean filter tests/data/tiling_mini_screen_masked.h5ad -o tests/data/tiling_mini_screen_annotated_wrong -s 0 -e 19 -w -b -ap 0.1 -sp 0.3 --translate --translate-gene LDLR"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(314)
def test_filter_tiling_screen_translate_fasta():
    cmd = "bean filter tests/data/tiling_mini_screen_masked.h5ad -o tests/data/tiling_mini_screen_annotated_wrong -s 0 -e 19 -w -b -ap 0.1 -sp 0.3 --translate --translate-fasta tests/data/ldlr_exons.fa"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(315)
def test_filter_tiling_screen_translate_genenames():
    cmd = "bean filter tests/data/tiling_mini_screen_masked.h5ad -o tests/data/tiling_mini_screen_alleleFiltered -s 0 -e 19 -w -b -t -ap 0.1 -sp 0.3 --translate --translate-genes-list tests/data/gene_symbols.txt"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


def test_translate_aa():
    abca1 = CDS.from_gene_name("ABCA1")
    allele_str = "chr9:104903664:0:+:G>A"
    assert str(abca1.get_aa_change(allele_str)) == "ABCA1:6:Q>*|", allele_str

    abca1 = CDS.from_gene_name("ABCA1")
    allele_str = "chr9:104903664:0:-:C>T"
    assert str(abca1.get_aa_change(allele_str)) == "ABCA1:6:Q>*|", allele_str


def test_pos_fasta():
    chrom, translated_seq, genomic_pos, strand = get_cds_seq_pos_from_fasta(
        "tests/data/abca1.fa"
    )
    chrom2, translated_seq2, genomic_pos2, strand2 = get_cds_seq_pos_from_gene_name(
        "ABCA1"
    )
    assert chrom == chrom2
    assert translated_seq == translated_seq2
    assert genomic_pos == genomic_pos2, genomic_pos[:10]
    assert strand == strand2


def test_translate_fasta():
    abca1 = CDS.from_fasta("tests/data/abca1.fa", "ABCA1")
    allele_str = "chr9:104903664:0:+:G>A"
    assert str(abca1.get_aa_change(allele_str)) == "ABCA1:6:Q>*|", allele_str

    abca1 = CDS.from_fasta("tests/data/abca1.fa", "ABCA1")
    allele_str = "chr9:104903664:0:-:C>T"
    assert str(abca1.get_aa_change(allele_str)) == "ABCA1:6:Q>*|", allele_str
