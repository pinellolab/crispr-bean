import pytest
import subprocess
from bean.annotate.translate_allele import CDS


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
