import pytest
import subprocess


@pytest.mark.order(11)
def test_filter_varscreen():
    cmd = "bean-filter tests/data/var_mini_screen_masked.h5ad -o tests/data/var_mini_screen_annotated -s 0 -e 19 -w -b -t -ap 0.1 -sp 0.3"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(12)
def test_filter_tiling_screen():
    cmd = "bean-filter tests/data/tiling_mini_screen_masked.h5ad -o tests/data/tiling_mini_screen_annotated -s 0 -e 19 -w -b -t -ap 0.1 -sp 0.3"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(13)
def test_filter_tiling_screen_translate_genename():
    cmd = "bean-filter tests/data/tiling_mini_screen_masked.h5ad -o tests/data/tiling_mini_screen_annotated_wrong -s 0 -e 19 -w -b -ap 0.1 -sp 0.3 --translate --translate-gene LDLR"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(14)
def test_filter_tiling_screen_translate_fasta():
    cmd = "bean-filter tests/data/tiling_mini_screen_masked.h5ad -o tests/data/tiling_mini_screen_annotated_wrong -s 0 -e 19 -w -b -ap 0.1 -sp 0.3 --translate --translate-fasta tests/data/ldlr_exons.fa"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(15)
def test_filter_tiling_screen_translate_genenames():
    cmd = "bean-filter tests/data/tiling_mini_screen_masked.h5ad -o tests/data/tiling_mini_screen_alleleFiltered -s 0 -e 19 -w -b -t -ap 0.1 -sp 0.3 --translate --translate-genes-list tests/data/gene_symbols.txt"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
