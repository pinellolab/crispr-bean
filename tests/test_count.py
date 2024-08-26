import pytest
import subprocess
from bean.mapping import GuideEditCounter
from bean.mapping.utils import _get_input_parser
from bean.mapping._supporting_fn import revcomp


@pytest.mark.order(102)
def test_count():
    cmd = "bean count --R1 tests/data/test_R1.fastq --R2 tests/data/test_R2.fastq -b A -f tests/data/test_guide_info.csv -o tests/test_res/ -r --guide-start-seq=GGAAAGGACGAAACACCG"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(103)
def test_guide_count():
    cmd = "bean count --R1 tests/data/test_R1.fastq --R2 tests/data/test_R2.fastq -b A -f tests/data/test_guide_info.csv -o tests/test_res/ -g --guide-start-seq=GGAAAGGACGAAACACCG"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(104)
def test_count_samples_match():
    cmd = "bean count-samples -i tests/data/sample_list.csv -b A -f tests/data/test_guide_info.csv -o tests/test_res/var/ -r --guide-start-seq=GGAAAGGACGAAACACCG --match-target-pos --rerun"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(105)
def test_count_samples_nofilter():
    cmd = "bean count-samples -i tests/data/sample_list.csv -b A -f tests/data/test_guide_info.csv -o tests/test_res/var/ -r --guide-start-seq=GGAAAGGACGAAACACCG --skip-filtering"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(106)
def test_count_samples_dual():
    cmd = "bean count-samples -i tests/data/sample_list.csv -b A,C -f tests/data/test_guide_info.csv -o tests/test_res/var/ -r --guide-start-seq=GGAAAGGACGAAACACCG"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(107)
def test_count_samples_bcstart():
    cmd = "bean count-samples -i tests/data/sample_list.csv -b A -f tests/data/test_guide_info.csv -o tests/test_res/var2/ -r --barcode-start-seq=GGAA"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


def test_barcode_start_idx():
    args = vars(
        _get_input_parser().parse_args(
            [
                "-b",
                "A",
                "-f",
                "tests/data/test_guide_info.csv",
            ]
        )
    )
    args["barcode_start_seq"] = "gtttgaattcgctagctaggtcttg"
    args["R1"] = "tests/data/test_R1.fastq"
    args["R2"] = "tests/data/test_R2.fastq"
    gc = GuideEditCounter(**args)
    sidx, bc = gc.get_barcode(
        "",
        revcomp(
            "agtggcaccgagtcggtgcttttttTAACAGTGTAATCTGGCGAGCCACTGTTCTTTGTACCAGAAgtttgaattcgctagctaggtcttgctgg".upper()
        ),
    )
    assert sidx == 29
    assert bc == "AGAA"


@pytest.mark.order(108)
def test_count_samples_tiling():
    cmd = "bean count-samples -i tests/data/sample_list_tiling.csv -b A -f tests/data/test_guide_info_tiling_chrom.csv -o tests/test_res/tiling/ -r --tiling"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(109)
def test_count_chroms():
    cmd = "bean count --R1 tests/data/test_tiling_R1.fastq --R2 tests/data/test_tiling_R2.fastq -b A -f tests/data/test_guide_info_tiling_chrom.csv -o tests/test_res/tiling_chrom/ -r --tiling"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
