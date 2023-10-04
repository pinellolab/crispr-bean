import pytest
import subprocess


@pytest.mark.order(2)
def test_count():
    cmd = "bean-count --R1 tests/data/test_R1.fastq --R2 tests/data/test_R2.fastq -b A -f tests/data/test_guide_info.csv -o tests/test_res/ -r --guide-start-seq=GGAAAGGACGAAACACCG"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(3)
def test_guide_count():
    cmd = "bean-count --R1 tests/data/test_R1.fastq --R2 tests/data/test_R2.fastq -b A -f tests/data/test_guide_info.csv -o tests/test_res/ -g --guide-start-seq=GGAAAGGACGAAACACCG"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(4)
def test_count_samples():
    cmd = "bean-count-samples --input tests/data/sample_list.csv -b A -f tests/data/test_guide_info.csv -o tests/test_res/ -r --guide-start-seq=GGAAAGGACGAAACACCG"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(5)
def test_count_chroms():
    cmd = "bean-count --R1 tests/data/test_tiling_R1.fastq --R2 tests/data/test_tiling_R2.fastq -b A -f tests/data/test_guide_info_tiling_chrom.csv -o tests/test_res/ -r"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
