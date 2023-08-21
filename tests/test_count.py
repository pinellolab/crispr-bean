import subprocess


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


def test_guide_count():
    cmd = "bean-count --R1 data/test_R1.fastq --R2 data/test_R2.fastq -b A -f data/test_guide_info.csv -o test_res/ -g --guide-start-seq=GGAAAGGACGAAACACCG"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


def test_count_samples():
    cmd = "python ../../bin/bean-count-samples --input data/sample_list.csv -b A -f data/test_guide_info.csv -o test_res/ -r --guide-start-seq=GGAAAGGACGAAACACCG"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
