import os


def test_count():
    os.system(
        "bean-count --R1 data/test_R1.fastq --R2 data/test_R2.fastq -b A -f data/test_guide_info.csv -o test_res/ -r"
    )
