import os


def test_count():
    os.system(
        "bean-count --R1 test_R1.fastq --R2 test_R2.fastq -b A -f test_guide_info.csv -o test_res/ -a -as"
    )
