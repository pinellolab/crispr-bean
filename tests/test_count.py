import os


def test_count():
    os.system(
        "bean-count --R1 data/test_R1.fastq --R2 data/test_R2.fastq -b A -f data/test_guide_info.csv -o test_res/ -r"
    )

def test_guide_count():
    os.system(
         "bean-count --R1 data/test_R1.fastq --R2 data/test_R2.fastq -b A -f data/test_guide_info.csv -o test_res/ -g"
            )

def test_count_samples():
    os.system(
        "python ../../bin/bean-count-samples --input data/sample_list.csv -b A -f data/test_guide_info.csv -o test_res/ -r --guide-start-seq=GGAAAGGACGAAACACCG"
    )
