import os

def test_filter():
    os.system(
        "bean-filter data/bean_count_sample_list.h5ad -o tmp -s 0 -e 19 -w -b -t -ap 0.1 -sp 0.3"
    )