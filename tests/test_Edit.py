import pytest
from bean import Edit


@pytest.mark.order(0)
def test_hash():
    e1 = Edit(1, "A", "G")
    d = {e1: 1}
    e2 = Edit(2, "A", "G")
    d[e2] = 2

    e3 = Edit(0, "A", "G", offset=1)
    assert e3 not in d


@pytest.mark.order(1)
def test_hash_chrom():
    e1 = Edit(1, "A", "G", "chr18")
    d = {e1: 1}
    e2 = Edit(2, "A", "G", "chr19")
    d[e2] = 2

    e3 = Edit(2, "A", "G", "chr19")
    assert e3 in d
