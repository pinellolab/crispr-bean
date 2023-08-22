import pytest
from bean import Edit


@pytest.mark.order(18)
def test_hash():
    e1 = Edit(1, "A", "G")
    d = {e1: 1}
    e2 = Edit(2, "A", "G")
    d[e2] = 2

    e3 = Edit(0, "A", "G", 1)
    assert e3 not in d
