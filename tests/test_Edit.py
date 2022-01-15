import crisprep

def test_hash():
	d = dict()
	e1 = Edit(1, "A", "G")
	d[e1] = 1

	e2 = Edit(2, "A", "G")
	d[e2] = 2

	e3 = crisprep.Edit(0, "A", "G", 1)
	assert e3 in d.keys()
	assert d[e3] == 1
