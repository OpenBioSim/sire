import sire as sr

import pytest


def test_stream_header(tmpdir, ala_mols):
    mols = ala_mols

    dir = tmpdir.mkdir("test_stream_header")

    s3file = str(dir.join("output.s3"))

    sr.stream.save(mols, s3file)

    m = sr.stream.load(s3file)

    assert mols.num_molecules() == m.num_molecules()

    h = sr.stream.get_data_header(s3file)

    assert not h.has_property("cat")

    assert h.property("cat") is None

    sr.stream.set_header_property("cat", [1, 2, 3, 4, 5, 6])
    sr.stream.set_header_property("mouse", "hello")
    sr.stream.set_header_property("fish", 3.141)

    sr.stream.save(mols, s3file)

    h = sr.stream.get_data_header(s3file)

    assert h.property("cat") == [1, 2, 3, 4, 5, 6]
    assert h.property("mouse") == "hello"
    assert h.property("fish") == 3.141

    assert h.has_property("cat")
    assert h.has_property("mouse")
    assert h.has_property("fish")

    m = sr.stream.load(s3file)

    assert mols.num_molecules() == m.num_molecules()

    sr.stream.set_header_property("cat", None)

    sr.stream.save(mols, s3file)

    h = sr.stream.get_data_header(s3file)

    assert h.property("cat") is None
    assert h.property("mouse") == "hello"
    assert h.property("fish") == 3.141

    assert not h.has_property("cat")
    assert h.has_property("mouse")
    assert h.has_property("fish")
