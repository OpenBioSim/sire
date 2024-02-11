import sire as sr

import pytest


def test_stream_lever(tmpdir):
    dir = tmpdir.mkdir("test_stream_lever")

    s3file = str(dir.join("lever.s3"))

    l = sr.cas.LambdaSchedule.standard_morph()

    l.add_force("bond")
    l.add_force("angle")

    l.add_lever("charge")
    l.add_lever("sigma")

    sr.stream.save(l, s3file)

    l2 = sr.stream.load(s3file)

    assert l == l2
