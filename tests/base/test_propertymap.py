import sire as sr


def test_propertymap():
    from sire.base import create_map

    m = create_map(
        {
            "cat": "meow",
            "dog": "woof",
            "number": 42,
            "space": sr.vol.Cartesian(),
        }
    )

    assert m.specified("cat")
    assert m.specified("dog")
    assert m.specified("number")
    assert m.specified("space")
    assert not m.specified("fish")

    assert m["cat"].has_source()
    assert not m["cat"].has_value()
    assert m["cat"].source() == "meow"

    assert m["dog"].has_source()
    assert not m["dog"].has_value()
    assert m["dog"].source() == "woof"

    assert not m["number"].has_source()
    assert m["number"].has_value()
    assert m["number"].value() == 42

    assert not m["space"].has_source()
    assert m["space"].has_value()
    assert m["space"].value() == sr.vol.Cartesian()


def test_propertymap_prefix():
    from sire.base import create_map

    m = create_map(
        {
            "cat": "meow",
            "cat0": "purr",
            "cat1": "hiss",
            "dog": "woof",
            "0dog": "bark",
            "1dog": "howl",
            "space": sr.vol.Cartesian(),
        }
    )

    props = ["cat", "space", "dog", "fish"]

    m0 = m.add_suffix("0", props)
    m1 = m.add_suffix("1", props)

    assert m0["cat"].source() == "purr"
    assert m1["cat"].source() == "hiss"

    assert m0["space"].value() == sr.vol.Cartesian()
    assert m0["space"].value() == sr.vol.Cartesian()

    assert m0["dog"].source() == "woof0"
    assert m1["dog"].source() == "woof1"

    assert m0["fish"].source() == "fish0"
    assert m1["fish"].source() == "fish1"

    m0 = m.add_prefix("0", props)
    m1 = m.add_prefix("1", props)

    print(m0)
    print(m1)

    assert m0["cat"].source() == "0meow"
    assert m1["cat"].source() == "1meow"

    assert m0["space"].value() == sr.vol.Cartesian()
    assert m0["space"].value() == sr.vol.Cartesian()

    assert m0["dog"].source() == "bark"
    assert m1["dog"].source() == "howl"

    assert m0["fish"].source() == "0fish"
    assert m1["fish"].source() == "1fish"
