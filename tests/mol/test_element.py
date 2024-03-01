import pytest
import sire as sr


def test_element():
    # Create an element
    el = sr.mol.Element("C")
    assert el.symbol() == "C"
    assert el.num_protons() == 6
    assert el.mass().to("g mol-1") == pytest.approx(12.01, abs=0.01)

    # Create an element from atomic number
    el = sr.mol.Element(6)
    assert el.symbol() == "C"

    # Create an element from name
    el = sr.mol.Element("Carbon")
    assert el.symbol() == "C"

    # Create an element from mass
    el = sr.mol.Element.element_with_mass(sr.u("12 g mol-1"))
    assert el.symbol() == "C"

    # Create a biological Iron
    el = sr.mol.Element.biological_element("Fe")

    assert el.symbol() == "Fe"
    assert el.num_protons() == 26
    assert el.mass().to("g mol-1") == pytest.approx(55.845)
    assert el.name() == "Iron"

    # remove this as a biological element
    sr.mol.Element.set_element_is_not_biological("Fe")

    el = sr.mol.Element.biological_element("Fe")

    assert el.symbol() != "Fe"
    assert el.num_protons() != 26

    # Add this back
    sr.mol.Element.set_element_is_biological("Fe")

    el = sr.mol.Element.biological_element("Fe")

    assert el.symbol() == "Fe"
    assert el.num_protons() == 26

    # get the list of biological elements
    elements = sr.mol.Element.get_biological_elements()

    assert sr.mol.Element("Fe") in elements
