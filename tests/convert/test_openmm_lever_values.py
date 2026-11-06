import sire as sr

import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_unchanged_levers_not_reported():
    """
    Levers that don't change, apart from round-off, aren't reported.
    """
    mols = sr.load_test_files("benzene.prm7", "benzene.rst")
    benzene = sr.morph.decouple(mols[0], as_new_molecule=False)

    p = benzene.perturbation().to_openmm()
    df = p.get_lever_values(schedule=benzene.property("schedule"))

    assert list(df.columns) == [
        "clj-charge-1",
        "clj-charge-7",
        "clj-epsilon-1",
        "clj-epsilon-7",
        "clj-alpha-1",
        "ghost/ghost-kappa-1",
    ]
