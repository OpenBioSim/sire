import math

import pytest
import sire as sr

pytestmark = pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)


def _convert(mols, platform, **options):
    m = {"cutoff_type": "PME", "cutoff": "9 A", "platform": platform}
    m.update(options)
    return sr.convert.to(mols, "openmm", map=m)


def _pme_parameters(omm):
    from openmm import NonbondedForce, unit

    for force in omm.getSystem().getForces():
        if isinstance(force, NonbondedForce):
            alpha, nx, ny, nz = force.getPMEParameters()
            return alpha.value_in_unit(unit.nanometer**-1), nx, ny, nz


def _box_lengths(omm):
    from openmm import unit

    return [
        math.sqrt(sum(x * x for x in v.value_in_unit(unit.nanometer)))
        for v in omm.getSystem().getDefaultPeriodicBoxVectors()
    ]


def test_pme_default(kigaki_mols, openmm_platform):
    omm = _convert(kigaki_mols, openmm_platform)

    alpha, _, _, _ = _pme_parameters(omm)

    assert alpha == 0


@pytest.mark.parametrize("grid", [32, [30, 32, 34]])
def test_pme_grid(kigaki_mols, openmm_platform, grid):
    omm = _convert(kigaki_mols, openmm_platform, pme_alpha=3.47, pme_grid=grid)

    alpha, nx, ny, nz = _pme_parameters(omm)

    expected = grid if isinstance(grid, list) else [grid] * 3

    assert alpha == pytest.approx(3.47)
    assert [nx, ny, nz] == expected


@pytest.mark.parametrize("fixture", ["kigaki_mols", "triclinic_protein"])
def test_pme_spacing(fixture, openmm_platform, request):
    mols = request.getfixturevalue(fixture)

    omm = _convert(mols, openmm_platform, pme_spacing="1.2 A")

    alpha, nx, ny, nz = _pme_parameters(omm)

    assert alpha == pytest.approx(math.sqrt(-math.log(2.0e-4)) / 0.9)
    assert [nx, ny, nz] == [math.ceil(x / 0.12) for x in _box_lengths(omm)]


@pytest.mark.parametrize(
    "options",
    [
        {"pme_alpha": 3.47},
        {"pme_grid": 32, "pme_spacing": "1.2 A"},
        {"pme_grid": 4},
        {"pme_grid": [32, 32]},
    ],
)
def test_pme_invalid(kigaki_mols, openmm_platform, options):
    with pytest.raises(RuntimeError, match="pme_"):
        _convert(kigaki_mols, openmm_platform, **options)
