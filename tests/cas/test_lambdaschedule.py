import sire as sr
import pytest
import random


def _assert_same_equation(x, eq1, eq2):
    for i in range(100):
        val = {x: random.uniform(0, 1)}

        assert eq1.evaluate(val) == pytest.approx(eq2.evaluate(val), 1e-5)


def test_charge_scale():
    l = sr.cas.LambdaSchedule.standard_morph()

    morph_equation = l.get_equation(stage="morph")

    l.add_charge_scale_stages()

    _assert_same_equation(l.lam(), l.get_equation(stage="morph"), morph_equation)

    assert l.get_stages() == ["decharge", "morph", "recharge"]

    assert l.get_levers() == ["charge"]

    assert l.get_constant("γ") == 0.2

    gamma = l.get_constant_symbol("γ")

    scaled_morph = gamma * morph_equation

    _assert_same_equation(
        l.lam(), l.get_equation(stage="morph", lever="charge"), scaled_morph
    )

    _assert_same_equation(
        gamma, l.get_equation(stage="morph", lever="charge"), scaled_morph
    )

    l.set_equation(stage="recharge", force="ghost/ghost", lever="charge", equation=0.5)

    l.set_equation(stage="*", force="ghost/ghost", lever="*", equation=1.5)

    _assert_same_equation(
        l.lam(),
        l.get_equation(stage="recharge", force="ghost/ghost", lever="charge"),
        sr.cas.Expression(0.5),
    )

    _assert_same_equation(
        l.lam(),
        l.get_equation(stage="recharge", lever="LJ"),
        (1.0 - ((1.0 - gamma) * (1.0 - l.lam()))) * l.final(),
    )

    _assert_same_equation(l.lam(), l.get_equation(stage="recharge"), l.final())

    _assert_same_equation(
        l.lam(),
        l.get_equation(stage="recharge", force="ghost/non-ghost", lever="LJ"),
        l.final(),
    )

    _assert_same_equation(
        l.lam(),
        l.get_equation(stage="recharge", force="ghost/ghost", lever="LJ"),
        sr.cas.Expression(1.5),
    )


def test_lambdaschedule():
    l = sr.cas.LambdaSchedule.standard_morph()

    morph_equation = (1 - l.lam()) * l.initial() + l.lam() * l.final()
    morph2_equation = (1 - l.lam() ** 2) * l.initial() + l.lam() ** 2 * l.final()
    morph3_equation = l.lam() * l.initial()
    morph4_equation = (1 - l.lam()) * l.final()
    morph5_equation = l.lam()

    assert l.get_stages() == ["morph"]

    assert len(l.get_levers()) == 0
    assert len(l.get_forces()) == 0

    _assert_same_equation(l.lam(), l.get_equation("morph"), morph_equation)

    l.set_equation(stage="morph", lever="charge", equation=morph2_equation)

    assert l.get_stages() == ["morph"]
    assert l.get_levers() == ["charge"]
    assert len(l.get_forces()) == 0

    _assert_same_equation(l.lam(), l.get_equation("morph", "charge"), morph2_equation)

    l.set_equation(stage="morph", force="CLJ", equation=morph3_equation)

    assert l.get_stages() == ["morph"]
    assert l.get_levers() == ["charge"]
    assert l.get_forces() == ["CLJ"]

    _assert_same_equation(
        l.lam(), l.get_equation(stage="morph", force="charge"), morph3_equation
    )

    l.set_equation(stage="morph", force="CLJ", lever="LJ", equation=morph4_equation)

    assert l.get_stages() == ["morph"]
    assert l.get_levers() == ["charge", "LJ"]
    assert l.get_forces() == ["CLJ"]

    _assert_same_equation(l.lam(), l.get_equation(stage="morph"), morph_equation)

    _assert_same_equation(
        l.lam(), l.get_equation(stage="morph", force="charge"), morph2_equation
    )

    _assert_same_equation(
        l.lam(), l.get_equation(stage="morph", force="charge"), morph3_equation
    )

    _assert_same_equation(
        l.lam(), l.get_equation(stage="morph", force="CLJ", lever="LJ"), morph4_equation
    )

    l.prepend_stage("scale_up", morph5_equation)

    assert l.get_stages() == ["scale_up", "morph"]

    assert l.get_stages() == ["scale_up", "morph"]
    assert l.get_levers() == ["charge", "LJ"]
    assert l.get_forces() == ["CLJ"]

    _assert_same_equation(l.lam(), l.get_equation(stage="morph"), morph_equation)

    _assert_same_equation(
        l.lam(), l.get_equation(stage="morph", force="charge"), morph2_equation
    )

    _assert_same_equation(
        l.lam(), l.get_equation(stage="morph", force="charge"), morph3_equation
    )

    _assert_same_equation(
        l.lam(), l.get_equation(stage="morph", force="CLJ", lever="LJ"), morph4_equation
    )

    _assert_same_equation(l.lam(), l.get_equation(stage="scale_up"), morph5_equation)


@pytest.mark.parametrize(
    "force, lever, contained",
    [("ghost-14", "kappa", True), ("ghost-14", "epsilon", True)],
)
def test_has_force_specific_equation(force, lever, contained):
    l = sr.cas.LambdaSchedule.standard_decouple()
    assert l.has_force_specific_equation("decouple", force, lever) == contained


def test_coupled_lever_default():
    """cmap::cmap_grid should follow torsion::torsion_k by default."""
    l = sr.cas.LambdaSchedule.standard_morph()
    morph_equation = (1 - l.lam()) * l.initial() + l.lam() * l.final()

    # With no custom equations, both should return the stage default.
    _assert_same_equation(
        l.lam(),
        l.get_equation(stage="morph", force="torsion", lever="torsion_k"),
        morph_equation,
    )
    _assert_same_equation(
        l.lam(),
        l.get_equation(stage="morph", force="cmap", lever="cmap_grid"),
        morph_equation,
    )


def test_coupled_lever_follows_torsion_k():
    """Setting a custom torsion_k equation should automatically apply to cmap_grid."""
    l = sr.cas.LambdaSchedule.standard_morph()
    custom_eq = l.lam() ** 2 * l.final() + (1 - l.lam() ** 2) * l.initial()

    l.set_equation(
        stage="morph", force="torsion", lever="torsion_k", equation=custom_eq
    )

    # cmap_grid should now follow the custom torsion_k equation via coupling.
    _assert_same_equation(
        l.lam(),
        l.get_equation(stage="morph", force="cmap", lever="cmap_grid"),
        custom_eq,
    )


def test_coupled_lever_explicit_override():
    """An explicit cmap_grid equation should take precedence over the coupling."""
    l = sr.cas.LambdaSchedule.standard_morph()
    torsion_eq = l.lam() ** 2 * l.final() + (1 - l.lam() ** 2) * l.initial()
    cmap_eq = l.initial()  # freeze CMAP at λ=0

    l.set_equation(
        stage="morph", force="torsion", lever="torsion_k", equation=torsion_eq
    )
    l.set_equation(stage="morph", force="cmap", lever="cmap_grid", equation=cmap_eq)

    # cmap_grid should use its own explicit equation, not torsion_k's.
    _assert_same_equation(
        l.lam(),
        l.get_equation(stage="morph", force="cmap", lever="cmap_grid"),
        cmap_eq,
    )
    # torsion_k should be unaffected.
    _assert_same_equation(
        l.lam(),
        l.get_equation(stage="morph", force="torsion", lever="torsion_k"),
        torsion_eq,
    )


def test_remove_coupled_lever():
    """Removing the coupling makes cmap_grid fall back to the stage default."""
    l = sr.cas.LambdaSchedule.standard_morph()
    morph_equation = (1 - l.lam()) * l.initial() + l.lam() * l.final()
    custom_eq = l.lam() ** 2 * l.final() + (1 - l.lam() ** 2) * l.initial()

    l.set_equation(
        stage="morph", force="torsion", lever="torsion_k", equation=custom_eq
    )
    l.remove_coupled_lever(force="cmap", lever="cmap_grid")

    # cmap_grid should now use the stage default, not follow torsion_k.
    _assert_same_equation(
        l.lam(),
        l.get_equation(stage="morph", force="cmap", lever="cmap_grid"),
        morph_equation,
    )


def test_couple_lever_custom():
    """coupleLever can set an arbitrary coupling between levers."""
    l = sr.cas.LambdaSchedule.standard_morph()
    custom_eq = l.lam() ** 2 * l.final() + (1 - l.lam() ** 2) * l.initial()

    # Couple bond_k to torsion_k (unusual, but should work).
    l.couple_lever(
        force="bond",
        lever="bond_k",
        fallback_force="torsion",
        fallback_lever="torsion_k",
    )
    l.set_equation(
        stage="morph", force="torsion", lever="torsion_k", equation=custom_eq
    )

    _assert_same_equation(
        l.lam(),
        l.get_equation(stage="morph", force="bond", lever="bond_k"),
        custom_eq,
    )
