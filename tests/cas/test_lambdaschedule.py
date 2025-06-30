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
