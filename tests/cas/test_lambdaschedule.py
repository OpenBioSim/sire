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


def test_stage_weights_default():
    """New stages should have weight 1.0 by default."""
    l = sr.cas.LambdaSchedule.standard_morph()

    assert l.get_stage_weights() == [1.0]
    assert l.get_stage_weight("morph") == pytest.approx(1.0)


def test_stage_weights_equal_split():
    """With equal weights, lambda space is split evenly (backward compat)."""
    l = sr.cas.LambdaSchedule.standard_morph()
    l.add_morph_stage("morph2")

    assert l.get_stage_weights() == [1.0, 1.0]

    # Boundary between stages is at 0.5
    assert l.get_stage(0.0) == "morph"
    assert l.get_stage(0.49) == "morph"
    assert l.get_stage(0.5) == "morph2"
    assert l.get_stage(1.0) == "morph2"

    # Local lambda at midpoint of each stage
    assert l.get_lambda_in_stage(0.25) == pytest.approx(0.5)
    assert l.get_lambda_in_stage(0.75) == pytest.approx(0.5)


def test_stage_weights_three_equal():
    """Three equal-weight stages each occupy one third of lambda space."""
    l = sr.cas.LambdaSchedule()
    l.add_stage("a", l.lam())
    l.add_stage("b", l.lam())
    l.add_stage("c", l.lam())

    assert l.get_stage_weights() == [1.0, 1.0, 1.0]

    assert l.get_stage(0.0) == "a"
    assert l.get_stage(1 / 3 - 0.001) == "a"
    assert l.get_stage(1 / 3) == "b"
    assert l.get_stage(2 / 3 - 0.001) == "b"
    assert l.get_stage(2 / 3) == "c"
    assert l.get_stage(1.0) == "c"

    # Local lambda at the start of each stage
    assert l.get_lambda_in_stage(0.0) == pytest.approx(0.0)
    assert l.get_lambda_in_stage(1 / 3) == pytest.approx(0.0, abs=1e-10)
    assert l.get_lambda_in_stage(2 / 3) == pytest.approx(0.0, abs=1e-10)


def test_stage_weights_unequal_two_stages():
    """A stage with weight 2 occupies twice the lambda range of weight 1."""
    l = sr.cas.LambdaSchedule()
    l.add_stage("small", l.lam(), weight=1.0)
    l.add_stage("large", l.lam(), weight=2.0)

    assert l.get_stage_weights() == [1.0, 2.0]

    # Stage boundary is at 1/3
    assert l.get_stage(0.0) == "small"
    assert l.get_stage(1 / 3 - 0.001) == "small"
    assert l.get_stage(1 / 3) == "large"
    assert l.get_stage(1.0) == "large"

    # Local lambda at midpoint of small stage (lambda=1/6)
    assert l.get_lambda_in_stage(1 / 6) == pytest.approx(0.5)

    # Local lambda at midpoint of large stage (lambda=1/3 + 1/3 = 2/3)
    assert l.get_lambda_in_stage(2 / 3) == pytest.approx(0.5)

    # lambda=0 is start of small stage
    assert l.get_lambda_in_stage(0.0) == pytest.approx(0.0)

    # lambda=1 is end of large stage
    assert l.get_lambda_in_stage(1.0) == pytest.approx(1.0)


def test_stage_weights_morph_value():
    """With weights [1, 2], morph produces correct interpolated values."""
    l = sr.cas.LambdaSchedule()
    l.add_morph_stage("a", weight=1.0)
    l.add_morph_stage("b", weight=2.0)

    # At global lambda=1/6 we are halfway through stage "a" (local lam=0.5)
    assert l.morph(initial=0.0, final=1.0, lambda_value=1 / 6) == pytest.approx(0.5)

    # At global lambda=1/3 we are at the start of stage "b" (local lam=0.0)
    assert l.morph(initial=0.0, final=1.0, lambda_value=1 / 3) == pytest.approx(0.0)

    # At global lambda=2/3 we are halfway through stage "b" (local lam=0.5)
    assert l.morph(initial=0.0, final=1.0, lambda_value=2 / 3) == pytest.approx(0.5)


def test_set_stage_weight():
    """set_stage_weight updates a weight after the stage is added."""
    l = sr.cas.LambdaSchedule.standard_morph()
    l.add_morph_stage("morph2")

    l.set_stage_weight("morph", 1.0)
    l.set_stage_weight("morph2", 3.0)

    assert l.get_stage_weight("morph") == pytest.approx(1.0)
    assert l.get_stage_weight("morph2") == pytest.approx(3.0)
    assert l.get_stage_weights() == pytest.approx([1.0, 3.0])

    # Boundary is now at 0.25 (1 out of 4 total weight units)
    assert l.get_stage(0.24) == "morph"
    assert l.get_stage(0.25) == "morph2"


def test_stage_weight_prepend():
    """prepend_stage respects the weight argument."""
    l = sr.cas.LambdaSchedule()
    l.add_stage("b", l.lam(), weight=1.0)
    l.prepend_stage("a", l.lam(), weight=3.0)

    assert l.get_stage_weights() == pytest.approx([3.0, 1.0])
    assert l.get_stages() == ["a", "b"]

    # Stage "a" occupies 3/4 of lambda space; boundary at 0.75
    assert l.get_stage(0.74) == "a"
    assert l.get_stage(0.75) == "b"


def test_stage_weight_insert():
    """insert_stage respects the weight argument."""
    l = sr.cas.LambdaSchedule()
    l.add_stage("first", l.lam(), weight=1.0)
    l.add_stage("last", l.lam(), weight=1.0)
    l.insert_stage(1, "middle", l.lam(), weight=2.0)

    assert l.get_stages() == ["first", "middle", "last"]
    assert l.get_stage_weights() == pytest.approx([1.0, 2.0, 1.0])

    # Total weight=4; boundaries at 0.25 and 0.75
    assert l.get_stage(0.24) == "first"
    assert l.get_stage(0.25) == "middle"
    assert l.get_stage(0.74) == "middle"
    assert l.get_stage(0.75) == "last"


def test_stage_weight_remove():
    """remove_stage also removes the corresponding weight."""
    l = sr.cas.LambdaSchedule()
    l.add_stage("a", l.lam(), weight=1.0)
    l.add_stage("b", l.lam(), weight=2.0)
    l.add_stage("c", l.lam(), weight=3.0)

    l.remove_stage("b")

    assert l.get_stages() == ["a", "c"]
    assert l.get_stage_weights() == pytest.approx([1.0, 3.0])


def test_stage_weight_clear():
    """clear() removes stage weights alongside stages."""
    l = sr.cas.LambdaSchedule()
    l.add_stage("a", l.lam(), weight=2.0)
    l.clear()

    assert l.get_stages() == []
    assert l.get_stage_weights() == []


def test_stage_weight_invalid():
    """A non-positive weight should raise an error."""
    l = sr.cas.LambdaSchedule()

    with pytest.raises(Exception):
        l.add_stage("a", l.lam(), weight=0.0)

    with pytest.raises(Exception):
        l.add_stage("b", l.lam(), weight=-1.0)

    l.add_stage("c", l.lam())

    with pytest.raises(Exception):
        l.set_stage_weight("c", 0.0)

    with pytest.raises(Exception):
        l.set_stage_weight("c", -0.5)


def test_stage_weight_add_morph_stage():
    """add_morph_stage accepts a weight argument."""
    l = sr.cas.LambdaSchedule()
    l.add_morph_stage("decharge", weight=1.0)
    l.add_morph_stage("morph", weight=2.0)
    l.add_morph_stage("recharge", weight=1.0)

    assert l.get_stage_weights() == pytest.approx([1.0, 2.0, 1.0])

    # "morph" occupies the middle half (0.25 to 0.75)
    assert l.get_stage(0.24) == "decharge"
    assert l.get_stage(0.25) == "morph"
    assert l.get_stage(0.74) == "morph"
    assert l.get_stage(0.75) == "recharge"


def test_stage_weight_add_decouple_stage():
    """add_decouple_stage accepts a weight argument."""
    l = sr.cas.LambdaSchedule()
    l.add_morph_stage("pre", weight=1.0)
    l.add_decouple_stage("decouple", weight=2.0)

    assert l.get_stage_weights() == pytest.approx([1.0, 2.0])
    assert l.get_stage(0.33) == "pre"
    assert l.get_stage(0.34) == "decouple"


def test_stage_weight_add_annihilate_stage():
    """add_annihilate_stage accepts a weight argument."""
    l = sr.cas.LambdaSchedule()
    l.add_annihilate_stage("annihilate", weight=3.0)
    l.add_morph_stage("morph", weight=1.0)

    assert l.get_stage_weights() == pytest.approx([3.0, 1.0])
    # Boundary at 0.75
    assert l.get_stage(0.74) == "annihilate"
    assert l.get_stage(0.75) == "morph"


def test_stage_weight_copy():
    """Copying a LambdaSchedule preserves weights."""
    import copy

    l = sr.cas.LambdaSchedule()
    l.add_morph_stage("a", weight=1.0)
    l.add_morph_stage("b", weight=3.0)

    l2 = copy.copy(l)
    assert l2.get_stage_weights() == pytest.approx([1.0, 3.0])

    # Modifying the copy should not affect the original
    l2.set_stage_weight("b", 1.0)
    assert l.get_stage_weight("b") == pytest.approx(3.0)


def test_stage_weight_equality():
    """Two schedules with different weights are not equal."""
    l1 = sr.cas.LambdaSchedule()
    l1.add_morph_stage("a", weight=1.0)
    l1.add_morph_stage("b", weight=2.0)

    l2 = sr.cas.LambdaSchedule()
    l2.add_morph_stage("a", weight=1.0)
    l2.add_morph_stage("b", weight=1.0)

    assert l1 != l2

    l2.set_stage_weight("b", 2.0)
    assert l1 == l2


def test_stage_weight_tostring():
    """toString shows weights when any stage has a non-default weight."""
    l = sr.cas.LambdaSchedule()
    l.add_morph_stage("a", weight=1.0)
    l.add_morph_stage("b", weight=2.0)

    s = str(l)
    assert "weight=1" in s
    assert "weight=2" in s

    # With all equal weights, no weight annotation
    l2 = sr.cas.LambdaSchedule.standard_morph()
    assert "weight" not in str(l2)


def test_stage_weight_get_stage_boundaries():
    """get_stage and get_lambda_in_stage agree at stage boundaries."""
    l = sr.cas.LambdaSchedule()
    l.add_morph_stage("a", weight=1.0)
    l.add_morph_stage("b", weight=3.0)

    # Total weight = 4; boundary at 0.25
    assert l.get_stage(0.0) == "a"
    assert l.get_stage(1.0) == "b"

    # Local lambda is 0.0 at start of each stage
    assert l.get_lambda_in_stage(0.0) == pytest.approx(0.0)
    assert l.get_lambda_in_stage(0.25) == pytest.approx(0.0, abs=1e-10)

    # Local lambda is 1.0 at end of last stage
    assert l.get_lambda_in_stage(1.0) == pytest.approx(1.0)

    # Midpoint of stage "a": lambda=0.125 → local lam=0.5
    assert l.get_lambda_in_stage(0.125) == pytest.approx(0.5)

    # Midpoint of stage "b": lambda=0.25 + 1.5/4 = 0.625 → local lam=0.5
    assert l.get_lambda_in_stage(0.625) == pytest.approx(0.5)
