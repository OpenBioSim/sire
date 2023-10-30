import pytest

import sire as sr


def test_options():
    class TestOptions(sr.options.Option):
        A = "a", "Option A"
        B = "b", "Option B"
        C = "c", "Option C"

        @staticmethod
        def create(option: str):
            return sr.options.Option._create(TestOptions, option)

        @staticmethod
        def options(include_docs: bool = False):
            return sr.options.Option._options(
                TestOptions, include_docs=include_docs
            )

    assert TestOptions.A == "a"
    assert TestOptions.B == "b"
    assert TestOptions.C == "c"

    assert TestOptions.A.value == "a"
    assert TestOptions.B.value == "b"
    assert TestOptions.C.value == "c"

    assert TestOptions.A.__doc__ == "Option A"
    assert TestOptions.B.__doc__ == "Option B"
    assert TestOptions.C.__doc__ == "Option C"

    assert TestOptions.create("a") == TestOptions.A
    assert TestOptions.create("b") == TestOptions.B
    assert TestOptions.create("c") == TestOptions.C

    assert TestOptions.options() == ["a", "b", "c"]
    assert TestOptions.options(include_docs=True) == [
        ("a", "Option A"),
        ("b", "Option B"),
        ("c", "Option C"),
    ]


@pytest.mark.parametrize(
    "cls, option, expected",
    [
        (sr.options.Platform, "CPU", "cpu"),
        (
            sr.options.Constraint,
            "  H-bonds-h-ANGLES   ",
            "h_bonds_h_angles",
        ),
        (sr.options.PerturbableConstraint, "None", "none"),
        (sr.options.Cutoff, "Particle Mesh Ewald", "pme"),
        (sr.options.Cutoff, "REACTION        FIELD", "rf"),
        (sr.options.Cutoff, "   No           CutoFF   ", "none"),
        (sr.options.Integrator, "Langevin", "langevin"),
        (sr.options.Integrator, "Langevin-Middle", "langevin_middle"),
    ],
)
def test_dynamics_options(cls, option, expected):
    assert cls.create(option) == expected
