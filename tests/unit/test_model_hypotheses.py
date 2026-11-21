"""`Hypothesis` and `HypothesisSet`: densities a map is computed for, without a new scan."""

import numpy as np
import pytest

from muoscene.geometry import BodyInfo, PathTable
from muoscene.model import Hypothesis, HypothesisSet

BODIES = (
    BodyInfo("soil", "ground", "ms-loam", 1.8),
    BodyInfo("cap", "structure", "ms-concrete", 2.3),
    BodyInfo("void", "anomaly", "ms-loam", 0.0, host="soil", drho_rel=-1.0),
    BodyInfo("air", "air", "ms-air", 0.0011),
)
UP = np.array([[0.0, 0.0, 1.0]])
PATHS = PathTable(BODIES, np.zeros((1, 3)), UP, np.array([[10.0, 2.0, 1.0, 0.3]]))  # m in every body


def test_constructors():
    assert Hypothesis.as_built() == Hypothesis("as built")
    assert Hypothesis.without("void").drho_rel == {"void": 0.0}
    assert Hypothesis.void(["a", "b"]).drho_rel == {"a": -1.0, "b": -1.0}
    assert Hypothesis.decompaction("void", -0.2).drho_rel == {"void": -0.2}
    assert (
        Hypothesis.without("void").name == "none" and Hypothesis.void("void", name="cavity").name == "cavity"
    )
    with pytest.raises(ValueError, match="at least one"):
        Hypothesis.without([])


def test_opacity_under_a_hypothesis_is_a_sum_over_the_same_table():
    """10 m of soil at 1.8, 2 m of concrete at 2.3 and 1 m of void: 1800 + 460 + 0 g/cm²."""
    assert Hypothesis.as_built().opacity_g_cm2(PATHS)[0] == pytest.approx(2260.0)
    assert Hypothesis.without("void").opacity_g_cm2(PATHS)[0] == pytest.approx(2440.0)
    assert Hypothesis.decompaction("void", -0.2).opacity_g_cm2(PATHS)[0] == pytest.approx(2260.0 + 144.0)
    drier = Hypothesis("drier", densities_g_cm3={"soil": 1.7}, drho_rel={"void": 0.0})
    by_material = drier.opacity_by_material_g_cm2(PATHS)
    assert by_material["ms-loam"][0] == pytest.approx(1870.0)  # the anomaly follows its host
    assert by_material["ms-concrete"][0] == pytest.approx(460.0)
    assert drier.to_dict() == {"densities_g_cm3": {"soil": 1.7}, "drho_rel": {"void": 0.0}}


def test_a_hypothesis_cannot_be_changed_and_checks_its_numbers():
    hypothesis = Hypothesis("drier", {"soil": 1.7}, {"void": -0.5})
    with pytest.raises(TypeError):
        hypothesis.densities_g_cm3["soil"] = 2.0
    with pytest.raises(AttributeError):
        hypothesis.name = "other"
    assert hash(hypothesis) == hash(Hypothesis("drier", {"soil": 1.7}, {"void": -0.5}))
    assert hypothesis != Hypothesis("drier", {"soil": 1.7})
    for arguments in (
        {"name": ""},
        {"name": "x", "densities_g_cm3": {"soil": -1.0}},
        {"name": "x", "densities_g_cm3": {"soil": float("nan")}},
        {"name": "x", "drho_rel": {"void": -1.5}},
        {"name": "x", "drho_rel": {"void": True}},
        {"name": "x", "drho_rel": {"": 0.0}},
    ):
        with pytest.raises(ValueError):
            Hypothesis(**arguments)
    with pytest.raises(TypeError, match="mapping"):
        Hypothesis("x", [("soil", 1.7)])


def test_a_set_has_an_order_and_different_names():
    none, void = Hypothesis.without("void"), Hypothesis.void("void")
    both = HypothesisSet([none, void])
    assert both.names == ("none", "void") and len(both) == 2
    assert list(both) == [none, void] and both["void"] is void and both[0] is none
    assert HypothesisSet(none).names == ("none",)
    assert both == HypothesisSet((none, void)) and hash(both) == hash(HypothesisSet((none, void)))
    with pytest.raises(ValueError, match="different names"):
        HypothesisSet([none, Hypothesis("none")])
    with pytest.raises(TypeError, match="at least one"):
        HypothesisSet([])
    with pytest.raises(TypeError):
        HypothesisSet(["none"])
    with pytest.raises(KeyError):
        both["other"]


def test_a_set_is_checked_against_the_bodies_of_a_scan():
    """A name the table does not have, a density for an anomaly, a change for a ground: refused by name."""
    HypothesisSet([Hypothesis.as_built(), Hypothesis.without("void")]).check(PATHS)
    for bad in (
        Hypothesis("a", drho_rel={"cavity": -1.0}),
        Hypothesis("b", densities_g_cm3={"void": 1.0}),
        Hypothesis("c", drho_rel={"soil": -0.1}),
        Hypothesis("d", densities_g_cm3={"air": 0.001}),
    ):
        with pytest.raises(ValueError, match=f"hypothesis '{bad.name}'"):
            HypothesisSet(bad).check(PATHS)
        with pytest.raises(ValueError, match=f"hypothesis '{bad.name}'"):
            bad.opacity_g_cm2(PATHS)
