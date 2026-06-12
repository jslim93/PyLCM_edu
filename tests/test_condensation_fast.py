import numpy as np
import pytest
from validation.golden_setup import run_condensation_only
from PyLCM.condensation_fast import drop_condensation_fast

GOLDEN = "tests/golden/condensation_golden.npz"


def test_fast_condensation_matches_golden():
    """The fast SoA/numba condensation path must reproduce the golden snapshot.

    Mirrors tests/test_golden_condensation.py but routes the run through
    drop_condensation_fast instead of the object-based reference.
    """
    ref = np.load(GOLDEN)
    out = run_condensation_only(condensation_fn=drop_condensation_fast)
    assert np.allclose(out["M"], ref["M"], rtol=1e-9, atol=0.0)
    assert np.allclose(out["A"], ref["A"], rtol=1e-9, atol=0.0)
    assert out["T"] == pytest.approx(float(ref["T"]), rel=1e-9)
    assert out["q"] == pytest.approx(float(ref["q"]), rel=1e-9)


def test_fast_condensation_is_deterministic():
    a = run_condensation_only(condensation_fn=drop_condensation_fast)
    b = run_condensation_only(condensation_fn=drop_condensation_fast)
    assert np.array_equal(a["M"], b["M"]) and a["T"] == b["T"]
