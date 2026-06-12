import numpy as np
import pytest
from validation.golden_setup import run_condensation_only

GOLDEN = "tests/golden/condensation_golden.npz"


def test_condensation_matches_golden():
    """The current condensation path must reproduce the committed golden snapshot.

    This anchors the SoA/numba rewrite: any fast condensation_fn passed to
    run_condensation_only must also pass this within rtol.
    """
    ref = np.load(GOLDEN)
    out = run_condensation_only()
    assert np.allclose(out["M"], ref["M"], rtol=1e-9, atol=0.0)
    assert np.allclose(out["A"], ref["A"], rtol=1e-9, atol=0.0)
    assert out["T"] == pytest.approx(float(ref["T"]), rel=1e-9)
    assert out["q"] == pytest.approx(float(ref["q"]), rel=1e-9)


def test_condensation_is_deterministic():
    a = run_condensation_only()
    b = run_condensation_only()
    assert np.array_equal(a["M"], b["M"]) and a["T"] == b["T"]
