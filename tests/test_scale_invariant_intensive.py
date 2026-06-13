import numpy as np
from validation.golden_setup import run_condensation_only

# Recorded from the pre-change run (air_mass=100, float A): condensation is
# scale-invariant, so T and q must be unchanged (to integer-rounding tolerance)
# after the PARCEL_AIR_MASS + integer-A change.
OLD_T = 291.302817
OLD_Q = 0.01358930


def test_condensation_T_q_are_scale_invariant():
    out = run_condensation_only()
    assert np.isclose(out["T"], OLD_T, rtol=1e-4)
    assert np.isclose(out["q"], OLD_Q, rtol=1e-3)
