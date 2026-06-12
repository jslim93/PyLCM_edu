import numpy as np
import matplotlib
matplotlib.use("Agg")
from pylcm_run import run_single_series


def test_mixing_off_matches_baseline():
    np.random.seed(0); base = run_single_series(n_ptcl=200, nt=100, collect_every=50)
    np.random.seed(0); same = run_single_series(n_ptcl=200, nt=100, collect_every=50, mixing=None)
    assert np.array_equal(base, same)              # mixing=None is a no-op


def test_inhomogeneous_mixing_reduces_cloud_water():
    from PyLCM.parcel import create_env_profiles
    from PyLCM.mixing import ParameterizedMixing
    qv, th, z_env = create_env_profiles(290.0, 0.010, 0.0, 95000.0, "Stable")
    mix = ParameterizedMixing(2e-3, 1.0, qv, th, z_env)
    np.random.seed(0); base = run_single_series(n_ptcl=200, nt=200, collect_every=200)
    np.random.seed(0); mixed = run_single_series(n_ptcl=200, nt=200, collect_every=200, mixing=mix)
    assert mixed[-1] <= base[-1] + 1e-9            # entrainment removes liquid water
