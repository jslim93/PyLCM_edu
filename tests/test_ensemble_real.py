import numpy as np
from ensemble import run_ensemble
from pylcm_run import run_single_series


def test_real_ensemble_runs_and_is_seed_deterministic():
    # Tiny run so the test stays fast; exercises the real condensation+collision path.
    kw = {"n_ptcl": 200, "nt": 100, "collect_every": 50, "diagnostic": "LWC"}
    m1, mean1, lo1, hi1 = run_ensemble(run_single_series, kw, n_members=2, n_jobs=1, base_seed=7)
    m2, *_ = run_ensemble(run_single_series, kw, n_members=2, n_jobs=1, base_seed=7)
    assert m1.shape == (2, 2)            # 2 members, nt/collect_every = 2 samples
    assert np.allclose(m1, m2)           # same seeds -> identical members
    assert np.all(lo1 <= mean1 + 1e-12) and np.all(mean1 - 1e-12 <= hi1)
    assert np.all(np.isfinite(m1))
