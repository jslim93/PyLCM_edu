import numpy as np
from PyLCM.ensemble import run_ensemble


def _toy_run(scale=1.0):
    # Stand-in for a PyLCM single run: stochastic series depending on the RNG,
    # so that fixing the seed makes the result reproducible.
    return scale * np.cumsum(np.random.random(20))


def test_same_seed_is_deterministic():
    m1, *_ = run_ensemble(_toy_run, {"scale": 1.0}, n_members=3, n_jobs=1, base_seed=42)
    m2, *_ = run_ensemble(_toy_run, {"scale": 1.0}, n_members=3, n_jobs=1, base_seed=42)
    assert np.allclose(m1, m2)


def test_envelope_brackets_mean():
    _, mean, lo, hi = run_ensemble(_toy_run, {"scale": 1.0}, n_members=20, n_jobs=1, base_seed=0)
    assert np.all(lo <= mean + 1e-9) and np.all(mean - 1e-9 <= hi)
