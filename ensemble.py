"""Parallel ensemble runner for PyLCM.

Each member is an independent stochastic realization (different RNG seed).
Members are embarrassingly parallel, so we fan them across CPU cores with joblib.
This gives both ensemble (mean + spread) capability and a large wall-clock win
for multi-run studies, without touching the validated single-run physics.
"""
import numpy as np
from joblib import Parallel, delayed


def run_member(seed, run_single, run_kwargs):
    """Run one ensemble member with a fixed seed.

    `run_single(**run_kwargs)` must return a 1-D array-like diagnostic time
    series (e.g. LWC or rain number concentration over time).
    """
    np.random.seed(seed)
    return np.asarray(run_single(**run_kwargs))


def run_ensemble(run_single, run_kwargs, n_members=10, n_jobs=-1, base_seed=0):
    """Run `n_members` members in parallel.

    Returns (members, mean, lo, hi) where members is (n_members, n_t) and
    lo/hi are the 10th/90th percentile envelope across members.
    """
    seeds = [base_seed + i for i in range(n_members)]
    results = Parallel(n_jobs=n_jobs)(
        delayed(run_member)(s, run_single, run_kwargs) for s in seeds
    )
    members = np.vstack(results)
    mean = members.mean(axis=0)
    lo = np.percentile(members, 10, axis=0)
    hi = np.percentile(members, 90, axis=0)
    return members, mean, lo, hi


def plot_envelope(time, mean, lo, hi, label="", ax=None):
    """Plot ensemble mean with a shaded 10-90 percentile envelope."""
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(time, mean, label=label)
    ax.fill_between(time, lo, hi, alpha=0.3)
    ax.legend()
    return ax
