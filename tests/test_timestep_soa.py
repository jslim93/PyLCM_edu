"""The persistent-SoA driver must stay physical (ghost-free, integer A, finite)
and reproduce the object path's liquid water within Monte-Carlo tolerance."""
import numpy as np
import matplotlib; matplotlib.use("Agg")

from PyLCM.timestep_soa import run_soa
from validation.phys_harness import run as run_obj


def test_soa_run_is_physical_and_integer():
    out, (M, A) = run_soa(seed=0, n_ptcl=500, nt=200, collisions=True)
    assert np.all(np.isfinite(M)) and np.all(np.isfinite(A))
    assert np.all(A == np.round(A))                 # integer multiplicity preserved
    assert np.sum((M <= 0) & (A > 0)) == 0          # no zero-radius ghosts
    last = out[max(out)]
    assert last["qc"] >= 0 and last["qr"] >= 0


def test_soa_liquid_water_matches_object_ensemble():
    # 3-seed ensemble mean of total water (qc+qr) — robust mass integral.
    def lwc(fn, **kw):
        vals = []
        for s in range(3):
            if fn is run_obj:
                d, _ = run_obj(seed=s, aerosol="maritime", n_ptcl=1000, nt=1000, collisions=True, collect=(1000,))
            else:
                d, _ = run_soa(seed=s, n_ptcl=1000, nt=1000, N_raw=(100., 20.), mu_um=(.08, .4),
                               sig=(1.6, 2.0), kappa=1.0, RH=0.88, collisions=True, collect=(1000,))
            vals.append(d[1000]["qc"] + d[1000]["qr"])
        return np.mean(vals)
    obj = lwc(run_obj)
    soa = lwc(run_soa)
    assert abs(soa - obj) / obj < 0.15, f"LWC object={obj:.3f} soa={soa:.3f}"
