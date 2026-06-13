"""Benchmark the persistent-SoA engine vs the validated object physics.

Two axes:
  SPEED   — wall time / throughput across particle counts (condensation+collision).
  PHYSICS — multi-seed ensemble agreement of q/N diagnostics (SoA must reproduce
            the object path, since it reuses our validated kernels).
"""
import warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import time
import numpy as np

from PyLCM.timestep_soa import run_soa
from validation.phys_harness import run as run_obj

MARITIME = dict(N_raw=(100., 20.), mu_um=(.08, .4), sig=(1.6, 2.0), kappa=1.0)


def speed():
    print("=" * 70)
    print("SPEED — object (validated) vs SoA, condensation+collision, nt=1500")
    print("=" * 70)
    print(f"{'n_ptcl':>8} {'object(s)':>10} {'SoA(s)':>9} {'speedup':>8} {'SoA Mpt/s':>11}")
    run_soa(seed=0, n_ptcl=300, nt=10)  # warm numba
    for n in (2000, 5000, 10000):
        t = time.time(); run_obj(seed=0, aerosol="maritime", n_ptcl=n, nt=1500,
                                 collisions=True, collect=(1500,)); to = time.time() - t
        t = time.time(); run_soa(seed=0, n_ptcl=n, nt=1500, RH=0.88, collisions=True,
                                 collect=(1500,), **MARITIME); ts = time.time() - t
        print(f"{n:>8} {to:>10.2f} {ts:>9.2f} {to/ts:>7.1f}x {n*1500/ts/1e6:>10.1f}")


def physics(seeds=5):
    print("\n" + "=" * 70)
    print(f"PHYSICS — {seeds}-seed ensemble mean at t=1500 (n=2000, maritime)")
    print("=" * 70)
    keys = ["qc+qr", "qr", "Nc", "Nr"]
    obj = {k: [] for k in keys}
    soa = {k: [] for k in keys}
    for s in range(seeds):
        od, _ = run_obj(seed=s, aerosol="maritime", n_ptcl=2000, nt=1500,
                        collisions=True, collect=(1500,))
        sd, _ = run_soa(seed=s, n_ptcl=2000, nt=1500, RH=0.88, collisions=True,
                        collect=(1500,), **MARITIME)
        o, p = od[1500], sd[1500]
        for d, r in ((obj, o), (soa, p)):
            d["qc+qr"].append(r["qc"] + r["qr"]); d["qr"].append(r["qr"])
            d["Nc"].append(r["NC"]); d["Nr"].append(r["NR"])
    print(f"{'diag':>7} {'object':>10} {'SoA':>10} {'rel diff':>10}")
    for k in keys:
        mo, ms = np.mean(obj[k]), np.mean(soa[k])
        print(f"{k:>7} {mo:>10.4f} {ms:>10.4f} {abs(ms-mo)/max(abs(mo),1e-9):>9.3f}")
    print("  (water integrals qc+qr/qr should match tightly; N counts are RNG-noisy)")


if __name__ == "__main__":
    speed()
    physics()
    print("\nDone.")
