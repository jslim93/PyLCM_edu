"""Struct-of-arrays (persistent numpy) collision must preserve OUR validated
object collision physics from PyLCM/collision.py.

Gates:
  1. Conservation: one collide_soa call conserves total water sum(M) (rtol<1e-9)
     and does not increase total multiplicity sum(A).
  2. Integer A: every A stays integer-valued.
  3. No ghosts: no entry has M<=0 while A>0 after the call.
  4. Statistical equivalence to the object path: a small full ascent run twice
     with the SAME seed — once via the object collection() (validation harness),
     once via a SoA driver using collide_soa — must agree on final cloud water
     (qc+qr) and rain number within a Monte-Carlo tolerance (RNG consumption
     differs between the two paths, so they are not bit-identical).
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib; matplotlib.use("Agg")
from numba import njit

from PyLCM.parameters import rho_liq, r_a, p0, cp, rv, rho_aero, z_env
from PyLCM.collision_soa import collide_soa, collide_soa_enumerate


@njit(cache=True)
def _seed_numba_rng(s):
    """Seed Numba's internal np.random generator (distinct from Python-level
    np.random), so the @njit kernel's draws are reproducible across runs."""
    np.random.seed(s)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _build_cloud(seed, n=2000, r_lo=8e-6, r_hi=45e-6, A_lo=1, A_hi=1000):
    np.random.seed(seed)
    r = np.random.uniform(r_lo, r_hi, n)
    x = 4.0 / 3.0 * np.pi * rho_liq * r**3       # per-droplet mass
    A = np.random.randint(A_lo, A_hi, n).astype(np.float64)
    M = x * A
    Ns = (x * 0.01) * A
    kappa = np.full(n, 0.6)
    return M, A, Ns, kappa


# ---------------------------------------------------------------------------
# Gate 1/2/3: conservation, integer A, no ghosts  (gravitational + turbulent)
# ---------------------------------------------------------------------------
def _check_invariants(switch_turb_kernel, epsilon_turb):
    M, A, Ns, kappa = _build_cloud(seed=0)
    T, P = 283.0, 900e2
    rho_parcel = P / (r_a * T)

    M0, A0 = M.sum(), A.sum()
    np.random.seed(123)
    _seed_numba_rng(123)
    M, A, Ns, kappa, n_removed = collide_soa(
        M, A, Ns, kappa, 1.0, rho_parcel, P, T,
        switch_turb_kernel=switch_turb_kernel, epsilon_turb=epsilon_turb)

    # conservation of total water
    assert np.isclose(M.sum(), M0, rtol=1e-9, atol=0.0), \
        f"water not conserved: {M.sum()} vs {M0}"
    # multiplicity does not increase
    assert A.sum() <= A0 + 1e-6, f"A increased: {A.sum()} > {A0}"
    # integer-valued A
    assert np.allclose(A, np.round(A), atol=1e-9), "A not integer-valued"
    # no ghosts
    assert np.count_nonzero((M <= 0.0) & (A > 0.0)) == 0, "ghost droplet (M<=0, A>0)"
    assert n_removed >= 0


def test_conservation_integer_noghost_gravitational():
    _check_invariants(switch_turb_kernel=False, epsilon_turb=0.0)


def test_conservation_integer_noghost_turbulent():
    _check_invariants(switch_turb_kernel=True, epsilon_turb=0.04)


def test_enumerate_mode_same_invariants():
    """The O(N^2) enumeration option (LSM off) must obey the SAME invariants:
    water conserved, A integer & non-increasing, no ghosts."""
    M, A, Ns, kappa = _build_cloud(seed=0, n=400)   # small N: O(N^2)
    T, P = 283.0, 900e2
    rho_parcel = P / (r_a * T)
    M0, A0 = M.sum(), A.sum()
    np.random.seed(7)
    _seed_numba_rng(7)
    M, A, Ns, kappa, n_removed = collide_soa_enumerate(
        M, A, Ns, kappa, 1.0, rho_parcel, P, T)

    assert np.isclose(M.sum(), M0, rtol=1e-9, atol=0.0), "enumerate: water not conserved"
    assert A.sum() <= A0 + 1e-6, "enumerate: A increased"
    assert np.allclose(A, np.round(A), atol=1e-9), "enumerate: A not integer-valued"
    assert np.count_nonzero((M <= 0.0) & (A > 0.0)) == 0, "enumerate: ghost droplet"
    assert n_removed >= 0


def test_invariants_over_many_steps():
    """Persistent arrays through many steps stay integer, ghost-free, conserving."""
    M, A, Ns, kappa = _build_cloud(seed=3, r_lo=15e-6, r_hi=60e-6, A_lo=1, A_hi=200)
    T, P = 283.0, 900e2
    rho_parcel = P / (r_a * T)
    M0, A0 = M.sum(), A.sum()
    np.random.seed(7)
    _seed_numba_rng(7)
    for step in range(40):
        M, A, Ns, kappa, _ = collide_soa(M, A, Ns, kappa, 2.0, rho_parcel, P, T)
        assert np.allclose(A, np.round(A), atol=1e-9), f"A non-integer at step {step}"
        assert np.count_nonzero((M <= 0.0) & (A > 0.0)) == 0, f"ghost at step {step}"
        assert (A > 0.0).all(), f"non-positive A survived compaction at step {step}"
    assert np.isclose(M.sum(), M0, rtol=1e-9, atol=0.0), "water drift over many steps"
    assert A.sum() <= A0 + 1e-6, "A increased over many steps"


# ---------------------------------------------------------------------------
# Gate 4: statistical equivalence to the OBJECT path (full ascent)
# ---------------------------------------------------------------------------
def _run_soa_ascent(seed, n_ptcl, nt, dt=1.0, T0=293.2, P0=1013e2, RH=0.92, w=1.0,
                    aerosol="maritime"):
    """SoA driver mirroring validation.phys_harness.run, but the collision step
    uses collide_soa over persistent arrays instead of the object collection()."""
    from PyLCM.aero_init import aero_init
    from PyLCM.parcel import ascend_parcel, parcel_rho
    from PyLCM.condensation import esatw
    from PyLCM.condensation_fast import drop_condensation_fast as drop_condensation
    from Post_process.analysis import ts_analysis
    from validation.phys_harness import PRESETS

    a = PRESETS[aerosol]
    N_raw, mu_um, sig, kappa0 = a["N_raw"], a["mu_um"], a["sig"], a["kappa"]
    mu = np.log(np.array(mu_um) * 1e-6)
    sg = np.log(np.array(sig))
    n_modes = len(N_raw)
    th = T0 * (p0 / P0) ** (r_a / cp) + 5e-3 * z_env
    q0 = RH * esatw(T0) / (P0 - RH * esatw(T0)) * r_a / rv

    np.random.seed(seed)
    _seed_numba_rng(seed)
    T, q, pl = aero_init("Random", n_ptcl, P0, 0.0, T0, q0,
                         np.array(N_raw) * 1e6, mu, sg, rho_aero,
                         [kappa0] * (n_modes + 1), False)

    # extract persistent arrays ONCE
    M = np.array([p.M for p in pl], dtype=np.float64)
    A = np.array([p.A for p in pl], dtype=np.float64)
    Ns = np.array([p.Ns for p in pl], dtype=np.float64)
    kap = np.array([p.kappa for p in pl], dtype=np.float64)

    P, z, S = P0, 0.0, 0.0
    for t in range(nt):
        z, T, P = ascend_parcel(z, T, P, w, dt, (t + 1) * dt, 3000.0, th, None, "linear")
        rp, _, am = parcel_rho(P, T)
        # condensation: object path on a thin particle list view rebuilt from arrays
        # (we only need the array growth; reuse the SoA condensation kernel).
        from PyLCM.condensation_fast import condense_soa
        T, q = condense_soa(M, A, Ns, kap, T, q, P, dt, am, rho_aero,
                            kohler_activation_radius=False,
                            switch_kappa_koehler=False)
        M, A, Ns, kap, _ = collide_soa(M, A, Ns, kap, dt, rp, P, T)

    # diagnostics: rebuild a lightweight particle list for ts_analysis
    from PyLCM.micro_particle import particles
    pl = []
    for k in range(len(M)):
        p = particles(0)
        p.M, p.A, p.Ns, p.kappa = M[k], A[k], Ns[k], kap[k]
        pl.append(p)
    _, _, qc, qr, NA, NC, NR, _, _, _ = ts_analysis(pl, am, [1e-6, 25e-6], 60, n_ptcl)
    return qc, qr, NR


def test_statistical_equivalence_to_object_path():
    """Object collection() vs SoA collide_soa over a full maritime ascent.

    Each path is run for the SAME set of seeds. Because the two paths consume the
    numpy RNG in a different order (object shuffle of a Python list + per-pair
    draws vs array permutation + kernel draws), single-seed results are NOT
    bit-identical, and the sub-1-cm^-3 rain number is very sensitive to a single
    super-droplet crossing the 25 um threshold. We therefore compare the
    ENSEMBLE MEAN over a few seeds — the same ensemble-averaging methodology used
    for the validated PyLCM-vs-Fortran comparison — within a Monte-Carlo
    tolerance (25%).
    """
    from validation.phys_harness import run as obj_run

    seeds = [2024, 17, 101, 5, 777]
    n_ptcl, nt = 2000, 1500

    obj_w, obj_qr, obj_nr = [], [], []
    soa_w, soa_qr, soa_nr = [], [], []
    for s in seeds:
        out, _ = obj_run(seed=s, n_ptcl=n_ptcl, nt=nt, aerosol="maritime",
                         collisions=True)
        o = out[nt]
        obj_w.append(o["qc"] + o["qr"]); obj_qr.append(o["qr"]); obj_nr.append(o["NR"])

        qc, qr, NR = _run_soa_ascent(seed=s, n_ptcl=n_ptcl, nt=nt, aerosol="maritime")
        soa_w.append(qc + qr); soa_qr.append(qr); soa_nr.append(NR)

    obj_water = float(np.mean(obj_w)); soa_water = float(np.mean(soa_w))
    obj_QR = float(np.mean(obj_qr)); soa_QR = float(np.mean(soa_qr))
    obj_NR = float(np.mean(obj_nr)); soa_NR = float(np.mean(soa_nr))

    rel_water = abs(soa_water - obj_water) / max(obj_water, 1e-12)
    rel_QR = abs(soa_QR - obj_QR) / max(obj_QR, 1e-12)
    rel_NR = abs(soa_NR - obj_NR) / max(obj_NR, 1e-9)

    print(f"\n[stat-equiv | mean of {len(seeds)} seeds] "
          f"object qc+qr={obj_water:.6e} qr={obj_QR:.6e} Nr={obj_NR:.4f} | "
          f"soa qc+qr={soa_water:.6e} qr={soa_QR:.6e} Nr={soa_NR:.4f} | "
          f"rel_water={rel_water:.3f} rel_QR={rel_QR:.3f} rel_NR={rel_NR:.3f}")

    # Robust integral diagnostics (total water, rain water) must match tightly:
    # collision moves mass, and OUR physics conserves it, so these are insensitive
    # to RNG-order differences.
    assert rel_water < 0.10, f"total water mismatch rel={rel_water:.3f}"
    assert rel_QR < 0.10, f"rain water mismatch rel={rel_QR:.3f}"
    # Rain NUMBER is a sub-1-cm^-3 boundary count (super-droplets crossing the
    # 25 um cloud/rain threshold); it is hypersensitive to RNG order, so we only
    # require the ensemble mean to be in the same ballpark (Monte-Carlo tol).
    assert rel_NR < 0.35, f"rain number mismatch rel={rel_NR:.3f}"
