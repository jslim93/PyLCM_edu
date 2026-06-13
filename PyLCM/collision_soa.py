"""Struct-of-arrays (persistent numpy) collision path.

Drop-in array version of OUR validated object collision in PyLCM/collision.py.
It replicates the SAME physics (integer multiplicities, zero-radius ghost guard,
Straub E_S09 coalescence, multi-collision p_crit with the LSM upscaling and the
NINT + p_crit_lim cap, gravitational and Wang-Ayala turbulent kernels) but
operates on persistent float64 arrays inside an @njit kernel so no per-step
object<->array conversion is needed.

The object reference in PyLCM/collision.py stays the source of truth; this module
must stay statistically equivalent to it (see tests/test_collision_soa.py).

REUSE of validated helpers
--------------------------
E_H80 and ws_drops_beard are already @njit in collision.py and are imported
directly. E_S09 and E_turb are pure-numeric Python functions there; we njit-wrap
their exact .py_func source (no edits). The Ayala (2008) turbulent helpers
phi_w / zhi_func / gck reference module globals and each other, so njit-wrapping
the originals in-place fails; we instead @njit them here as VERBATIM copies of the
collision.py source (mechanical port, identical arithmetic) so gck can close over
the njit phi_w / zhi_func.
"""
import math
import numpy as np
from numba import njit

from PyLCM.parameters import (
    rho_liq, pi, g, muelq, PARCEL_AIR_MASS, seperation_radius_ts,
)
# Already-@njit, validated helpers — import and reuse directly.
from PyLCM.collision import E_H80, ws_drops_beard
from PyLCM.condensation import sigma_air_liq

# E_S09, E_turb and ws_drops_stokes are plain (pure-numeric) Python in
# collision.py; njit-wrap the EXACT source. E_S09 calls sigma_air_liq (already
# @njit). This is the same code object, just compiled — no physics change.
import PyLCM.collision as _coll
E_S09 = njit(cache=True)(_coll.E_S09)
E_turb = njit(cache=True)(_coll.E_turb)
ws_drops_stokes = njit(cache=True)(_coll.ws_drops_stokes)


# ---------------------------------------------------------------------------
# Wang-Ayala turbulent helpers — VERBATIM @njit copies of collision.py sources.
# (njit-wrapping the originals in-place fails because gck references the plain
# Python phi_w / zhi_func globals; copying them here lets gck close over the
# njit versions. Arithmetic is byte-for-byte identical to collision.py.)
# ---------------------------------------------------------------------------
@njit(cache=True)
def phi_w(a, b, vsett, tau0):
    aa1 = 1.0 / tau0 + 1.0 / a + vsett / b
    return 1.0 / aa1 - 0.5 * vsett / b / aa1**2


@njit(cache=True)
def zhi_func(a, b, vsett1, tau1, vsett2, tau2):
    aa1 = vsett2 / b - 1.0 / tau2 - 1.0 / a
    aa2 = vsett1 / b + 1.0 / tau1 + 1.0 / a
    aa3 = (vsett1 - vsett2) / b + 1.0 / tau1 + 1.0 / tau2
    aa4 = (vsett2 / b)**2 - (1.0 / tau2 + 1.0 / a)**2
    aa5 = vsett2 / b + 1.0 / tau2 + 1.0 / a
    aa6 = 1.0 / tau1 - 1.0 / a + (1.0 / tau2 + 1.0 / a) * vsett1 / vsett2

    result = ((1.0 / aa1 - 1.0 / aa2) * (vsett1 - vsett2) * 0.5 /
              b / aa3**2 +
              (4.0 / aa4 - 1.0 / aa5**2 - 1.0 / aa1**2) *
              vsett2 * 0.5 / b / aa6 +
              (2.0 * (b / aa2 - b / aa1) -
               vsett1 / aa2**2 + vsett2 / aa1**2) * 0.5 / b / aa3)
    return result


@njit(cache=True)
def gck(r1, r2, ws1, ws2, epsilon, tke):
    rho_dummy = 1.2  # reference air density for kinematic viscosity

    urms = math.sqrt(2.0 / 3.0 * tke)

    vis_kin = muelq / rho_dummy  # kinematic viscosity

    lam = urms * math.sqrt(15.0 * vis_kin / epsilon)       # Taylor microscale
    lambda_re = urms**2 * math.sqrt(15.0 / epsilon / vis_kin)  # Taylor-Re
    tl = urms**2 / epsilon
    lf = 0.5 * urms**3 / epsilon
    tauk = math.sqrt(vis_kin / epsilon)
    eta = (vis_kin**3 / epsilon)**0.25
    vk = eta / tauk

    ao = (11.0 + 7.0 * lambda_re) / (205.0 + lambda_re)
    tt = math.sqrt(2.0 * lambda_re / (math.sqrt(15.0) * ao)) * tauk

    tau1 = ws1 / g   # inertial time scale
    st1 = tau1 / tauk  # Stokes number
    tau2 = ws2 / g
    st2 = tau2 / tauk

    # Average radial relative velocity at contact (wrfin)
    z = tt / tl
    be = math.sqrt(2.0) * lam / lf
    bbb = math.sqrt(1.0 - 2.0 * be**2)
    d1 = (1.0 + bbb) / (2.0 * bbb)
    e1 = lf * (1.0 + bbb) * 0.5
    d2 = (1.0 - bbb) * 0.5 / bbb
    e2 = lf * (1.0 - bbb) * 0.5
    ccc = math.sqrt(1.0 - 2.0 * z**2)
    b1 = (1.0 + ccc) * 0.5 / ccc
    c1 = tl * (1.0 + ccc) * 0.5
    b2 = (1.0 - ccc) * 0.5 / ccc
    c2 = tl * (1.0 - ccc) * 0.5

    v1 = ws1
    t1 = tau1
    v2 = ws2
    t2 = tau2
    rrp = r1 + r2

    v1xysq = (b1 * d1 * phi_w(c1, e1, v1, t1) - b1 * d2 * phi_w(c1, e2, v1, t1)
              - b2 * d1 * phi_w(c2, e1, v1, t1) + b2 * d2 * phi_w(c2, e2, v1, t1))
    v1xysq = v1xysq * urms**2 / t1
    vrms1xy = math.sqrt(v1xysq)

    v2xysq = (b1 * d1 * phi_w(c1, e1, v2, t2) - b1 * d2 * phi_w(c1, e2, v2, t2)
              - b2 * d1 * phi_w(c2, e1, v2, t2) + b2 * d2 * phi_w(c2, e2, v2, t2))
    v2xysq = v2xysq * urms**2 / t2
    vrms2xy = math.sqrt(v2xysq)

    # Sort so v1 >= v2 for the cross-correlation term
    if ws1 >= ws2:
        v1, t1 = ws1, tau1
        v2, t2 = ws2, tau2
    else:
        v1, t1 = ws2, tau2
        v2, t2 = ws1, tau1

    v1v2xy = (b1 * d1 * zhi_func(c1, e1, v1, t1, v2, t2)
              - b1 * d2 * zhi_func(c1, e2, v1, t1, v2, t2)
              - b2 * d1 * zhi_func(c2, e1, v1, t1, v2, t2)
              + b2 * d2 * zhi_func(c2, e2, v1, t1, v2, t2))
    fr = d1 * math.exp(-rrp / e1) - d2 * math.exp(-rrp / e2)
    v1v2xy = v1v2xy * fr * urms**2 / tau1 / tau2
    wrtur2xy = vrms1xy**2 + vrms2xy**2 - 2.0 * v1v2xy
    if wrtur2xy < 0.0:
        wrtur2xy = 0.0
    wrgrav2 = pi / 8.0 * (ws2 - ws1)**2
    wrfin = math.sqrt((2.0 / pi) * (wrtur2xy + wrgrav2))

    # Radial distribution function (grfin)
    sst = max(st1, st2)

    xx = -0.1988 * sst**4 + 1.5275 * sst**3 - 4.2942 * sst**2 + 5.3406 * sst
    if xx < 0.0:
        xx = 0.0
    yy = 0.1886 * math.exp(20.306 / lambda_re)

    c1_gr = xx / (g / vk * tauk)**yy

    ao_gr = ao + (pi / 8.0) * (g / vk * tauk)**2
    fao_gr = 20.115 * math.sqrt(ao_gr / lambda_re)
    rc = math.sqrt(fao_gr * abs(st2 - st1)) * eta

    grfin = ((eta**2 + rc**2) / (rrp**2 + rc**2))**(c1_gr * 0.5)
    if grfin < 1.0:
        grfin = 1.0

    return 2.0 * pi * rrp**2 * wrfin * grfin


# Critical mass for accretion/autoconversion classification (matches collision.py).
_MASS_CRIT = (seperation_radius_ts ** 3) * 4.0 / 3.0 * pi * rho_liq
# Larger drop must exceed 10 um to collide (matches collision.py skip guard).
_MASS_10UM = (10.0e-6 ** 3) * 4.0 / 3.0 * pi * rho_liq


@njit(cache=True)
def _collision_kernel(M, A, Ns, kappa, idx, half, nptcl,
                      dt, rho_parcel, p_env, T_parcel,
                      switch_E_constant, switch_vt_simple,
                      switch_turb_kernel, epsilon_turb):
    """Serial in-place LSM collision over the shuffled index array `idx`.

    Replicates determine_collision() + same_weights_update()/liquid_update_collection()
    from PyLCM/collision.py EXACTLY, on arrays. Pairs idx[k] with idx[k+half].
    Returns n_collisions (number of accepted pairs).
    """
    V_parcel = PARCEL_AIR_MASS / rho_parcel
    n_coll = 0

    for k in range(half):
        i = idx[k]
        j = idx[k + half]

        Ai = A[i]
        Aj = A[j]
        Mi = M[i]
        Mj = M[j]

        # --- skip guards (collection() in collision.py) ---
        # at least one real particle to collect
        if min(Ai, Aj) <= 0.0:
            continue
        # a single-droplet super-droplet has nothing to give away
        if max(Ai, Aj) <= 1.0:
            continue

        # ghost guard (determine_collision): degenerate super-droplet, zero radius
        if Mi <= 0.0 or Mj <= 0.0 or Ai <= 0.0 or Aj <= 0.0:
            continue

        # individual (per-droplet) masses
        xi = Mi / Ai
        xj = Mj / Aj

        # larger droplet must exceed 10 um
        if max(xi, xj) < _MASS_10UM:
            continue

        # radii from per-droplet mass
        R_n = (xi / (4.0 / 3.0 * pi * rho_liq)) ** 0.33333333333
        R_m = (xj / (4.0 / 3.0 * pi * rho_liq)) ** 0.33333333333

        # terminal velocities: Beard (1976) or simplified Stokes
        if switch_vt_simple:
            v_r1 = ws_drops_stokes(R_n, rho_parcel, rho_liq)
            v_r2 = ws_drops_stokes(R_m, rho_parcel, rho_liq)
        else:
            v_r1 = ws_drops_beard(R_n, rho_parcel, rho_liq, p_env, T_parcel)
            v_r2 = ws_drops_beard(R_m, rho_parcel, rho_liq, p_env, T_parcel)

        v_r = abs(v_r1 - v_r2)

        # collision efficiency: Hall (1980) or constant E=1
        if switch_E_constant:
            E_coll = 1.0
        else:
            E_coll = E_H80(R_m, R_n)

        # collection kernel: gravitational or Wang-Ayala turbulent
        if switch_turb_kernel and epsilon_turb > 1.0e-10:
            urms_est = 2.02 * (epsilon_turb / 0.04)**(1.0 / 3.0)
            tke_est = 1.5 * urms_est**2
            K = (E_coll * gck(R_n, R_m, v_r1, v_r2, epsilon_turb, tke_est)
                 * E_turb(R_n, R_m, epsilon_turb)
                 * E_S09(R_m, R_n, v_r, rho_liq, T_parcel))
        else:
            K = (pi * (R_m + R_n) ** 2 * v_r * E_coll
                 * E_S09(R_m, R_n, v_r, rho_liq, T_parcel))

        p_crit = max(Ai, Aj) * K / V_parcel * dt
        p_crit = p_crit * nptcl * (nptcl - 1) / (half * 2)

        x_rand = np.random.random()

        if p_crit > x_rand:
            # resolve number of collisions (multi-collision, Shima et al. 2009)
            if p_crit <= 1.0:
                p_int = 1
            else:
                p_int = max(int(round(p_crit)), 1)  # NINT equivalent
                A_max = max(Ai, Aj)
                A_min = min(Ai, Aj)
                p_crit_lim = max(int((A_max - 1) / A_min), 1)
                p_int = min(p_int, p_crit_lim)

            # ---- update ----
            if Ai == Aj:
                # same weighting factor: SAM floor split
                _same_weights_update(M, A, Ns, kappa, i, j)
            else:
                _liquid_update_collection(M, A, Ns, kappa, i, j, p_int)

            n_coll += 1

    return n_coll


@njit(cache=True)
def _liquid_update_collection(M, A, Ns, kappa, i, j, p_crit):
    """Asymmetric update: smaller-A gains, larger-A loses (collision.py).

    Reconstructs M = x_int * A_new for the loser (no subtraction) and keeps A
    integer-valued. Operates in place on entries i, j.
    """
    # int1 = smaller A (gains), int2 = larger A (loses)
    if A[i] < A[j]:
        a = i
        b = j
    else:
        a = j
        b = i

    x_int = M[b] / A[b]
    xs_int = Ns[b] / A[b]

    v_a = M[a] / A[a] / rho_liq
    v_b = M[b] / A[b] / rho_liq

    A_a = A[a]
    # smaller-A super-droplet gains p_crit collisions worth of mass
    M[a] = M[a] + A_a * x_int * p_crit
    Ns[a] = Ns[a] + A_a * xs_int * p_crit
    kappa[a] = (v_a * kappa[a] + v_b * kappa[b]) / (v_a + v_b)

    # larger-A super-droplet loses number; reconstruct M, Ns proportional to A
    A[b] = A[b] - A_a * p_crit
    M[b] = x_int * A[b]
    Ns[b] = xs_int * A[b]


@njit(cache=True)
def _same_weights_update(M, A, Ns, kappa, i, j):
    """Equal-A update: SAM integer floor split A//2 with <1 dissolve branch."""
    xn = M[i] / A[i]
    xm = M[j] / A[j]
    xsn = Ns[i] / A[i]
    xsm = Ns[j] / A[j]

    v_i = M[i] / A[i] / rho_liq
    v_j = M[j] / A[j] / rho_liq

    A_total = A[i]
    A_half = float(int(A_total) // 2)

    new_kappa = (v_i * kappa[i] + v_j * kappa[j]) / (v_i + v_j)

    if A_half < 1.0:
        # too few to split: merge all into i, dissolve j
        A[i] = A_total
        M[i] = (xn + xm) * A_total
        Ns[i] = (xsn + xsm) * A_total
        A[j] = 0.0
        M[j] = 0.0
        Ns[j] = 0.0
    else:
        A[i] = A_half
        A[j] = A_total - A_half
        M[i] = (xn + xm) * A[i]
        M[j] = (xn + xm) * A[j]
        Ns[i] = (xsn + xsm) * A[i]
        Ns[j] = (xsn + xsm) * A[j]

    kappa[i] = new_kappa
    kappa[j] = new_kappa


def collide_soa(M, A, Ns, kappa, dt, rho_parcel, p_env, T_parcel,
                switch_E_constant=False, switch_vt_simple=False,
                switch_turb_kernel=False, epsilon_turb=0.0):
    """Array (SoA) collision-coalescence over persistent float64 arrays.

    Mutates M, A, Ns, kappa in place (A holds integer values as float64), then
    compacts out the A<=0 entries and RETURNS the compacted arrays plus the
    number of removed super-droplets.

    Parameters mirror determine_collision()/collection() in PyLCM/collision.py.

    Returns
    -------
    M, A, Ns, kappa : np.ndarray
        Compacted arrays (A>0 entries only).
    n_removed : int
        Number of super-droplets dropped (A<=0) during compaction.
    """
    n = len(M)
    if n < 2:
        return M, A, Ns, kappa, 0

    # LSM: shuffle, pair k with k+half (half = n//2). The upscaling
    # nptcl*(nptcl-1)/(half*2) is applied inside the kernel, exactly as in
    # determine_collision.
    idx = np.random.permutation(n)
    half = n // 2

    _collision_kernel(M, A, Ns, kappa, idx, half, n,
                      float(dt), float(rho_parcel), float(p_env), float(T_parcel),
                      bool(switch_E_constant), bool(switch_vt_simple),
                      bool(switch_turb_kernel), float(epsilon_turb))

    # compact: drop A<=0 (degenerate / dissolved super-droplets)
    keep = A > 0.0
    n_removed = int(n - np.count_nonzero(keep))
    if n_removed > 0:
        M = np.ascontiguousarray(M[keep])
        A = np.ascontiguousarray(A[keep])
        Ns = np.ascontiguousarray(Ns[keep])
        kappa = np.ascontiguousarray(kappa[keep])

    return M, A, Ns, kappa, n_removed
