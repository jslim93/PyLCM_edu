"""Tier-2 guard: the np.random.permutation shuffle keeps the full
condensation+collision run physically valid and reproducible-under-seed.

The shuffle change is stochastic (it alters the random sequence), so it cannot
be bit-exact against the old random.shuffle path. Instead we prove the
properties that MUST hold:

  (a) a small full run produces no NaN,
  (b) the collision step family conserves total water mass (sum of particle.M)
      to a tight relative tolerance every step (collisions only redistribute
      mass; with sedimentation removal off, no mass leaves the system), and
  (c) the run is reproducible: two runs under the same np.random.seed give
      identical particle state.
"""
import numpy as np

from PyLCM.parameters import r_a, cp, rv, rho_aero, rho_liq, p0
from PyLCM.parameters import z_env
from PyLCM.aero_init import aero_init
from PyLCM.parcel import ascend_parcel, parcel_rho
from PyLCM.condensation import esatw
from PyLCM.condensation_fast import drop_condensation_fast as drop_condensation
from PyLCM.collision import collection


def _run(seed, n_ptcl=500, nt=200, check_collision_mass=False):
    np.random.seed(seed)
    T0 = 293.2; P0 = 1013.0e2; RH = 0.92; w = 1.0; z0 = 0.0; zmax = 3000.0; dt = 1.0
    N_raw = np.array([118.0, 11.0, 0.72])
    mu_aero = np.log(np.array([0.019, 0.056, 0.46]) * 1e-6)
    sigma_aero = np.log(np.array([3.3, 1.6, 2.2]))
    kappa_arr = [1.6] * 4
    N_aero = N_raw * 1e6
    theta_init = T0 * (p0 / P0) ** (r_a / cp)
    theta_profiles = theta_init + 5e-3 * z_env
    e_s = esatw(T0)
    q0 = RH * e_s / (P0 - RH * e_s) * r_a / rv

    T_parcel, q_parcel, particles_list = aero_init(
        "Random", n_ptcl, P0, z0, T0, q0,
        N_aero, mu_aero, sigma_aero, rho_aero, kappa_arr, False)
    P_parcel = P0; z_parcel = z0; S_lst = 0.0

    max_collision_rel_err = 0.0
    for t in range(nt):
        tt = (t + 1) * dt
        z_parcel, T_parcel, P_parcel = ascend_parcel(
            z_parcel, T_parcel, P_parcel, w, dt, tt, zmax, theta_profiles, None, "linear")
        rho_p, _, air_mass = parcel_rho(P_parcel, T_parcel)
        ct = at = et = da = 0.0
        particles_list, T_parcel, q_parcel, S_lst, ct, at, et, da = drop_condensation(
            particles_list, T_parcel, q_parcel, P_parcel, nt, dt, air_mass, S_lst,
            rho_aero, False, ct, at, et, da, False)

        if check_collision_mass:
            m_before = sum(p.M for p in particles_list if p.A > 0)

        ac = au = pr = 0.0
        # sedi_removal defaults to False -> collisions only redistribute mass.
        particles_list, ac, au, pr = collection(
            dt, particles_list, rho_p, rho_liq, P_parcel, T_parcel, ac, au, pr,
            False, z_parcel, zmax, w)

        if check_collision_mass:
            m_after = sum(p.M for p in particles_list if p.A > 0)
            if m_before > 0:
                rel = abs(m_after - m_before) / m_before
                max_collision_rel_err = max(max_collision_rel_err, rel)

    # Snapshot of final state for reproducibility comparison.
    state = np.array([[p.M, p.A, p.Ns, p.kappa] for p in particles_list])
    return state, max_collision_rel_err


def test_full_run_no_nan():
    state, _ = _run(seed=12345)
    assert state.size > 0
    assert np.all(np.isfinite(state)), "full condensation+collision run produced NaN/Inf"


def test_collision_step_conserves_total_water_mass():
    _, max_rel_err = _run(seed=2024, check_collision_mass=True)
    # Collisions redistribute mass among super-droplets; with sedimentation
    # removal off, the per-step collision mass change must be ~0 (floating-point).
    assert max_rel_err < 1e-10, f"collision step changed total M by rel {max_rel_err:.2e}"


def test_full_run_reproducible_under_seed():
    state_a, _ = _run(seed=777)
    state_b, _ = _run(seed=777)
    assert state_a.shape == state_b.shape, "non-reproducible particle count under same seed"
    assert np.array_equal(state_a, state_b), "run not bit-identical under the same np.random.seed"
