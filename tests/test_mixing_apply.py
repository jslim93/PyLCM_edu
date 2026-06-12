import numpy as np
import matplotlib
matplotlib.use("Agg")
from PyLCM.micro_particle import particles
from PyLCM.parcel import create_env_profiles
from PyLCM.mixing import ParameterizedMixing


def _setup():
    qv, th, z_env = create_env_profiles(290.0, 0.010, 0.0, 95000.0, "Stable")
    ps = []
    for _ in range(100):
        p = particles(1); p.M, p.A, p.Ns, p.kappa = 1.0e-9, 50.0, 1e-18, 0.5
        ps.append(p)
    return qv, th, z_env, ps


def test_redistribution_adds_no_net_water():
    # Entrainment legitimately exchanges water with the environment; the IHMD
    # redistribution (liquid -> vapor) must itself be conservative. So the TOTAL
    # water change must equal EXACTLY the bulk entrainment exchange and nothing more.
    qv, th, z_env, ps = _setup()
    air_mass = 100.0
    lam, w, dt = 5e-4, 1.0, 1.0
    mix = ParameterizedMixing(lambda_ent=lam, ihmd=0.5, qv_profiles=qv,
                              theta_profiles=th, z_env=z_env)
    T, q = 288.0, 0.009
    total0 = q * air_mass + sum(p.M for p in ps)
    q_env = float(np.interp(500.0, z_env, qv))
    expected_change = (lam * w * dt) * (q_env - q) * air_mass   # bulk exchange only
    ps, T1, q1 = mix.apply(ps, T, q, 90000.0, 500.0, dt=dt, w=w, air_mass=air_mass)
    total1 = q1 * air_mass + sum(p.M for p in ps)
    assert np.isclose(total1 - total0, expected_change, atol=1e-12)
    assert np.isfinite(T1) and np.isfinite(q1)


def test_disabled_is_noop():
    qv, th, z_env, ps = _setup()
    mix = ParameterizedMixing(lambda_ent=0.0, ihmd=0.5, qv_profiles=qv,
                              theta_profiles=th, z_env=z_env)
    M_before = [p.M for p in ps]
    ps, T1, q1 = mix.apply(ps, 288.0, 0.009, 90000.0, 500.0, dt=1.0, w=1.0, air_mass=100.0)
    assert [p.M for p in ps] == M_before and (T1, q1) == (288.0, 0.009)


def test_apply_is_deterministic():
    qv, th, z_env, ps1 = _setup()
    _, _, _, ps2 = _setup()
    m1 = ParameterizedMixing(8e-4, 0.7, qv, th, z_env)
    m2 = ParameterizedMixing(8e-4, 0.7, qv, th, z_env)
    r1 = m1.apply(ps1, 288.0, 0.009, 90000.0, 500.0, 1.0, 1.0, 100.0)
    r2 = m2.apply(ps2, 288.0, 0.009, 90000.0, 500.0, 1.0, 1.0, 100.0)
    assert [p.M for p in r1[0]] == [p.M for p in r2[0]] and r1[1:] == r2[1:]
