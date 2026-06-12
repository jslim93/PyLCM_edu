import numpy as np
from PyLCM.micro_particle import particles
from PyLCM.mixing import redistribute_droplets


def _cloud(n=50, M=2.0e-9, A=100.0):
    out = []
    for _ in range(n):
        p = particles(1); p.M, p.A, p.Ns, p.kappa = M, A, 1e-18, 0.5
        out.append(p)
    return out


def test_nq_power_law_holds_exactly():
    for ihmd in (0.0, 0.5, 1.0):
        ps = _cloud()
        N0 = sum(p.A for p in ps); q0 = sum(p.M for p in ps)
        evap = redistribute_droplets(ps, ihmd=ihmd, frac=0.3)
        N1 = sum(p.A for p in ps); q1 = sum(p.M for p in ps)
        assert abs(evap - 0.3 * q0) < 1e-18                      # evaporated mass
        assert np.isclose(N1 / N0, (q1 / q0) ** ihmd, rtol=1e-9)  # the IHMD law


def test_homogeneous_conserves_number():
    ps = _cloud(); N0 = sum(p.A for p in ps)
    redistribute_droplets(ps, ihmd=0.0, frac=0.4)
    assert abs(sum(p.A for p in ps) - N0) < 1e-9                 # N unchanged
    assert all(np.isclose(p.M / p.A, (2.0e-9 / 100.0) * 0.6) for p in ps)  # all shrink


def test_inhomogeneous_preserves_droplet_size():
    ps = _cloud(); N0 = sum(p.A for p in ps)
    redistribute_droplets(ps, ihmd=1.0, frac=0.4)
    assert sum(p.A for p in ps) < N0                            # number drops
    assert all(np.isclose(p.M / p.A, 2.0e-9 / 100.0) for p in ps)  # size unchanged
