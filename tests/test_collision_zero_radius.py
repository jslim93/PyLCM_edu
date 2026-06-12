from PyLCM.collision import determine_collision
from PyLCM.micro_particle import particles


def _p(M, A):
    p = particles(1); p.M, p.A, p.Ns, p.kappa = M, A, 1e-18, 0.5
    return p


def test_fully_evaporated_droplet_does_not_crash_collision():
    # A super-droplet that fully evaporated (M=0) has zero radius. Previously this
    # caused a ZeroDivisionError in E_H80 (min(r1/r2, r2/r1)); it must now be skipped.
    p1 = _p(0.0, 100.0)            # fully evaporated
    p2 = _p(5.0e-9, 10.0)
    res = determine_collision(1.0, p1, p2, 1.0, 1000.0, 9.0e4, 285.0, 50, 100)
    assert res[4] == 0            # p_crit == 0: no collision, no crash


def test_zero_multiplicity_droplet_does_not_crash_collision():
    p1 = _p(2.0e-9, 0.0)          # zero multiplicity
    p2 = _p(5.0e-9, 10.0)
    res = determine_collision(1.0, p1, p2, 1.0, 1000.0, 9.0e4, 285.0, 50, 100)
    assert res[4] == 0
