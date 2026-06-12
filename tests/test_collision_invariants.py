from PyLCM.collision import liquid_update_collection, same_weights_update
from tests.conftest import make_particle


def test_collision_conserves_total_water_mass():
    p1 = make_particle(M=2.0e-9, A=100.0)
    p2 = make_particle(M=5.0e-9, A=10.0)
    total_before = p1.M + p2.M
    a1, a2, acc, aut = liquid_update_collection(p1, p2, 0.0, 0.0, p_crit=1)
    assert abs((a1.M + a2.M) - total_before) < 1e-20


def test_collision_decreases_total_number():
    p1 = make_particle(M=2.0e-9, A=100.0)
    p2 = make_particle(M=5.0e-9, A=10.0)
    a_before = p1.A + p2.A
    a1, a2, _, _ = liquid_update_collection(p1, p2, 0.0, 0.0, p_crit=1)
    assert (a1.A + a2.A) < a_before


def test_same_weight_collision_conserves_mass():
    p1 = make_particle(M=3.0e-9, A=50.0)
    p2 = make_particle(M=4.0e-9, A=50.0)
    total_before = p1.M + p2.M
    a1, a2, _, _ = same_weights_update(p1, p2, 0.0, 0.0)
    assert abs((a1.M + a2.M) - total_before) < 1e-20
