import numpy as np
import matplotlib; matplotlib.use("Agg")
from validation.phys_harness import run


def test_turbulent_collision_makes_no_ghosts():
    # The config that previously crashed at t~971 with a zero-radius ghost.
    _, pl = run(seed=0, aerosol="maritime", n_ptcl=2000, nt=1800,
                collisions=True, switch_turb=True, eps=0.04)
    # No super-droplet may have zero/negative liquid mass while keeping A>0,
    # and every surviving A must be a positive integer.
    for p in pl:
        assert p.A == round(p.A)
        if p.A > 0:
            assert p.M > 0.0
