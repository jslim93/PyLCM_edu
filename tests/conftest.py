import pytest
from PyLCM.micro_particle import particles


def make_particle(M, A, Ns=1e-18, kappa=0.5):
    p = particles(1)
    p.M, p.A, p.Ns, p.kappa = M, A, Ns, kappa
    return p
