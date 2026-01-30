import numpy as np
from matplotlib import pyplot as plt
import random

from PyLCM.parameters import *
from PyLCM.parameters import *
from PyLCM.parcel import *
from PyLCM.condensation import *

class particles:
    def __init__(self,n):
        self.id     = 0
        self.M      = 1.0 # mass
        self.A      = 1.0 # weighting factor
        self.Ns     = 1.0 # Aerosol mass
        self.kappa  = 0.5 # kappa parameter
        self.z      = 0.0 # particle vertical location
        # Ice-related properties
        self.micro_type = 1   # +1: liquid, +2: liquid with IN, -1: ice (hom), -2: ice (het/IN)
        self.aaxis  = 0.0     # a-axis half-diameter (m), ice only
        self.caxis  = 0.0     # c-axis half-diameter (m), ice only
        self.dns    = rho_liq # apparent density (kg/m3)
        self.IN_sfc = 0.0     # INP surface area (m^2), for ABIFM
        self.wsedi  = 0.0     # terminal velocity (m/s)

    def shuffle(particles_list):
        random.shuffle(particles_list)

    def is_ice(self):
        """Check if particle is ice"""
        return self.micro_type < 0

    def is_liquid(self):
        """Check if particle is liquid"""
        return self.micro_type > 0

    def has_IN(self):
        """Check if particle has ice nucleating particle"""
        return abs(self.micro_type) == 2