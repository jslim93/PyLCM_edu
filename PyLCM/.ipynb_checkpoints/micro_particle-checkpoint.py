import numpy as np
from matplotlib import pyplot as plt
import random

from PyLCM.parameters import *
from PyLCM.parameters import *
from PyLCM.parcel import *
from PyLCM.condensation import *

class particles:
    def __init__(self,n):
        self.id     = n
        self.M      = 1.0 # mass
        self.A      = 1.0 # weighting factor
        self.Ns     = 1.0 # Aerosol mass
    def shuffle(particles_list):
        random.shuffle(particles_list)