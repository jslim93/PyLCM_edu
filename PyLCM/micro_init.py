import numpy as np
from matplotlib import pyplot as plt
from PyLCM.parameters import *
import random

class particles:
    def __init__(self,n):
        self.id     = n
        self.M      = 1.0 # mass
        self.A      = 1.0 # weighting factor
    def shuffle(particles_list):
        random.shuffle(particles_list)