import numpy as np
from matplotlib import pyplot as plt
from parameters import *
import random

class particles:
    def __init__(self,n):
        self.id     = n
        self.M      = 1.0
        self.A      = 1.0 #weighting factor
        self.r      = 1.0
    def __shuffle__(self):
        random.shuffle(self.particles)
        