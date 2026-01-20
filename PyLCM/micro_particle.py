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
    def shuffle(particles_list):
        random.shuffle(particles_list)


# ============================================================================
# NEW: Array-Based Particle Representation for Performance
# ============================================================================

class ParticleArrays:
    """
    High-performance array-based particle representation.
    
    This class stores particle properties in numpy arrays instead of Python objects,
    enabling vectorized operations and 3-5x speedup on large particle ensembles.
    
    Attributes:
        n (int): Number of particles
        M (np.ndarray): Mass of each particle [kg]
        A (np.ndarray): Weighting factor (multiplicity)
        Ns (np.ndarray): Aerosol mass [kg]
        kappa (np.ndarray): Hygroscopicity parameter
        z (np.ndarray): Particle vertical location [m]
        id (np.ndarray): Particle ID
    """
    
    def __init__(self, n_particles):
        """
        Initialize particle arrays.
        
        Args:
            n_particles (int): Number of particles to create
        """
        self.n = n_particles
        self.M = np.ones(n_particles, dtype=np.float64)
        self.A = np.ones(n_particles, dtype=np.float64)
        self.Ns = np.ones(n_particles, dtype=np.float64)
        self.kappa = np.full(n_particles, 0.5, dtype=np.float64)
        self.z = np.zeros(n_particles, dtype=np.float64)
        self.id = np.arange(n_particles, dtype=np.int32)
    
    @classmethod
    def from_particle_list(cls, particles_list):
        """
        Convert a list of particle objects to ParticleArrays.
        
        Args:
            particles_list (list): List of particle objects
            
        Returns:
            ParticleArrays: Array-based representation
        """
        n = len(particles_list)
        arrays = cls(n)
        
        for i, p in enumerate(particles_list):
            arrays.M[i] = p.M
            arrays.A[i] = p.A
            arrays.Ns[i] = p.Ns
            arrays.kappa[i] = p.kappa
            arrays.z[i] = p.z
            arrays.id[i] = p.id
        
        return arrays
    
    def to_particle_list(self):
        """
        Convert ParticleArrays back to list of particle objects.
        Useful for backward compatibility.
        
        Returns:
            list: List of particle objects
        """
        particles_list = []
        for i in range(self.n):
            p = particles(i)
            p.M = self.M[i]
            p.A = self.A[i]
            p.Ns = self.Ns[i]
            p.kappa = self.kappa[i]
            p.z = self.z[i]
            p.id = self.id[i]
            particles_list.append(p)
        
        return particles_list
    
    def shuffle(self):
        """Shuffle particles (in-place)."""
        indices = np.random.permutation(self.n)
        self.M = self.M[indices]
        self.A = self.A[indices]
        self.Ns = self.Ns[indices]
        self.kappa = self.kappa[indices]
        self.z = self.z[indices]
        self.id = self.id[indices]
    
    def copy(self):
        """Create a deep copy of the particle arrays."""
        new = ParticleArrays(self.n)
        new.M = self.M.copy()
        new.A = self.A.copy()
        new.Ns = self.Ns.copy()
        new.kappa = self.kappa.copy()
        new.z = self.z.copy()
        new.id = self.id.copy()
        return new
    
    def get_radii(self):
        """
        Calculate particle radii from mass.
        
        Returns:
            np.ndarray: Particle radii [m]
        """
        # r = (3M / (4π ρ A))^(1/3)
        return (3.0 * self.M / (4.0 * np.pi * rho_liq * self.A)) ** (1.0/3.0)
    
    def set_radii(self, radii):
        """
        Set particle mass from radii.
        
        Args:
            radii (np.ndarray): Particle radii [m]
        """
        # M = (4/3) π r³ ρ A
        self.M = (4.0/3.0) * np.pi * radii**3 * rho_liq * self.A
    
    def __len__(self):
        """Return number of particles."""
        return self.n

    def __repr__(self):
        return f"ParticleArrays(n={self.n}, total_mass={np.sum(self.M):.3e} kg)"

    def compact(self):
        """
        Remove particles with A <= 0 (in-place).

        This method compacts the arrays by removing invalid particles,
        which is essential after collision events that may set A to 0.

        Returns:
            int: Number of particles removed
        """
        valid = self.A > 0
        n_valid = np.sum(valid)
        n_removed = self.n - n_valid

        if n_removed > 0:
            self.M = self.M[valid]
            self.A = self.A[valid]
            self.Ns = self.Ns[valid]
            self.kappa = self.kappa[valid]
            self.z = self.z[valid]
            self.id = self.id[valid]
            self.n = n_valid

        return n_removed

    def remove_particles(self, mask):
        """
        Remove particles indicated by boolean mask (in-place).

        Args:
            mask (np.ndarray): Boolean array where True indicates particles to REMOVE

        Returns:
            int: Number of particles removed
        """
        keep = ~mask
        n_keep = np.sum(keep)
        n_removed = self.n - n_keep

        if n_removed > 0:
            self.M = self.M[keep]
            self.A = self.A[keep]
            self.Ns = self.Ns[keep]
            self.kappa = self.kappa[keep]
            self.z = self.z[keep]
            self.id = self.id[keep]
            self.n = n_keep

        return n_removed

    def add_particles(self, other):
        """
        Add particles from another ParticleArrays instance.

        Args:
            other (ParticleArrays): Particles to add
        """
        self.M = np.concatenate([self.M, other.M])
        self.A = np.concatenate([self.A, other.A])
        self.Ns = np.concatenate([self.Ns, other.Ns])
        self.kappa = np.concatenate([self.kappa, other.kappa])
        self.z = np.concatenate([self.z, other.z])
        self.id = np.concatenate([self.id, other.id])
        self.n = len(self.M)

    def get_total_mass(self):
        """Return total liquid water mass (summed over all particles)."""
        return np.sum(self.M)

    def get_total_number(self):
        """Return total number concentration (sum of all multiplicities)."""
        return np.sum(self.A)

    def get_aerosol_radii(self, rho_aero=2170.0):
        """
        Calculate aerosol core radii from aerosol mass.

        Args:
            rho_aero (float): Aerosol density [kg/m³]

        Returns:
            np.ndarray: Aerosol core radii [m]
        """
        return (self.Ns / (self.A * 4.0/3.0 * np.pi * rho_aero)) ** (1.0/3.0)