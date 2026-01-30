"""
Profiling script for PyLCM_edu to identify performance bottlenecks.
Run with: python profile_model.py
"""

import cProfile
import pstats
import io
from pstats import SortKey

# Mock widget class
class MockWidget:
    def __init__(self, value):
        self.value = value

# Import after class definition
from PyLCM.parameters import *
from PyLCM.timestep_routine import timesteps_function
from PyLCM.parcel import create_env_profiles
from PyLCM.condensation import esatw

# Configuration for profiling
n_particles = 10000  # Test with 10,000 particles
P_parcel = 101300.0
RH_parcel = 0.88
T_parcel = 293.2
w_parcel = 1.0
nt = 100  # Reduced timesteps for faster profiling
dt = 1.0
ascending_mode = 'linear'
z_parcel = 0.0
max_z = 3000.0
do_condensation = True
do_collision = False
aero_mode = 'Weighting_factor'
use_koehler = False
use_kappa = False
do_sedi_removal = False
do_entrainment = False

# Aerosol parameters (simplified single mode)
class MockGridWidget:
    def __init__(self):
        self.widgets = {}
        def create_mock(val):
            return type('obj', (object,), {'value': val})
        
        # Mode 1: Active
        self.widgets[(1, 0)] = create_mock(118.0)    # N
        self.widgets[(2, 0)] = create_mock(0.019)    # r (µm)
        self.widgets[(3, 0)] = create_mock(3.3)      # sigma
        self.widgets[(4, 0)] = create_mock(1.6)      # kappa
        
        # Modes 2-4: Inactive
        for i in range(1, 4):
            self.widgets[(1, i)] = create_mock(0.0)
            self.widgets[(2, i)] = create_mock(0.0)
            self.widgets[(3, i)] = create_mock(0.0)
            self.widgets[(4, i)] = create_mock(1.6)
    
    def __getitem__(self, key):
        return self.widgets.get(key, type('obj', (object,), {'value': 0.0}))

gridwidget = MockGridWidget()

# Create environment profiles
qv_init = RH_parcel * esatw(T_parcel) / (P_parcel - RH_parcel * esatw(T_parcel)) * r_a / rv
qv_profiles, theta_profiles, z_env = create_env_profiles(T_parcel, qv_init, z_parcel, P_parcel, 'Stable', plot_profiles=False)

print("Starting profiling...")
print(f"Configuration: {n_particles} particles, {nt} timesteps, condensation={do_condensation}, collision={do_collision}")

# Create profiler
profiler = cProfile.Profile()

# Profile the main function
profiler.enable()

try:
    results = timesteps_function(
        MockWidget(n_particles), MockWidget(P_parcel), MockWidget(RH_parcel), 
        MockWidget(T_parcel), MockWidget(w_parcel), MockWidget(nt), 
        MockWidget(dt), rm_spec, MockWidget(ascending_mode), 
        MockWidget('text_fast'),
        MockWidget(z_parcel), MockWidget(max_z), 
        MockWidget(do_condensation), MockWidget(do_collision), 
        MockWidget(aero_mode), gridwidget, 
        use_koehler, use_kappa, MockWidget(do_sedi_removal),
        0.0, do_entrainment, qv_profiles, theta_profiles, 0, 0
    )
    print("\nSimulation completed successfully!")
    
except Exception as e:
    print(f"\nError during simulation: {e}")

profiler.disable()

# Analyze results
print("\n" + "="*80)
print("PROFILING RESULTS")
print("="*80)

s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(profiler, stream=s).sort_stats(sortby)

print("\nTop 20 functions by cumulative time:")
print("-" * 80)
ps.print_stats(20)
print(s.getvalue())

# Save detailed stats to file
with open('profile_results.txt', 'w') as f:
    ps = pstats.Stats(profiler, stream=f).sort_stats(sortby)
    f.write("="*80 + "\n")
    f.write("DETAILED PROFILING RESULTS\n")
    f.write("="*80 + "\n\n")
    f.write(f"Configuration: {n_particles} particles, {nt} timesteps\n")
    f.write(f"Condensation: {do_condensation}, Collision: {do_collision}\n\n")
    ps.print_stats()

print("\nDetailed results saved to profile_results.txt")
