import streamlit as st
import numpy as np
import pandas as pd
import time
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Page config must be the first Streamlit command
st.set_page_config(
    page_title="PyLCM Educational Model",
    page_icon="☁️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Import model modules
# Try to use optimized modules if available
try:
    from PyLCM import collision_optimized as collision
    from PyLCM import condensation_optimized as condensation
    from PyLCM.timestep_routine_arrays import timesteps_function_arrays as timesteps_function
    st.sidebar.success("✨ Using OPTIMIZED array-based physics (15-30x faster!)")
except ImportError:
    try:
        from PyLCM import collision_optimized as collision
        from PyLCM import condensation_optimized as condensation
        from PyLCM.timestep_routine import timesteps_function
        st.sidebar.warning("Using optimized physics modules (Numba JIT), but not array-based.")
    except ImportError:
        from PyLCM import collision, condensation
        from PyLCM.timestep_routine import timesteps_function
        st.sidebar.warning("Using standard physics modules (Optimizations not found)")

from PyLCM.parameters import *
from PyLCM.micro_particle import *
from PyLCM.aero_init import *
from PyLCM.parcel import *
from PyLCM.parcel import create_env_profiles

# Title and Description
st.title("☁️ PyLCM Educational Microphysics Model")
st.markdown("""
This is an interactive interface for the Lagrangian Cloud Model (LCM) parcel model. 
Adjust the parameters in the sidebar to simulate cloud microphysics processes.
""")

# --- Sidebar Configuration ---
st.sidebar.header("Model Configuration")

# Preset Scenarios
st.sidebar.subheader("Quick Start Presets")
preset = st.sidebar.selectbox(
    "Load a preset scenario",
    options=["Custom", "Cumulus Cloud", "Stratus/Fog", "Clean Maritime", "Polluted Continental"],
    help="Select a preset to auto-fill parameters for common cloud types"
)

# Define preset configurations
PRESETS = {
    "Cumulus Cloud": {
        "T": 293.2, "P": 101300.0, "RH": 0.88, "w": 2.0, "z": 0.0, "max_z": 4000.0,
        "n_particles": 2000, "ascending_mode": "linear",
        "aero": [
            {'n': 300.0, 'r': 0.02, 'sigma': 2.5, 'kappa': 0.6},
            {'n': 50.0, 'r': 0.08, 'sigma': 1.8, 'kappa': 0.6},
            {'n': 1.0, 'r': 0.5, 'sigma': 2.0, 'kappa': 0.6},
            {'n': 0.0, 'r': 0.0, 'sigma': 0.0, 'kappa': 0.6}
        ],
        "desc": "Strong updraft with continental aerosol - typical cumulus development"
    },
    "Stratus/Fog": {
        "T": 283.0, "P": 101300.0, "RH": 0.95, "w": 0.3, "z": 0.0, "max_z": 1000.0,
        "n_particles": 1000, "ascending_mode": "linear",
        "aero": [
            {'n': 100.0, 'r': 0.03, 'sigma': 2.0, 'kappa': 1.0},
            {'n': 20.0, 'r': 0.1, 'sigma': 1.5, 'kappa': 1.0},
            {'n': 0.5, 'r': 0.4, 'sigma': 2.0, 'kappa': 1.0},
            {'n': 0.0, 'r': 0.0, 'sigma': 0.0, 'kappa': 1.0}
        ],
        "desc": "Weak updraft, high RH - stratus or fog formation"
    },
    "Clean Maritime": {
        "T": 295.0, "P": 101300.0, "RH": 0.85, "w": 1.0, "z": 0.0, "max_z": 3000.0,
        "n_particles": 1000, "ascending_mode": "linear",
        "aero": [
            {'n': 40.0, 'r': 0.04, 'sigma': 2.0, 'kappa': 1.2},
            {'n': 10.0, 'r': 0.15, 'sigma': 1.8, 'kappa': 1.2},
            {'n': 0.3, 'r': 0.8, 'sigma': 2.2, 'kappa': 1.2},
            {'n': 0.0, 'r': 0.0, 'sigma': 0.0, 'kappa': 1.2}
        ],
        "desc": "Low aerosol concentration - fewer but larger droplets, faster rain"
    },
    "Polluted Continental": {
        "T": 293.2, "P": 101300.0, "RH": 0.80, "w": 1.5, "z": 0.0, "max_z": 3500.0,
        "n_particles": 3000, "ascending_mode": "linear",
        "aero": [
            {'n': 1000.0, 'r': 0.015, 'sigma': 2.5, 'kappa': 0.3},
            {'n': 200.0, 'r': 0.05, 'sigma': 2.0, 'kappa': 0.3},
            {'n': 5.0, 'r': 0.3, 'sigma': 2.0, 'kappa': 0.3},
            {'n': 0.0, 'r': 0.0, 'sigma': 0.0, 'kappa': 0.3}
        ],
        "desc": "High aerosol load - many small droplets, rain suppression"
    }
}

if preset != "Custom":
    st.sidebar.success(f"📋 {PRESETS[preset]['desc']}")

# 1. Model Steering
st.sidebar.subheader("Time Integration")
dt = st.sidebar.number_input(
    "Time step (dt) [s]",
    min_value=0.001, max_value=5.0, value=1.0, step=0.1,
    help="Integration time step. Smaller values are more accurate but slower. 1s works well for most cases."
)
nt = st.sidebar.number_input(
    "Number of steps (nt)",
    min_value=100, max_value=10000, value=3600, step=100,
    help="Total number of time steps. Total simulation time = dt × nt"
)
total_time = dt * nt
st.sidebar.info(f"Total simulation time: {total_time:.1f} s ({total_time/60:.1f} min)")

# 2. Physics Toggles
st.sidebar.subheader("Physics Processes")
do_condensation = st.sidebar.checkbox(
    "Condensation", value=True,
    help="Enable diffusional growth/evaporation of droplets. Essential for cloud formation."
)
do_collision = st.sidebar.checkbox(
    "Collision-Coalescence", value=False,
    help="Enable droplet collisions and coalescence. Required for rain formation. Computationally expensive."
)
do_sedi_removal = st.sidebar.checkbox(
    "Sedimentation Removal", value=False,
    help="Remove particles that fall below the parcel. Simulates precipitation loss."
)

# 3. Particle Settings
st.sidebar.subheader("Particles")
default_n_particles = PRESETS[preset]["n_particles"] if preset != "Custom" else 1000
n_particles = st.sidebar.number_input(
    "Number of Superdroplets",
    min_value=100, max_value=10000,
    value=default_n_particles, step=100,
    help="Each superdroplet represents many real particles. More = better statistics but slower. 1000-2000 recommended."
)

# 4. Parcel Initial Conditions
st.sidebar.subheader("Parcel Initial Conditions")

# Get default values from preset or use standard defaults
default_T = PRESETS[preset]["T"] if preset != "Custom" else 293.2
default_P = PRESETS[preset]["P"] if preset != "Custom" else 101300.0
default_RH = PRESETS[preset]["RH"] if preset != "Custom" else 0.88
default_w = PRESETS[preset]["w"] if preset != "Custom" else 1.0
default_z = PRESETS[preset]["z"] if preset != "Custom" else 0.0
default_max_z = PRESETS[preset]["max_z"] if preset != "Custom" else 3000.0
default_ascending_mode = PRESETS[preset]["ascending_mode"] if preset != "Custom" else "linear"

T_parcel = st.sidebar.number_input(
    "Temperature (T) [K]",
    min_value=200.0, max_value=320.0, value=default_T, step=0.1,
    help="Initial air parcel temperature. Typical surface values: 280-300 K"
)
P_parcel = st.sidebar.number_input(
    "Pressure (P) [Pa]",
    min_value=10000.0, max_value=105000.0, value=default_P, step=100.0,
    help="Initial air pressure. Sea level ≈ 101325 Pa. Decreases with altitude."
)
RH_parcel = st.sidebar.slider(
    "Relative Humidity (RH)",
    min_value=0.01, max_value=1.5, value=default_RH, step=0.01,
    help="Initial relative humidity (1.0 = 100% = saturation). Higher RH → earlier cloud formation."
)
w_parcel = st.sidebar.slider(
    "Updraft Velocity (w) [m/s]",
    min_value=0.1, max_value=20.0, value=default_w, step=0.1,
    help="Vertical velocity of the parcel. Stronger updrafts → higher supersaturation → more activation."
)
z_parcel = st.sidebar.number_input(
    "Initial Height (z) [m]",
    min_value=0.0, max_value=5000.0, value=default_z, step=100.0,
    help="Starting altitude of the parcel above ground level."
)
max_z = st.sidebar.number_input(
    "Max Height (z_max) [m]",
    min_value=1000.0, max_value=20000.0, value=default_max_z, step=100.0,
    help="Maximum altitude the parcel can reach before stopping ascent."
)

ascending_mode = st.sidebar.selectbox(
    "Ascending Mode",
    options=['linear', 'sine', 'in_cloud_oscillation'],
    index=['linear', 'sine', 'in_cloud_oscillation'].index(default_ascending_mode),
    help="How the parcel ascends: linear (constant w), sine (smooth acceleration), in_cloud_oscillation (up and down)"
)

# 5. Aerosol Settings
st.sidebar.subheader("Aerosol Initialization")
aero_mode = st.sidebar.radio(
    "Initialization Mode",
    options=['Random', 'Weighting_factor'], index=0,
    help="'Random' samples from distribution randomly; 'Weighting_factor' uses deterministic sampling with weights."
)
use_kappa = st.sidebar.checkbox(
    "Use Hygroscopicity (Kappa)", value=False,
    help="Enable κ-Köhler theory for activation. κ represents aerosol hygroscopicity (sea salt ≈ 1.2, sulfate ≈ 0.6)."
)
use_koehler = st.sidebar.checkbox(
    "Use Koehler Activation Radius", value=False,
    help="Initialize droplets at their Köhler equilibrium radius instead of dry aerosol radius."
)

# Aerosol Modes
st.sidebar.markdown("### Aerosol Size Distribution Modes")
st.sidebar.caption("Multi-modal lognormal distribution representing different aerosol populations")
aero_modes_data = {}

with st.sidebar.expander("Configure Aerosol Modes", expanded=True):
    tabs = st.tabs(["Mode 1", "Mode 2", "Mode 3", "Mode 4"])

    # Get defaults from preset or use standard defaults
    if preset != "Custom":
        defaults = PRESETS[preset]["aero"]
    else:
        defaults = [
            {'n': 118.0, 'r': 0.019, 'sigma': 3.3, 'kappa': 1.6}, # Mode 1 - Aitken mode
            {'n': 11.0,  'r': 0.056, 'sigma': 1.6, 'kappa': 1.6}, # Mode 2 - Accumulation mode
            {'n': 0.72,  'r': 0.46,  'sigma': 2.2, 'kappa': 1.6}, # Mode 3 - Coarse mode
            {'n': 0.0,   'r': 0.0,   'sigma': 0.0, 'kappa': 1.6}  # Mode 4 - Optional
        ]

    for i, tab in enumerate(tabs):
        with tab:
            mode_names = ["Aitken/nucleation", "Accumulation", "Coarse", "Optional"]
            st.markdown(f"**Mode {i+1}** ({mode_names[i]})")
            n = st.number_input(
                f"N [cm⁻³] #{i+1}", value=defaults[i]['n'], step=1.0, key=f"n_{i}",
                help="Number concentration of aerosol particles"
            )
            r = st.number_input(
                f"Mean radius [µm] #{i+1}", value=defaults[i]['r'], format="%.3f", step=0.001, key=f"r_{i}",
                help="Geometric mean radius of lognormal distribution"
            )
            sigma = st.number_input(
                f"Std. dev. (σ) #{i+1}", value=defaults[i]['sigma'], step=0.1, key=f"s_{i}",
                help="Geometric standard deviation (width of distribution)"
            )
            kappa = st.number_input(
                f"Kappa (κ) #{i+1}", value=defaults[i]['kappa'], step=0.1, key=f"k_{i}",
                help="Hygroscopicity parameter: sea salt≈1.2, sulfate≈0.6, organic≈0.1-0.3"
            )

            aero_modes_data[i] = {'n': n, 'r': r, 'sigma': sigma, 'kappa': kappa}

# Prepare aerosol grid widget structure (mocking the widget object expected by backend)
class MockGridWidget:
    def __init__(self, modes_data):
        self.widgets = {}
        
        # Helper to create a mock widget with a value attribute
        def create_mock(val):
            return type('obj', (object,), {'value': val})

        for i in range(4):
            # Row indices in widget.py seem to be 1-based for data: 
            # Row 1: N, Row 2: mu, Row 3: sigma, Row 4: kappa
            # Col indices: 0 for Mode 1, 1 for Mode 2, etc.
            
            data = modes_data[i]
            col = i
            
            self.widgets[(1, col)] = create_mock(data['n'])       # N_aero in cm⁻³
            self.widgets[(2, col)] = create_mock(data['r'])       # mu in µm (NOT meters!)
            self.widgets[(3, col)] = create_mock(data['sigma'])   # sigma (dimensionless)
            self.widgets[(4, col)] = create_mock(data['kappa'])   # hygroscopicity (dimensionless)

    def __getitem__(self, key):
        # key is expected to be a tuple (row, col)
        return self.widgets.get(key, type('obj', (object,), {'value': 0.0}))

gridwidget = MockGridWidget(aero_modes_data)

# 6. Entrainment (hidden - module under development)
# Note: Entrainment functionality is disabled in this educational version
# The underlying code is preserved for future development
do_entrainment = False
entrainment_rate = 0.0
stability = 'Stable'
entr_start = 0
entr_end = 0

# --- Main Execution ---

if st.button("Run Simulation", type="primary"):
    
    # Progress bar
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    # Create environment profiles (Required for pressure calculation even if entrainment is off)
    # Need to import or reimplement create_env_profiles logic if it relies on widgets
    # Assuming we can pass values directly if we modify the call or mock widgets
    qv_init = RH_parcel * esatw(T_parcel) / (P_parcel - RH_parcel * esatw(T_parcel)) * r_a / rv
    
    # Use selected stability if entrainment is on, otherwise default to Stable (or whatever makes sense)
    current_stability = stability if do_entrainment else 'Stable'
    
    qv_profiles, theta_profiles, z_env = create_env_profiles(T_parcel, qv_init, z_parcel, P_parcel, current_stability, plot_profiles=False)

    # Run Simulation
    status_text.text("Initializing model...")
    
    # We need to wrap the execution to capture the output
    # The original function runs the whole loop. We might want to modify it to yield progress
    # For now, we run it as a block.
    
    start_time = time.time()
    
    # Output capture for real-time logs
    log_placeholder = st.empty()
    
    # Progress callback for Streamlit
    def update_progress(percent, message):
        progress_bar.progress(percent)
        status_text.text(message)
    
    # Call the main routine (array-based version uses direct parameters)
    try:
        results = timesteps_function(
            n_particles=n_particles,
            P_parcel=P_parcel,
            RH_parcel=RH_parcel,
            T_parcel=T_parcel,
            w_parcel=w_parcel,
            nt=nt,
            dt=dt,
            rm_spec=rm_spec,
            ascending_mode=ascending_mode,
            z_parcel=z_parcel,
            max_z=max_z,
            do_condensation=do_condensation,
            do_collision=do_collision,
            mode_aero_init=aero_mode,
            N_aero=[mode_data['n'] * 1e6 for mode_data in aero_modes_data.values()],  # Convert cm⁻³ to m⁻³
            mu_aero=[mode_data['r'] * 1e-6 for mode_data in aero_modes_data.values()],  # Convert µm to m
            sigma_aero=[mode_data['sigma'] for mode_data in aero_modes_data.values()],
            k_aero=[mode_data['kappa'] for mode_data in aero_modes_data.values()],
            kohler_activation_radius=use_koehler,
            switch_kappa_koehler=use_kappa,
            do_sedi_removal=do_sedi_removal,
            entrainment_rate=entrainment_rate,
            switch_entrainment=do_entrainment,
            qv_profiles=qv_profiles,
            theta_profiles=theta_profiles,
            entrainment_start=entr_start,
            entrainment_end=entr_end,
            output_interval=max(1, nt // 100),  # Compute diagnostics every ~1% of simulation
            verbose=False,  # Disable verbose mode for cleaner interface
            progress_callback=update_progress  # Pass progress callback for real-time updates
        )
        
        # Unpack results
        (nt_out, dt_out, time_array, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array, 
         qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, spectra_arr, con_ts, act_ts, evp_ts, dea_ts, 
         acc_ts, aut_ts, precip_ts, particles_array, rc_liq_avg_array, rc_liq_std_array) = results

        end_time = time.time()
        st.success(f"✅ Simulation completed in {end_time - start_time:.2f} seconds!")
        progress_bar.progress(100)
        
        # --- Visualization Dashboard ---
        st.markdown("---")
        st.header("📊 Simulation Results")
        
        # Key Statistics Summary at top
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Final Height", f"{z_parcel_array[-1]:.0f} m", f"+{z_parcel_array[-1]-z_parcel_array[0]:.0f} m")
        with col2:
            st.metric("Final RH", f"{RH_parcel_array[-1]*100:.1f}%", f"{(RH_parcel_array[-1]-RH_parcel_array[0])*100:+.1f}%")
        with col3:
            st.metric("Max Cloud Water", f"{np.max(qc_ts):.2f} g/kg")
        with col4:
            st.metric("Max Rain Water", f"{np.max(qr_ts):.2f} g/kg")
        
        # Create comprehensive tabs
        tab1, tab2, tab3, tab4, tab5 = st.tabs([
            "🌡️ Thermodynamics", 
            "💧 Microphysics", 
            "📈 Spectrum Evolution",
            "⚙️ Process Rates",
            "📋 Summary"
        ])
        
        # === TAB 1: Thermodynamics ===
        with tab1:
            st.subheader("Parcel Thermodynamics")
            
            col1, col2 = st.columns(2)
            
            with col1:
                # Height and RH evolution
                fig_thermo = make_subplots(specs=[[{"secondary_y": True}]])
                fig_thermo.add_trace(
                    go.Scatter(x=time_array, y=z_parcel_array, name="Height", 
                              line=dict(color='blue', width=2)), 
                    secondary_y=False
                )
                fig_thermo.add_trace(
                    go.Scatter(x=time_array, y=RH_parcel_array*100, name="RH", 
                              line=dict(color='red', width=2, dash='dot')), 
                    secondary_y=True
                )
                fig_thermo.update_xaxes(title_text="Time (s)")
                fig_thermo.update_yaxes(title_text="Height (m)", secondary_y=False, color='blue')
                fig_thermo.update_yaxes(title_text="Relative Humidity (%)", secondary_y=True, color='red')
                fig_thermo.update_layout(
                    title="Parcel Ascent and Saturation",
                    hovermode='x unified',
                    height=400
                )
                st.plotly_chart(fig_thermo, use_container_width=True)
                
                # Vertical profiles
                fig_vert = make_subplots(rows=1, cols=2, 
                                        subplot_titles=("Temperature Profile", "RH Profile"))
                fig_vert.add_trace(
                    go.Scatter(x=T_parcel_array, y=z_parcel_array, 
                              line=dict(color='orange', width=2)),
                    row=1, col=1
                )
                fig_vert.add_trace(
                    go.Scatter(x=RH_parcel_array*100, y=z_parcel_array,
                              line=dict(color='purple', width=2)),
                    row=1, col=2
                )
                fig_vert.update_xaxes(title_text="Temperature (K)", row=1, col=1)
                fig_vert.update_xaxes(title_text="RH (%)", row=1, col=2)
                fig_vert.update_yaxes(title_text="Height (m)", row=1, col=1)
                fig_vert.update_layout(height=400, showlegend=False)
                st.plotly_chart(fig_vert, use_container_width=True)
            
            with col2:
                # Temperature time series
                fig_temp = go.Figure()
                fig_temp.add_trace(go.Scatter(
                    x=time_array, y=T_parcel_array,
                    fill='tozeroy', fillcolor='rgba(255,100,100,0.2)',
                    line=dict(color='darkred', width=2)
                ))
                fig_temp.update_layout(
                    title="Temperature Evolution",
                    xaxis_title="Time (s)", 
                    yaxis_title="Temperature (K)",
                    height=400
                )
                st.plotly_chart(fig_temp, use_container_width=True)
                
                # Vapor mixing ratio
                fig_q = go.Figure()
                fig_q.add_trace(go.Scatter(
                    x=time_array, y=q_parcel_array*1000,
                    fill='tozeroy', fillcolor='rgba(100,100,255,0.2)',
                    line=dict(color='darkblue', width=2)
                ))
                fig_q.update_layout(
                    title="Water Vapor Mixing Ratio",
                    xaxis_title="Time (s)",
                    yaxis_title="qv (g/kg)",
                    height=400
                )
                st.plotly_chart(fig_q, use_container_width=True)
        
        # === TAB 2: Microphysics ===
        with tab2:
            st.subheader("Cloud and Rain Microphysics")
            
            col1, col2 = st.columns(2)
            
            with col1:
                # Mixing ratios
                fig_mix = go.Figure()
                fig_mix.add_trace(go.Scatter(
                    x=time_array, y=qc_ts, name="Cloud Water (qc)",
                    fill='tozeroy', line=dict(color='cyan', width=2)
                ))
                fig_mix.add_trace(go.Scatter(
                    x=time_array, y=qr_ts, name="Rain Water (qr)",
                    fill='tozeroy', line=dict(color='blue', width=2)
                ))
                fig_mix.update_layout(
                    title="Liquid Water Mixing Ratios",
                    xaxis_title="Time (s)",
                    yaxis_title="Mixing Ratio (g/kg)",
                    height=400,
                    hovermode='x unified'
                )
                st.plotly_chart(fig_mix, use_container_width=True)
                
                # Mean radius evolution
                fig_radius = go.Figure()
                fig_radius.add_trace(go.Scatter(
                    x=time_array, y=rc_liq_avg_array*1e6,
                    error_y=dict(type='data', array=rc_liq_std_array*1e6, visible=True),
                    name="Mean cloud droplet radius",
                    line=dict(color='green', width=2)
                ))
                fig_radius.update_layout(
                    title="Mean Cloud Droplet Radius (±1σ)",
                    xaxis_title="Time (s)",
                    yaxis_title="Radius (µm)",
                    height=400
                )
                st.plotly_chart(fig_radius, use_container_width=True)
            
            with col2:
                # Number concentrations
                fig_conc = go.Figure()
                fig_conc.add_trace(go.Scatter(
                    x=time_array, y=nc_ts, name="Cloud Droplets (Nc)",
                    line=dict(color='lightblue', width=2)
                ))
                fig_conc.add_trace(go.Scatter(
                    x=time_array, y=nr_ts, name="Rain Drops (Nr)",
                    line=dict(color='darkblue', width=2)
                ))
                fig_conc.update_layout(
                    title="Number Concentrations",
                    xaxis_title="Time (s)",
                    yaxis_title="Concentration (#/mg)",
                    yaxis_type="log",
                    height=400,
                    hovermode='x unified'
                )
                st.plotly_chart(fig_conc, use_container_width=True)
                
                # Vertical profiles of liquid water
                fig_vert_lw = go.Figure()
                fig_vert_lw.add_trace(go.Scatter(
                    x=qc_ts, y=z_parcel_array, name="Cloud Water",
                    line=dict(color='cyan', width=2)
                ))
                fig_vert_lw.add_trace(go.Scatter(
                    x=qr_ts, y=z_parcel_array, name="Rain Water",
                    line=dict(color='blue', width=2)
                ))
                fig_vert_lw.update_layout(
                    title="Vertical Liquid Water Profile",
                    xaxis_title="Mixing Ratio (g/kg)",
                    yaxis_title="Height (m)",
                    height=400
                )
                st.plotly_chart(fig_vert_lw, use_container_width=True)
        
        # === TAB 3: Spectrum Evolution ===
        with tab3:
            st.subheader("Droplet Size Distribution Evolution")
            
            # Create animated spectrum plot
            if len(spectra_arr.shape) > 1 and spectra_arr.shape[0] > 1:
                st.info("🎬 Use the slider below to see how the droplet size distribution evolves over time!")
                
                # Time slider
                time_idx = st.slider(
                    "Select Time Step",
                    min_value=0,
                    max_value=len(time_array)-1,
                    value=len(time_array)-1,
                    format=f"t = %.0f s"
                )
                
                # Get spectrum at selected time
                spectrum_t = spectra_arr[time_idx]
                
                # Create spectrum plot
                fig_spec = go.Figure()
                fig_spec.add_trace(go.Scatter(
                    x=rm_spec*1e6,
                    y=spectrum_t,
                    mode='lines+markers',
                    line=dict(color='darkgreen', width=3),
                    marker=dict(size=6),
                    name=f"t = {time_array[time_idx]:.1f} s"
                ))
                
                # Add initial spectrum for comparison
                fig_spec.add_trace(go.Scatter(
                    x=rm_spec*1e6,
                    y=spectra_arr[0],
                    mode='lines',
                    line=dict(color='gray', width=1, dash='dash'),
                    name="Initial (t=0)"
                ))
                
                fig_spec.update_xaxes(type="log", title_text="Radius (µm)")
                fig_spec.update_yaxes(type="log", title_text="dN/dlogr (#/m³)")
                fig_spec.update_layout(
                    title=f"Droplet Size Distribution at t = {time_array[time_idx]:.1f} s (z = {z_parcel_array[time_idx]:.0f} m)",
                    height=500,
                    hovermode='x unified'
                )
                st.plotly_chart(fig_spec, use_container_width=True)
                
                # Spectrum evolution heatmap
                st.subheader("Spectrum Evolution Heatmap")
                
                # Sample every N timesteps for plotting
                sample_rate = max(1, len(time_array) // 100)
                time_sample = time_array[::sample_rate]
                spectra_sample = spectra_arr[::sample_rate, :]
                
                fig_heatmap = go.Figure(data=go.Heatmap(
                    x=rm_spec*1e6,
                    y=time_sample,
                    z=np.log10(spectra_sample + 1),  # log scale for better visualization
                    colorscale='Viridis',
                    colorbar=dict(title="log₁₀(dN/dlogr)")
                ))
                fig_heatmap.update_xaxes(type="log", title_text="Radius (µm)")
                fig_heatmap.update_yaxes(title_text="Time (s)")
                fig_heatmap.update_layout(
                    title="Droplet Spectrum Evolution Over Time",
                    height=400
                )
                st.plotly_chart(fig_heatmap, use_container_width=True)
            else:
                st.warning("Spectrum data not available or insufficient timesteps for animation.")
        
        # === TAB 4: Process Rates ===
        with tab4:
            st.subheader("Microphysical Process Rates")
            
            col1, col2 = st.columns(2)
            
            with col1:
                # Condensation processes
                fig_cond = go.Figure()
                fig_cond.add_trace(go.Scatter(
                    x=time_array, y=con_ts, name="Condensation",
                    line=dict(color='green', width=2)
                ))
                fig_cond.add_trace(go.Scatter(
                    x=time_array, y=evp_ts, name="Evaporation",
                    line=dict(color='red', width=2)
                ))
                fig_cond.add_trace(go.Scatter(
                    x=time_array, y=act_ts, name="Activation",
                    line=dict(color='blue', width=2, dash='dot')
                ))
                fig_cond.add_trace(go.Scatter(
                    x=time_array, y=dea_ts, name="Deactivation",
                    line=dict(color='orange', width=2, dash='dot')
                ))
                fig_cond.update_layout(
                    title="Condensation/Evaporation Rates",
                    xaxis_title="Time (s)",
                    yaxis_title="Rate (g/kg/s)",
                    height=400,
                    hovermode='x unified'
                )
                st.plotly_chart(fig_cond, use_container_width=True)
            
            with col2:
                # Collision processes (if enabled)
                fig_coll = go.Figure()
                fig_coll.add_trace(go.Scatter(
                    x=time_array, y=acc_ts, name="Accretion",
                    line=dict(color='purple', width=2)
                ))
                fig_coll.add_trace(go.Scatter(
                    x=time_array, y=aut_ts, name="Autoconversion",
                    line=dict(color='brown', width=2)
                ))
                fig_coll.update_layout(
                    title="Collision/Coalescence Rates",
                    xaxis_title="Time (s)",
                    yaxis_title="Rate (g/kg/s)",
                    height=400,
                    hovermode='x unified'
                )
                st.plotly_chart(fig_coll, use_container_width=True)
            
            # Summary statistics
            st.subheader("Process Statistics")
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Total Condensation", f"{np.sum(con_ts)*dt:.2f} g/kg")
            with col2:
                st.metric("Total Evaporation", f"{np.sum(evp_ts)*dt:.2f} g/kg")
            with col3:
                st.metric("Total Accretion", f"{np.sum(acc_ts)*dt:.2f} g/kg")
            with col4:
                st.metric("Total Autoconv.", f"{np.sum(aut_ts)*dt:.2f} g/kg")
        
        # === TAB 5: Summary ===
        with tab5:
            st.subheader("Simulation Summary")
            
            # Simulation parameters
            st.markdown("### 📋 Simulation Parameters")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.markdown(f"""
                **Time Integration**
                - Steps: {nt}
                - Time step: {dt} s
                - Total time: {nt*dt:.0f} s
                - Runtime: {end_time - start_time:.2f} s
                """)
            with col2:
                # Note: T_widget, P_widget, RH_widget, w_widget are not defined in this scope.
                # Assuming they refer to the initial values T_parcel, P_parcel, RH_parcel, w_parcel
                st.markdown(f"""
                **Parcel Initial Conditions**
                - Temperature: {T_parcel} K
                - Pressure: {P_parcel/100:.0f} hPa
                - RH: {RH_parcel*100:.0f}%
                - Updraft: {w_parcel} m/s
                """)
            with col3:
                st.markdown(f"""
                **Physics**
                - Particles: {n_particles}
                - Condensation: {"✅" if do_condensation else "❌"}
                - Collision: {"✅" if do_collision else "❌"}
                - Sedimentation: {"✅" if do_sedi_removal else "❌"}
                """)
            
            # Key results
            st.markdown("### 🎯 Key Results")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown(f"""
                **Cloud Properties**
                - Max cloud water: {np.max(qc_ts):.3f} g/kg
                - Max cloud droplets: {np.max(nc_ts):.0f} #/mg
                - Mean droplet radius: {rc_liq_avg_array[-1]*1e6:.2f} ± {rc_liq_std_array[-1]*1e6:.2f} µm
                - Final height: {z_parcel_array[-1]:.0f} m
                """)
            with col2:
                st.markdown(f"""
                **Precipitation**
                - Max rain water: {np.max(qr_ts):.3f} g/kg
                - Max rain drops: {np.max(nr_ts):.0f} #/mg
                - Total condensed: {np.sum(con_ts)*dt:.3f} g/kg
                - Total precipitated: {np.sum(precip_ts)*dt:.3f} g/kg
                """)
            
            # Educational notes
            st.markdown("### 📚 Educational Notes")
            st.info("""
            **Understanding the Results:**
            - 🌡️ **Thermodynamics**: Watch how temperature decreases and RH increases as the parcel rises and expands adiabatically
            - 💧 **Microphysics**: Cloud droplets form when RH > 100%. Rain forms through collision-coalescence
            - 📈 **Spectrum**: The droplet size distribution shows how particles grow from aerosols to cloud droplets to rain
            - ⚙️ **Processes**: Condensation dominates cloud formation; collision creates rain
            
            **Try These Experiments:**
            - Increase updraft velocity → stronger supersaturation → more activation
            - Increase aerosol concentration → more, smaller droplets → less rain
            - Enable collision → watch autoconversion create rain drops
            """)

    except ValueError as e:
        st.error(f"⚠️ Invalid parameter value: {str(e)}")
        st.info("💡 Try adjusting your input parameters. Common issues:\n"
                "- RH too high (>1.5) or too low (<0.01)\n"
                "- Temperature outside realistic range\n"
                "- Zero or negative aerosol concentrations")
    except MemoryError:
        st.error("💾 Out of memory! Try reducing the number of particles or time steps.")
    except Exception as e:
        st.error(f"❌ Simulation failed: {str(e)}")
        with st.expander("Show technical details"):
            st.exception(e)
        st.info("💡 If this persists, try:\n"
                "- Reducing particle count\n"
                "- Using shorter simulation time\n"
                "- Disabling collision (most computationally expensive)")

else:
    st.info("💡 Configure parameters in the sidebar and click 'Run Simulation' to start!")
    st.markdown("""
    ### Quick Start Guide:
    1. **Set time parameters** - Try 100 steps × 1s for quick results
    2. **Choose physics** - Enable condensation (essential), optionally collision
    3. **Configure parcel** - Default values work well for cumulus clouds
    4. **Set aerosol** - Use defaults or explore different modes
    5. **Click Run!** - Watch the real-time progress bar
    
    ### What to Expect:
    - 1,000 particles: ~10-20 seconds
    - 5,000 particles: ~30-60 seconds  
    - 10,000 particles: ~1-2 minutes
    """)
