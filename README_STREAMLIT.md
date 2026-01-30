# PyLCM Streamlit GUI

A modern, interactive web interface for the PyLCM Educational Microphysics Model.

## Features

- **Interactive Sidebar**: Adjust all model parameters (time step, parcel properties, aerosols, etc.) in real-time.
- **Optimized Physics**: Automatically uses Numba-accelerated modules if available (5-10x faster).
- **Visualization**: Interactive Plotly charts for:
  - Time series (Height, RH, Mixing Ratios, Number Concentrations)
  - Vertical profiles (Temperature, Liquid Water)
  - Droplet size distributions
- **No Installation**: Runs locally in your browser.

## How to Run

1. **Install Dependencies** (if not already installed):
   ```bash
   pip install streamlit plotly
   # or
   conda install streamlit plotly
   ```

2. **Run the App**:
   ```bash
   streamlit run app_streamlit.py
   ```

3. **Open in Browser**:
   The app will automatically open in your default web browser (usually at `http://localhost:8501`).

## Usage Tips

- **Optimizations**: The app checks for `PyLCM.collision_optimized` and `PyLCM.condensation_optimized`. If found, it enables them for faster simulations.
- **Entrainment**: Enable entrainment in the sidebar to simulate mixing with environmental air.
- **Ascending Modes**: Choose between 'linear', 'sine', or 'in_cloud_oscillation' to simulate different updraft dynamics.

## Troubleshooting

- **"ModuleNotFoundError"**: Ensure you are in the `PyLCM_edu` directory and your environment has `numpy`, `scipy`, `pandas`, `plotly`, and `streamlit` installed.
- **Slow Performance**: Reduce `nt` (number of steps) or `n_particles` if the simulation is too slow, or ensure Numba optimizations are working.
