# PyLCM Educational GUI (Streamlit)

An interactive browser GUI over the validated PyLCM SoA engine
(`PyLCM.timestep_soa.run_soa`). It adds no physics — every result comes from the
engine. Because results are cached on every parameter, moving any slider
re-runs the model essentially instantly.

The visual design deliberately mirrors the author's own notebooks: the 3×2
time-series panel of `PyLCM/animation.py` and the DSD contour + viridis spectra
of `Post_process/print_plot.py`.

## Run

    pip install -r requirements.txt        # includes streamlit, plotly
    streamlit run app_streamlit.py

Opens at http://localhost:8501.

## Modes and levels

- **Mode**
  - *Guided experiment* — pick one of the five lessons below.
  - *Free exploration* — adjust the full parameter set freely; Advanced level
    here unlocks the multi-mode aerosol editor.
- **Level**
  - *Beginner* — plain-language explanation + core controls (T₀, P₀, RH,
    updraft w, ascending mode, an aerosol **preset**, collisions, entrainment).
  - *Advanced* — adds the governing equation (`st.latex`), numerics
    (super-droplet count, steps nt, dt), the turbulent collision kernel + ε,
    the per-mode aerosol editor, and a **Diagnostics** tab.

## Full parameter set

- **Parcel** — initial temperature T₀, pressure P₀, relative humidity RH,
  updraft w, and **ascending mode** (`linear` / `sine` / `in_cloud_oscillation`).
  Advanced exposes the step count nt, step length dt, and super-droplet count.
- **Aerosol** — a quick **preset** (default / maritime / continental / arctic),
  or, in Advanced + Free exploration, a full **multi-mode editor**: choose 1–4
  lognormal modes, each with N (cm⁻³), mean radius (µm), σ, and κ. A **GCCN**
  checkbox appends a coarse giant-CCN mode (N≈0.01 cm⁻³, r≈2 µm, κ≈1.2) that
  seeds early raindrops.
- **Physics** — collision–coalescence on/off, the Wang & Ayala turbulent kernel
  with dissipation ε (Advanced), and entrainment (rate λ + inhomogeneous mixing
  degree IHMD). Every control carries a `help=` tooltip.

## The five guided lessons

1. **Aerosol activation (Köhler)** — vary RH / updraft / aerosol; watch the
   cloud-droplet number Nc respond.
2. **Condensational growth** — collisions off; the DSD narrows with height.
3. **Collision → rain** — collisions on with the maritime preset; the DSD
   broadens and a rain mode (Nr, qr) appears.
4. **Maritime vs continental (A/B)** — runs *both* presets and overlays their
   six-panel time series, final DSD, particle populations, and profiles.
5. **Entrainment mixing (IHMD)** — with `lambda_ent ≈ 5e-4`, runs IHMD ∈
   {0, 0.5, 1} against a no-entrainment baseline. Homogeneous mixing (IHMD=0)
   keeps Nc while droplets shrink; inhomogeneous mixing (IHMD=1) drops Nc. Law:
   `N_c/N_{c,0} = (q_c/q_{c,0})^IHMD` (Lim & Hoffmann, 2023).

## Result tabs

- **Time series** — the six-panel `make_subplots(3×2)` matching `animation.py`:
  RH (%), vapour q_v (g/kg), height z (m), temperature T (K), the three mixing
  ratios qa/qc/qr (g/kg) as three lines, and the three number concentrations
  na/nc/nr (cm⁻³) on a log y-axis. Overlays multiple runs in A/B and sweeps.
- **DSD** — the `Post_process` design: a filled time×radius **contour** of dN
  (log colour) with the mean-radius `rv` line overlaid, plus **spectra lines**
  at incremental times coloured by a viridis time colormap. Multi-run mode
  overlays the final spectra instead.
- **Particle distribution** — the final super-droplet population: each marker is
  one super-droplet at its radius `r = (M/(A·4/3·π·ρ_liq))^(1/3)`, with marker
  size/colour scaling by multiplicity A, plus a multiplicity-weighted histogram
  per radius bin. Shows the spread of sizes directly.
- **Vertical profiles** — temperature vs height and LWC (qc+qr) vs height.
- **Diagnostics** (Advanced only) — final Nc / Nr / Na, supersaturation, LWC,
  and mean radius as `st.metric`.

## Implementation notes

- Dense time sampling (~120 collected times) feeds the time series and the DSD
  contour; `dsd_n` is stacked into a 2-D `(time × radius)` array for the
  heatmap.
- The cache key includes every parameter (tuples passed as tuples for
  hashability), so identical configurations return instantly.
- Tests (`tests/test_gui_smoke.py`) never launch a Streamlit server; they
  exercise the engine adapter and verify the app module spec resolves.
