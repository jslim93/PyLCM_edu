# PyLCM Educational GUI (Streamlit)

An interactive browser GUI over the validated PyLCM SoA engine
(`PyLCM.timestep_soa.run_soa`). It adds no physics — every result comes from the
engine. Because a 2000×1500 ascent runs in ~0.3 s and results are cached, moving
any slider re-runs the model essentially instantly.

## Run

    pip install -r requirements.txt        # includes streamlit, plotly
    streamlit run app_streamlit.py

Opens at http://localhost:8501.

## Modes and levels

- **Mode**
  - *Guided experiment* — pick one of the five lessons below.
  - *Free exploration* — adjust any parameter freely.
- **Level**
  - *Beginner* — plain-language explanation + core sliders (RH, updraft w,
    aerosol preset, collisions on/off; plus an IHMD slider in the entrainment
    lesson).
  - *Advanced* — adds the governing equation (`st.latex`), extra numerics
    (super-droplet count, steps, turbulent collision kernel + ε), and a
    **Diagnostics** tab with final Nc / Nr / Na.

## The five guided lessons

1. **Aerosol activation (Köhler)** — vary RH / updraft / aerosol; watch the
   cloud-droplet number Nc respond.
2. **Condensational growth** — collisions off; the DSD narrows with height.
3. **Collision → rain** — collisions on with the maritime preset; the DSD
   broadens and a rain mode (Nr, qr) appears.
4. **Maritime vs continental (A/B)** — runs *both* presets and overlays their
   final DSD and qr-vs-height.
5. **Entrainment mixing (IHMD)** — with `lambda_ent ≈ 5e-4` and collisions off,
   runs IHMD ∈ {0, 0.5, 1} against a no-entrainment baseline. Homogeneous
   mixing (IHMD=0) keeps Nc while droplets shrink; inhomogeneous mixing
   (IHMD=1) drops Nc. Law: `N_c/N_{c,0} = (q_c/q_{c,0})^IHMD`
   (Lim & Hoffmann, 2023).

## Result tabs

- **DSD** (central plot) — overlaid log-x lines of `dsd_r` [µm] vs `dsd_n`
  [cm⁻³] at the collected heights.
- **Time series** — qc and qr vs height.
- **Vertical profile** — temperature vs height.
- **Diagnostics** (Advanced only) — final Nc / Nr / Na as metrics.

## Notes

- The optional ascent animation from the design spec is **deferred** — a future
  addition (e.g. a Plotly frames animation or a slider over collected times).
  The core app does not depend on it.
