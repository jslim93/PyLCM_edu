"""PyLCM Educational Streamlit GUI.

A thin browser presentation layer over the validated SoA engine
(`PyLCM.timestep_soa.run_soa`). No new physics — all results come from the
engine. Guided experiments + free exploration, beginner/advanced level toggle,
and instant re-run via `@st.cache_data`.

Run:  streamlit run app_streamlit.py
"""

import numpy as np
import streamlit as st
import plotly.graph_objects as go

from PyLCM.timestep_soa import run_soa

st.set_page_config(page_title="PyLCM Educational Model", page_icon="☁️",
                   layout="wide")

# Aerosol presets: (N_raw cm^-3 per mode, mode mean diameter µm, geom. std, kappa)
PRESETS = {
    "default":     dict(N_raw=(118., 11., .72), mu_um=(.019, .056, .46),
                        sig=(3.3, 1.6, 2.2), kappa=1.6),
    "maritime":    dict(N_raw=(100., 20.),      mu_um=(.08, .4),
                        sig=(1.6, 2.0),         kappa=1.0),
    "continental": dict(N_raw=(3200., 2900.),   mu_um=(.012, .04),
                        sig=(1.7, 2.0),         kappa=0.3),
    "arctic":      dict(N_raw=(15., 5.),        mu_um=(.05, .2),
                        sig=(1.6, 2.0),         kappa=0.5),
}


@st.cache_data(show_spinner=False)
def cached_run(seed, n_ptcl, nt, RH, w, preset, collisions, switch_turb, eps,
               lambda_ent, ihmd):
    """Cached wrapper around run_soa keyed on every parameter, so an unchanged
    configuration returns instantly. Collects six evenly spaced heights."""
    p = PRESETS[preset]
    collect = tuple(range(nt // 6, nt + 1, nt // 6))
    out, _ = run_soa(seed=seed, n_ptcl=n_ptcl, nt=nt, T0=293.2, P0=1013e2,
                     RH=RH, w=w, N_raw=p["N_raw"], mu_um=p["mu_um"],
                     sig=p["sig"], kappa=p["kappa"], collisions=collisions,
                     switch_turb=switch_turb, eps=eps, lambda_ent=lambda_ent,
                     ihmd=ihmd, collect=collect)
    return out


# --- plot helpers (pure: dict-of-time -> Plotly figure) ------------------------

def dsd_figure(runs):
    """runs: list of (label, out-dict). Overlaid log-x DSD lines at each height."""
    fig = go.Figure()
    for label, out in runs:
        times = sorted(out)
        for i, t in enumerate(times):
            d = out[t]
            name = f"z={d['z']:.0f}m" if len(runs) == 1 else \
                   f"{label} z={d['z']:.0f}m"
            # emphasize the final (top) height for multi-run overlays
            width = 3 if (len(runs) > 1 and i == len(times) - 1) else 1.5
            fig.add_scatter(x=d["dsd_r"] * 1e6, y=d["dsd_n"], mode="lines",
                            name=name, line=dict(width=width))
    fig.update_layout(xaxis_title="radius (µm)", yaxis_title="dN (cm⁻³)",
                      xaxis_type="log", legend_title="height")
    return fig


def timeseries_figure(runs):
    """qc and qr vs height for each run."""
    fig = go.Figure()
    dash = ["solid", "dash", "dot", "dashdot"]
    for j, (label, out) in enumerate(runs):
        times = sorted(out)
        z = [out[t]["z"] for t in times]
        for k, lab in [("qc", "cloud water"), ("qr", "rain water")]:
            nm = lab if len(runs) == 1 else f"{label} {lab}"
            fig.add_scatter(x=z, y=[out[t][k] for t in times],
                            mode="lines+markers", name=nm,
                            line=dict(dash=dash[j % len(dash)]))
    fig.update_layout(xaxis_title="height (m)",
                      yaxis_title="mixing ratio (g/kg)")
    return fig


def profile_figure(runs):
    """Temperature vs height."""
    fig = go.Figure()
    for label, out in runs:
        times = sorted(out)
        z = [out[t]["z"] for t in times]
        nm = "T" if len(runs) == 1 else label
        fig.add_scatter(x=[out[t]["T"] for t in times], y=z,
                        mode="lines+markers", name=nm)
    fig.update_layout(xaxis_title="T (°C)", yaxis_title="height (m)")
    return fig


def render_tabs(runs, level):
    """Render the four result tabs for one or more runs."""
    labels = ["DSD", "Time series", "Vertical profile"]
    if level == "Advanced":
        labels.append("Diagnostics")
    tabs = st.tabs(labels)
    with tabs[0]:
        st.caption("Droplet size distribution — the central diagnostic. "
                   "Growth narrows it; collisions broaden it toward rain sizes.")
        st.plotly_chart(dsd_figure(runs), use_container_width=True)
    with tabs[1]:
        st.plotly_chart(timeseries_figure(runs), use_container_width=True)
    with tabs[2]:
        st.plotly_chart(profile_figure(runs), use_container_width=True)
    if level == "Advanced":
        with tabs[3]:
            cols = st.columns(len(runs))
            for c, (label, out) in zip(cols, runs):
                last = out[sorted(out)[-1]]
                c.markdown(f"**{label}**" if len(runs) > 1 else "**Final state**")
                c.metric("Nc (cm⁻³)", f"{last['NC']:.2f}")
                c.metric("Nr (cm⁻³)", f"{last['NR']:.3f}")
                c.metric("Na residual (cm⁻³)", f"{last['NA']:.2f}")


# --- sidebar -------------------------------------------------------------------

st.sidebar.title("☁️ PyLCM")
mode = st.sidebar.radio("Mode", ["Guided experiment", "Free exploration"],
                        help="Guided walks through five lessons; Free exposes "
                             "all parameters.")
level = st.sidebar.radio("Level", ["Beginner", "Advanced"], horizontal=True,
                         help="Advanced adds equations, extra parameters and a "
                              "diagnostics tab.")

LESSONS = ["1 · Aerosol activation", "2 · Condensational growth",
           "3 · Collision → rain", "4 · Maritime vs continental",
           "5 · Entrainment mixing (IHMD)"]
lesson = (st.sidebar.selectbox("Experiment", LESSONS)
          if mode == "Guided experiment" else None)

st.sidebar.subheader("Parcel")
RH = st.sidebar.slider("Initial RH", 0.80, 1.0, 0.92, 0.01,
                       help="Higher RH → cloud forms sooner.")
w = st.sidebar.slider("Updraft w (m/s)", 0.1, 5.0, 1.0, 0.1,
                      help="Stronger updraft → higher supersaturation → more "
                           "droplets activate.")
preset = st.sidebar.selectbox(
    "Aerosol", list(PRESETS),
    help="Maritime: few large CCN → fast rain. Continental: many small CCN → "
         "suppressed rain.")
collisions = st.sidebar.checkbox("Collision–coalescence", True,
                                 help="Required for rain to form.")

# IHMD only matters for the entrainment lesson; expose it at both levels there.
ihmd = 0.0
lambda_ent = 0.0
if lesson == LESSONS[4]:
    ihmd = st.sidebar.slider(
        "IHMD (mixing degree)", 0.0, 1.0, 0.5, 0.1,
        help="0 = homogeneous (all droplets shrink, number kept); "
             "1 = inhomogeneous (some droplets evaporate entirely, number drops).")
    lambda_ent = 5e-4

if level == "Advanced":
    st.sidebar.subheader("Numerics")
    n_ptcl = int(st.sidebar.number_input("Super-droplets", 500, 20000, 2000, 500))
    nt = int(st.sidebar.number_input("Steps (nt)", 300, 3600, 1500, 300))
    switch_turb = st.sidebar.checkbox("Turbulent collision kernel", False,
                                      help="Wang & Ayala enhancement.")
    eps = (st.sidebar.number_input("ε (m²/s³)", 0.0, 0.1, 0.01, 0.01)
           if switch_turb else 0.0)
else:
    n_ptcl, nt, switch_turb, eps = 2000, 1500, False, 0.0


# --- lesson copy ---------------------------------------------------------------

EXPLAIN = {
    LESSONS[0]: (
        "**Aerosol activation (Köhler theory).** As the parcel rises, it cools "
        "and supersaturation grows. Once it crosses each particle's critical "
        "value, that haze particle *activates* into a cloud droplet. Raise the "
        "updraft or RH and watch the cloud-droplet number **Nc** rise.",
        r"S_{crit} \propto \sqrt{A^3 / B} \quad\text{(Köhler curve peak)}"),
    LESSONS[1]: (
        "**Condensational growth.** With collisions OFF, every droplet grows "
        "only by vapour diffusion. Small droplets grow faster than large ones "
        "(dr/dt ∝ 1/r), so the size distribution **narrows** with height.",
        r"\frac{dr}{dt} = \frac{G\,S}{r}"),
    LESSONS[2]: (
        "**Collision → rain.** With collisions ON (maritime aerosol), larger "
        "droplets fall faster and collect smaller ones. The distribution "
        "**broadens** and a rain mode (Nr, qr) develops.",
        r"K(r_1,r_2) = \pi (r_1+r_2)^2\, |v_1-v_2|\, E"),
    LESSONS[3]: (
        "**Maritime vs continental (A/B).** Both clouds get the same updraft, "
        "but maritime air has *few large* CCN that rain quickly, while "
        "continental air has *many small* CCN that make many tiny droplets and "
        "**suppress** rain at the same water content.",
        r"N_c \uparrow \;\Rightarrow\; \bar r \downarrow \;\Rightarrow\; "
        r"\text{rain suppressed}"),
    LESSONS[4]: (
        "**Entrainment mixing (IHMD).** Dry environmental air is entrained and "
        "evaporates cloud water. *Homogeneous* mixing (IHMD=0) shrinks every "
        "droplet but keeps the number; *inhomogeneous* mixing (IHMD=1) "
        "evaporates whole droplets but keeps the survivors' size. We run "
        "IHMD ∈ {0, 0.5, 1} and overlay the final spectra. "
        "(Lim & Hoffmann, 2023.)",
        r"N_c / N_{c,0} = \left(q_c / q_{c,0}\right)^{\mathrm{IHMD}}"),
}


# --- main panel ----------------------------------------------------------------

if lesson:
    txt, eqn = EXPLAIN[lesson]
    st.header(lesson)
    st.markdown(txt)
    if level == "Advanced":
        st.latex(eqn)
    # lesson-driven parameter overrides
    if lesson == LESSONS[1]:
        collisions = False
    elif lesson == LESSONS[2]:
        collisions = True
        preset = "maritime"
else:
    st.header("Free exploration")
    st.markdown("Adjust any parameter in the sidebar; the model re-runs "
                "instantly. The DSD tab is the central diagnostic.")


if lesson == LESSONS[3]:
    # A/B: run both presets, overlay.
    with st.spinner("Running maritime and continental parcels…"):
        out_m = cached_run(0, n_ptcl, nt, RH, w, "maritime", True,
                           switch_turb, eps, 0.0, 0.0)
        out_c = cached_run(0, n_ptcl, nt, RH, w, "continental", True,
                           switch_turb, eps, 0.0, 0.0)
    render_tabs([("maritime", out_m), ("continental", out_c)], level)

elif lesson == LESSONS[4]:
    # IHMD sweep: baseline (no entrainment) + IHMD in {0, 0.5, 1}.
    with st.spinner("Running entrainment-mixing sweep…"):
        runs = [("no entrainment",
                 cached_run(0, n_ptcl, nt, RH, w, preset, False,
                            switch_turb, eps, 0.0, 0.0))]
        for val in (0.0, 0.5, 1.0):
            runs.append((f"IHMD={val:g}",
                         cached_run(0, n_ptcl, nt, RH, w, preset, False,
                                    switch_turb, eps, 5e-4, val)))
    st.info("Homogeneous (IHMD=0) keeps Nc while droplets shrink; "
            "inhomogeneous (IHMD=1) drops Nc. Compare the Diagnostics tab "
            "(Advanced) or the leftward shift / height of the final DSD.")
    render_tabs(runs, level)

else:
    # single-run lessons (1, 2, 3) and free exploration
    with st.spinner("Running parcel ascent…"):
        out = cached_run(0, n_ptcl, nt, RH, w, preset, collisions,
                         switch_turb, eps, lambda_ent, ihmd)
    render_tabs([(preset, out)], level)
