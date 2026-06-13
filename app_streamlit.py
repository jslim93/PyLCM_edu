"""PyLCM Educational Streamlit GUI (rich edition).

A browser presentation layer over the validated SoA engine
(`PyLCM.timestep_soa.run_soa`). No new physics — every number comes from the
engine. The visual design deliberately mirrors the author's own notebooks:

  * the 3×2 time-series panel of `PyLCM/animation.py` (RH, q_v, z, T, the three
    mixing ratios qa/qc/qr, and the three number concentrations na/nc/nr on a
    log axis), and
  * the DSD views of `Post_process/print_plot.py` — a filled time×radius
    contour of dN with the mean-radius line overlaid, plus spectra lines at
    incremental heights coloured by a viridis time colormap.

It also exposes the full parameter set: a per-mode aerosol editor (1–4 modes
with N, mean radius, σ, κ each) with an optional GCCN coarse mode, the three
ascending modes, collisions / turbulent kernel, and entrainment (λ + IHMD).

Run:  streamlit run app_streamlit.py
"""

import numpy as np
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from PyLCM.timestep_soa import run_soa
from PyLCM.parameters import rho_liq, pi

st.set_page_config(page_title="PyLCM Educational Model", page_icon="☁️",
                   layout="wide")

# Aerosol presets: per-mode (N cm^-3, mean radius µm, geom. std) + scalar kappa
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

# A coarse "giant CCN" mode appended when the GCCN checkbox is on.
GCCN_MODE = dict(N=0.01, r=2.0, sig=1.5, kappa=1.2)

VIRIDIS = "viridis"


# --- engine wrapper (cached, hashable args) -----------------------------------

@st.cache_data(show_spinner=False)
def cached_run(seed, n_ptcl, nt, dt, T0, P0, RH, w, ascending_mode,
               N_raw, mu_um, sig, kappa, collisions, switch_turb, eps,
               lambda_ent, ihmd, n_collect=120):
    """Cached wrapper around run_soa. All tuple args are passed as tuples
    (hashable) so identical configurations return instantly. Returns BOTH the
    dense time-sampled diagnostics dict and the final (M, A) arrays."""
    step = max(1, nt // n_collect)
    collect = tuple(range(step, nt + 1, step))
    out, (M, A) = run_soa(
        seed=seed, n_ptcl=n_ptcl, nt=nt, dt=dt, T0=T0, P0=P0, RH=RH, w=w,
        N_raw=N_raw, mu_um=mu_um, sig=sig, kappa=kappa,
        ascending_mode=ascending_mode, collisions=collisions,
        switch_turb=switch_turb, eps=eps, lambda_ent=lambda_ent, ihmd=ihmd,
        collect=collect)
    return out, np.asarray(M), np.asarray(A)


# --- data shaping helpers ------------------------------------------------------

def _series(out, key):
    ts = sorted(out)
    return np.array(ts), np.array([out[t][key] for t in ts])


def _dsd_stack(out):
    """Stack dsd_n over collected times into a 2-D array (n_time, n_bins).
    Returns (times, radius_µm, dN_array, mean_radius_µm)."""
    ts = sorted(out)
    radii = out[ts[0]]["dsd_r"] * 1e6          # µm, bin centers
    stack = np.array([out[t]["dsd_n"] for t in ts])
    rv = np.array([out[t]["rv"] for t in ts])
    return np.array(ts), radii, stack, rv


# --- figure builders (pure: out-dict[s] -> Plotly figure) ----------------------

def timeseries_figure(runs, dt):
    """The author's 3×2 panel from animation.py, generalised to overlay runs:
       (1,1) RH %    (1,2) q_v g/kg    (2,1) z m    (2,2) T K
       (3,1) qa/qc/qr g/kg (3 lines)   (3,2) na/nc/nr cm^-3 (3 lines, log y)."""
    fig = make_subplots(
        rows=3, cols=2,
        subplot_titles=("Relative Humidity RH (%)", "Vapour q<sub>v</sub> (g/kg)",
                        "Height z (m)", "Temperature T (K)",
                        "Mixing ratios q<sub>x</sub> (g/kg)",
                        "Number conc. n<sub>x</sub> (cm⁻³)"))
    multi = len(runs) > 1
    dashes = ["solid", "dash", "dot", "dashdot"]
    for j, (label, out, *_rest) in enumerate(runs):
        dash = dashes[j % len(dashes)]
        pre = f"{label} " if multi else ""
        ts, _ = _series(out, "RH")
        tsec = ts * dt
        def add(key, row, col, color, name, scale=1.0):
            _, y = _series(out, key)
            fig.add_trace(go.Scatter(
                x=tsec, y=y * scale, mode="lines", name=pre + name,
                line=dict(color=None if multi else color, dash=dash),
                legendgroup=label, showlegend=(row == 3)),
                row=row, col=col)
        add("RH", 1, 1, "lightblue", "RH", 100.0)
        add("qv", 1, 2, "green", "q_v")
        add("z", 2, 1, "black", "z")
        add("T_K", 2, 2, "red", "T")
        add("qa", 3, 1, "blue", "q_a (aerosol)")
        add("qc", 3, 1, "orange", "q_c (cloud)")
        add("qr", 3, 1, "green", "q_r (rain)")
        add("NA", 3, 2, "blue", "n_a (aerosol)")
        add("NC", 3, 2, "orange", "n_c (cloud)")
        add("NR", 3, 2, "green", "n_r (rain)")
    fig.update_yaxes(type="log", row=3, col=2)
    fig.update_xaxes(title_text="Time (s)", row=3, col=1)
    fig.update_xaxes(title_text="Time (s)", row=3, col=2)
    fig.update_layout(height=820, legend=dict(orientation="h", y=-0.08),
                      margin=dict(t=40))
    return fig


def dsd_contour_figure(out, dt):
    """Filled time×radius contour of dN (log color) with the mean-radius line
    overlaid — the print_plot.spec_plot design, in Plotly."""
    ts, radii, stack, rv = _dsd_stack(out)
    tsec = ts * dt
    z = stack.T.copy()                          # (n_radius, n_time)
    z[z <= 0] = np.nan
    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        x=tsec, y=radii, z=np.log10(z), colorscale=VIRIDIS,
        colorbar=dict(title="log₁₀ dN<br>(cm⁻³)"),
        hovertemplate="t=%{x:.0f}s<br>r=%{y:.2f}µm<br>log₁₀dN=%{z:.2f}<extra></extra>"))
    fig.add_trace(go.Scatter(
        x=tsec, y=rv, mode="lines", line=dict(color="black", width=2.5),
        name="mean radius r̄ (µm)"))
    fig.update_yaxes(type="log", title_text="Radius r (µm)")
    fig.update_xaxes(title_text="Time (s)")
    fig.update_layout(height=520, title="DSD time evolution",
                      legend=dict(orientation="h", y=1.05))
    return fig


def dsd_spectra_figure(out, dt):
    """Spectra lines at incremental heights, coloured by a viridis time colormap
    — the print_plot subplot_array_function lower-right panel."""
    ts, radii, stack, _ = _dsd_stack(out)
    z = stack.copy()
    z[z <= 0] = np.nan
    n = len(ts)
    fig = go.Figure()
    # subsample to ~12 spectra to keep the plot readable
    idx = np.linspace(0, n - 1, min(12, n)).astype(int)
    for i in idx:
        frac = i / max(1, n - 1)
        # sample viridis via Plotly's built-in colorscale sampler
        color = _viridis_at(frac)
        fig.add_trace(go.Scatter(
            x=radii, y=z[i], mode="lines",
            line=dict(color=color),
            name=f"t={ts[i] * dt:.0f}s", showlegend=False,
            hovertemplate="r=%{x:.2f}µm<br>dN=%{y:.3g}<extra></extra>"))
    # a dummy colorbar trace for the time legend
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode="markers",
        marker=dict(colorscale=VIRIDIS, cmin=0, cmax=ts[-1] * dt,
                    color=[0], colorbar=dict(title="Time (s)")),
        showlegend=False))
    fig.update_xaxes(type="log", title_text="Radius r (µm)")
    fig.update_yaxes(type="log", title_text="dN (cm⁻³)")
    fig.update_layout(height=520, title="DSD spectra coloured by time")
    return fig


def _viridis_at(frac):
    import plotly.colors as pc
    return pc.sample_colorscale(VIRIDIS, [float(np.clip(frac, 0, 1))])[0]


def particle_figure(M, A):
    """Final super-droplet population by radius. Each super-droplet is one
    point; marker size / opacity ∝ its multiplicity A. A weighted histogram
    (number concentration per log-radius bin) is overlaid for context."""
    m = A > 0
    r = np.zeros_like(M)
    r[m] = (M[m] / (A[m] * 4.0 / 3.0 * pi * rho_liq)) ** (1.0 / 3.0)
    r_um = r[m] * 1e6
    Aw = A[m]
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    # scatter: each super-droplet, size by multiplicity
    if r_um.size:
        s = 4 + 16 * (np.log10(Aw + 1) / max(1.0, np.log10(Aw.max() + 1)))
        fig.add_trace(go.Scattergl(
            x=r_um, y=Aw, mode="markers",
            marker=dict(size=s, color=np.log10(Aw + 1), colorscale=VIRIDIS,
                        opacity=0.5, colorbar=dict(title="log₁₀ A")),
            name="super-droplets",
            hovertemplate="r=%{x:.3f}µm<br>A=%{y:.3g}<extra></extra>"),
            secondary_y=False)
    # weighted histogram of multiplicity vs radius
    if r_um.size:
        edges = np.logspace(np.log10(max(1e-3, r_um.min())),
                            np.log10(r_um.max() + 1e-9), 40)
        hist, e = np.histogram(r_um, bins=edges, weights=Aw)
        centers = np.sqrt(e[:-1] * e[1:])
        fig.add_trace(go.Bar(
            x=centers, y=hist, name="Σ A per bin", opacity=0.35,
            marker_color="steelblue"), secondary_y=True)
    fig.update_xaxes(type="log", title_text="Droplet radius r (µm)")
    fig.update_yaxes(type="log", title_text="multiplicity A (per super-droplet)",
                     secondary_y=False)
    fig.update_yaxes(title_text="Σ A in bin", secondary_y=True)
    fig.update_layout(height=520, title="Final super-droplet distribution",
                      legend=dict(orientation="h", y=1.05))
    return fig


def vertical_profile_figure(runs):
    """T vs z and LWC (qc+qr) vs z."""
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("Temperature T (°C) vs height",
                                        "LWC q_c+q_r (g/kg) vs height"))
    multi = len(runs) > 1
    dashes = ["solid", "dash", "dot", "dashdot"]
    for j, (label, out, *_rest) in enumerate(runs):
        dash = dashes[j % len(dashes)]
        _, z = _series(out, "z")
        _, T = _series(out, "T")
        _, qc = _series(out, "qc")
        _, qr = _series(out, "qr")
        nm = label if multi else None
        fig.add_trace(go.Scatter(x=T, y=z, mode="lines", line=dict(dash=dash),
                                 name=(nm or "T"), legendgroup=label,
                                 showlegend=multi), row=1, col=1)
        fig.add_trace(go.Scatter(x=qc + qr, y=z, mode="lines",
                                 line=dict(dash=dash), name=(nm or "LWC"),
                                 legendgroup=label, showlegend=False),
                      row=1, col=2)
    fig.update_yaxes(title_text="Height z (m)", row=1, col=1)
    fig.update_xaxes(title_text="T (°C)", row=1, col=1)
    fig.update_xaxes(title_text="q_c + q_r (g/kg)", row=1, col=2)
    fig.update_layout(height=520)
    return fig


# --- tab rendering -------------------------------------------------------------

def render_tabs(runs, level, dt):
    """runs: list of (label, out, M, A). Multi-run overlays share the time
    series / profiles; single-run also gets DSD contour + particle view."""
    single = len(runs) == 1
    labels = ["Time series", "DSD", "Particle distribution", "Vertical profiles"]
    if level == "Advanced":
        labels.append("Diagnostics")
    tabs = st.tabs(labels)

    with tabs[0]:
        st.caption("Six-panel parcel evolution (after animation.py): RH, vapour, "
                   "height, temperature, the three mixing ratios, and the three "
                   "number concentrations (log axis).")
        st.plotly_chart(timeseries_figure(runs, dt), use_container_width=True)

    with tabs[1]:
        st.caption("Droplet size distribution. Left: time×radius contour of dN "
                   "with the mean-radius line. Right: spectra at incremental "
                   "heights, coloured by time. Growth narrows the spectrum; "
                   "collisions broaden it toward rain sizes.")
        if single:
            label, out = runs[0][0], runs[0][1]
            c1, c2 = st.columns(2)
            with c1:
                st.plotly_chart(dsd_contour_figure(out, dt),
                                use_container_width=True)
            with c2:
                st.plotly_chart(dsd_spectra_figure(out, dt),
                                use_container_width=True)
        else:
            # overlay final spectra of each run
            fig = go.Figure()
            for label, out, *_ in runs:
                ts, radii, stack, _ = _dsd_stack(out)
                y = stack[-1].copy(); y[y <= 0] = np.nan
                fig.add_trace(go.Scatter(x=radii, y=y, mode="lines",
                                         name=label))
            fig.update_xaxes(type="log", title_text="Radius r (µm)")
            fig.update_yaxes(type="log", title_text="dN (cm⁻³)")
            fig.update_layout(height=520, title="Final DSD overlay")
            st.plotly_chart(fig, use_container_width=True)

    with tabs[2]:
        st.caption("Final super-droplet population. Each marker is one "
                   "super-droplet at its radius; size/colour scale with its "
                   "multiplicity A. The bars sum multiplicity per radius bin.")
        if single:
            _, _, M, A = runs[0]
            st.plotly_chart(particle_figure(M, A), use_container_width=True)
        else:
            cols = st.columns(len(runs))
            for col, (label, _, M, A) in zip(cols, runs):
                with col:
                    st.markdown(f"**{label}**")
                    st.plotly_chart(particle_figure(M, A),
                                    use_container_width=True)

    with tabs[3]:
        st.caption("Vertical structure of the ascent.")
        st.plotly_chart(vertical_profile_figure(runs), use_container_width=True)

    if level == "Advanced":
        with tabs[4]:
            cols = st.columns(len(runs))
            for c, (label, out, M, A) in zip(cols, runs):
                last = out[sorted(out)[-1]]
                c.markdown(f"**{label}**" if len(runs) > 1 else "**Final state**")
                c.metric("Nc (cm⁻³)", f"{last['NC']:.2f}")
                c.metric("Nr (cm⁻³)", f"{last['NR']:.3f}")
                c.metric("Na residual (cm⁻³)", f"{last['NA']:.2f}")
                c.metric("Supersaturation (%)", f"{(last['RH'] - 1) * 100:+.3f}")
                c.metric("LWC qc+qr (g/kg)", f"{last['qc'] + last['qr']:.3f}")
                c.metric("Mean radius (µm)", f"{last['rv']:.2f}")


# --- sidebar -------------------------------------------------------------------

st.sidebar.title("☁️ PyLCM")
mode = st.sidebar.radio("Mode", ["Guided experiment", "Free exploration"],
                        help="Guided walks through five lessons with explanations; "
                             "Free exposes the full parameter set.")
level = st.sidebar.radio("Level", ["Beginner", "Advanced"], horizontal=True,
                         help="Advanced adds equations, numerics controls, the "
                              "turbulent kernel, and a diagnostics tab.")

LESSONS = ["1 · Aerosol activation", "2 · Condensational growth",
           "3 · Collision → rain", "4 · Maritime vs continental",
           "5 · Entrainment mixing (IHMD)"]
lesson = (st.sidebar.selectbox("Experiment", LESSONS)
          if mode == "Guided experiment" else None)

# --- Parcel ---
st.sidebar.subheader("Parcel")
T0 = st.sidebar.slider("Initial T₀ (K)", 273.0, 303.0, 293.2, 0.1,
                       help="Cloud-base temperature of the rising parcel.")
P0 = st.sidebar.slider("Initial P₀ (hPa)", 700.0, 1030.0, 1013.0, 1.0,
                       help="Cloud-base pressure.") * 100.0
RH = st.sidebar.slider("Initial RH", 0.80, 1.0, 0.92, 0.01,
                       help="Higher RH → cloud forms sooner.")
w = st.sidebar.slider("Updraft w (m/s)", 0.1, 5.0, 1.0, 0.1,
                      help="Stronger updraft → higher supersaturation → more "
                           "droplets activate.")
ascending_mode = st.sidebar.selectbox(
    "Ascending mode", ["linear", "sine", "in_cloud_oscillation"],
    help="linear: constant w. sine: smooth acceleration/deceleration. "
         "in_cloud_oscillation: rises then oscillates up/down inside the cloud.")

if level == "Advanced":
    nt = int(st.sidebar.number_input("Steps (nt)", 300, 3600, 1500, 100,
                                     help="Number of time steps in the ascent."))
    dt = float(st.sidebar.number_input("dt (s)", 0.25, 5.0, 1.0, 0.25,
                                       help="Time-step length."))
    n_ptcl = int(st.sidebar.number_input("Super-droplets", 500, 20000, 2000, 500,
                                         help="Number of Lagrangian "
                                              "super-droplets."))
else:
    nt, dt, n_ptcl = 1500, 1.0, 2000

# --- Aerosol ---
st.sidebar.subheader("Aerosol")
gccn = st.sidebar.checkbox(
    "Add GCCN (giant CCN)", False,
    help="Append a coarse mode (N≈0.01 cm⁻³, r≈2 µm, κ≈1.2). A few giant, "
         "highly hygroscopic nuclei seed early raindrops.")

if level == "Advanced" and mode == "Free exploration":
    # full multi-mode editor
    n_modes = int(st.sidebar.number_input("Number of modes", 1, 4, 3, 1,
                                          help="Lognormal aerosol modes."))
    base = PRESETS["default"]
    defaults = [
        dict(N=118., r=.019, sig=3.3, kappa=1.6),
        dict(N=11., r=.056, sig=1.6, kappa=1.6),
        dict(N=.72, r=.46, sig=2.2, kappa=1.6),
        dict(N=.1, r=1.0, sig=1.5, kappa=1.2),
    ]
    N_list, mu_list, sig_list, k_list = [], [], [], []
    with st.sidebar.expander("Configure aerosol modes", expanded=True):
        mode_tabs = st.tabs([f"Mode {i + 1}" for i in range(n_modes)])
        names = ["Aitken", "Accumulation", "Coarse", "Extra"]
        for i, mtab in enumerate(mode_tabs):
            with mtab:
                st.markdown(f"**Mode {i + 1}** ({names[i]})")
                N = st.number_input("N (cm⁻³)", value=defaults[i]["N"],
                                    step=1.0, key=f"N{i}",
                                    help="Number concentration of this mode.")
                r = st.number_input("Mean radius (µm)", value=defaults[i]["r"],
                                    format="%.3f", step=0.001, key=f"r{i}",
                                    help="Geometric mean radius.")
                sg = st.number_input("σ (geom. std)", value=defaults[i]["sig"],
                                     step=0.1, key=f"s{i}",
                                     help="Geometric standard deviation (width).")
                kp = st.number_input("κ (hygroscopicity)",
                                     value=defaults[i]["kappa"], step=0.1,
                                     key=f"k{i}",
                                     help="κ: sea salt≈1.2, sulfate≈0.6, "
                                          "organic≈0.1–0.3.")
                N_list.append(N); mu_list.append(r)
                sig_list.append(sg); k_list.append(kp)
    N_raw = tuple(N_list); mu_um = tuple(mu_list)
    sig = tuple(sig_list); kappa = tuple(k_list)
    preset = "custom"
else:
    preset = st.sidebar.selectbox(
        "Preset", list(PRESETS),
        help="default: 3-mode continental-ish. maritime: few large CCN → fast "
             "rain. continental: many small CCN → suppressed rain. arctic: "
             "ultra-clean.")
    p = PRESETS[preset]
    N_raw, mu_um, sig, kappa = p["N_raw"], p["mu_um"], p["sig"], p["kappa"]

# Append GCCN coarse mode (kappa becomes a per-mode tuple if it was scalar)
if gccn:
    N_raw = tuple(N_raw) + (GCCN_MODE["N"],)
    mu_um = tuple(mu_um) + (GCCN_MODE["r"],)
    sig = tuple(sig) + (GCCN_MODE["sig"],)
    if np.isscalar(kappa):
        kappa = (kappa,) * (len(N_raw) - 1) + (GCCN_MODE["kappa"],)
    else:
        kappa = tuple(kappa) + (GCCN_MODE["kappa"],)

# --- Physics ---
st.sidebar.subheader("Physics")
collisions = st.sidebar.checkbox("Collision–coalescence", True,
                                 help="Required for rain to form.")
if level == "Advanced":
    switch_turb = st.sidebar.checkbox("Turbulent collision kernel", False,
                                      help="Wang & Ayala turbulence enhancement.")
    eps = (st.sidebar.number_input("ε (m²/s³)", 0.0, 0.5, 0.01, 0.01,
                                   help="Turbulent dissipation rate.")
           if switch_turb else 0.0)
else:
    switch_turb, eps = False, 0.0

lambda_ent = st.sidebar.slider(
    "Entrainment λ", 0.0, 2e-3, 0.0, 1e-4, format="%.4f",
    help="Entrainment rate. >0 mixes in dry environmental air, evaporating "
         "cloud water (used in lesson 5 and free mode).")
ihmd = st.sidebar.slider(
    "IHMD (mixing degree)", 0.0, 1.0, 0.5, 0.1,
    help="0 = homogeneous (all droplets shrink, number kept); "
         "1 = inhomogeneous (whole droplets evaporate, survivors keep size). "
         "N_c/N_c0 = (q_c/q_c0)^IHMD.")


# --- lesson copy ---------------------------------------------------------------

EXPLAIN = {
    LESSONS[0]: (
        "**Aerosol activation (Köhler theory).** As the parcel rises it cools "
        "and supersaturation grows. Once it crosses each particle's critical "
        "value, that haze particle *activates* into a cloud droplet. Raise the "
        "updraft or RH and watch the cloud-droplet number **Nc** climb.",
        r"S_{crit} \propto \sqrt{A^3 / B} \quad\text{(Köhler curve peak)}"),
    LESSONS[1]: (
        "**Condensational growth.** With collisions OFF, droplets grow only by "
        "vapour diffusion. Small droplets grow faster than large ones "
        "(dr/dt ∝ 1/r), so the spectrum **narrows** with height.",
        r"\frac{dr}{dt} = \frac{G\,S}{r}"),
    LESSONS[2]: (
        "**Collision → rain.** With collisions ON (maritime aerosol), larger "
        "droplets fall faster and collect smaller ones. The distribution "
        "**broadens** and a rain mode (Nr, qr) develops.",
        r"K(r_1,r_2) = \pi (r_1+r_2)^2\, |v_1-v_2|\, E"),
    LESSONS[3]: (
        "**Maritime vs continental (A/B).** Same updraft, different aerosol: "
        "maritime air has *few large* CCN that rain quickly, while continental "
        "air has *many small* CCN that make many tiny droplets and **suppress** "
        "rain at the same water content.",
        r"N_c \uparrow \;\Rightarrow\; \bar r \downarrow \;\Rightarrow\; "
        r"\text{rain suppressed}"),
    LESSONS[4]: (
        "**Entrainment mixing (IHMD).** Dry environmental air is entrained and "
        "evaporates cloud water. *Homogeneous* mixing (IHMD=0) shrinks every "
        "droplet but keeps the number; *inhomogeneous* mixing (IHMD=1) "
        "evaporates whole droplets but keeps the survivors' size. We run "
        "IHMD ∈ {0, 0.5, 1} and overlay the spectra. (Lim & Hoffmann, 2023.)",
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
        p = PRESETS["maritime"]
        N_raw, mu_um, sig, kappa = p["N_raw"], p["mu_um"], p["sig"], p["kappa"]
        if gccn:
            N_raw += (GCCN_MODE["N"],); mu_um += (GCCN_MODE["r"],)
            sig += (GCCN_MODE["sig"],)
            kappa = (kappa,) * len(p["N_raw"]) + (GCCN_MODE["kappa"],)
else:
    st.header("Free exploration")
    st.markdown("Adjust any sidebar parameter; the model re-runs instantly "
                "(cached). Time series + DSD are the central diagnostics; the "
                "particle tab shows the raw super-droplet population.")


if lesson == LESSONS[3]:
    # A/B: run both presets, overlay.
    with st.spinner("Running maritime and continental parcels…"):
        pm, pc = PRESETS["maritime"], PRESETS["continental"]
        om, Mm, Am = cached_run(0, n_ptcl, nt, dt, T0, P0, RH, w, ascending_mode,
                                pm["N_raw"], pm["mu_um"], pm["sig"], pm["kappa"],
                                True, switch_turb, eps, 0.0, 0.0)
        oc, Mc, Ac = cached_run(0, n_ptcl, nt, dt, T0, P0, RH, w, ascending_mode,
                                pc["N_raw"], pc["mu_um"], pc["sig"], pc["kappa"],
                                True, switch_turb, eps, 0.0, 0.0)
    render_tabs([("maritime", om, Mm, Am), ("continental", oc, Mc, Ac)],
                level, dt)

elif lesson == LESSONS[4]:
    # IHMD sweep: baseline + IHMD in {0, 0.5, 1}.
    with st.spinner("Running entrainment-mixing sweep…"):
        runs = []
        ob, Mb, Ab = cached_run(0, n_ptcl, nt, dt, T0, P0, RH, w, ascending_mode,
                                N_raw, mu_um, sig, kappa, False, switch_turb,
                                eps, 0.0, 0.0)
        runs.append(("no entrainment", ob, Mb, Ab))
        for val in (0.0, 0.5, 1.0):
            ov, Mv, Av = cached_run(0, n_ptcl, nt, dt, T0, P0, RH, w, ascending_mode,
                                    N_raw, mu_um, sig, kappa, False, switch_turb,
                                    eps, 5e-4, val)
            runs.append((f"IHMD={val:g}", ov, Mv, Av))
    st.info("Homogeneous (IHMD=0) keeps Nc while droplets shrink; inhomogeneous "
            "(IHMD=1) drops Nc. Compare the Diagnostics tab (Advanced) and the "
            "final-DSD overlay.")
    render_tabs(runs, level, dt)

else:
    # single-run lessons (1, 2, 3) and free exploration
    n_eff = nt
    with st.spinner("Running parcel ascent…"):
        out, M, A = cached_run(0, n_eff, nt, dt, T0, P0, RH, w, ascending_mode,
                               tuple(N_raw), tuple(mu_um), tuple(sig),
                               kappa if np.isscalar(kappa) else tuple(kappa),
                               collisions, switch_turb, eps, lambda_ent, ihmd)
    render_tabs([(preset, out, M, A)], level, dt)
