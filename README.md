# PyLCM — Educational Lagrangian Cloud Model

![CI](https://github.com/jslim93/PyLCM_edu/actions/workflows/ci.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A teaching parcel-model simulator for warm-cloud microphysics using Lagrangian
super-droplets. It integrates an ascending air parcel and resolves droplet
**condensation/evaporation** (Köhler theory) and **collision–coalescence**
(stochastic super-droplet method), so students can watch a cloud droplet
spectrum form and rain develop from first principles.

## Install (one command)

With conda (recommended):

```bash
conda env create -f environment.yml && conda activate PyLCM && pip install -e .
```

Or with pip into an existing Python ≥3.11 environment:

```bash
pip install -e . && pip install -r requirements.txt
```

> **Note on versions:** `numpy` and `numba` are version-coupled (numba 0.65+
> supports numpy up to 2.4). The pins in `environment.yml` / `requirements.txt`
> keep this pairing consistent — don't loosen them independently.

## Quickstart

```bash
jupyter notebook PyLCM_edu.ipynb
```

Guided teaching notebooks:

- `PyLCM_Part1_Foundations.ipynb` — the physics, module by module.
- `PyLCM_Part2_Experiments.ipynb` — hands-on numerical experiments.

## Scientific basis

- **Köhler theory** for activation and condensational growth.
- **Hall (1980)** gravitational collision-efficiency table.
- **Straub et al. (2009)** coalescence efficiency `E_S09`.
- **Beard (1976)** droplet terminal velocity.
- **Shima et al. (2009)** Linear Sampling Method for stochastic collisions.
- Optional **Ayala et al. (2008) / Wang & Grabowski (2009)** turbulent collision kernel.

## Validation

`validation/collision_validation.ipynb` answers the common question — *are the
collision results correct versus the original Fortran / SAM6-LCM?* In short:
**yes.** PyLCM shares its core collision algorithm (Linear Sampling, Hall 1980,
Beard 1976, multi-collision limiter) with SAM6-LCM and the Fortran box model. The
one substantive difference is that PyLCM applies the Straub (2009) coalescence
efficiency `E_S09` in **both** the gravitational and turbulent kernels, whereas
SAM6-LCM omits it. PyLCM therefore **agrees with the Fortran box-model reference**
and is *more* physically complete than SAM6-LCM in the turbulent regime. Physics
invariants (mass conservation, positivity, numerical stability) are checked in CI
under `tests/`.

## Ensemble runs

`ensemble.py` runs `N` independent stochastic members **in parallel across CPU
cores** (joblib) and returns the mean plus a 10–90 percentile envelope — useful
for separating physical signal from Monte-Carlo noise. See its module docstring.

## Performance

The single-run hot path uses numba JIT on the inner numeric kernels (`esatw`,
`sigma_air_liq`, `radius_liquid_euler`). A baseline profile lives in
`validation/PROFILE.md`; further single-run speedup requires a struct-of-arrays
refactor of the per-particle loop, planned for a future release. For multi-run
studies the parallel ensemble is the main lever today.

## Known limitations

- **Entrainment is EXPERIMENTAL and not physically validated.** It is hard-gated:
  calling `PyLCM.entrainment.basic_entrainment(...)` raises unless you pass
  `experimental=True`. A validated entrainment scheme is planned for **v1.1**.

## How to cite

> J. Lim, *PyLCM v1.0* (2026). Educational Lagrangian Cloud Model.

## Contact

J.lim@physik.uni-muenchen.de

## License

MIT — see [`LICENSE`](LICENSE).
