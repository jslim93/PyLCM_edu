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

PyLCM is benchmarked against the original Fortran / SAM6-LCM reference — bulk
thermodynamics agree within ~0.24 °C and liquid water within ~1.5%. Details are in
`validation/collision_validation.ipynb`, and physics invariants (mass conservation,
positivity, no-NaN, integer-multiplicity collisions) run in CI under `tests/`.

## Command-line runs (no notebook)

Edit the parameters in **`input.yaml`** and run the model from the terminal — the
same idea as the Fortran model's `input.nml`:

```bash
python run.py                      # single run, reads input.yaml
python run.py myconfig.yaml        # use a different config
python run.py input.yaml --ensemble 20 --jobs -1   # 20-member parallel ensemble
```

- **Single run** writes a multi-column time-series CSV (`time, z, T, qc, qr, Nc, Nr,
  Na, LWC`) to `output.file` and prints a summary.
- **`--ensemble N`** runs `N` stochastic members and writes the mean + 10–90
  percentile envelope of the chosen `output.diagnostic`.
- **`--jobs J`** sets how many CPU cores the ensemble uses (joblib): `-1` (default)
  uses **all** cores, `1` runs serially, `4` uses four. Only matters with `--ensemble`.

`input.yaml` exposes the parcel state, aerosol preset, collision/turbulence switches,
and entrainment (λ, IHMD). The notebooks remain the primary interactive interface.

The underlying library functions live in the package: `PyLCM.pylcm_run.run_single_series`
(one ascent → diagnostic series) and `PyLCM.ensemble.run_ensemble` (parallel members).

## Performance

The condensation hot loop runs on a struct-of-arrays numba kernel
(`PyLCM/condensation_fast.py`), **5.4× faster than the object loop and
bit-for-bit identical** to it (locked by a golden regression test). The collision
helpers `E_H80` and `ws_drops_beard` are JIT-compiled, and the particle shuffle
uses `np.random.permutation`. Together these give **~3× faster full runs**
(condensation + collision). `ensemble.py` then runs members **in parallel across
cores**, so a K-member ensemble is roughly 3K× faster than serial object-based
runs. Profile and benchmarks: `validation/PROFILE.md`.

## Entrainment mixing (IHMD)

Warm-cloud entrainment mixing is parameterized by two controls (`PyLCM/mixing.py`):
an entrainment rate **λ** and the **Inhomogeneous Mixing Degree (IHMD)** of
Lim & Hoffmann (2023). Mixing runs before condensation each step, diluting the
parcel toward the environment and redistributing cloud liquid so that

> **N_c / N_{c,0} = (q_c / q_{c,0})^IHMD**

IHMD=0 is homogeneous (every droplet shrinks, number conserved); IHMD=1 is
inhomogeneous (a subset evaporates entirely, survivors keep size). The notebook
`validation/entrainment_mixing.ipynb` sweeps IHMD and renders the DSD contrast and
the N–r_v mixing diagram. A `MixingModel` interface leaves a seam for a future
**LEM** backend (`LEMMixing`, Phase 3b); ice/mixed-phase mixing is Phase 4.

## Known limitations

- The legacy `PyLCM.entrainment.basic_entrainment(...)` is superseded by
  `ParameterizedMixing` (the homogeneous limit, IHMD=0) and remains hard-gated
  behind `experimental=True`.
- `LEMMixing` is an interface stub (raises `NotImplementedError`) pending Phase 3b.

## How to cite

> J. Lim, *PyLCM v1.0* (2026). Educational Lagrangian Cloud Model.

## Contact

J.lim@physik.uni-muenchen.de

## License

MIT — see [`LICENSE`](LICENSE).
