#!/usr/bin/env python3
"""Command-line runner for PyLCM, driven by a YAML config (input.yaml).

The Jupyter notebooks remain the primary interface; this is for scripted/batch runs.

Usage:
    python run_pylcm.py [config.yaml]            # single run -> time-series CSV
    python run_pylcm.py input.yaml --ensemble 10 # N stochastic members -> envelope CSV
    python run_pylcm.py input.yaml --out run1.csv --jobs -1

Single run  -> multi-column CSV: time, z, T, qc, qr, Nc, Nr, Na, LWC (like the
               Fortran time_series.dat).
Ensemble    -> mean + 10-90 percentile envelope of the chosen `diagnostic`.
"""
import argparse
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")  # create_env_profiles() calls plt.show(); keep it headless

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required for the CLI:  pip install pyyaml")

from validation.phys_harness import run, PRESETS
from PyLCM.parameters import r_a, rv
from PyLCM.condensation import esatw
from PyLCM.parcel import create_env_profiles
from PyLCM.mixing import ParameterizedMixing


DEFAULTS = {
    "parcel": {"T0": 293.2, "P0": 101300.0, "RH": 0.92, "w": 1.0, "nt": 3000, "dt": 1.0},
    "aerosol": {"n_ptcl": 2000, "preset": "default"},
    "physics": {"collisions": True, "turbulent": False, "epsilon": 0.0},
    "entrainment": {"enabled": False, "lambda": 0.0, "ihmd": 0.0},
    "output": {"collect_every": 100, "diagnostic": "LWC", "file": "output.csv"},
}


def load_config(path):
    try:
        with open(path) as f:
            user = yaml.safe_load(f) or {}
    except FileNotFoundError:
        sys.exit(f"config file not found: {path}")
    merged = {sect: dict(DEFAULTS[sect]) for sect in DEFAULTS}
    for sect, vals in user.items():
        merged.setdefault(sect, {})
        merged[sect].update(vals or {})
    return merged


def build_mixing(cfg):
    ent = cfg["entrainment"]
    if not ent.get("enabled", False):
        return None
    p = cfg["parcel"]
    es = esatw(p["T0"])
    q0 = p["RH"] * es / (p["P0"] - p["RH"] * es) * r_a / rv
    qv, th, z_env = create_env_profiles(p["T0"], q0, 0.0, p["P0"], "Stable")
    return ParameterizedMixing(ent.get("lambda", 0.0), ent.get("ihmd", 0.0), qv, th, z_env)


def run_kwargs(cfg, mixing, collect, seed=0):
    p, a, phys = cfg["parcel"], cfg["aerosol"], cfg["physics"]
    if a["preset"] not in PRESETS:
        sys.exit(f"unknown aerosol preset '{a['preset']}'; choose from {list(PRESETS)}")
    return dict(seed=seed, n_ptcl=a["n_ptcl"], nt=p["nt"], dt=p["dt"],
                T0=p["T0"], P0=p["P0"], RH=p["RH"], w=p["w"],
                aerosol=a["preset"], collisions=phys["collisions"],
                switch_turb=phys["turbulent"], eps=phys["epsilon"],
                mixing=mixing, collect=collect)


def get_diag(d, name):
    table = {"LWC": d["qc"] + d["qr"], "qc": d["qc"], "qr": d["qr"],
             "Nc": d["NC"], "Nr": d["NR"], "Na": d["NA"]}
    if name not in table:
        sys.exit(f"unknown diagnostic '{name}'; choose from {list(table)}")
    return table[name]


def write_timeseries_csv(diag, path):
    cols = ["time_s", "z_m", "T_C", "qc_gkg", "qr_gkg", "Nc_cm3", "Nr_cm3", "Na_cm3", "LWC_gkg"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for t in sorted(diag):
            d = diag[t]
            w.writerow([t, f"{d['z']:.1f}", f"{d['T']:.3f}", f"{d['qc']:.5f}",
                        f"{d['qr']:.5f}", f"{d['NC']:.3f}", f"{d['NR']:.3f}",
                        f"{d['NA']:.3f}", f"{d['qc'] + d['qr']:.5f}"])


def main():
    ap = argparse.ArgumentParser(description="Run PyLCM from a YAML config.")
    ap.add_argument("config", nargs="?", default="input.yaml",
                    help="YAML config file (default: input.yaml)")
    ap.add_argument("--ensemble", type=int, default=0, metavar="N",
                    help="run N stochastic ensemble members")
    ap.add_argument("--jobs", type=int, default=-1,
                    help="parallel jobs for the ensemble (-1 = all cores)")
    ap.add_argument("--out", default=None, help="output CSV (overrides config output.file)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    nt = cfg["parcel"]["nt"]
    every = cfg["output"]["collect_every"]
    collect = tuple(range(every, nt + 1, every))
    mixing = build_mixing(cfg)
    out_file = args.out or cfg["output"]["file"]
    diagnostic = cfg["output"]["diagnostic"]

    if args.ensemble > 0:
        from joblib import Parallel, delayed
        times = list(collect)

        def member(i):
            diag, _ = run(**run_kwargs(cfg, mixing, collect, seed=i))
            return np.array([get_diag(diag[t], diagnostic) for t in sorted(diag)])

        members = np.vstack(Parallel(n_jobs=args.jobs)(
            delayed(member)(i) for i in range(args.ensemble)))
        mean = members.mean(0)
        lo = np.percentile(members, 10, 0)
        hi = np.percentile(members, 90, 0)
        ens_file = out_file if out_file.endswith(".csv") else "ensemble.csv"
        with open(ens_file, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["time_s", f"{diagnostic}_mean", f"{diagnostic}_p10", f"{diagnostic}_p90"])
            for k, t in enumerate(times):
                w.writerow([t, f"{mean[k]:.5f}", f"{lo[k]:.5f}", f"{hi[k]:.5f}"])
        print(f"Ensemble: {args.ensemble} members, diagnostic={diagnostic}, jobs={args.jobs}")
        print(f"  final {diagnostic}: mean={mean[-1]:.4f}  [p10={lo[-1]:.4f}, p90={hi[-1]:.4f}]")
        print(f"  envelope saved -> {ens_file}")
    else:
        diag, _ = run(**run_kwargs(cfg, mixing, collect))
        write_timeseries_csv(diag, out_file)
        last = diag[max(diag)]
        a = cfg["aerosol"]
        print(f"Single run: n_ptcl={a['n_ptcl']}, nt={nt}, preset={a['preset']}, "
              f"collisions={cfg['physics']['collisions']}, entrainment={cfg['entrainment']['enabled']}")
        print(f"  final: T={last['T']:.2f}C  qc={last['qc']:.4f}  qr={last['qr']:.4f}  "
              f"Nc={last['NC']:.2f}  Nr={last['NR']:.3f}  LWC={last['qc'] + last['qr']:.4f} g/kg")
        print(f"  time series ({len(diag)} rows) saved -> {out_file}")


if __name__ == "__main__":
    main()
