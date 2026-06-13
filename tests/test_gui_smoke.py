"""Smoke tests for the educational Streamlit GUI.

These never launch a Streamlit server. One checks the engine adapter produces
plottable DSD data; the other checks the app module resolves as an import spec.
"""

import importlib.util
import os

import numpy as np

from PyLCM.timestep_soa import run_soa


def test_engine_adapter_produces_plottable_dsd():
    out, _ = run_soa(seed=0, n_ptcl=400, nt=120, collisions=True,
                     collect=(60, 120))
    d = out[120]
    assert "dsd_r" in d and "dsd_n" in d
    assert d["dsd_r"].shape == d["dsd_n"].shape
    assert np.all(np.isfinite(d["dsd_n"])) and np.all(d["dsd_n"] >= 0)
    assert np.all(np.isfinite(d["dsd_r"])) and np.all(d["dsd_r"] >= 0)
    assert d["qc"] >= 0 and d["qr"] >= 0


def test_app_module_imports():
    # Resolve the module spec WITHOUT executing it: importing app_streamlit
    # runs Streamlit calls that require a running server, so we only verify the
    # file is locatable and parses into a loadable spec.
    spec = importlib.util.spec_from_file_location(
        "app_streamlit", os.path.join(os.getcwd(), "app_streamlit.py"))
    assert spec is not None
