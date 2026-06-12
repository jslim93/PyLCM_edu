"""Generate the condensation golden snapshot from the CURRENT object-based code.

Run once (and only re-run intentionally if the reference physics legitimately
changes). Writes tests/golden/condensation_golden.npz.
"""
import os
import numpy as np
from validation.golden_setup import run_condensation_only

out = run_condensation_only()
os.makedirs("tests/golden", exist_ok=True)
np.savez("tests/golden/condensation_golden.npz",
         M=out["M"], A=out["A"], T=out["T"], q=out["q"])
print(f"Golden saved: {len(out['M'])} particles, T={out['T']:.6f}, q={out['q']:.8f}")
