#!/usr/bin/env python3
"""Ensemble comparison: 10 runs × 100,000 particles, condensation + collision."""
import numpy as np
import time, sys, warnings, json
warnings.filterwarnings('ignore')
sys.path.insert(0, '.')
from PyLCM.parameters import *
from PyLCM.micro_particle import *
from PyLCM.aero_init import *
from PyLCM.parcel import *
from PyLCM.condensation import *
from PyLCM.collision import *
from Post_process.analysis import *

def run_single(seed, n_ptcl=100000, nt=3000):
    T0=293.2; P0=1013.0e2; RH=0.88; w=1.0; z0=0.0; zmax=3000.0; dt=1.0
    N_raw=np.array([118.0,11.0,0.72])
    mu_aero=np.log(np.array([0.019,0.056,0.46])*1e-6)
    sigma_aero=np.log(np.array([3.3,1.6,2.2]))
    kappa_arr=[1.6]*4
    N_aero=N_raw*1e6
    theta_init=T0*(p0/P0)**(r_a/cp)
    theta_profiles=theta_init+5e-3*z_env
    e_s=esatw(T0)
    q0=RH*e_s/(P0-RH*e_s)*r_a/rv
    rm_spec=[1e-6, 25e-6]

    np.random.seed(seed)
    T_parcel,q_parcel,particles_list=aero_init('Random',n_ptcl,P0,z0,T0,q0,
        N_aero,mu_aero,sigma_aero,rho_aero,kappa_arr,False)
    P_parcel=P0; z_parcel=z0; S_lst=0.0

    # Storage for time series at key times
    results = {}
    t0 = time.time()

    for t in range(nt):
        tt=(t+1)*dt
        z_parcel,T_parcel,P_parcel=ascend_parcel(z_parcel,T_parcel,P_parcel,w,dt,tt,
            zmax,theta_profiles,None,'linear')
        rho_p,V_p,air_mass=parcel_rho(P_parcel,T_parcel)
        ct=at=et=da=0.0
        particles_list,T_parcel,q_parcel,S_lst,ct,at,et,da=drop_condensation(
            particles_list,T_parcel,q_parcel,P_parcel,nt,dt,air_mass,S_lst,
            rho_aero,False,ct,at,et,da,False)
        ac=au=pr=0.0
        particles_list,ac,au,pr=collection(dt,particles_list,rho_p,rho_liq,
            P_parcel,T_parcel,ac,au,pr,False,z_parcel,zmax,w)

        if (t+1) in [500,1000,1500,2000,2500,3000]:
            sp,qa,qc,qr,na,nc,nr,_,ra,rs=ts_analysis(particles_list,air_mass,rm_spec,60,n_ptcl)
            rho_p2,_,_=parcel_rho(P_parcel,T_parcel)
            results[t+1] = {
                'T': T_parcel-273.15, 'z': z_parcel, 'P': P_parcel,
                'rho': rho_p2, 'qc': qc, 'qr': qr,
                'Nc': nc, 'Nr': nr, 'Na': na,
                'LWC': (qc+qr)*rho_p2
            }

    elapsed = time.time()-t0
    return results, elapsed

if __name__ == '__main__':
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 42
    member = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    print(f'[Member {member}] Starting with seed={seed}, n=100000, nt=3000')
    results, elapsed = run_single(seed)
    print(f'[Member {member}] Done in {elapsed:.1f}s ({elapsed/3000*1e3:.0f}ms/step)')

    # Print results
    for t in sorted(results.keys()):
        r = results[t]
        print(f'  t={t:4d} T={r["T"]:6.2f}C qc={r["qc"]:7.4f} qr={r["qr"]:7.4f} '
              f'Nc={r["Nc"]:7.2f} Nr={r["Nr"]:7.2f} Na={r["Na"]:7.2f} LWC={r["LWC"]:7.3f}')

    # Save to JSON
    outfile = f'/tmp/pylcm_ensemble_{member:02d}.json'
    # Convert for JSON
    out = {}
    for t, r in results.items():
        out[str(t)] = {k: float(v) for k, v in r.items()}
    out['elapsed'] = elapsed
    out['seed'] = seed
    with open(outfile, 'w') as f:
        json.dump(out, f)
    print(f'  Saved to {outfile}')
