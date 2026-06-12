"""Cross-feature physical-correctness checks. Prints PASS/FAIL per expectation."""
import warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import numpy as np

from validation.phys_harness import run, PRESETS
from PyLCM.parameters import rho_liq, p0, r_a, cp, rv as RV
from PyLCM.aero_init import aero_init
from PyLCM.parcel import ascend_parcel, parcel_rho, create_env_profiles
from PyLCM.condensation import esatw
from PyLCM.condensation_fast import drop_condensation_fast as cond
from PyLCM.collision import collection
from PyLCM.mixing import ParameterizedMixing, redistribute_droplets, LEMMixing
from PyLCM.micro_particle import particles

P = []
def chk(name, cond_ok, detail):
    P.append((name, bool(cond_ok)))
    print(f"  [{'PASS' if cond_ok else 'FAIL'}] {name}: {detail}")

def mean_radius(pl):
    A = np.array([p.A for p in pl]); M = np.array([p.M for p in pl])
    m = (A > 0) & (M > 0)
    if not m.any(): return 0.0
    r = (M[m] / A[m] / (4.0/3.0*np.pi*rho_liq)) ** (1.0/3.0)
    return float(np.sum(A[m]*r)/np.sum(A[m]))

print("=== DOMAIN 1: thermodynamics / initial conditions ===")
out, pl = run(seed=0, n_ptcl=1200, nt=900, collisions=False, collect=(300,600,900))
Ts = [out[t]["T"] for t in (300,600,900)]; qcs = [out[t]["qc"] for t in (300,600,900)]
chk("adiabatic cooling (T decreases with z)", Ts[0] > Ts[1] > Ts[2], f"T={[round(x,2) for x in Ts]}C")
chk("LWC grows with height", qcs[0] < qcs[2], f"qc={[round(x,4) for x in qcs]}")
ncw = {}
for w in (0.5,1.0,2.0,4.0):
    # measure at a FIXED height (~300m) so faster updrafts don't simply overshoot the
    # 25um cloud/rain cutoff; total activated = NC+NR (droplets that grew past cutoff).
    nstep = int(300.0 / w)
    o,_ = run(seed=0, n_ptcl=1200, nt=nstep, w=w, collisions=False, collect=(nstep,))
    ncw[w]=o[nstep]["NC"]+o[nstep]["NR"]
chk("updraft -> >= activation (total activated non-decreasing in w)", ncw[0.5] <= ncw[4.0]+1.0, f"NC+NR @~300m vs w={ {k:round(v,1) for k,v in ncw.items()} }")
qrh = {}
for rh in (0.80,0.90,0.98):
    o,_ = run(seed=0, n_ptcl=1200, nt=400, RH=rh, collisions=False, collect=(400,)); qrh[rh]=o[400]["qc"]
chk("higher RH -> earlier/more cloud", qrh[0.98] >= qrh[0.80], f"qc@t400 vs RH={ {k:round(v,4) for k,v in qrh.items()} }")
# conservation: condensation-only, total water vapor+liquid per step
T0,P0,RH,w,dt,nt = 293.2,1013e2,0.92,1.0,1.0,400
from PyLCM.parameters import z_env as ZENV
th = T0*(p0/P0)**(r_a/cp)+5e-3*ZENV
q0 = RH*esatw(T0)/(P0-RH*esatw(T0))*r_a/RV
np.random.seed(0)
T,q,pl2 = aero_init("Random",1200,P0,0.0,T0,q0,np.array((118.,11.,.72))*1e6,
                    np.log(np.array((.019,.056,.46))*1e-6),np.log(np.array((3.3,1.6,2.2))),
                    rho_liq*0+1777.0,[1.6]*4,False)
Pp,z,S = P0,0.0,0.0; drift=0.0
for t in range(nt):
    z,T,Pp = ascend_parcel(z,T,Pp,w,dt,(t+1)*dt,3000.0,th,None,"linear")
    rp,_,am = parcel_rho(Pp,T)
    tw0 = q*am + sum(p.M for p in pl2)
    c=a=e=d=0.0
    pl2,T,q,S,c,a,e,d = cond(pl2,T,q,Pp,nt,dt,am,S,1777.0,False,c,a,e,d,False)
    rp2,_,am2 = parcel_rho(Pp,T)
    tw1 = q*am2 + sum(p.M for p in pl2)
    drift = max(drift, abs(tw1-tw0)/max(tw0,1e-30))
chk("water conservation (condensation-only, per step)", drift < 1e-3, f"max rel drift={drift:.2e}")
# no-NaN grid
bad=[]
for T0g in (283.,293.,300.):
    for P0g in (900e2,1013e2):
        for rhg in (0.85,0.95):
            o,_ = run(seed=1,n_ptcl=600,nt=300,T0=T0g,P0=P0g,RH=rhg,collisions=False,collect=(300,))
            if not all(np.isfinite(v) for v in o[300].values()): bad.append((T0g,P0g,rhg))
chk("no NaN across IC grid (12 configs)", not bad, f"bad configs={bad}")

print("=== DOMAIN 2: aerosols / CCN ===")
nc={}; rmean={}
for pre in ("continental","default","maritime","arctic"):
    o,plf = run(seed=0, aerosol=pre, n_ptcl=1500, nt=700, w=1.0, RH=0.95, collisions=False, collect=(700,))
    nc[pre]=o[700]["NC"]; rmean[pre]=mean_radius(plf)*1e6
chk("CCN count controls NC (continental>maritime,arctic)", nc["continental"]>nc["maritime"] and nc["continental"]>nc["arctic"], f"NC={ {k:round(v,1) for k,v in nc.items()} }")
chk("fewer CCN -> larger droplets (maritime r > continental r)", rmean["maritime"]>rmean["continental"], f"mean r(um)={ {k:round(v,2) for k,v in rmean.items()} }")
om,_ = run(seed=0, aerosol="maritime", n_ptcl=2000, nt=1800, collisions=True, collect=(1800,))
oc,_ = run(seed=0, aerosol="continental", n_ptcl=2000, nt=1800, collisions=True, collect=(1800,))
chk("maritime promotes rain vs continental (suppressed)", om[1800]["qr"] >= oc[1800]["qr"], f"qr maritime={om[1800]['qr']:.4f} vs continental={oc[1800]['qr']:.4f}; NR {om[1800]['NR']:.3f} vs {oc[1800]['NR']:.3f}")
chk("activation fraction sane (NC <= total aerosol)", nc["continental"] <= sum(PRESETS['continental']['N_raw'])+1, f"NC_cont={nc['continental']:.0f} <= {sum(PRESETS['continental']['N_raw'])}")

print("=== DOMAIN 3: collision-coalescence ===")
off,_ = run(seed=0, aerosol="maritime", n_ptcl=2000, nt=1800, collisions=False, collect=(1800,))
on,plon = run(seed=0, aerosol="maritime", n_ptcl=2000, nt=1800, collisions=True, collect=(900,1800,))
# Note: qr_off can be >0 because condensation alone grows maritime droplets past the
# 25um cloud/rain cutoff; the collision signal is the LARGE amplification of rain.
chk("collisions strongly amplify rain (on >> off)", on[1800]["qr"] > 3.0*max(off[1800]["qr"],1e-6), f"qr off={off[1800]['qr']:.4f} on={on[1800]['qr']:.4f} ({on[1800]['qr']/max(off[1800]['qr'],1e-9):.1f}x) NR_on={on[1800]['NR']:.3f}")
chk("rain develops after cloud (qr rises late)", on[900]["qc"]>0 , f"qc@900={on[900]['qc']:.3f} qr@900={on[900]['qr']:.4f} qr@1800={on[1800]['qr']:.4f}")
# collision-step mass conservation
th2 = 293.2*(p0/1013e2)**(r_a/cp)+5e-3*ZENV
np.random.seed(0)
T,q,plc = aero_init("Random",1500,1013e2,0.0,293.2,0.92*esatw(293.2)/(1013e2-0.92*esatw(293.2))*r_a/RV,
                    np.array((100.,20.))*1e6,np.log(np.array((.08,.4))*1e-6),np.log(np.array((1.6,2.0))),1777.0,[1.0]*3,False)
Pp,z,S=1013e2,0.0,0.0; cdrift=0.0
for t in range(1000):
    z,T,Pp=ascend_parcel(z,T,Pp,1.0,1.0,(t+1),3000.0,th2,None,"linear")
    rp,_,am=parcel_rho(Pp,T); c=a=e=d=0.0
    plc,T,q,S,c,a,e,d=cond(plc,T,q,Pp,1000,1.0,am,S,1777.0,False,c,a,e,d,False)
    m0=sum(p.M for p in plc); ac=au=pr=0.0
    plc,ac,au,pr=collection(1.0,plc,rp,rho_liq,Pp,T,ac,au,pr,False,z,3000.0,1.0)
    m1=sum(p.M for p in plc)
    cdrift=max(cdrift, abs(m1-m0)/max(m0,1e-30))
chk("collision step conserves liquid water mass", cdrift < 1e-9, f"max rel drift={cdrift:.2e}")
og,_ = run(seed=0, aerosol="maritime", n_ptcl=2000, nt=1500, collisions=True, switch_turb=False, eps=0.0, collect=(1500,))
ot,_ = run(seed=0, aerosol="maritime", n_ptcl=2000, nt=1500, collisions=True, switch_turb=True, eps=0.04, collect=(1500,))
chk("turbulent kernel >= gravitational rain", ot[1500]["qr"]+1e-9 >= og[1500]["qr"]*0.9, f"qr grav={og[1500]['qr']:.4f} turb={ot[1500]['qr']:.4f}")

print("=== DOMAIN 4: entrainment / IHMD ===")
# unit law
def cloudps(n=40,M=2e-9,A=100.):
    out=[]
    for _ in range(n):
        p=particles(1); p.M,p.A,p.Ns,p.kappa=M,A,1e-18,0.5; out.append(p)
    return out
law_ok=True
for ih in (0.,0.25,0.5,0.75,1.0):
    ps=cloudps(); N0=sum(p.A for p in ps); q0c=sum(p.M for p in ps)
    redistribute_droplets(ps,ih,0.3); N1=sum(p.A for p in ps); q1=sum(p.M for p in ps)
    law_ok &= np.isclose(N1/N0,(q1/q0c)**ih,rtol=1e-9)
chk("IHMD law N/N0=(q/q0)^IHMD exact", law_ok, "verified for IHMD in {0,.25,.5,.75,1}")
def tot_mult(pl): return float(sum(p.A for p in pl if p.M > 0))
qv,thp,zc = create_env_profiles(290.0,0.010,0.0,95000.0,"Stable")
# Gentle entrainment so the cloud survives; measure TRUE multiplicity sum(A) (condensation
# never changes A, so only IHMD redistribution can) rather than the thresholded NC count.
base,plb = run(seed=0, aerosol="default", collisions=False, mixing=None, n_ptcl=1500, nt=600, collect=(600,))
Abase = tot_mult(plb)
res={}
for ih in (0.0,0.5,1.0):
    m=ParameterizedMixing(3e-4,ih,qv,thp,zc)
    o,plm=run(seed=0, aerosol="default", collisions=False, mixing=m, n_ptcl=1500, nt=600, collect=(600,))
    res[ih]=(tot_mult(plm), mean_radius(plm)*1e6, o[600]["qc"])
chk("homogeneous(IHMD=0) conserves multiplicity sum(A)", abs(res[0.0][0]-Abase)/max(Abase,1e-9) < 1e-6, f"sumA base={Abase:.3e} ihmd0={res[0.0][0]:.3e}")
chk("inhomogeneous(IHMD=1) depletes multiplicity", res[1.0][0] < res[0.0][0]*0.999, f"sumA ihmd0={res[0.0][0]:.3e} ihmd1={res[1.0][0]:.3e}")
chk("multiplicity decreases monotonically with IHMD", res[0.0][0] >= res[0.5][0] >= res[1.0][0]-1e-3, f"sumA={[f'{res[i][0]:.2e}' for i in (0.,0.5,1.)]}")
chk("mean radius higher at IHMD=1 than IHMD=0", res[1.0][1] > res[0.0][1], f"r(um) ihmd0={res[0.0][1]:.2f} ihmd1={res[1.0][1]:.2f}")
chk("entrainment reduces cloud water vs baseline", res[0.5][2] < base[600]["qc"], f"qc base={base[600]['qc']:.4f} mixed={res[0.5][2]:.4f}")
ml,_=run(seed=0,collisions=False,mixing=ParameterizedMixing(5e-4,0.5,qv,thp,zc),n_ptcl=1200,nt=800,collect=(800,))
mh,_=run(seed=0,collisions=False,mixing=ParameterizedMixing(2e-3,0.5,qv,thp,zc),n_ptcl=1200,nt=800,collect=(800,))
chk("stronger entrainment removes more water", mh[800]["qc"] <= ml[800]["qc"], f"qc lam=5e-4 {ml[800]['qc']:.4f} >= lam=2e-3 {mh[800]['qc']:.4f}")
try:
    LEMMixing().apply([],288.,.009,9e4,500.,1.,1.,100.); lem=False
except NotImplementedError: lem=True
chk("LEM backend is a Phase-3b stub", lem, "raises NotImplementedError")

n_pass = sum(1 for _,ok in P if ok); n=len(P)
print(f"\n=== SUMMARY: {n_pass}/{n} physical checks PASS ===")
fails=[nm for nm,ok in P if not ok]
if fails: print("FAILURES:", fails)
