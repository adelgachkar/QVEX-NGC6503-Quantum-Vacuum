#!/usr/bin/env python3
"""Euler–Tolman / circular-flow consistency check on SPARC NGC 6503.

SPARC columns: R[kpc], Vobs[km/s], e_Vobs[km/s], Vgas[km/s],
Vdisk[km/s], Vbul[km/s].  The source is RC_MRT.mrt, whose published table
contains 31 NGC6503 points. If the MRT file is not alongside this script,
the exact rows copied from it below are used.

IMPORTANT: SPARC gives kinematics and baryonic velocity contributions, not
rho(r), p_r(r), p_t(r), or Q_r. Therefore the script validates equations,
units, grid/data wiring, and zero-pressure algebraic limit; it does NOT claim
an independent pressure-profile or gravity-model fit.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np

# MRT row tuple order required by task: r, v_obs, v_gas, v_disk, v_bul, err_v
SPARC_NGC6503 = np.array([(0.76, 77.0, 4.04, 64.484, 0.0, 10.0), (1.52, 100.0, 7.62, 112.309, 0.0, 5.0), (2.28, 110.0, 9.72, 130.062, 0.0, 2.0), (3.04, 118.0, 12.42, 122.567, 0.0, 2.0), (3.79, 121.0, 13.35, 111.215, 0.0, 2.0), (4.55, 121.0, 14.78, 101.696, 0.0, 2.0), (5.31, 117.0, 17.55, 94.127, 0.0, 2.0), (6.07, 116.0, 19.2, 88.128, 0.0, 2.0), (6.83, 116.0, 20.33, 83.851, 0.0, 2.0), (7.59, 116.0, 21.66, 79.894, 0.0, 2.0), (8.35, 115.0, 22.79, 75.909, 0.0, 1.0), (9.11, 115.0, 23.92, 72.563, 0.0, 1.0), (9.86, 116.0, 24.84, 70.037, 0.0, 1.0), (10.64, 117.0, 25.46, 67.639, 0.0, 1.0), (11.38, 117.0, 26.38, 65.195, 0.0, 1.0), (12.12, 118.0, 27.2, 63.08, 0.0, 1.0), (12.86, 117.0, 28.13, 61.129, 0.0, 1.0), (13.7, 116.0, 28.74, 59.123, 0.0, 1.0), (14.44, 115.0, 28.64, 57.49, 0.0, 1.0), (15.18, 117.0, 28.33, 55.977, 0.0, 1.0), (15.91, 116.0, 28.13, 54.6, 0.0, 1.0), (16.65, 115.0, 28.33, 53.296, 0.0, 2.0), (17.49, 115.0, 28.23, 51.91, 0.0, 2.0), (18.23, 115.0, 27.82, 50.798, 0.0, 2.0), (18.97, 114.0, 27.2, 49.74, 0.0, 2.0), (19.71, 115.0, 26.9, 48.746, 0.0, 2.0), (20.45, 114.0, 26.69, 47.816, 0.0, 2.0), (21.29, 116.0, 26.49, 46.822, 0.0, 2.0), (22.03, 115.0, 26.18, 45.992, 0.0, 3.0), (22.76, 118.0, 25.87, 45.235, 0.0, 4.0), (23.5, 115.0, 25.46, 44.488, 0.0, 10.0)], dtype=float)


@dataclass
class EulerTolmanHydrodynamics:
    """Spherical anisotropic-fluid balance, SI units.

    rho [kg m^-3], p_r,p_t [Pa], r [m], dlnB_dr [m^-1],
    Q_r [N m^-3]. Covariant static balance:
      dp_r/dr + .5*(rho*c^2+p_r)*dlnB/dr + 2*(p_r-p_t)/r = Q_r.
    The weak-field static limit uses dlnB/dr=2*g/c^2 and gives
      g = (Q_r - dp_r/dr - 2*(p_r-p_t)/r)/rho.

    For a circularly streaming fluid (Newtonian radial Euler equation),
      v_phi^2/r = g_N + [dp_r/dr + 2*(p_r-p_t)/r - Q_r]/rho.
    This rotational term is not present in the static TOV/Tolman equation;
    it is added explicitly in the weak-field circular-flow relation.
    """
    c: float = 299792458.0

    @staticmethod
    def _arrays(r, rho, p_r, p_t, dp_r_dr, Q_r):
        r,rho,p_r,p_t,dp_r_dr,Q_r=np.broadcast_arrays(*[np.asarray(x,dtype=float) for x in (r,rho,p_r,p_t,dp_r_dr,Q_r)])
        if np.any(~np.isfinite(r+rho+p_r+p_t+dp_r_dr+Q_r)) or np.any(r<=0) or np.any(rho<=0):
            raise ValueError("all fields must be finite; r and rho must be positive")
        return r,rho,p_r,p_t,dp_r_dr,Q_r

    def tolman_residual(self,r,rho,p_r,p_t,dp_r_dr,dlnB_dr,Q_r=0.0):
        r,rho,p_r,p_t,dp_r_dr,Q_r=self._arrays(r,rho,p_r,p_t,dp_r_dr,Q_r)
        return dp_r_dr+0.5*(rho*self.c**2+p_r)*np.asarray(dlnB_dr,dtype=float)+2*(p_r-p_t)/r-Q_r

    def dlnB_dr_from_circular_speed(self,r,v_circ):
        r,v=np.broadcast_arrays(np.asarray(r,dtype=float),np.asarray(v_circ,dtype=float))
        if np.any(r<=0) or np.any(v<0) or np.any(~np.isfinite(r+v)): raise ValueError("invalid r or speed")
        return 2*v**2/(r*self.c**2)

    def pressure_acceleration(self,r,rho,dp_r_dr,p_r,p_t,Q_r=0.0):
        r,rho,p_r,p_t,dp_r_dr,Q_r=self._arrays(r,rho,p_r,p_t,dp_r_dr,Q_r)
        return (dp_r_dr+2*(p_r-p_t)/r-Q_r)/rho

    def circular_speed_squared(self,r,g_N,rho,dp_r_dr,p_r,p_t,Q_r=0.0):
        """Weak-field circular-flow balance; g_N is inward gravity magnitude."""
        return np.asarray(r)*np.asarray(g_N)+np.asarray(r)*self.pressure_acceleration(r,rho,dp_r_dr,p_r,p_t,Q_r)

    def static_gravity_from_pressure(self,r,rho,dp_r_dr,p_r,p_t,Q_r=0.0):
        return -self.pressure_acceleration(r,rho,dp_r_dr,p_r,p_t,Q_r)

def load_mrt_or_builtin():
    path=Path(__file__).with_name("RC_MRT.mrt")
    if path.exists():
        out=[]
        for line in path.read_text(encoding="ascii",errors="ignore").splitlines():
            fields=line.split()
            if fields and fields[0]=="NGC6503":
                # MRT: ID,R,Vobs,eVobs,Vgas,Vdisk,...,Vbul
                out.append([float(fields[i]) for i in (1,2,4,5,9,3)])
        if len(out)!=31: raise AssertionError(f"expected 31 MRT rows, found {len(out)}")
        return np.asarray(out),"RC_MRT.mrt"
    return SPARC_NGC6503,"embedded MRT transcription"

def verify():
    d,source=load_mrt_or_builtin()
    assert d.shape==(31,6), d.shape
    r_kpc,vobs,vgas,vdisk,vbul,err=d.T
    assert np.all(np.diff(r_kpc)>0) and np.all(err>0)
    assert np.all(np.isfinite(d))
    # Confirm exact field ordering/units range for this SPARC table.
    assert np.isclose(r_kpc[0],0.76) and np.isclose(r_kpc[-1],23.50)
    assert np.all((vobs>0)&(vobs<300))
    kpc=3.085677581491367e19; kms=1000.0
    r=r_kpc*kpc; v=vobs*kms
    model=EulerTolmanHydrodynamics()
    # Weak-field circular kinematics: metric gradient equals 2 v^2/(r c^2).
    dlnB=model.dlnB_dr_from_circular_speed(r,v)
    g=v*v/r
    assert np.allclose(0.5*model.c**2*dlnB,g,rtol=2e-15,atol=0)
    # Zero-pressure test: circular balance reduces exactly to v^2/r=g_N.
    zero=np.zeros_like(r); rho=np.ones_like(r)
    v2=model.circular_speed_squared(r,g,rho,zero,zero,zero)
    assert np.allclose(v2,v*v,rtol=2e-15,atol=0)
    # Non-zero anisotropic pressure/source test and identity with full residual.
    pr=np.linspace(1e-12,2e-12,len(r)); pt=0.8*pr
    dp=np.gradient(pr,r,edge_order=2); Q=np.full_like(r,1e-30)
    acc=model.pressure_acceleration(r,rho,dp,pr,pt,Q)
    v2p=model.circular_speed_squared(r,g,rho,dp,pr,pt,Q)
    assert np.allclose(v2p/r-g,acc,rtol=2e-15,atol=1e-30)
    # Static conservation is separately checked with dlnB solved from its RHS.
    static_dlnB=2*(Q-dp-2*(pr-pt)/r)/(rho*model.c**2+pr)
    residual=model.tolman_residual(r,rho,pr,pt,dp,static_dlnB,Q)
    assert np.max(np.abs(residual))<1e-27
    print(f"PASS: {len(r)} SPARC points from {source}; R={r_kpc[0]:.2f}..{r_kpc[-1]:.2f} kpc")
    print(f"Vobs={vobs.min():.1f}..{vobs.max():.1f} km/s; max err={err.max():.1f} km/s")
    print(f"Tolman static residual max={np.max(np.abs(residual)):.3e} SI")
    print(f"Weak-field parameter max(v^2/c^2)={np.max(v*v/model.c**2):.3e}")
    print("PASS: metric-gradient identity, pressure-free circular limit, anisotropic/source algebra.")
    print("CAVEAT: SPARC has no density/pressure/shear/source profiles; this is equation/data consistency, not a pressure-driven fit.")

if __name__=="__main__": verify()
