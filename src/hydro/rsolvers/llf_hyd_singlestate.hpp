#ifndef HYDRO_RSOLVERS_LLF_HYD_SINGLESTATE_HPP_
#define HYDRO_RSOLVERS_LLF_HYD_SINGLESTATE_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file llf_hyd_singlestate.hpp
//! \brief various Local Lax Friedrichs (LLF) Riemann solvers, also known as Rusanov's
//! method, for NR/SR/GR hydrodynamics.  This flux is very diffusive, even more diffusive
//! than HLLE, and so it is not recommended for use in applications.  However, it is
//! useful for testing, or for problems where other Riemann solvers fail.
//!
//! Each solver in this file works on a single L/R state
//!
//! REFERENCES:
//! - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//!   Springer-Verlag, Berlin, (1999) chpt. 10.

namespace hydro {
//----------------------------------------------------------------------------------------
//! \fn void SingleStateLLF_HYD
//  \brief The LLF Riemann solver for hydrodynamics for a single L/R state

KOKKOS_INLINE_FUNCTION
void SingleStateLLF_Hyd(const HydPrim1D &wl, const HydPrim1D &wr, const EOS_Data &eos,
                        HydCons1D &flux) {
  Real qa = wl.d*wl.vx;
  Real qb = wr.d*wr.vx;

  // Compute sum of L/R fluxes
  HydCons1D fsum;
  fsum.d  = qa        + qb;
  fsum.mx = qa*wl.vx + qb*wr.vx;
  fsum.my = qa*wl.vy + qb*wr.vy;
  fsum.mz = qa*wl.vz + qb*wr.vz;

  Real el,er,pl,pr;
  if (eos.is_ideal) {
    pl = eos.IdealGasPressure(wl.e);
    pr = eos.IdealGasPressure(wr.e);
    el = wl.e + 0.5*wl.d*(SQR(wl.vx) + SQR(wl.vy) + SQR(wl.vz));
    er = wr.e + 0.5*wr.d*(SQR(wr.vx) + SQR(wr.vy) + SQR(wr.vz));
    fsum.mx += (pl + pr);
    fsum.e  = (el + pl)*wl.vx + (er + pr)*wr.vx;
  } else {
    fsum.mx += SQR(eos.iso_cs)*(wl.d + wr.d);
  }

  // Compute max wave speed in L,R states (see Toro eq. 10.43)
  if (eos.is_ideal) {
    qa = eos.IdealHydroSoundSpeed(wl.d, pl);
    qb = eos.IdealHydroSoundSpeed(wr.d, pr);
  } else {
    qa = eos.iso_cs;
    qb = eos.iso_cs;
  }
  Real a = fmax( (fabs(wl.vx) + qa), (fabs(wr.vx) + qb) );

  // Compute difference in L/R states dU, multiplied by max wave speed
  HydCons1D du;
  du.d  = a*(wr.d       - wl.d);
  du.mx = a*(wr.d*wr.vx - wl.d*wl.vx);
  du.my = a*(wr.d*wr.vy - wl.d*wl.vy);
  du.mz = a*(wr.d*wr.vz - wl.d*wl.vz);
  if (eos.is_ideal) du.e = a*(er - el);

  // Compute the LLF flux at interface (see Toro eq. 10.42).
  flux.d  = 0.5*(fsum.d  - du.d );
  flux.mx = 0.5*(fsum.mx - du.mx);
  flux.my = 0.5*(fsum.my - du.my);
  flux.mz = 0.5*(fsum.mz - du.mz);
  if (eos.is_ideal) {flux.e = 0.5*(fsum.e - du.e);}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SingleStateLLF_SRHyd
//  \brief The LLF Riemann solver for SR hydrodynamics for a single L/R state

KOKKOS_INLINE_FUNCTION
void SingleStateLLF_SRHyd(const HydPrim1D &wl, const HydPrim1D &wr, const EOS_Data &eos,
                        HydCons1D &flux) {
  // Recall in SR the primitive variables are (\rho, u^i, P_g), where
  //  \rho is the mass density in the comoving/fluid frame,
  //  u^i = \gamma v^i are the spatial components of the 4-velocity (v^i is the 3-vel),
  //  P_g is the pressure.

  Real u0l  = sqrt(1.0 + SQR(wl.vz) + SQR(wl.vy) + SQR(wl.vx)); // Lorentz fact in L
  Real u0r  = sqrt(1.0 + SQR(wr.vz) + SQR(wr.vy) + SQR(wr.vx)); // Lorentz fact in R
  // FIXME ERM: Ideal fluid for now
  Real wgas_l = wl.d + eos.gamma * wl.e;  // total enthalpy in L-state
  Real wgas_r = wr.d + eos.gamma * wr.e;  // total enthalpy in R-state

  // Compute wave speeds in L,R states (see Toro eq. 10.43)
  Real pl = eos.IdealGasPressure(wl.e);
  Real lp_l, lm_l;
  eos.IdealSRHydroSoundSpeeds(wl.d, pl, wl.vx, u0l, lp_l, lm_l);

  Real pr = eos.IdealGasPressure(wr.e);
  Real lp_r, lm_r;
  eos.IdealSRHydroSoundSpeeds(wr.d, pr, wr.vx, u0r, lp_r, lm_r);

  Real qa = fmax(-fmin(lm_l,lm_r), 0.0);
  Real a = fmax(fmax(lp_l,lp_r), qa);

  // Compute sum of L/R fluxes
  qa = wgas_l * wl.vx;
  Real qb = wgas_r * wr.vx;

  HydCons1D fsum;
  fsum.d  = wl.d * wl.vx + wr.d * wr.vx;
  fsum.mx = qa*wl.vx + qb*wr.vx + (pl + pr);
  fsum.my = qa*wl.vy + qb*wr.vy;
  fsum.mz = qa*wl.vz + qb*wr.vz;
  fsum.e  = qa*u0l + qb*u0r;

  // Compute difference dU = U_R - U_L multiplied by max wave speed
  HydCons1D du;
  qa = wgas_r*u0r;
  qb = wgas_l*u0l;
  Real er = qa*u0r - pr;
  Real el = qb*u0l - pl;
  du.d  = a*(u0r*wr.d  - u0l*wl.d);
  du.mx = a*( qa*wr.vx -  qb*wl.vx);
  du.my = a*( qa*wr.vy -  qb*wl.vy);
  du.mz = a*( qa*wr.vz -  qb*wl.vz);
  du.e  = a*(er - el);

  // Compute the LLF flux at the interface
  flux.d  = 0.5*(fsum.d  - du.d );
  flux.mx = 0.5*(fsum.mx - du.mx);
  flux.my = 0.5*(fsum.my - du.my);
  flux.mz = 0.5*(fsum.mz - du.mz);
  flux.e  = 0.5*(fsum.e  - du.e );

  // We evolve tau = E - D
  flux.e -= flux.d;

  return;
}

} // namespace hydro
#endif // HYDRO_RSOLVERS_LLF_HYD_SINGLESTATE_HPP_
