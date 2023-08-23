#ifndef MHD_RSOLVERS_LLF_MHD_SINGLESTATE_HPP_
#define MHD_RSOLVERS_LLF_MHD_SINGLESTATE_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file llf_mhd_singlestate.hpp
//! \brief various Local Lax Friedrichs (LLF) Riemann solvers, also known as Rusanov's
//! method, for NR/SR/GR MHD. This flux is very diffusive, even more diffusive than HLLE,
//! and so it is not recommended for use in applications.  However, it is useful for
//! testing, or for problems where other Riemann solvers fail.
//!
//! Each solver in this file works on a single L/R state
//!
//! REFERENCES:
//! - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//!   Springer-Verlag, Berlin, (1999) chpt. 10.

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn void SingleStateLLF_MHD
//! \brief The LLF Riemann solver for MHD for a single L/R state

KOKKOS_INLINE_FUNCTION
void SingleStateLLF_MHD(const MHDPrim1D &wl, const MHDPrim1D &wr, const Real &bxi,
                        const EOS_Data &eos, MHDCons1D &flux) {
  // Compute sum of L/R fluxes
  Real qa = wl.d*wl.vx;
  Real qb = wr.d*wr.vx;
  Real qc = 0.5*(SQR(wl.by) + SQR(wl.bz) - SQR(bxi));
  Real qd = 0.5*(SQR(wr.by) + SQR(wr.bz) - SQR(bxi));

  MHDCons1D fsum;
  fsum.d  = qa       + qb;
  fsum.mx = qa*wl.vx + qb*wr.vx + qc + qd;
  fsum.my = qa*wl.vy + qb*wr.vy - bxi*(wl.by + wr.by);
  fsum.mz = qa*wl.vz + qb*wr.vz - bxi*(wl.bz + wr.bz);
  fsum.by = wl.by*wl.vx + wr.by*wr.vx - bxi*(wl.vy + wr.vy);
  fsum.bz = wl.bz*wl.vx + wr.bz*wr.vx - bxi*(wl.vz + wr.vz);

  Real el,er,pl,pr;
  if (eos.is_ideal) {
    pl = eos.IdealGasPressure(wl.e);
    pr = eos.IdealGasPressure(wr.e);
    el = wl.e + 0.5*wl.d*(SQR(wl.vx)+SQR(wl.vy)+SQR(wl.vz)) + qc + SQR(bxi);
    er = wr.e + 0.5*wr.d*(SQR(wr.vx)+SQR(wr.vy)+SQR(wr.vz)) + qd + SQR(bxi);
    fsum.mx += (pl + pr);
    fsum.e  = (el + pl + qc)*wl.vx + (er + pr + qd)*wr.vx;
    fsum.e  -= bxi*(wl.by*wl.vy + wl.bz*wl.vz);
    fsum.e  -= bxi*(wr.by*wr.vy + wr.bz*wr.vz);
  } else {
    fsum.mx += SQR(eos.iso_cs)*(wl.d + wr.d);
  }

  // Compute max wave speed in L,R states (see Toro eq. 10.43)
  if (eos.is_ideal) {
    qa = eos.IdealMHDFastSpeed(wl.d, pl, bxi, wl.by, wl.bz);
    qb = eos.IdealMHDFastSpeed(wr.d, pr, bxi, wr.by, wr.bz);
  } else {
    qa = eos.IdealMHDFastSpeed(wl.d, bxi, wl.by, wl.bz);
    qb = eos.IdealMHDFastSpeed(wr.d, bxi, wr.by, wr.bz);
  }
  Real a = fmax( (fabs(wl.vx) + qa), (fabs(wr.vx) + qb) );

  // Compute difference in L/R states dU, multiplied by max wave speed
  MHDCons1D du;
  du.d  = a*(wr.d       - wl.d);
  du.mx = a*(wr.d*wr.vx - wl.d*wl.vx);
  du.my = a*(wr.d*wr.vy - wl.d*wl.vy);
  du.mz = a*(wr.d*wr.vz - wl.d*wl.vz);
  if (eos.is_ideal) du.e = a*(er - el);
  du.by = a*(wr.by - wl.by);
  du.bz = a*(wr.bz - wl.bz);

  // Compute the LLF flux at interface (see Toro eq. 10.42).
  flux.d  = 0.5*(fsum.d  - du.d);
  flux.mx = 0.5*(fsum.mx - du.mx);
  flux.my = 0.5*(fsum.my - du.my);
  flux.mz = 0.5*(fsum.mz - du.mz);
  if (eos.is_ideal) {flux.e = 0.5*(fsum.e  - du.e);}
  flux.by = -0.5*(fsum.by - du.by);
  flux.bz =  0.5*(fsum.bz - du.bz);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SingleStateLLF_SRMHD
//! \brief The LLF Riemann solver for SR MHD for a single L/R state

KOKKOS_INLINE_FUNCTION
void SingleStateLLF_SRMHD(const MHDPrim1D &wl, const MHDPrim1D &wr, const Real &bxi,
                          const EOS_Data &eos, MHDCons1D &flux) {
  // Calculate 4-magnetic field in left state
  Real gam_l = sqrt(1.0 + SQR(wl.vx) + SQR(wl.vy) + SQR(wl.vz));
  Real b_l[4];
  b_l[0] = bxi*wl.vx + wl.by*wl.vy + wl.bz*wl.vz;
  b_l[1] = (bxi   + b_l[0] * wl.vx) / gam_l;
  b_l[2] = (wl.by + b_l[0] * wl.vy) / gam_l;
  b_l[3] = (wl.bz + b_l[0] * wl.vz) / gam_l;
  Real b_sq_l = -SQR(b_l[0]) + SQR(b_l[1]) + SQR(b_l[2]) + SQR(b_l[3]);

  // Calculate 4-magnetic field in right state
  Real gam_r = sqrt(1.0 + SQR(wr.vx) + SQR(wr.vy) + SQR(wr.vz));
  Real b_r[4];
  b_r[0] = bxi*wr.vx + wr.by*wr.vy + wr.bz*wr.vz;
  b_r[1] = (bxi   + b_r[0] * wr.vx) / gam_r;
  b_r[2] = (wr.by + b_r[0] * wr.vy) / gam_r;
  b_r[3] = (wr.bz + b_r[0] * wr.vz) / gam_r;
  Real b_sq_r = -SQR(b_r[0]) + SQR(b_r[1]) + SQR(b_r[2]) + SQR(b_r[3]);

  // Calculate left wavespeeds
  Real pl = eos.IdealGasPressure(wl.e);
  Real lm_l, lp_l;
  eos.IdealSRMHDFastSpeeds(wl.d, pl, wl.vx, gam_l, b_sq_l, lp_l, lm_l);

  // Calculate right wavespeeds
  Real pr = eos.IdealGasPressure(wr.e);
  Real lm_r, lp_r;
  eos.IdealSRMHDFastSpeeds(wr.d, pr, wr.vx, gam_r, b_sq_r, lp_r, lm_r);

  // Calculate extremal wavespeeds
  Real lambda_l = fmin(lm_l, lm_r);  // (MB 55)
  Real lambda_r = fmax(lp_l, lp_r);  // (MB 55)
  Real lambda = fmax(lambda_r, -lambda_l);

  // Calculate conserved quantities in L region (MUB 8)
  MHDCons1D consl;
  Real wgas_l = wl.d + eos.gamma * wl.e;
  Real wtot_l = wgas_l + b_sq_l;
  Real ptot_l = pl + 0.5*b_sq_l;
  consl.d  = wl.d * gam_l;
  consl.e  = wtot_l * gam_l * gam_l - b_l[0] * b_l[0] - ptot_l;
  consl.mx = wtot_l * wl.vx * gam_l - b_l[1] * b_l[0];
  consl.my = wtot_l * wl.vy * gam_l - b_l[2] * b_l[0];
  consl.mz = wtot_l * wl.vz * gam_l - b_l[3] * b_l[0];
  consl.by = b_l[2] * gam_l - b_l[0] * wl.vy;
  consl.bz = b_l[3] * gam_l - b_l[0] * wl.vz;

  // Calculate fluxes in L region (MUB 15)
  MHDCons1D fl;
  fl.d  = wl.d * wl.vx;
  fl.e  = wtot_l * gam_l * wl.vx - b_l[0] * b_l[1];
  fl.mx = wtot_l * wl.vx * wl.vx - b_l[1] * b_l[1] + ptot_l;
  fl.my = wtot_l * wl.vy * wl.vx - b_l[2] * b_l[1];
  fl.mz = wtot_l * wl.vz * wl.vx - b_l[3] * b_l[1];
  fl.by = b_l[2] * wl.vx - b_l[1] * wl.vy;
  fl.bz = b_l[3] * wl.vx - b_l[1] * wl.vz;

  // Calculate conserved quantities in R region (MUB 8)
  MHDCons1D consr;
  Real wgas_r = wr.d + eos.gamma * wr.e;
  Real wtot_r = wgas_r + b_sq_r;
  Real ptot_r = pr + 0.5*b_sq_r;
  consr.d  = wr.d * gam_r;
  consr.e  = wtot_r * gam_r * gam_r - b_r[0] * b_r[0] - ptot_r;
  consr.mx = wtot_r * wr.vx * gam_r - b_r[1] * b_r[0];
  consr.my = wtot_r * wr.vy * gam_r - b_r[2] * b_r[0];
  consr.mz = wtot_r * wr.vz * gam_r - b_r[3] * b_r[0];
  consr.by = b_r[2] * gam_r - b_r[0] * wr.vy;
  consr.bz = b_r[3] * gam_r - b_r[0] * wr.vz;

  // Calculate fluxes in R region (MUB 15)
  MHDCons1D fr;
  fr.d  = wr.d * wr.vx;
  fr.e  = wtot_r * gam_r * wr.vx - b_r[0] * b_r[1];
  fr.mx = wtot_r * wr.vx * wr.vx - b_r[1] * b_r[1] + ptot_r;
  fr.my = wtot_r * wr.vy * wr.vx - b_r[2] * b_r[1];
  fr.mz = wtot_r * wr.vz * wr.vx - b_r[3] * b_r[1];
  fr.by = b_r[2] * wr.vx - b_r[1] * wr.vy;
  fr.bz = b_r[3] * wr.vx - b_r[1] * wr.vz;

  // Compute the LLF flux at the interface
  flux.d  = 0.5 * (fl.d  + fr.d  - lambda * (consr.d  - consl.d ));
  flux.e  = 0.5 * (fl.e  + fr.e  - lambda * (consr.e  - consl.e ));
  flux.mx = 0.5 * (fl.mx + fr.mx - lambda * (consr.mx - consl.mx));
  flux.my = 0.5 * (fl.my + fr.my - lambda * (consr.my - consl.my));
  flux.mz = 0.5 * (fl.mz + fr.mz - lambda * (consr.mz - consl.mz));
  flux.by = -0.5 * (fl.by + fr.by - lambda * (consr.by - consl.by));
  flux.bz =  0.5 * (fl.bz + fr.bz - lambda * (consr.bz - consl.bz));

  // We evolve tau = E - D
  flux.e -= flux.d;

  return;
}

} // namespace mhd
#endif // MHD_RSOLVERS_LLF_MHD_SINGLESTATE_HPP_
