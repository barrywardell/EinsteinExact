/*  File produced by Kranc */

#define KRANC_C

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"
#include "Differencing.h"
#include "cctk_Loop.h"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define INV(x) ((CCTK_REAL)1.0 / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * SQR(x))
#define QAD(x) (SQR(SQR(x)))

static void ModifiedSchwarzschildBL_initial_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  ptrdiff_t const di = 1;
  ptrdiff_t const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const cdi = sizeof(CCTK_REAL) * di;
  ptrdiff_t const cdj = sizeof(CCTK_REAL) * dj;
  ptrdiff_t const cdk = sizeof(CCTK_REAL) * dk;
  CCTK_REAL const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL const dz = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL const dt = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL const t = ToReal(cctk_time);
  CCTK_REAL const dxi = INV(dx);
  CCTK_REAL const dyi = INV(dy);
  CCTK_REAL const dzi = INV(dz);
  CCTK_REAL const khalf = 0.5;
  CCTK_REAL const kthird = 1/3.0;
  CCTK_REAL const ktwothird = 2.0/3.0;
  CCTK_REAL const kfourthird = 4.0/3.0;
  CCTK_REAL const keightthird = 8.0/3.0;
  CCTK_REAL const hdxi = 0.5 * dxi;
  CCTK_REAL const hdyi = 0.5 * dyi;
  CCTK_REAL const hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  CCTK_LOOP3(ModifiedSchwarzschildBL_initial,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL xL = x[index];
    CCTK_REAL yL = y[index];
    CCTK_REAL zL = z[index];
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL xx1 = xL;
    
    CCTK_REAL xx2 = yL;
    
    CCTK_REAL xx3 = zL;
    
    CCTK_REAL position1 = ToReal(positionx);
    
    CCTK_REAL position2 = ToReal(positiony);
    
    CCTK_REAL position3 = ToReal(positionz);
    
    CCTK_REAL shiftadd1 = ToReal(shiftaddx);
    
    CCTK_REAL shiftadd2 = ToReal(shiftaddy);
    
    CCTK_REAL shiftadd3 = ToReal(shiftaddz);
    
    CCTK_REAL csetemp0 = cos(ToReal(phi));
    
    CCTK_REAL csetemp1 = cos(ToReal(psi));
    
    CCTK_REAL csetemp2 = cos(ToReal(theta));
    
    CCTK_REAL csetemp3 = sin(ToReal(phi));
    
    CCTK_REAL csetemp4 = sin(ToReal(psi));
    
    CCTK_REAL Jac11 = csetemp0*csetemp1 - csetemp2*csetemp3*csetemp4;
    
    CCTK_REAL Jac12 = csetemp1*csetemp3 + csetemp0*csetemp2*csetemp4;
    
    CCTK_REAL csetemp5 = sin(ToReal(theta));
    
    CCTK_REAL Jac13 = csetemp4*csetemp5;
    
    CCTK_REAL Jac21 = -(csetemp1*csetemp2*csetemp3) - csetemp0*csetemp4;
    
    CCTK_REAL Jac22 = csetemp0*csetemp1*csetemp2 - csetemp3*csetemp4;
    
    CCTK_REAL Jac23 = csetemp1*csetemp5;
    
    CCTK_REAL Jac31 = csetemp3*csetemp5;
    
    CCTK_REAL Jac32 = -(csetemp0*csetemp5);
    
    CCTK_REAL Jac33 = csetemp2;
    
    CCTK_REAL InvJac11 = Jac11;
    
    CCTK_REAL InvJac12 = Jac21;
    
    CCTK_REAL InvJac13 = Jac31;
    
    CCTK_REAL InvJac21 = Jac12;
    
    CCTK_REAL InvJac22 = Jac22;
    
    CCTK_REAL InvJac23 = Jac32;
    
    CCTK_REAL InvJac31 = Jac13;
    
    CCTK_REAL InvJac32 = Jac23;
    
    CCTK_REAL InvJac33 = Jac33;
    
    CCTK_REAL T = ToReal(lapsefactor)*(t - ToReal(positiont));
    
    CCTK_REAL csetemp6 = -(shiftadd1*T);
    
    CCTK_REAL csetemp7 = -(shiftadd2*T);
    
    CCTK_REAL csetemp8 = -(shiftadd3*T);
    
    CCTK_REAL XX1 = Jac11*(csetemp6 - position1 + xx1) + Jac12*(csetemp7 - 
      position2 + xx2) + Jac13*(csetemp8 - position3 + xx3);
    
    CCTK_REAL XX2 = Jac21*(csetemp6 - position1 + xx1) + Jac22*(csetemp7 - 
      position2 + xx2) + Jac23*(csetemp8 - position3 + xx3);
    
    CCTK_REAL XX3 = Jac31*(csetemp6 - position1 + xx1) + Jac32*(csetemp7 - 
      position2 + xx2) + Jac33*(csetemp8 - position3 + xx3);
    
    CCTK_REAL X = XX1;
    
    CCTK_REAL Y = XX2;
    
    CCTK_REAL Z = XX3;
    
    CCTK_REAL rXYZ = sqrt(SQR(X) + SQR(Y) + SQR(Z));
    
    CCTK_REAL csetemp9 = 0.5*ToReal(M);
    
    CCTK_REAL csetemp10 = INV(rXYZ);
    
    CCTK_REAL csetemp11 = INV(ToReal(M));
    
    CCTK_REAL csetemp12 = SQR(ToReal(M));
    
    CCTK_REAL csetemp13 = INV(SQR(Pi));
    
    CCTK_REAL G11 = 0.00390625*QAD(4 + csetemp10*(8*fmin(csetemp9,rXYZ) + 
      csetemp11*(csetemp12*(csetemp13 - 
      csetemp13*cos(4*csetemp11*Pi*fmin(csetemp9,rXYZ))) - 
      8*SQR(fmin(csetemp9,rXYZ)))));
    
    CCTK_REAL G21 = 0;
    
    CCTK_REAL G31 = 0;
    
    CCTK_REAL G22 = 0.00390625*QAD(4 + csetemp10*(8*fmin(csetemp9,rXYZ) + 
      csetemp11*(csetemp12*(csetemp13 - 
      csetemp13*cos(4*csetemp11*Pi*fmin(csetemp9,rXYZ))) - 
      8*SQR(fmin(csetemp9,rXYZ)))));
    
    CCTK_REAL G32 = 0;
    
    CCTK_REAL G33 = 0.00390625*QAD(4 + csetemp10*(8*fmin(csetemp9,rXYZ) + 
      csetemp11*(csetemp12*(csetemp13 - 
      csetemp13*cos(4*csetemp11*Pi*fmin(csetemp9,rXYZ))) - 
      8*SQR(fmin(csetemp9,rXYZ)))));
    
    CCTK_REAL K11 = 0;
    
    CCTK_REAL K21 = 0;
    
    CCTK_REAL K31 = 0;
    
    CCTK_REAL K22 = 0;
    
    CCTK_REAL K32 = 0;
    
    CCTK_REAL K33 = 0;
    
    CCTK_REAL alpp = 1;
    
    CCTK_REAL dtalpp = 0;
    
    CCTK_REAL betap1 = 0;
    
    CCTK_REAL betap2 = 0;
    
    CCTK_REAL betap3 = 0;
    
    CCTK_REAL dtbetap1 = 0;
    
    CCTK_REAL dtbetap2 = 0;
    
    CCTK_REAL dtbetap3 = 0;
    
    CCTK_REAL csetemp14 = SQR(Jac11);
    
    CCTK_REAL csetemp15 = SQR(Jac21);
    
    CCTK_REAL csetemp16 = SQR(Jac31);
    
    CCTK_REAL gxxL = csetemp14*G11 + csetemp15*G22 + csetemp16*G33 + 
      2*(G32*Jac21*Jac31 + Jac11*(G21*Jac21 + G31*Jac31));
    
    CCTK_REAL gxyL = Jac12*(G11*Jac11 + G21*Jac21 + G31*Jac31) + 
      Jac22*(G21*Jac11 + G22*Jac21 + G32*Jac31) + (G31*Jac11 + G32*Jac21 + 
      G33*Jac31)*Jac32;
    
    CCTK_REAL gxzL = Jac13*(G11*Jac11 + G21*Jac21 + G31*Jac31) + 
      Jac23*(G21*Jac11 + G22*Jac21 + G32*Jac31) + (G31*Jac11 + G32*Jac21 + 
      G33*Jac31)*Jac33;
    
    CCTK_REAL csetemp17 = SQR(Jac12);
    
    CCTK_REAL csetemp18 = SQR(Jac22);
    
    CCTK_REAL csetemp19 = SQR(Jac32);
    
    CCTK_REAL gyyL = csetemp17*G11 + csetemp18*G22 + csetemp19*G33 + 
      2*(G32*Jac22*Jac32 + Jac12*(G21*Jac22 + G31*Jac32));
    
    CCTK_REAL gyzL = Jac13*(G11*Jac12 + G21*Jac22 + G31*Jac32) + 
      Jac23*(G21*Jac12 + G22*Jac22 + G32*Jac32) + (G31*Jac12 + G32*Jac22 + 
      G33*Jac32)*Jac33;
    
    CCTK_REAL csetemp20 = SQR(Jac13);
    
    CCTK_REAL csetemp21 = SQR(Jac23);
    
    CCTK_REAL csetemp22 = SQR(Jac33);
    
    CCTK_REAL gzzL = csetemp20*G11 + csetemp21*G22 + csetemp22*G33 + 
      2*(G32*Jac23*Jac33 + Jac13*(G21*Jac23 + G31*Jac33));
    
    CCTK_REAL kxxL = csetemp14*K11 + csetemp15*K22 + 2*(Jac11*(Jac21*K21 + 
      Jac31*K31) + Jac21*Jac31*K32) + csetemp16*K33;
    
    CCTK_REAL kxyL = Jac11*(Jac12*K11 + Jac22*K21 + Jac32*K31) + 
      Jac21*(Jac12*K21 + Jac22*K22 + Jac32*K32) + Jac31*(Jac12*K31 + 
      Jac22*K32 + Jac32*K33);
    
    CCTK_REAL kxzL = Jac11*(Jac13*K11 + Jac23*K21 + Jac33*K31) + 
      Jac21*(Jac13*K21 + Jac23*K22 + Jac33*K32) + Jac31*(Jac13*K31 + 
      Jac23*K32 + Jac33*K33);
    
    CCTK_REAL kyyL = csetemp17*K11 + csetemp18*K22 + 2*(Jac12*(Jac22*K21 + 
      Jac32*K31) + Jac22*Jac32*K32) + csetemp19*K33;
    
    CCTK_REAL kyzL = Jac12*(Jac13*K11 + Jac23*K21 + Jac33*K31) + 
      Jac22*(Jac13*K21 + Jac23*K22 + Jac33*K32) + Jac32*(Jac13*K31 + 
      Jac23*K32 + Jac33*K33);
    
    CCTK_REAL kzzL = csetemp20*K11 + csetemp21*K22 + 2*(Jac13*(Jac23*K21 + 
      Jac33*K31) + Jac23*Jac33*K32) + csetemp22*K33;
    
    CCTK_REAL alpL = alpp*ToReal(lapsefactor);
    
    CCTK_REAL dtalpL = dtalpp*ToReal(lapsefactor);
    
    CCTK_REAL betaxL = betap1*InvJac11 + betap2*InvJac12 + betap3*InvJac13 
      + shiftadd1;
    
    CCTK_REAL betayL = betap1*InvJac21 + betap2*InvJac22 + betap3*InvJac23 
      + shiftadd2;
    
    CCTK_REAL betazL = betap1*InvJac31 + betap2*InvJac32 + betap3*InvJac33 
      + shiftadd3;
    
    CCTK_REAL dtbetaxL = dtbetap1*InvJac11 + dtbetap2*InvJac12 + 
      dtbetap3*InvJac13;
    
    CCTK_REAL dtbetayL = dtbetap1*InvJac21 + dtbetap2*InvJac22 + 
      dtbetap3*InvJac23;
    
    CCTK_REAL dtbetazL = dtbetap1*InvJac31 + dtbetap2*InvJac32 + 
      dtbetap3*InvJac33;
    
    /* Copy local copies back to grid functions */
    alp[index] = alpL;
    betax[index] = betaxL;
    betay[index] = betayL;
    betaz[index] = betazL;
    dtalp[index] = dtalpL;
    dtbetax[index] = dtbetaxL;
    dtbetay[index] = dtbetayL;
    dtbetaz[index] = dtbetazL;
    gxx[index] = gxxL;
    gxy[index] = gxyL;
    gxz[index] = gxzL;
    gyy[index] = gyyL;
    gyz[index] = gyzL;
    gzz[index] = gzzL;
    kxx[index] = kxxL;
    kxy[index] = kxyL;
    kxz[index] = kxzL;
    kyy[index] = kyyL;
    kyz[index] = kyzL;
    kzz[index] = kzzL;
  }
  CCTK_ENDLOOP3(ModifiedSchwarzschildBL_initial);
}

extern "C" void ModifiedSchwarzschildBL_initial(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ModifiedSchwarzschildBL_initial_Body");
  }
  
  if (cctk_iteration % ModifiedSchwarzschildBL_initial_calc_every != ModifiedSchwarzschildBL_initial_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "admbase::curv",
    "admbase::dtlapse",
    "admbase::dtshift",
    "admbase::lapse",
    "admbase::metric",
    "admbase::shift",
    "grid::coordinates"};
  GenericFD_AssertGroupStorage(cctkGH, "ModifiedSchwarzschildBL_initial", 7, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, ModifiedSchwarzschildBL_initial_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ModifiedSchwarzschildBL_initial_Body");
  }
}
