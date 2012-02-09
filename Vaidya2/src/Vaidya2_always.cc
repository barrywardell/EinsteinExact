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
#define QAD(x) (SQR(SQR(x)))
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))

static void Vaidya2_always_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_LOOP3(Vaidya2_always,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL = alp[index];
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
    
    CCTK_REAL Jac11 = cos(ToReal(phi))*cos(ToReal(psi)) - 
      cos(ToReal(theta))*sin(ToReal(phi))*sin(ToReal(psi));
    
    CCTK_REAL Jac12 = cos(ToReal(psi))*sin(ToReal(phi)) + 
      cos(ToReal(phi))*cos(ToReal(theta))*sin(ToReal(psi));
    
    CCTK_REAL Jac13 = sin(ToReal(psi))*sin(ToReal(theta));
    
    CCTK_REAL Jac21 = 
      -(cos(ToReal(psi))*cos(ToReal(theta))*sin(ToReal(phi))) - 
      cos(ToReal(phi))*sin(ToReal(psi));
    
    CCTK_REAL Jac22 = cos(ToReal(phi))*cos(ToReal(psi))*cos(ToReal(theta)) 
      - sin(ToReal(phi))*sin(ToReal(psi));
    
    CCTK_REAL Jac23 = cos(ToReal(psi))*sin(ToReal(theta));
    
    CCTK_REAL Jac31 = sin(ToReal(phi))*sin(ToReal(theta));
    
    CCTK_REAL Jac32 = -(cos(ToReal(phi))*sin(ToReal(theta)));
    
    CCTK_REAL Jac33 = cos(ToReal(theta));
    
    CCTK_REAL InvJac11 = Jac11;
    
    CCTK_REAL InvJac12 = Jac21;
    
    CCTK_REAL InvJac13 = Jac31;
    
    CCTK_REAL InvJac21 = Jac12;
    
    CCTK_REAL InvJac22 = Jac22;
    
    CCTK_REAL InvJac23 = Jac32;
    
    CCTK_REAL InvJac31 = Jac13;
    
    CCTK_REAL InvJac32 = Jac23;
    
    CCTK_REAL InvJac33 = Jac33;
    
    CCTK_REAL T = t - ToReal(positiont);
    
    CCTK_REAL XX1 = -(Jac11*(position1 + shiftadd1*T - xx1)) - 
      Jac12*(position2 + shiftadd2*T - xx2) - Jac13*(position3 + shiftadd3*T 
      - xx3);
    
    CCTK_REAL XX2 = -(Jac21*(position1 + shiftadd1*T - xx1)) - 
      Jac22*(position2 + shiftadd2*T - xx2) - Jac23*(position3 + shiftadd3*T 
      - xx3);
    
    CCTK_REAL XX3 = -(Jac31*(position1 + shiftadd1*T - xx1)) - 
      Jac32*(position2 + shiftadd2*T - xx2) - Jac33*(position3 + shiftadd3*T 
      - xx3);
    
    CCTK_REAL X = XX1;
    
    CCTK_REAL Y = XX2;
    
    CCTK_REAL Z = XX3;
    
    CCTK_REAL rXYZ = sqrt(SQR(X) + SQR(Y) + SQR(Z));
    
    CCTK_REAL mTXYZ = (1 + SQR(tanh(INV(ToReal(M))*(T + sqrt(SQR(X) + 
      SQR(Y) + SQR(Z)))*ToReal(dM))))*ToReal(M);
    
    alpL = INV(sqrt(INV(CUB(rXYZ))*(CUB(rXYZ) + 2*mTXYZ*(SQR(X) + SQR(Y) + 
      SQR(Z)))));
    
    CCTK_REAL dtalpL = -2*INV(CUB(rXYZ))*pow(INV(CUB(rXYZ))*(CUB(rXYZ) + 
      2*mTXYZ*(SQR(X) + SQR(Y) + SQR(Z))),-1.5)*(SQR(X) + SQR(Y) + 
      SQR(Z))*SQR(INV(cosh((rXYZ + T)*INV(ToReal(M))*ToReal(dM))))*tanh((rXYZ 
      + T)*INV(ToReal(M))*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL G11 = 1 + 2*mTXYZ*INV(CUB(rXYZ))*SQR(X);
    
    CCTK_REAL G21 = 2*mTXYZ*X*Y*INV(CUB(rXYZ));
    
    CCTK_REAL G31 = 2*mTXYZ*X*Z*INV(CUB(rXYZ));
    
    CCTK_REAL G22 = 1 + 2*mTXYZ*INV(CUB(rXYZ))*SQR(Y);
    
    CCTK_REAL G32 = 2*mTXYZ*Y*Z*INV(CUB(rXYZ));
    
    CCTK_REAL G33 = 1 + 2*mTXYZ*INV(CUB(rXYZ))*SQR(Z);
    
    CCTK_REAL K11 = -2*INV(alpL*QAD(rXYZ)*(CUB(rXYZ) + 2*mTXYZ*(SQR(X) + 
      SQR(Y) + SQR(Z))))*(mTXYZ*(-pow(rXYZ,5) + SQR(X)*(2*CUB(rXYZ) + 
      mTXYZ*(SQR(X) + SQR(Y) + SQR(Z)))) - QAD(rXYZ)*SQR(X*INV(cosh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))))*tanh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL K21 = -2*X*Y*INV(alpL*QAD(rXYZ)*(CUB(rXYZ) + 
      2*mTXYZ*(SQR(X) + SQR(Y) + SQR(Z))))*(mTXYZ*(2*CUB(rXYZ) + 
      mTXYZ*(SQR(X) + SQR(Y) + SQR(Z))) - QAD(rXYZ)*SQR(INV(cosh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))))*tanh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL K31 = -2*X*Z*INV(alpL*QAD(rXYZ)*(CUB(rXYZ) + 
      2*mTXYZ*(SQR(X) + SQR(Y) + SQR(Z))))*(mTXYZ*(2*CUB(rXYZ) + 
      mTXYZ*(SQR(X) + SQR(Y) + SQR(Z))) - QAD(rXYZ)*SQR(INV(cosh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))))*tanh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL K22 = -2*INV(alpL*QAD(rXYZ)*(CUB(rXYZ) + 2*mTXYZ*(SQR(X) + 
      SQR(Y) + SQR(Z))))*(mTXYZ*(-pow(rXYZ,5) + SQR(Y)*(2*CUB(rXYZ) + 
      mTXYZ*(SQR(X) + SQR(Y) + SQR(Z)))) - QAD(rXYZ)*SQR(Y*INV(cosh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))))*tanh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL K32 = -2*Y*Z*INV(alpL*QAD(rXYZ)*(CUB(rXYZ) + 
      2*mTXYZ*(SQR(X) + SQR(Y) + SQR(Z))))*(mTXYZ*(2*CUB(rXYZ) + 
      mTXYZ*(SQR(X) + SQR(Y) + SQR(Z))) - QAD(rXYZ)*SQR(INV(cosh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))))*tanh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL K33 = -2*INV(alpL*QAD(rXYZ)*(CUB(rXYZ) + 2*mTXYZ*(SQR(X) + 
      SQR(Y) + SQR(Z))))*(mTXYZ*(-pow(rXYZ,5) + SQR(Z)*(2*CUB(rXYZ) + 
      mTXYZ*(SQR(X) + SQR(Y) + SQR(Z)))) - QAD(rXYZ)*SQR(Z*INV(cosh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))))*tanh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL betap1 = 2*mTXYZ*rXYZ*X*INV(CUB(rXYZ) + 2*mTXYZ*(SQR(X) + 
      SQR(Y) + SQR(Z)));
    
    CCTK_REAL betap2 = 2*mTXYZ*rXYZ*Y*INV(CUB(rXYZ) + 2*mTXYZ*(SQR(X) + 
      SQR(Y) + SQR(Z)));
    
    CCTK_REAL betap3 = 2*mTXYZ*rXYZ*Z*INV(CUB(rXYZ) + 2*mTXYZ*(SQR(X) + 
      SQR(Y) + SQR(Z)));
    
    CCTK_REAL dtbetap1 = 4*X*INV(SQR(CUB(rXYZ) + 2*mTXYZ*(SQR(X) + SQR(Y) 
      + SQR(Z))))*QAD(rXYZ)*SQR(INV(cosh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))))*tanh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL dtbetap2 = 4*Y*INV(SQR(CUB(rXYZ) + 2*mTXYZ*(SQR(X) + SQR(Y) 
      + SQR(Z))))*QAD(rXYZ)*SQR(INV(cosh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))))*tanh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL dtbetap3 = 4*Z*INV(SQR(CUB(rXYZ) + 2*mTXYZ*(SQR(X) + SQR(Y) 
      + SQR(Z))))*QAD(rXYZ)*SQR(INV(cosh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))))*tanh((rXYZ + 
      T)*INV(ToReal(M))*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL gxxL = 2*(G32*Jac21*Jac31 + Jac11*(G21*Jac21 + G31*Jac31)) + 
      G11*SQR(Jac11) + G22*SQR(Jac21) + G33*SQR(Jac31);
    
    CCTK_REAL gxyL = Jac12*(G11*Jac11 + G21*Jac21 + G31*Jac31) + 
      Jac22*(G21*Jac11 + G22*Jac21 + G32*Jac31) + (G31*Jac11 + G32*Jac21 + 
      G33*Jac31)*Jac32;
    
    CCTK_REAL gxzL = Jac13*(G11*Jac11 + G21*Jac21 + G31*Jac31) + 
      Jac23*(G21*Jac11 + G22*Jac21 + G32*Jac31) + (G31*Jac11 + G32*Jac21 + 
      G33*Jac31)*Jac33;
    
    CCTK_REAL gyyL = 2*(G32*Jac22*Jac32 + Jac12*(G21*Jac22 + G31*Jac32)) + 
      G11*SQR(Jac12) + G22*SQR(Jac22) + G33*SQR(Jac32);
    
    CCTK_REAL gyzL = Jac13*(G11*Jac12 + G21*Jac22 + G31*Jac32) + 
      Jac23*(G21*Jac12 + G22*Jac22 + G32*Jac32) + (G31*Jac12 + G32*Jac22 + 
      G33*Jac32)*Jac33;
    
    CCTK_REAL gzzL = 2*(G32*Jac23*Jac33 + Jac13*(G21*Jac23 + G31*Jac33)) + 
      G11*SQR(Jac13) + G22*SQR(Jac23) + G33*SQR(Jac33);
    
    CCTK_REAL kxxL = 2*(Jac11*(Jac21*K21 + Jac31*K31) + Jac21*Jac31*K32) + 
      K11*SQR(Jac11) + K22*SQR(Jac21) + K33*SQR(Jac31);
    
    CCTK_REAL kxyL = Jac11*(Jac12*K11 + Jac22*K21 + Jac32*K31) + 
      Jac21*(Jac12*K21 + Jac22*K22 + Jac32*K32) + Jac31*(Jac12*K31 + 
      Jac22*K32 + Jac32*K33);
    
    CCTK_REAL kxzL = Jac11*(Jac13*K11 + Jac23*K21 + Jac33*K31) + 
      Jac21*(Jac13*K21 + Jac23*K22 + Jac33*K32) + Jac31*(Jac13*K31 + 
      Jac23*K32 + Jac33*K33);
    
    CCTK_REAL kyyL = 2*(Jac12*(Jac22*K21 + Jac32*K31) + Jac22*Jac32*K32) + 
      K11*SQR(Jac12) + K22*SQR(Jac22) + K33*SQR(Jac32);
    
    CCTK_REAL kyzL = Jac12*(Jac13*K11 + Jac23*K21 + Jac33*K31) + 
      Jac22*(Jac13*K21 + Jac23*K22 + Jac33*K32) + Jac32*(Jac13*K31 + 
      Jac23*K32 + Jac33*K33);
    
    CCTK_REAL kzzL = 2*(Jac13*(Jac23*K21 + Jac33*K31) + Jac23*Jac33*K32) + 
      K11*SQR(Jac13) + K22*SQR(Jac23) + K33*SQR(Jac33);
    
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
  CCTK_ENDLOOP3(Vaidya2_always);
}

extern "C" void Vaidya2_always(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering Vaidya2_always_Body");
  }
  
  if (cctk_iteration % Vaidya2_always_calc_every != Vaidya2_always_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"admbase::curv","admbase::dtlapse","admbase::dtshift","admbase::lapse","admbase::metric","admbase::shift","grid::coordinates"};
  GenericFD_AssertGroupStorage(cctkGH, "Vaidya2_always", 7, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, Vaidya2_always_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving Vaidya2_always_Body");
  }
}
