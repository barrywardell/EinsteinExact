/*  File produced by Kranc */

#define KRANC_C

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"
#include "Differencing.h"

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

void gauge_wave_always_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare the variables used for looping over grid points */
  CCTK_INT i, j, k;
  // CCTK_INT index = INITVALUE;
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering gauge_wave_always_Body");
  }
  
  if (cctk_iteration % gauge_wave_always_calc_every != gauge_wave_always_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"admbase::curv","admbase::lapse","admbase::metric","admbase::shift","grid::coordinates"};
  GenericFD_AssertGroupStorage(cctkGH, "gauge_wave_always", 5, groups);
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  CCTK_REAL const dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL const dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL const dz = CCTK_DELTA_SPACE(2);
  int const di = 1;
  int const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  int const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  CCTK_REAL const dxi = 1.0 / dx;
  CCTK_REAL const dyi = 1.0 / dy;
  CCTK_REAL const dzi = 1.0 / dz;
  CCTK_REAL const khalf = 0.5;
  CCTK_REAL const kthird = 1/3.0;
  CCTK_REAL const ktwothird = 2.0/3.0;
  CCTK_REAL const kfourthird = 4.0/3.0;
  CCTK_REAL const keightthird = 8.0/3.0;
  CCTK_REAL const hdxi = 0.5 * dxi;
  CCTK_REAL const hdyi = 0.5 * dyi;
  CCTK_REAL const hdzi = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  
  /* Loop over the grid points */
  for (k = min[2]; k < max[2]; k++)
  {
    for (j = min[1]; j < max[1]; j++)
    {
      for (i = min[0]; i < max[0]; i++)
      {
         int  const  index  =  CCTK_GFINDEX3D(cctkGH,i,j,k) ;
        /* Declare derivatives */
        
        /* Assign local copies of grid functions */
        CCTK_REAL  xL = x[index];
        CCTK_REAL  yL = y[index];
        CCTK_REAL  zL = z[index];
        
        /* Include user supplied include files */
        
        /* Precompute derivatives */
        
        /* Calculate temporaries and grid functions */
        CCTK_REAL xx1 = xL;
        
        CCTK_REAL xx2 = yL;
        
        CCTK_REAL xx3 = zL;
        
        CCTK_REAL Jac11 = Cos(phi)*Cos(psi) - Cos(theta)*Sin(phi)*Sin(psi);
        
        CCTK_REAL Jac12 = Cos(psi)*Sin(phi) + Cos(phi)*Cos(theta)*Sin(psi);
        
        CCTK_REAL Jac13 = Sin(psi)*Sin(theta);
        
        CCTK_REAL Jac21 = -(Cos(psi)*Cos(theta)*Sin(phi)) - Cos(phi)*Sin(psi);
        
        CCTK_REAL Jac22 = Cos(phi)*Cos(psi)*Cos(theta) - Sin(phi)*Sin(psi);
        
        CCTK_REAL Jac23 = Cos(psi)*Sin(theta);
        
        CCTK_REAL Jac31 = Sin(phi)*Sin(theta);
        
        CCTK_REAL Jac32 = -(Cos(phi)*Sin(theta));
        
        CCTK_REAL Jac33 = Cos(theta);
        
        CCTK_REAL XX1 = Jac11*xx1 + Jac21*xx2 + Jac31*xx3;
        
        CCTK_REAL X = XX1;
        
        CCTK_REAL G11 = 1;
        
        CCTK_REAL G12 = 0;
        
        CCTK_REAL G13 = 0;
        
        CCTK_REAL G21 = 0;
        
        CCTK_REAL G22 = 1;
        
        CCTK_REAL G23 = 0;
        
        CCTK_REAL G31 = 0;
        
        CCTK_REAL G32 = 0;
        
        CCTK_REAL G33 = 1;
        
        G11 = 1 - amp*Sin(2*Pi*(X - cctk_time)*INV(period));
        
        CCTK_REAL K11 = 0;
        
        CCTK_REAL K12 = 0;
        
        CCTK_REAL K13 = 0;
        
        CCTK_REAL K21 = 0;
        
        CCTK_REAL K22 = 0;
        
        CCTK_REAL K23 = 0;
        
        CCTK_REAL K31 = 0;
        
        CCTK_REAL K32 = 0;
        
        CCTK_REAL K33 = 0;
        
        K11 = -(amp*Pi*Cos(2*Pi*(X - cctk_time)*INV(period))*INV(period)*pow(1 - 
          amp*Sin(2*Pi*(X - cctk_time)*INV(period)),-khalf));
        
        CCTK_REAL alpL = pow(1 - amp*Sin(2*Pi*(X - cctk_time)*INV(period)),0.5);
        
        CCTK_REAL betaxL = 0;
        
        CCTK_REAL betayL = 0;
        
        CCTK_REAL betazL = 0;
        
        CCTK_REAL gxxL = (G31*Jac11 + (G23 + G32)*Jac12)*Jac13 + Jac11*((G12 + 
          G21)*Jac12 + G13*Jac13) + G11*SQR(Jac11) + G22*SQR(Jac12) + 
          G33*SQR(Jac13);
        
        CCTK_REAL gxyL = (G11*Jac11 + G21*Jac12 + G31*Jac13)*Jac21 + 
          (G12*Jac11 + G22*Jac12 + G32*Jac13)*Jac22 + (G13*Jac11 + G23*Jac12 + 
          G33*Jac13)*Jac23;
        
        CCTK_REAL gxzL = (G11*Jac11 + G21*Jac12 + G31*Jac13)*Jac31 + 
          (G12*Jac11 + G22*Jac12 + G32*Jac13)*Jac32 + (G13*Jac11 + G23*Jac12 + 
          G33*Jac13)*Jac33;
        
        CCTK_REAL gyyL = (G31*Jac21 + (G23 + G32)*Jac22)*Jac23 + Jac21*((G12 + 
          G21)*Jac22 + G13*Jac23) + G11*SQR(Jac21) + G22*SQR(Jac22) + 
          G33*SQR(Jac23);
        
        CCTK_REAL gyzL = (G11*Jac21 + G21*Jac22 + G31*Jac23)*Jac31 + 
          (G12*Jac21 + G22*Jac22 + G32*Jac23)*Jac32 + (G13*Jac21 + G23*Jac22 + 
          G33*Jac23)*Jac33;
        
        CCTK_REAL gzzL = (G31*Jac31 + (G23 + G32)*Jac32)*Jac33 + Jac31*((G12 + 
          G21)*Jac32 + G13*Jac33) + G11*SQR(Jac31) + G22*SQR(Jac32) + 
          G33*SQR(Jac33);
        
        CCTK_REAL kxxL = Jac11*(Jac12*(K12 + K21) + Jac13*(K13 + K31)) + 
          Jac12*Jac13*(K23 + K32) + K11*SQR(Jac11) + K22*SQR(Jac12) + 
          K33*SQR(Jac13);
        
        CCTK_REAL kxyL = Jac11*(Jac21*K11 + Jac22*K12 + Jac23*K13) + 
          Jac12*(Jac21*K21 + Jac22*K22 + Jac23*K23) + Jac13*(Jac21*K31 + 
          Jac22*K32 + Jac23*K33);
        
        CCTK_REAL kxzL = Jac11*(Jac31*K11 + Jac32*K12 + Jac33*K13) + 
          Jac12*(Jac31*K21 + Jac32*K22 + Jac33*K23) + Jac13*(Jac31*K31 + 
          Jac32*K32 + Jac33*K33);
        
        CCTK_REAL kyyL = Jac21*(Jac22*(K12 + K21) + Jac23*(K13 + K31)) + 
          Jac22*Jac23*(K23 + K32) + K11*SQR(Jac21) + K22*SQR(Jac22) + 
          K33*SQR(Jac23);
        
        CCTK_REAL kyzL = Jac21*(Jac31*K11 + Jac32*K12 + Jac33*K13) + 
          Jac22*(Jac31*K21 + Jac32*K22 + Jac33*K23) + Jac23*(Jac31*K31 + 
          Jac32*K32 + Jac33*K33);
        
        CCTK_REAL kzzL = Jac31*(Jac32*(K12 + K21) + Jac33*(K13 + K31)) + 
          Jac32*Jac33*(K23 + K32) + K11*SQR(Jac31) + K22*SQR(Jac32) + 
          K33*SQR(Jac33);
        
        
        /* Copy local copies back to grid functions */
        alp[index] = alpL;
        betax[index] = betaxL;
        betay[index] = betayL;
        betaz[index] = betazL;
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
    }
  }
}

void gauge_wave_always(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &gauge_wave_always_Body);
}
