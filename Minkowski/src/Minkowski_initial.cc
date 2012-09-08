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

static void Minkowski_initial_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_LOOP3(Minkowski_initial,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL = alp[index];
    CCTK_REAL betaxL = betax[index];
    CCTK_REAL betayL = betay[index];
    CCTK_REAL betazL = betaz[index];
    CCTK_REAL dtbetaxL = dtbetax[index];
    CCTK_REAL dtbetayL = dtbetay[index];
    CCTK_REAL dtbetazL = dtbetaz[index];
    CCTK_REAL gxxL = gxx[index];
    CCTK_REAL gxyL = gxy[index];
    CCTK_REAL gxzL = gxz[index];
    CCTK_REAL gyyL = gyy[index];
    CCTK_REAL gyzL = gyz[index];
    CCTK_REAL gzzL = gzz[index];
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL csetemp0 = SQR(ToReal(boostx));
    
    CCTK_REAL csetemp1 = SQR(ToReal(boosty));
    
    CCTK_REAL csetemp2 = SQR(ToReal(boostz));
    
    CCTK_REAL invXform1L00 = -(INV((-1 + csetemp0 + csetemp1 + 
      csetemp2)*(1 + sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2))*ToReal(lapsefactor))*(1 - csetemp0 - csetemp1 - csetemp2 + 
      sqrt(1 - csetemp0 - csetemp1 - csetemp2))*(1 + 
      ToReal(boostx)*ToReal(shiftaddx) + ToReal(boosty)*ToReal(shiftaddy) + 
      ToReal(boostz)*ToReal(shiftaddz)));
    
    CCTK_REAL invXform1L01 = INV((-1 + csetemp0 + csetemp1 + csetemp2)*(1 
      + sqrt(1 - csetemp0 - csetemp1 - csetemp2))*ToReal(lapsefactor))*(1 - 
      csetemp0 - csetemp1 - csetemp2 + sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2))*ToReal(boostx);
    
    CCTK_REAL invXform1L02 = INV((-1 + csetemp0 + csetemp1 + csetemp2)*(1 
      + sqrt(1 - csetemp0 - csetemp1 - csetemp2))*ToReal(lapsefactor))*(1 - 
      csetemp0 - csetemp1 - csetemp2 + sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2))*ToReal(boosty);
    
    CCTK_REAL invXform1L03 = INV((-1 + csetemp0 + csetemp1 + csetemp2)*(1 
      + sqrt(1 - csetemp0 - csetemp1 - csetemp2))*ToReal(lapsefactor))*(1 - 
      csetemp0 - csetemp1 - csetemp2 + sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2))*ToReal(boostz);
    
    CCTK_REAL invXform1L10 = INV((-1 + csetemp0 + csetemp1 + csetemp2)*(1 
      + sqrt(1 - csetemp0 - csetemp1 - csetemp2)))*(-((csetemp0 + (-1 + 
      csetemp1 + csetemp2)*(1 + sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2)))*ToReal(shiftaddx)) + ToReal(boostx)*(1 - csetemp0 - 
      csetemp1 - csetemp2 + sqrt(1 - csetemp0 - csetemp1 - csetemp2) + sqrt(1 
      - csetemp0 - csetemp1 - csetemp2)*ToReal(boosty)*ToReal(shiftaddy) + 
      sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2)*ToReal(boostz)*ToReal(shiftaddz)));
    
    CCTK_REAL invXform1L11 = INV((-1 + csetemp0 + csetemp1 + csetemp2)*(1 
      + sqrt(1 - csetemp0 - csetemp1 - csetemp2)))*(csetemp0 + (-1 + csetemp1 
      + csetemp2)*(1 + sqrt(1 - csetemp0 - csetemp1 - csetemp2)));
    
    CCTK_REAL invXform1L12 = -(INV(-1 + csetemp0 + csetemp1 + csetemp2 - 
      sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2))*ToReal(boostx)*ToReal(boosty));
    
    CCTK_REAL invXform1L13 = -(INV(-1 + csetemp0 + csetemp1 + csetemp2 - 
      sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2))*ToReal(boostx)*ToReal(boostz));
    
    CCTK_REAL invXform1L20 = INV((-1 + csetemp0 + csetemp1 + csetemp2)*(1 
      + sqrt(1 - csetemp0 - csetemp1 - csetemp2)))*(-((csetemp1 + (-1 + 
      csetemp0 + csetemp2)*(1 + sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2)))*ToReal(shiftaddy)) + ToReal(boosty)*(1 - csetemp0 - 
      csetemp1 - csetemp2 + sqrt(1 - csetemp0 - csetemp1 - csetemp2) + sqrt(1 
      - csetemp0 - csetemp1 - csetemp2)*ToReal(boostx)*ToReal(shiftaddx) + 
      sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2)*ToReal(boostz)*ToReal(shiftaddz)));
    
    CCTK_REAL invXform1L21 = -(INV(-1 + csetemp0 + csetemp1 + csetemp2 - 
      sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2))*ToReal(boostx)*ToReal(boosty));
    
    CCTK_REAL invXform1L22 = INV((-1 + csetemp0 + csetemp1 + csetemp2)*(1 
      + sqrt(1 - csetemp0 - csetemp1 - csetemp2)))*(csetemp1 + (-1 + csetemp0 
      + csetemp2)*(1 + sqrt(1 - csetemp0 - csetemp1 - csetemp2)));
    
    CCTK_REAL invXform1L23 = -(INV(-1 + csetemp0 + csetemp1 + csetemp2 - 
      sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2))*ToReal(boosty)*ToReal(boostz));
    
    CCTK_REAL invXform1L30 = INV((-1 + csetemp0 + csetemp1 + csetemp2)*(1 
      + sqrt(1 - csetemp0 - csetemp1 - csetemp2)))*(ToReal(boostz)*(1 - 
      csetemp0 - csetemp1 - csetemp2 + sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2) + sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2)*ToReal(boostx)*ToReal(shiftaddx) + sqrt(1 - csetemp0 - 
      csetemp1 - csetemp2)*ToReal(boosty)*ToReal(shiftaddy)) - (-1 + csetemp0 
      + csetemp1 + csetemp2 + (-1 + csetemp0 + csetemp1)*sqrt(1 - csetemp0 - 
      csetemp1 - csetemp2))*ToReal(shiftaddz));
    
    CCTK_REAL invXform1L31 = -(INV(-1 + csetemp0 + csetemp1 + csetemp2 - 
      sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2))*ToReal(boostx)*ToReal(boostz));
    
    CCTK_REAL invXform1L32 = -(INV(-1 + csetemp0 + csetemp1 + csetemp2 - 
      sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2))*ToReal(boosty)*ToReal(boostz));
    
    CCTK_REAL invXform1L33 = INV((-1 + csetemp0 + csetemp1 + csetemp2)*(1 
      + sqrt(1 - csetemp0 - csetemp1 - csetemp2)))*(-1 + csetemp0 + csetemp1 
      + csetemp2 + (-1 + csetemp0 + csetemp1)*sqrt(1 - csetemp0 - csetemp1 - 
      csetemp2));
    
    CCTK_REAL invXform2L00 = 1;
    
    CCTK_REAL invXform2L01 = 0;
    
    CCTK_REAL invXform2L02 = 0;
    
    CCTK_REAL invXform2L03 = 0;
    
    CCTK_REAL invXform2L10 = 0;
    
    CCTK_REAL csetemp3 = cos(ToReal(phi));
    
    CCTK_REAL csetemp4 = cos(ToReal(psi));
    
    CCTK_REAL csetemp5 = cos(ToReal(theta));
    
    CCTK_REAL csetemp6 = sin(ToReal(phi));
    
    CCTK_REAL csetemp7 = sin(ToReal(psi));
    
    CCTK_REAL csetemp8 = sin(ToReal(theta));
    
    CCTK_REAL invXform2L11 = INV((SQR(csetemp3) + 
      SQR(csetemp6))*(SQR(csetemp4) + SQR(csetemp7))*(SQR(csetemp5) + 
      SQR(csetemp8)))*(-(csetemp5*csetemp6*csetemp7) + 
      csetemp3*csetemp4*(SQR(csetemp5) + SQR(csetemp8)));
    
    CCTK_REAL invXform2L12 = INV((SQR(csetemp3) + 
      SQR(csetemp6))*(SQR(csetemp4) + SQR(csetemp7))*(SQR(csetemp5) + 
      SQR(csetemp8)))*(csetemp3*csetemp5*csetemp7 + 
      csetemp4*csetemp6*(SQR(csetemp5) + SQR(csetemp8)));
    
    CCTK_REAL invXform2L13 = csetemp7*csetemp8*INV((SQR(csetemp4) + 
      SQR(csetemp7))*(SQR(csetemp5) + SQR(csetemp8)));
    
    CCTK_REAL invXform2L20 = 0;
    
    CCTK_REAL invXform2L21 = -(INV((SQR(csetemp3) + 
      SQR(csetemp6))*(SQR(csetemp4) + SQR(csetemp7))*(SQR(csetemp5) + 
      SQR(csetemp8)))*(csetemp4*csetemp5*csetemp6 + 
      csetemp3*csetemp7*(SQR(csetemp5) + SQR(csetemp8))));
    
    CCTK_REAL invXform2L22 = INV((SQR(csetemp3) + 
      SQR(csetemp6))*(SQR(csetemp4) + SQR(csetemp7))*(SQR(csetemp5) + 
      SQR(csetemp8)))*(csetemp3*csetemp4*csetemp5 - 
      csetemp6*csetemp7*(SQR(csetemp5) + SQR(csetemp8)));
    
    CCTK_REAL invXform2L23 = csetemp4*csetemp8*INV((SQR(csetemp4) + 
      SQR(csetemp7))*(SQR(csetemp5) + SQR(csetemp8)));
    
    CCTK_REAL invXform2L30 = 0;
    
    CCTK_REAL invXform2L31 = csetemp6*csetemp8*INV((SQR(csetemp3) + 
      SQR(csetemp6))*(SQR(csetemp5) + SQR(csetemp8)));
    
    CCTK_REAL invXform2L32 = -(csetemp3*csetemp8*INV((SQR(csetemp3) + 
      SQR(csetemp6))*(SQR(csetemp5) + SQR(csetemp8))));
    
    CCTK_REAL invXform2L33 = csetemp5*INV(SQR(csetemp5) + SQR(csetemp8));
    
    CCTK_REAL invXformL00 = invXform1L00*invXform2L00 + 
      invXform1L10*invXform2L01 + invXform1L20*invXform2L02 + 
      invXform1L30*invXform2L03;
    
    CCTK_REAL invXformL01 = invXform1L01*invXform2L00 + 
      invXform1L11*invXform2L01 + invXform1L21*invXform2L02 + 
      invXform1L31*invXform2L03;
    
    CCTK_REAL invXformL02 = invXform1L02*invXform2L00 + 
      invXform1L12*invXform2L01 + invXform1L22*invXform2L02 + 
      invXform1L32*invXform2L03;
    
    CCTK_REAL invXformL03 = invXform1L03*invXform2L00 + 
      invXform1L13*invXform2L01 + invXform1L23*invXform2L02 + 
      invXform1L33*invXform2L03;
    
    CCTK_REAL invXformL10 = invXform1L00*invXform2L10 + 
      invXform1L10*invXform2L11 + invXform1L20*invXform2L12 + 
      invXform1L30*invXform2L13;
    
    CCTK_REAL invXformL11 = invXform1L01*invXform2L10 + 
      invXform1L11*invXform2L11 + invXform1L21*invXform2L12 + 
      invXform1L31*invXform2L13;
    
    CCTK_REAL invXformL12 = invXform1L02*invXform2L10 + 
      invXform1L12*invXform2L11 + invXform1L22*invXform2L12 + 
      invXform1L32*invXform2L13;
    
    CCTK_REAL invXformL13 = invXform1L03*invXform2L10 + 
      invXform1L13*invXform2L11 + invXform1L23*invXform2L12 + 
      invXform1L33*invXform2L13;
    
    CCTK_REAL invXformL20 = invXform1L00*invXform2L20 + 
      invXform1L10*invXform2L21 + invXform1L20*invXform2L22 + 
      invXform1L30*invXform2L23;
    
    CCTK_REAL invXformL21 = invXform1L01*invXform2L20 + 
      invXform1L11*invXform2L21 + invXform1L21*invXform2L22 + 
      invXform1L31*invXform2L23;
    
    CCTK_REAL invXformL22 = invXform1L02*invXform2L20 + 
      invXform1L12*invXform2L21 + invXform1L22*invXform2L22 + 
      invXform1L32*invXform2L23;
    
    CCTK_REAL invXformL23 = invXform1L03*invXform2L20 + 
      invXform1L13*invXform2L21 + invXform1L23*invXform2L22 + 
      invXform1L33*invXform2L23;
    
    CCTK_REAL invXformL30 = invXform1L00*invXform2L30 + 
      invXform1L10*invXform2L31 + invXform1L20*invXform2L32 + 
      invXform1L30*invXform2L33;
    
    CCTK_REAL invXformL31 = invXform1L01*invXform2L30 + 
      invXform1L11*invXform2L31 + invXform1L21*invXform2L32 + 
      invXform1L31*invXform2L33;
    
    CCTK_REAL invXformL32 = invXform1L02*invXform2L30 + 
      invXform1L12*invXform2L31 + invXform1L22*invXform2L32 + 
      invXform1L32*invXform2L33;
    
    CCTK_REAL invXformL33 = invXform1L03*invXform2L30 + 
      invXform1L13*invXform2L31 + invXform1L23*invXform2L32 + 
      invXform1L33*invXform2L33;
    
    CCTK_REAL tg400 = -1;
    
    CCTK_REAL tg401 = 0;
    
    CCTK_REAL tg402 = 0;
    
    CCTK_REAL tg403 = 0;
    
    CCTK_REAL tg411 = 1;
    
    CCTK_REAL tg412 = 0;
    
    CCTK_REAL tg413 = 0;
    
    CCTK_REAL tg422 = 1;
    
    CCTK_REAL tg423 = 0;
    
    CCTK_REAL tg433 = 1;
    
    CCTK_REAL tdg4000 = 0;
    
    CCTK_REAL tdg4001 = 0;
    
    CCTK_REAL tdg4002 = 0;
    
    CCTK_REAL tdg4003 = 0;
    
    CCTK_REAL tdg4010 = 0;
    
    CCTK_REAL tdg4011 = 0;
    
    CCTK_REAL tdg4012 = 0;
    
    CCTK_REAL tdg4013 = 0;
    
    CCTK_REAL tdg4020 = 0;
    
    CCTK_REAL tdg4021 = 0;
    
    CCTK_REAL tdg4022 = 0;
    
    CCTK_REAL tdg4023 = 0;
    
    CCTK_REAL tdg4030 = 0;
    
    CCTK_REAL tdg4031 = 0;
    
    CCTK_REAL tdg4032 = 0;
    
    CCTK_REAL tdg4033 = 0;
    
    CCTK_REAL tdg4110 = 0;
    
    CCTK_REAL tdg4111 = 0;
    
    CCTK_REAL tdg4112 = 0;
    
    CCTK_REAL tdg4113 = 0;
    
    CCTK_REAL tdg4120 = 0;
    
    CCTK_REAL tdg4121 = 0;
    
    CCTK_REAL tdg4122 = 0;
    
    CCTK_REAL tdg4123 = 0;
    
    CCTK_REAL tdg4130 = 0;
    
    CCTK_REAL tdg4131 = 0;
    
    CCTK_REAL tdg4132 = 0;
    
    CCTK_REAL tdg4133 = 0;
    
    CCTK_REAL tdg4220 = 0;
    
    CCTK_REAL tdg4221 = 0;
    
    CCTK_REAL tdg4222 = 0;
    
    CCTK_REAL tdg4223 = 0;
    
    CCTK_REAL tdg4230 = 0;
    
    CCTK_REAL tdg4231 = 0;
    
    CCTK_REAL tdg4232 = 0;
    
    CCTK_REAL tdg4233 = 0;
    
    CCTK_REAL tdg4330 = 0;
    
    CCTK_REAL tdg4331 = 0;
    
    CCTK_REAL tdg4332 = 0;
    
    CCTK_REAL tdg4333 = 0;
    
    CCTK_REAL g400 = 2*(invXformL00*(invXformL10*tg401 + invXformL20*tg402 
      + invXformL30*tg403) + invXformL10*(invXformL20*tg412 + 
      invXformL30*tg413) + invXformL20*invXformL30*tg423) + 
      tg400*SQR(invXformL00) + tg411*SQR(invXformL10) + 
      tg422*SQR(invXformL20) + tg433*SQR(invXformL30);
    
    CCTK_REAL csetemp9 = invXformL01*tg400;
    
    CCTK_REAL csetemp10 = invXformL11*tg401;
    
    CCTK_REAL csetemp11 = invXformL21*tg402;
    
    CCTK_REAL csetemp12 = invXformL31*tg403;
    
    CCTK_REAL csetemp13 = invXformL01*tg401;
    
    CCTK_REAL csetemp14 = invXformL11*tg411;
    
    CCTK_REAL csetemp15 = invXformL21*tg412;
    
    CCTK_REAL csetemp16 = invXformL31*tg413;
    
    CCTK_REAL csetemp17 = invXformL01*tg402;
    
    CCTK_REAL csetemp18 = invXformL11*tg412;
    
    CCTK_REAL csetemp19 = invXformL21*tg422;
    
    CCTK_REAL csetemp20 = invXformL31*tg423;
    
    CCTK_REAL csetemp21 = invXformL01*tg403;
    
    CCTK_REAL csetemp22 = invXformL11*tg413;
    
    CCTK_REAL csetemp23 = invXformL21*tg423;
    
    CCTK_REAL csetemp24 = invXformL31*tg433;
    
    CCTK_REAL g401 = (csetemp10 + csetemp11 + csetemp12 + 
      csetemp9)*invXformL00 + (csetemp13 + csetemp14 + csetemp15 + 
      csetemp16)*invXformL10 + (csetemp17 + csetemp18 + csetemp19 + 
      csetemp20)*invXformL20 + (csetemp21 + csetemp22 + csetemp23 + 
      csetemp24)*invXformL30;
    
    CCTK_REAL csetemp25 = invXformL02*tg400;
    
    CCTK_REAL csetemp26 = invXformL12*tg401;
    
    CCTK_REAL csetemp27 = invXformL22*tg402;
    
    CCTK_REAL csetemp28 = invXformL32*tg403;
    
    CCTK_REAL csetemp29 = invXformL02*tg401;
    
    CCTK_REAL csetemp30 = invXformL12*tg411;
    
    CCTK_REAL csetemp31 = invXformL22*tg412;
    
    CCTK_REAL csetemp32 = invXformL32*tg413;
    
    CCTK_REAL csetemp33 = invXformL02*tg402;
    
    CCTK_REAL csetemp34 = invXformL12*tg412;
    
    CCTK_REAL csetemp35 = invXformL22*tg422;
    
    CCTK_REAL csetemp36 = invXformL32*tg423;
    
    CCTK_REAL csetemp37 = invXformL02*tg403;
    
    CCTK_REAL csetemp38 = invXformL12*tg413;
    
    CCTK_REAL csetemp39 = invXformL22*tg423;
    
    CCTK_REAL csetemp40 = invXformL32*tg433;
    
    CCTK_REAL g402 = (csetemp25 + csetemp26 + csetemp27 + 
      csetemp28)*invXformL00 + (csetemp29 + csetemp30 + csetemp31 + 
      csetemp32)*invXformL10 + (csetemp33 + csetemp34 + csetemp35 + 
      csetemp36)*invXformL20 + (csetemp37 + csetemp38 + csetemp39 + 
      csetemp40)*invXformL30;
    
    CCTK_REAL csetemp41 = invXformL03*tg400;
    
    CCTK_REAL csetemp42 = invXformL13*tg401;
    
    CCTK_REAL csetemp43 = invXformL23*tg402;
    
    CCTK_REAL csetemp44 = invXformL33*tg403;
    
    CCTK_REAL csetemp45 = invXformL03*tg401;
    
    CCTK_REAL csetemp46 = invXformL13*tg411;
    
    CCTK_REAL csetemp47 = invXformL23*tg412;
    
    CCTK_REAL csetemp48 = invXformL33*tg413;
    
    CCTK_REAL csetemp49 = invXformL03*tg402;
    
    CCTK_REAL csetemp50 = invXformL13*tg412;
    
    CCTK_REAL csetemp51 = invXformL23*tg422;
    
    CCTK_REAL csetemp52 = invXformL33*tg423;
    
    CCTK_REAL csetemp53 = invXformL03*tg403;
    
    CCTK_REAL csetemp54 = invXformL13*tg413;
    
    CCTK_REAL csetemp55 = invXformL23*tg423;
    
    CCTK_REAL csetemp56 = invXformL33*tg433;
    
    CCTK_REAL g403 = (csetemp41 + csetemp42 + csetemp43 + 
      csetemp44)*invXformL00 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*invXformL10 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*invXformL20 + (csetemp53 + csetemp54 + csetemp55 + 
      csetemp56)*invXformL30;
    
    CCTK_REAL g411 = (csetemp10 + csetemp11 + csetemp12 + 
      csetemp9)*invXformL01 + (csetemp13 + csetemp14 + csetemp15 + 
      csetemp16)*invXformL11 + (csetemp17 + csetemp18 + csetemp19 + 
      csetemp20)*invXformL21 + (csetemp21 + csetemp22 + csetemp23 + 
      csetemp24)*invXformL31;
    
    CCTK_REAL g412 = (csetemp25 + csetemp26 + csetemp27 + 
      csetemp28)*invXformL01 + (csetemp29 + csetemp30 + csetemp31 + 
      csetemp32)*invXformL11 + (csetemp33 + csetemp34 + csetemp35 + 
      csetemp36)*invXformL21 + (csetemp37 + csetemp38 + csetemp39 + 
      csetemp40)*invXformL31;
    
    CCTK_REAL g413 = (csetemp41 + csetemp42 + csetemp43 + 
      csetemp44)*invXformL01 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*invXformL11 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*invXformL21 + (csetemp53 + csetemp54 + csetemp55 + 
      csetemp56)*invXformL31;
    
    CCTK_REAL g422 = (csetemp25 + csetemp26 + csetemp27 + 
      csetemp28)*invXformL02 + (csetemp29 + csetemp30 + csetemp31 + 
      csetemp32)*invXformL12 + (csetemp33 + csetemp34 + csetemp35 + 
      csetemp36)*invXformL22 + (csetemp37 + csetemp38 + csetemp39 + 
      csetemp40)*invXformL32;
    
    CCTK_REAL g423 = (csetemp41 + csetemp42 + csetemp43 + 
      csetemp44)*invXformL02 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*invXformL12 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*invXformL22 + (csetemp53 + csetemp54 + csetemp55 + 
      csetemp56)*invXformL32;
    
    CCTK_REAL g433 = (csetemp41 + csetemp42 + csetemp43 + 
      csetemp44)*invXformL03 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*invXformL13 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*invXformL23 + (csetemp53 + csetemp54 + csetemp55 + 
      csetemp56)*invXformL33;
    
    CCTK_REAL csetemp57 = invXformL00*tdg4000;
    
    CCTK_REAL csetemp58 = invXformL10*tdg4001;
    
    CCTK_REAL csetemp59 = invXformL20*tdg4002;
    
    CCTK_REAL csetemp60 = invXformL30*tdg4003;
    
    CCTK_REAL csetemp61 = invXformL00*tdg4010;
    
    CCTK_REAL csetemp62 = invXformL10*tdg4011;
    
    CCTK_REAL csetemp63 = invXformL20*tdg4012;
    
    CCTK_REAL csetemp64 = invXformL30*tdg4013;
    
    CCTK_REAL csetemp65 = invXformL00*tdg4020;
    
    CCTK_REAL csetemp66 = invXformL10*tdg4021;
    
    CCTK_REAL csetemp67 = invXformL20*tdg4022;
    
    CCTK_REAL csetemp68 = invXformL30*tdg4023;
    
    CCTK_REAL csetemp69 = invXformL00*tdg4030;
    
    CCTK_REAL csetemp70 = invXformL10*tdg4031;
    
    CCTK_REAL csetemp71 = invXformL20*tdg4032;
    
    CCTK_REAL csetemp72 = invXformL30*tdg4033;
    
    CCTK_REAL csetemp73 = invXformL00*tdg4110;
    
    CCTK_REAL csetemp74 = invXformL10*tdg4111;
    
    CCTK_REAL csetemp75 = invXformL20*tdg4112;
    
    CCTK_REAL csetemp76 = invXformL30*tdg4113;
    
    CCTK_REAL csetemp77 = invXformL00*tdg4120;
    
    CCTK_REAL csetemp78 = invXformL10*tdg4121;
    
    CCTK_REAL csetemp79 = invXformL20*tdg4122;
    
    CCTK_REAL csetemp80 = invXformL30*tdg4123;
    
    CCTK_REAL csetemp81 = invXformL00*tdg4130;
    
    CCTK_REAL csetemp82 = invXformL10*tdg4131;
    
    CCTK_REAL csetemp83 = invXformL20*tdg4132;
    
    CCTK_REAL csetemp84 = invXformL30*tdg4133;
    
    CCTK_REAL csetemp85 = invXformL00*tdg4220;
    
    CCTK_REAL csetemp86 = invXformL10*tdg4221;
    
    CCTK_REAL csetemp87 = invXformL20*tdg4222;
    
    CCTK_REAL csetemp88 = invXformL30*tdg4223;
    
    CCTK_REAL csetemp89 = invXformL00*tdg4230;
    
    CCTK_REAL csetemp90 = invXformL10*tdg4231;
    
    CCTK_REAL csetemp91 = invXformL20*tdg4232;
    
    CCTK_REAL csetemp92 = invXformL30*tdg4233;
    
    CCTK_REAL csetemp93 = invXformL00*tdg4330;
    
    CCTK_REAL csetemp94 = invXformL10*tdg4331;
    
    CCTK_REAL csetemp95 = invXformL20*tdg4332;
    
    CCTK_REAL csetemp96 = invXformL30*tdg4333;
    
    CCTK_REAL dg4000 = invXformL00*((csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*invXformL00 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*invXformL10 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*invXformL20 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*invXformL30) + invXformL10*((csetemp61 + csetemp62 + 
      csetemp63 + csetemp64)*invXformL00 + (csetemp73 + csetemp74 + csetemp75 
      + csetemp76)*invXformL10 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*invXformL20 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*invXformL30) + invXformL20*((csetemp65 + csetemp66 + 
      csetemp67 + csetemp68)*invXformL00 + (csetemp77 + csetemp78 + csetemp79 
      + csetemp80)*invXformL10 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*invXformL20 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL30) + invXformL30*((csetemp69 + csetemp70 + 
      csetemp71 + csetemp72)*invXformL00 + (csetemp81 + csetemp82 + csetemp83 
      + csetemp84)*invXformL10 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL20 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*invXformL30);
    
    CCTK_REAL dg4010 = invXformL00*((csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*invXformL01 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*invXformL11 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*invXformL21 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*invXformL31) + invXformL10*((csetemp61 + csetemp62 + 
      csetemp63 + csetemp64)*invXformL01 + (csetemp73 + csetemp74 + csetemp75 
      + csetemp76)*invXformL11 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*invXformL21 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*invXformL31) + invXformL20*((csetemp65 + csetemp66 + 
      csetemp67 + csetemp68)*invXformL01 + (csetemp77 + csetemp78 + csetemp79 
      + csetemp80)*invXformL11 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*invXformL21 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL31) + invXformL30*((csetemp69 + csetemp70 + 
      csetemp71 + csetemp72)*invXformL01 + (csetemp81 + csetemp82 + csetemp83 
      + csetemp84)*invXformL11 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL21 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*invXformL31);
    
    CCTK_REAL csetemp97 = invXformL01*tdg4000;
    
    CCTK_REAL csetemp98 = invXformL11*tdg4001;
    
    CCTK_REAL csetemp99 = invXformL21*tdg4002;
    
    CCTK_REAL csetemp100 = invXformL31*tdg4003;
    
    CCTK_REAL csetemp101 = invXformL01*tdg4010;
    
    CCTK_REAL csetemp102 = invXformL11*tdg4011;
    
    CCTK_REAL csetemp103 = invXformL21*tdg4012;
    
    CCTK_REAL csetemp104 = invXformL31*tdg4013;
    
    CCTK_REAL csetemp105 = invXformL01*tdg4020;
    
    CCTK_REAL csetemp106 = invXformL11*tdg4021;
    
    CCTK_REAL csetemp107 = invXformL21*tdg4022;
    
    CCTK_REAL csetemp108 = invXformL31*tdg4023;
    
    CCTK_REAL csetemp109 = invXformL01*tdg4030;
    
    CCTK_REAL csetemp110 = invXformL11*tdg4031;
    
    CCTK_REAL csetemp111 = invXformL21*tdg4032;
    
    CCTK_REAL csetemp112 = invXformL31*tdg4033;
    
    CCTK_REAL csetemp113 = invXformL01*tdg4110;
    
    CCTK_REAL csetemp114 = invXformL11*tdg4111;
    
    CCTK_REAL csetemp115 = invXformL21*tdg4112;
    
    CCTK_REAL csetemp116 = invXformL31*tdg4113;
    
    CCTK_REAL csetemp117 = invXformL01*tdg4120;
    
    CCTK_REAL csetemp118 = invXformL11*tdg4121;
    
    CCTK_REAL csetemp119 = invXformL21*tdg4122;
    
    CCTK_REAL csetemp120 = invXformL31*tdg4123;
    
    CCTK_REAL csetemp121 = invXformL01*tdg4130;
    
    CCTK_REAL csetemp122 = invXformL11*tdg4131;
    
    CCTK_REAL csetemp123 = invXformL21*tdg4132;
    
    CCTK_REAL csetemp124 = invXformL31*tdg4133;
    
    CCTK_REAL csetemp125 = invXformL01*tdg4220;
    
    CCTK_REAL csetemp126 = invXformL11*tdg4221;
    
    CCTK_REAL csetemp127 = invXformL21*tdg4222;
    
    CCTK_REAL csetemp128 = invXformL31*tdg4223;
    
    CCTK_REAL csetemp129 = invXformL01*tdg4230;
    
    CCTK_REAL csetemp130 = invXformL11*tdg4231;
    
    CCTK_REAL csetemp131 = invXformL21*tdg4232;
    
    CCTK_REAL csetemp132 = invXformL31*tdg4233;
    
    CCTK_REAL csetemp133 = invXformL01*tdg4330;
    
    CCTK_REAL csetemp134 = invXformL11*tdg4331;
    
    CCTK_REAL csetemp135 = invXformL21*tdg4332;
    
    CCTK_REAL csetemp136 = invXformL31*tdg4333;
    
    CCTK_REAL dg4011 = invXformL00*((csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL01 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*invXformL11 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*invXformL21 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL31) + invXformL10*((csetemp101 + csetemp102 + 
      csetemp103 + csetemp104)*invXformL01 + (csetemp113 + csetemp114 + 
      csetemp115 + csetemp116)*invXformL11 + (csetemp117 + csetemp118 + 
      csetemp119 + csetemp120)*invXformL21 + (csetemp121 + csetemp122 + 
      csetemp123 + csetemp124)*invXformL31) + invXformL20*((csetemp105 + 
      csetemp106 + csetemp107 + csetemp108)*invXformL01 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*invXformL11 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*invXformL21 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*invXformL31) + 
      invXformL30*((csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL01 + (csetemp121 + csetemp122 + csetemp123 + 
      csetemp124)*invXformL11 + (csetemp129 + csetemp130 + csetemp131 + 
      csetemp132)*invXformL21 + (csetemp133 + csetemp134 + csetemp135 + 
      csetemp136)*invXformL31);
    
    CCTK_REAL csetemp137 = invXformL02*tdg4000;
    
    CCTK_REAL csetemp138 = invXformL12*tdg4001;
    
    CCTK_REAL csetemp139 = invXformL22*tdg4002;
    
    CCTK_REAL csetemp140 = invXformL32*tdg4003;
    
    CCTK_REAL csetemp141 = invXformL02*tdg4010;
    
    CCTK_REAL csetemp142 = invXformL12*tdg4011;
    
    CCTK_REAL csetemp143 = invXformL22*tdg4012;
    
    CCTK_REAL csetemp144 = invXformL32*tdg4013;
    
    CCTK_REAL csetemp145 = invXformL02*tdg4020;
    
    CCTK_REAL csetemp146 = invXformL12*tdg4021;
    
    CCTK_REAL csetemp147 = invXformL22*tdg4022;
    
    CCTK_REAL csetemp148 = invXformL32*tdg4023;
    
    CCTK_REAL csetemp149 = invXformL02*tdg4030;
    
    CCTK_REAL csetemp150 = invXformL12*tdg4031;
    
    CCTK_REAL csetemp151 = invXformL22*tdg4032;
    
    CCTK_REAL csetemp152 = invXformL32*tdg4033;
    
    CCTK_REAL csetemp153 = invXformL02*tdg4110;
    
    CCTK_REAL csetemp154 = invXformL12*tdg4111;
    
    CCTK_REAL csetemp155 = invXformL22*tdg4112;
    
    CCTK_REAL csetemp156 = invXformL32*tdg4113;
    
    CCTK_REAL csetemp157 = invXformL02*tdg4120;
    
    CCTK_REAL csetemp158 = invXformL12*tdg4121;
    
    CCTK_REAL csetemp159 = invXformL22*tdg4122;
    
    CCTK_REAL csetemp160 = invXformL32*tdg4123;
    
    CCTK_REAL csetemp161 = invXformL02*tdg4130;
    
    CCTK_REAL csetemp162 = invXformL12*tdg4131;
    
    CCTK_REAL csetemp163 = invXformL22*tdg4132;
    
    CCTK_REAL csetemp164 = invXformL32*tdg4133;
    
    CCTK_REAL csetemp165 = invXformL02*tdg4220;
    
    CCTK_REAL csetemp166 = invXformL12*tdg4221;
    
    CCTK_REAL csetemp167 = invXformL22*tdg4222;
    
    CCTK_REAL csetemp168 = invXformL32*tdg4223;
    
    CCTK_REAL csetemp169 = invXformL02*tdg4230;
    
    CCTK_REAL csetemp170 = invXformL12*tdg4231;
    
    CCTK_REAL csetemp171 = invXformL22*tdg4232;
    
    CCTK_REAL csetemp172 = invXformL32*tdg4233;
    
    CCTK_REAL csetemp173 = invXformL02*tdg4330;
    
    CCTK_REAL csetemp174 = invXformL12*tdg4331;
    
    CCTK_REAL csetemp175 = invXformL22*tdg4332;
    
    CCTK_REAL csetemp176 = invXformL32*tdg4333;
    
    CCTK_REAL dg4012 = invXformL00*((csetemp137 + csetemp138 + csetemp139 
      + csetemp140)*invXformL01 + (csetemp141 + csetemp142 + csetemp143 + 
      csetemp144)*invXformL11 + (csetemp145 + csetemp146 + csetemp147 + 
      csetemp148)*invXformL21 + (csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL31) + invXformL10*((csetemp141 + csetemp142 + 
      csetemp143 + csetemp144)*invXformL01 + (csetemp153 + csetemp154 + 
      csetemp155 + csetemp156)*invXformL11 + (csetemp157 + csetemp158 + 
      csetemp159 + csetemp160)*invXformL21 + (csetemp161 + csetemp162 + 
      csetemp163 + csetemp164)*invXformL31) + invXformL20*((csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*invXformL01 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*invXformL11 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*invXformL21 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*invXformL31) + 
      invXformL30*((csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL01 + (csetemp161 + csetemp162 + csetemp163 + 
      csetemp164)*invXformL11 + (csetemp169 + csetemp170 + csetemp171 + 
      csetemp172)*invXformL21 + (csetemp173 + csetemp174 + csetemp175 + 
      csetemp176)*invXformL31);
    
    CCTK_REAL csetemp177 = invXformL03*tdg4000;
    
    CCTK_REAL csetemp178 = invXformL13*tdg4001;
    
    CCTK_REAL csetemp179 = invXformL23*tdg4002;
    
    CCTK_REAL csetemp180 = invXformL33*tdg4003;
    
    CCTK_REAL csetemp181 = invXformL03*tdg4010;
    
    CCTK_REAL csetemp182 = invXformL13*tdg4011;
    
    CCTK_REAL csetemp183 = invXformL23*tdg4012;
    
    CCTK_REAL csetemp184 = invXformL33*tdg4013;
    
    CCTK_REAL csetemp185 = invXformL03*tdg4020;
    
    CCTK_REAL csetemp186 = invXformL13*tdg4021;
    
    CCTK_REAL csetemp187 = invXformL23*tdg4022;
    
    CCTK_REAL csetemp188 = invXformL33*tdg4023;
    
    CCTK_REAL csetemp189 = invXformL03*tdg4030;
    
    CCTK_REAL csetemp190 = invXformL13*tdg4031;
    
    CCTK_REAL csetemp191 = invXformL23*tdg4032;
    
    CCTK_REAL csetemp192 = invXformL33*tdg4033;
    
    CCTK_REAL csetemp193 = invXformL03*tdg4110;
    
    CCTK_REAL csetemp194 = invXformL13*tdg4111;
    
    CCTK_REAL csetemp195 = invXformL23*tdg4112;
    
    CCTK_REAL csetemp196 = invXformL33*tdg4113;
    
    CCTK_REAL csetemp197 = invXformL03*tdg4120;
    
    CCTK_REAL csetemp198 = invXformL13*tdg4121;
    
    CCTK_REAL csetemp199 = invXformL23*tdg4122;
    
    CCTK_REAL csetemp200 = invXformL33*tdg4123;
    
    CCTK_REAL csetemp201 = invXformL03*tdg4130;
    
    CCTK_REAL csetemp202 = invXformL13*tdg4131;
    
    CCTK_REAL csetemp203 = invXformL23*tdg4132;
    
    CCTK_REAL csetemp204 = invXformL33*tdg4133;
    
    CCTK_REAL csetemp205 = invXformL03*tdg4220;
    
    CCTK_REAL csetemp206 = invXformL13*tdg4221;
    
    CCTK_REAL csetemp207 = invXformL23*tdg4222;
    
    CCTK_REAL csetemp208 = invXformL33*tdg4223;
    
    CCTK_REAL csetemp209 = invXformL03*tdg4230;
    
    CCTK_REAL csetemp210 = invXformL13*tdg4231;
    
    CCTK_REAL csetemp211 = invXformL23*tdg4232;
    
    CCTK_REAL csetemp212 = invXformL33*tdg4233;
    
    CCTK_REAL csetemp213 = invXformL03*tdg4330;
    
    CCTK_REAL csetemp214 = invXformL13*tdg4331;
    
    CCTK_REAL csetemp215 = invXformL23*tdg4332;
    
    CCTK_REAL csetemp216 = invXformL33*tdg4333;
    
    CCTK_REAL dg4013 = invXformL00*((csetemp177 + csetemp178 + csetemp179 
      + csetemp180)*invXformL01 + (csetemp181 + csetemp182 + csetemp183 + 
      csetemp184)*invXformL11 + (csetemp185 + csetemp186 + csetemp187 + 
      csetemp188)*invXformL21 + (csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL31) + invXformL10*((csetemp181 + csetemp182 + 
      csetemp183 + csetemp184)*invXformL01 + (csetemp193 + csetemp194 + 
      csetemp195 + csetemp196)*invXformL11 + (csetemp197 + csetemp198 + 
      csetemp199 + csetemp200)*invXformL21 + (csetemp201 + csetemp202 + 
      csetemp203 + csetemp204)*invXformL31) + invXformL20*((csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*invXformL01 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*invXformL11 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*invXformL21 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*invXformL31) + 
      invXformL30*((csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL01 + (csetemp201 + csetemp202 + csetemp203 + 
      csetemp204)*invXformL11 + (csetemp209 + csetemp210 + csetemp211 + 
      csetemp212)*invXformL21 + (csetemp213 + csetemp214 + csetemp215 + 
      csetemp216)*invXformL31);
    
    CCTK_REAL dg4020 = invXformL00*((csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*invXformL02 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*invXformL12 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*invXformL22 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*invXformL32) + invXformL10*((csetemp61 + csetemp62 + 
      csetemp63 + csetemp64)*invXformL02 + (csetemp73 + csetemp74 + csetemp75 
      + csetemp76)*invXformL12 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*invXformL22 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*invXformL32) + invXformL20*((csetemp65 + csetemp66 + 
      csetemp67 + csetemp68)*invXformL02 + (csetemp77 + csetemp78 + csetemp79 
      + csetemp80)*invXformL12 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*invXformL22 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL32) + invXformL30*((csetemp69 + csetemp70 + 
      csetemp71 + csetemp72)*invXformL02 + (csetemp81 + csetemp82 + csetemp83 
      + csetemp84)*invXformL12 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL22 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*invXformL32);
    
    CCTK_REAL dg4021 = invXformL00*((csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL02 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*invXformL12 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*invXformL22 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL32) + invXformL10*((csetemp101 + csetemp102 + 
      csetemp103 + csetemp104)*invXformL02 + (csetemp113 + csetemp114 + 
      csetemp115 + csetemp116)*invXformL12 + (csetemp117 + csetemp118 + 
      csetemp119 + csetemp120)*invXformL22 + (csetemp121 + csetemp122 + 
      csetemp123 + csetemp124)*invXformL32) + invXformL20*((csetemp105 + 
      csetemp106 + csetemp107 + csetemp108)*invXformL02 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*invXformL12 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*invXformL22 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*invXformL32) + 
      invXformL30*((csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL02 + (csetemp121 + csetemp122 + csetemp123 + 
      csetemp124)*invXformL12 + (csetemp129 + csetemp130 + csetemp131 + 
      csetemp132)*invXformL22 + (csetemp133 + csetemp134 + csetemp135 + 
      csetemp136)*invXformL32);
    
    CCTK_REAL dg4022 = invXformL00*((csetemp137 + csetemp138 + csetemp139 
      + csetemp140)*invXformL02 + (csetemp141 + csetemp142 + csetemp143 + 
      csetemp144)*invXformL12 + (csetemp145 + csetemp146 + csetemp147 + 
      csetemp148)*invXformL22 + (csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL32) + invXformL10*((csetemp141 + csetemp142 + 
      csetemp143 + csetemp144)*invXformL02 + (csetemp153 + csetemp154 + 
      csetemp155 + csetemp156)*invXformL12 + (csetemp157 + csetemp158 + 
      csetemp159 + csetemp160)*invXformL22 + (csetemp161 + csetemp162 + 
      csetemp163 + csetemp164)*invXformL32) + invXformL20*((csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*invXformL02 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*invXformL12 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*invXformL22 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*invXformL32) + 
      invXformL30*((csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL02 + (csetemp161 + csetemp162 + csetemp163 + 
      csetemp164)*invXformL12 + (csetemp169 + csetemp170 + csetemp171 + 
      csetemp172)*invXformL22 + (csetemp173 + csetemp174 + csetemp175 + 
      csetemp176)*invXformL32);
    
    CCTK_REAL dg4023 = invXformL00*((csetemp177 + csetemp178 + csetemp179 
      + csetemp180)*invXformL02 + (csetemp181 + csetemp182 + csetemp183 + 
      csetemp184)*invXformL12 + (csetemp185 + csetemp186 + csetemp187 + 
      csetemp188)*invXformL22 + (csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL32) + invXformL10*((csetemp181 + csetemp182 + 
      csetemp183 + csetemp184)*invXformL02 + (csetemp193 + csetemp194 + 
      csetemp195 + csetemp196)*invXformL12 + (csetemp197 + csetemp198 + 
      csetemp199 + csetemp200)*invXformL22 + (csetemp201 + csetemp202 + 
      csetemp203 + csetemp204)*invXformL32) + invXformL20*((csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*invXformL02 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*invXformL12 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*invXformL22 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*invXformL32) + 
      invXformL30*((csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL02 + (csetemp201 + csetemp202 + csetemp203 + 
      csetemp204)*invXformL12 + (csetemp209 + csetemp210 + csetemp211 + 
      csetemp212)*invXformL22 + (csetemp213 + csetemp214 + csetemp215 + 
      csetemp216)*invXformL32);
    
    CCTK_REAL dg4030 = invXformL00*((csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*invXformL03 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*invXformL13 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*invXformL23 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*invXformL33) + invXformL10*((csetemp61 + csetemp62 + 
      csetemp63 + csetemp64)*invXformL03 + (csetemp73 + csetemp74 + csetemp75 
      + csetemp76)*invXformL13 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*invXformL23 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*invXformL33) + invXformL20*((csetemp65 + csetemp66 + 
      csetemp67 + csetemp68)*invXformL03 + (csetemp77 + csetemp78 + csetemp79 
      + csetemp80)*invXformL13 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*invXformL23 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL33) + invXformL30*((csetemp69 + csetemp70 + 
      csetemp71 + csetemp72)*invXformL03 + (csetemp81 + csetemp82 + csetemp83 
      + csetemp84)*invXformL13 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL23 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*invXformL33);
    
    CCTK_REAL dg4031 = invXformL00*((csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL03 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*invXformL13 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*invXformL23 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL33) + invXformL10*((csetemp101 + csetemp102 + 
      csetemp103 + csetemp104)*invXformL03 + (csetemp113 + csetemp114 + 
      csetemp115 + csetemp116)*invXformL13 + (csetemp117 + csetemp118 + 
      csetemp119 + csetemp120)*invXformL23 + (csetemp121 + csetemp122 + 
      csetemp123 + csetemp124)*invXformL33) + invXformL20*((csetemp105 + 
      csetemp106 + csetemp107 + csetemp108)*invXformL03 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*invXformL13 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*invXformL23 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*invXformL33) + 
      invXformL30*((csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL03 + (csetemp121 + csetemp122 + csetemp123 + 
      csetemp124)*invXformL13 + (csetemp129 + csetemp130 + csetemp131 + 
      csetemp132)*invXformL23 + (csetemp133 + csetemp134 + csetemp135 + 
      csetemp136)*invXformL33);
    
    CCTK_REAL dg4032 = invXformL00*((csetemp137 + csetemp138 + csetemp139 
      + csetemp140)*invXformL03 + (csetemp141 + csetemp142 + csetemp143 + 
      csetemp144)*invXformL13 + (csetemp145 + csetemp146 + csetemp147 + 
      csetemp148)*invXformL23 + (csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL33) + invXformL10*((csetemp141 + csetemp142 + 
      csetemp143 + csetemp144)*invXformL03 + (csetemp153 + csetemp154 + 
      csetemp155 + csetemp156)*invXformL13 + (csetemp157 + csetemp158 + 
      csetemp159 + csetemp160)*invXformL23 + (csetemp161 + csetemp162 + 
      csetemp163 + csetemp164)*invXformL33) + invXformL20*((csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*invXformL03 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*invXformL13 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*invXformL23 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*invXformL33) + 
      invXformL30*((csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL03 + (csetemp161 + csetemp162 + csetemp163 + 
      csetemp164)*invXformL13 + (csetemp169 + csetemp170 + csetemp171 + 
      csetemp172)*invXformL23 + (csetemp173 + csetemp174 + csetemp175 + 
      csetemp176)*invXformL33);
    
    CCTK_REAL dg4033 = invXformL00*((csetemp177 + csetemp178 + csetemp179 
      + csetemp180)*invXformL03 + (csetemp181 + csetemp182 + csetemp183 + 
      csetemp184)*invXformL13 + (csetemp185 + csetemp186 + csetemp187 + 
      csetemp188)*invXformL23 + (csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL33) + invXformL10*((csetemp181 + csetemp182 + 
      csetemp183 + csetemp184)*invXformL03 + (csetemp193 + csetemp194 + 
      csetemp195 + csetemp196)*invXformL13 + (csetemp197 + csetemp198 + 
      csetemp199 + csetemp200)*invXformL23 + (csetemp201 + csetemp202 + 
      csetemp203 + csetemp204)*invXformL33) + invXformL20*((csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*invXformL03 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*invXformL13 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*invXformL23 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*invXformL33) + 
      invXformL30*((csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL03 + (csetemp201 + csetemp202 + csetemp203 + 
      csetemp204)*invXformL13 + (csetemp209 + csetemp210 + csetemp211 + 
      csetemp212)*invXformL23 + (csetemp213 + csetemp214 + csetemp215 + 
      csetemp216)*invXformL33);
    
    CCTK_REAL dg4110 = invXformL01*((csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*invXformL01 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*invXformL11 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*invXformL21 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*invXformL31) + invXformL11*((csetemp61 + csetemp62 + 
      csetemp63 + csetemp64)*invXformL01 + (csetemp73 + csetemp74 + csetemp75 
      + csetemp76)*invXformL11 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*invXformL21 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*invXformL31) + invXformL21*((csetemp65 + csetemp66 + 
      csetemp67 + csetemp68)*invXformL01 + (csetemp77 + csetemp78 + csetemp79 
      + csetemp80)*invXformL11 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*invXformL21 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL31) + invXformL31*((csetemp69 + csetemp70 + 
      csetemp71 + csetemp72)*invXformL01 + (csetemp81 + csetemp82 + csetemp83 
      + csetemp84)*invXformL11 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL21 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*invXformL31);
    
    CCTK_REAL dg4111 = invXformL01*((csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL01 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*invXformL11 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*invXformL21 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL31) + invXformL11*((csetemp101 + csetemp102 + 
      csetemp103 + csetemp104)*invXformL01 + (csetemp113 + csetemp114 + 
      csetemp115 + csetemp116)*invXformL11 + (csetemp117 + csetemp118 + 
      csetemp119 + csetemp120)*invXformL21 + (csetemp121 + csetemp122 + 
      csetemp123 + csetemp124)*invXformL31) + invXformL21*((csetemp105 + 
      csetemp106 + csetemp107 + csetemp108)*invXformL01 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*invXformL11 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*invXformL21 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*invXformL31) + 
      invXformL31*((csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL01 + (csetemp121 + csetemp122 + csetemp123 + 
      csetemp124)*invXformL11 + (csetemp129 + csetemp130 + csetemp131 + 
      csetemp132)*invXformL21 + (csetemp133 + csetemp134 + csetemp135 + 
      csetemp136)*invXformL31);
    
    CCTK_REAL dg4112 = invXformL01*((csetemp137 + csetemp138 + csetemp139 
      + csetemp140)*invXformL01 + (csetemp141 + csetemp142 + csetemp143 + 
      csetemp144)*invXformL11 + (csetemp145 + csetemp146 + csetemp147 + 
      csetemp148)*invXformL21 + (csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL31) + invXformL11*((csetemp141 + csetemp142 + 
      csetemp143 + csetemp144)*invXformL01 + (csetemp153 + csetemp154 + 
      csetemp155 + csetemp156)*invXformL11 + (csetemp157 + csetemp158 + 
      csetemp159 + csetemp160)*invXformL21 + (csetemp161 + csetemp162 + 
      csetemp163 + csetemp164)*invXformL31) + invXformL21*((csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*invXformL01 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*invXformL11 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*invXformL21 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*invXformL31) + 
      invXformL31*((csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL01 + (csetemp161 + csetemp162 + csetemp163 + 
      csetemp164)*invXformL11 + (csetemp169 + csetemp170 + csetemp171 + 
      csetemp172)*invXformL21 + (csetemp173 + csetemp174 + csetemp175 + 
      csetemp176)*invXformL31);
    
    CCTK_REAL dg4113 = invXformL01*((csetemp177 + csetemp178 + csetemp179 
      + csetemp180)*invXformL01 + (csetemp181 + csetemp182 + csetemp183 + 
      csetemp184)*invXformL11 + (csetemp185 + csetemp186 + csetemp187 + 
      csetemp188)*invXformL21 + (csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL31) + invXformL11*((csetemp181 + csetemp182 + 
      csetemp183 + csetemp184)*invXformL01 + (csetemp193 + csetemp194 + 
      csetemp195 + csetemp196)*invXformL11 + (csetemp197 + csetemp198 + 
      csetemp199 + csetemp200)*invXformL21 + (csetemp201 + csetemp202 + 
      csetemp203 + csetemp204)*invXformL31) + invXformL21*((csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*invXformL01 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*invXformL11 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*invXformL21 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*invXformL31) + 
      invXformL31*((csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL01 + (csetemp201 + csetemp202 + csetemp203 + 
      csetemp204)*invXformL11 + (csetemp209 + csetemp210 + csetemp211 + 
      csetemp212)*invXformL21 + (csetemp213 + csetemp214 + csetemp215 + 
      csetemp216)*invXformL31);
    
    CCTK_REAL dg4120 = invXformL01*((csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*invXformL02 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*invXformL12 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*invXformL22 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*invXformL32) + invXformL11*((csetemp61 + csetemp62 + 
      csetemp63 + csetemp64)*invXformL02 + (csetemp73 + csetemp74 + csetemp75 
      + csetemp76)*invXformL12 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*invXformL22 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*invXformL32) + invXformL21*((csetemp65 + csetemp66 + 
      csetemp67 + csetemp68)*invXformL02 + (csetemp77 + csetemp78 + csetemp79 
      + csetemp80)*invXformL12 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*invXformL22 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL32) + invXformL31*((csetemp69 + csetemp70 + 
      csetemp71 + csetemp72)*invXformL02 + (csetemp81 + csetemp82 + csetemp83 
      + csetemp84)*invXformL12 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL22 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*invXformL32);
    
    CCTK_REAL dg4121 = invXformL01*((csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL02 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*invXformL12 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*invXformL22 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL32) + invXformL11*((csetemp101 + csetemp102 + 
      csetemp103 + csetemp104)*invXformL02 + (csetemp113 + csetemp114 + 
      csetemp115 + csetemp116)*invXformL12 + (csetemp117 + csetemp118 + 
      csetemp119 + csetemp120)*invXformL22 + (csetemp121 + csetemp122 + 
      csetemp123 + csetemp124)*invXformL32) + invXformL21*((csetemp105 + 
      csetemp106 + csetemp107 + csetemp108)*invXformL02 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*invXformL12 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*invXformL22 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*invXformL32) + 
      invXformL31*((csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL02 + (csetemp121 + csetemp122 + csetemp123 + 
      csetemp124)*invXformL12 + (csetemp129 + csetemp130 + csetemp131 + 
      csetemp132)*invXformL22 + (csetemp133 + csetemp134 + csetemp135 + 
      csetemp136)*invXformL32);
    
    CCTK_REAL dg4122 = invXformL01*((csetemp137 + csetemp138 + csetemp139 
      + csetemp140)*invXformL02 + (csetemp141 + csetemp142 + csetemp143 + 
      csetemp144)*invXformL12 + (csetemp145 + csetemp146 + csetemp147 + 
      csetemp148)*invXformL22 + (csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL32) + invXformL11*((csetemp141 + csetemp142 + 
      csetemp143 + csetemp144)*invXformL02 + (csetemp153 + csetemp154 + 
      csetemp155 + csetemp156)*invXformL12 + (csetemp157 + csetemp158 + 
      csetemp159 + csetemp160)*invXformL22 + (csetemp161 + csetemp162 + 
      csetemp163 + csetemp164)*invXformL32) + invXformL21*((csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*invXformL02 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*invXformL12 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*invXformL22 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*invXformL32) + 
      invXformL31*((csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL02 + (csetemp161 + csetemp162 + csetemp163 + 
      csetemp164)*invXformL12 + (csetemp169 + csetemp170 + csetemp171 + 
      csetemp172)*invXformL22 + (csetemp173 + csetemp174 + csetemp175 + 
      csetemp176)*invXformL32);
    
    CCTK_REAL dg4123 = invXformL01*((csetemp177 + csetemp178 + csetemp179 
      + csetemp180)*invXformL02 + (csetemp181 + csetemp182 + csetemp183 + 
      csetemp184)*invXformL12 + (csetemp185 + csetemp186 + csetemp187 + 
      csetemp188)*invXformL22 + (csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL32) + invXformL11*((csetemp181 + csetemp182 + 
      csetemp183 + csetemp184)*invXformL02 + (csetemp193 + csetemp194 + 
      csetemp195 + csetemp196)*invXformL12 + (csetemp197 + csetemp198 + 
      csetemp199 + csetemp200)*invXformL22 + (csetemp201 + csetemp202 + 
      csetemp203 + csetemp204)*invXformL32) + invXformL21*((csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*invXformL02 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*invXformL12 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*invXformL22 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*invXformL32) + 
      invXformL31*((csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL02 + (csetemp201 + csetemp202 + csetemp203 + 
      csetemp204)*invXformL12 + (csetemp209 + csetemp210 + csetemp211 + 
      csetemp212)*invXformL22 + (csetemp213 + csetemp214 + csetemp215 + 
      csetemp216)*invXformL32);
    
    CCTK_REAL dg4130 = invXformL01*((csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*invXformL03 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*invXformL13 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*invXformL23 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*invXformL33) + invXformL11*((csetemp61 + csetemp62 + 
      csetemp63 + csetemp64)*invXformL03 + (csetemp73 + csetemp74 + csetemp75 
      + csetemp76)*invXformL13 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*invXformL23 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*invXformL33) + invXformL21*((csetemp65 + csetemp66 + 
      csetemp67 + csetemp68)*invXformL03 + (csetemp77 + csetemp78 + csetemp79 
      + csetemp80)*invXformL13 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*invXformL23 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL33) + invXformL31*((csetemp69 + csetemp70 + 
      csetemp71 + csetemp72)*invXformL03 + (csetemp81 + csetemp82 + csetemp83 
      + csetemp84)*invXformL13 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL23 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*invXformL33);
    
    CCTK_REAL dg4131 = invXformL01*((csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL03 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*invXformL13 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*invXformL23 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL33) + invXformL11*((csetemp101 + csetemp102 + 
      csetemp103 + csetemp104)*invXformL03 + (csetemp113 + csetemp114 + 
      csetemp115 + csetemp116)*invXformL13 + (csetemp117 + csetemp118 + 
      csetemp119 + csetemp120)*invXformL23 + (csetemp121 + csetemp122 + 
      csetemp123 + csetemp124)*invXformL33) + invXformL21*((csetemp105 + 
      csetemp106 + csetemp107 + csetemp108)*invXformL03 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*invXformL13 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*invXformL23 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*invXformL33) + 
      invXformL31*((csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL03 + (csetemp121 + csetemp122 + csetemp123 + 
      csetemp124)*invXformL13 + (csetemp129 + csetemp130 + csetemp131 + 
      csetemp132)*invXformL23 + (csetemp133 + csetemp134 + csetemp135 + 
      csetemp136)*invXformL33);
    
    CCTK_REAL dg4132 = invXformL01*((csetemp137 + csetemp138 + csetemp139 
      + csetemp140)*invXformL03 + (csetemp141 + csetemp142 + csetemp143 + 
      csetemp144)*invXformL13 + (csetemp145 + csetemp146 + csetemp147 + 
      csetemp148)*invXformL23 + (csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL33) + invXformL11*((csetemp141 + csetemp142 + 
      csetemp143 + csetemp144)*invXformL03 + (csetemp153 + csetemp154 + 
      csetemp155 + csetemp156)*invXformL13 + (csetemp157 + csetemp158 + 
      csetemp159 + csetemp160)*invXformL23 + (csetemp161 + csetemp162 + 
      csetemp163 + csetemp164)*invXformL33) + invXformL21*((csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*invXformL03 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*invXformL13 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*invXformL23 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*invXformL33) + 
      invXformL31*((csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL03 + (csetemp161 + csetemp162 + csetemp163 + 
      csetemp164)*invXformL13 + (csetemp169 + csetemp170 + csetemp171 + 
      csetemp172)*invXformL23 + (csetemp173 + csetemp174 + csetemp175 + 
      csetemp176)*invXformL33);
    
    CCTK_REAL dg4133 = invXformL01*((csetemp177 + csetemp178 + csetemp179 
      + csetemp180)*invXformL03 + (csetemp181 + csetemp182 + csetemp183 + 
      csetemp184)*invXformL13 + (csetemp185 + csetemp186 + csetemp187 + 
      csetemp188)*invXformL23 + (csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL33) + invXformL11*((csetemp181 + csetemp182 + 
      csetemp183 + csetemp184)*invXformL03 + (csetemp193 + csetemp194 + 
      csetemp195 + csetemp196)*invXformL13 + (csetemp197 + csetemp198 + 
      csetemp199 + csetemp200)*invXformL23 + (csetemp201 + csetemp202 + 
      csetemp203 + csetemp204)*invXformL33) + invXformL21*((csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*invXformL03 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*invXformL13 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*invXformL23 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*invXformL33) + 
      invXformL31*((csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL03 + (csetemp201 + csetemp202 + csetemp203 + 
      csetemp204)*invXformL13 + (csetemp209 + csetemp210 + csetemp211 + 
      csetemp212)*invXformL23 + (csetemp213 + csetemp214 + csetemp215 + 
      csetemp216)*invXformL33);
    
    CCTK_REAL dg4220 = invXformL02*((csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*invXformL02 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*invXformL12 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*invXformL22 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*invXformL32) + invXformL12*((csetemp61 + csetemp62 + 
      csetemp63 + csetemp64)*invXformL02 + (csetemp73 + csetemp74 + csetemp75 
      + csetemp76)*invXformL12 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*invXformL22 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*invXformL32) + invXformL22*((csetemp65 + csetemp66 + 
      csetemp67 + csetemp68)*invXformL02 + (csetemp77 + csetemp78 + csetemp79 
      + csetemp80)*invXformL12 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*invXformL22 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL32) + invXformL32*((csetemp69 + csetemp70 + 
      csetemp71 + csetemp72)*invXformL02 + (csetemp81 + csetemp82 + csetemp83 
      + csetemp84)*invXformL12 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL22 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*invXformL32);
    
    CCTK_REAL dg4221 = invXformL02*((csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL02 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*invXformL12 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*invXformL22 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL32) + invXformL12*((csetemp101 + csetemp102 + 
      csetemp103 + csetemp104)*invXformL02 + (csetemp113 + csetemp114 + 
      csetemp115 + csetemp116)*invXformL12 + (csetemp117 + csetemp118 + 
      csetemp119 + csetemp120)*invXformL22 + (csetemp121 + csetemp122 + 
      csetemp123 + csetemp124)*invXformL32) + invXformL22*((csetemp105 + 
      csetemp106 + csetemp107 + csetemp108)*invXformL02 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*invXformL12 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*invXformL22 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*invXformL32) + 
      invXformL32*((csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL02 + (csetemp121 + csetemp122 + csetemp123 + 
      csetemp124)*invXformL12 + (csetemp129 + csetemp130 + csetemp131 + 
      csetemp132)*invXformL22 + (csetemp133 + csetemp134 + csetemp135 + 
      csetemp136)*invXformL32);
    
    CCTK_REAL dg4222 = invXformL02*((csetemp137 + csetemp138 + csetemp139 
      + csetemp140)*invXformL02 + (csetemp141 + csetemp142 + csetemp143 + 
      csetemp144)*invXformL12 + (csetemp145 + csetemp146 + csetemp147 + 
      csetemp148)*invXformL22 + (csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL32) + invXformL12*((csetemp141 + csetemp142 + 
      csetemp143 + csetemp144)*invXformL02 + (csetemp153 + csetemp154 + 
      csetemp155 + csetemp156)*invXformL12 + (csetemp157 + csetemp158 + 
      csetemp159 + csetemp160)*invXformL22 + (csetemp161 + csetemp162 + 
      csetemp163 + csetemp164)*invXformL32) + invXformL22*((csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*invXformL02 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*invXformL12 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*invXformL22 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*invXformL32) + 
      invXformL32*((csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL02 + (csetemp161 + csetemp162 + csetemp163 + 
      csetemp164)*invXformL12 + (csetemp169 + csetemp170 + csetemp171 + 
      csetemp172)*invXformL22 + (csetemp173 + csetemp174 + csetemp175 + 
      csetemp176)*invXformL32);
    
    CCTK_REAL dg4223 = invXformL02*((csetemp177 + csetemp178 + csetemp179 
      + csetemp180)*invXformL02 + (csetemp181 + csetemp182 + csetemp183 + 
      csetemp184)*invXformL12 + (csetemp185 + csetemp186 + csetemp187 + 
      csetemp188)*invXformL22 + (csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL32) + invXformL12*((csetemp181 + csetemp182 + 
      csetemp183 + csetemp184)*invXformL02 + (csetemp193 + csetemp194 + 
      csetemp195 + csetemp196)*invXformL12 + (csetemp197 + csetemp198 + 
      csetemp199 + csetemp200)*invXformL22 + (csetemp201 + csetemp202 + 
      csetemp203 + csetemp204)*invXformL32) + invXformL22*((csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*invXformL02 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*invXformL12 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*invXformL22 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*invXformL32) + 
      invXformL32*((csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL02 + (csetemp201 + csetemp202 + csetemp203 + 
      csetemp204)*invXformL12 + (csetemp209 + csetemp210 + csetemp211 + 
      csetemp212)*invXformL22 + (csetemp213 + csetemp214 + csetemp215 + 
      csetemp216)*invXformL32);
    
    CCTK_REAL dg4230 = invXformL02*((csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*invXformL03 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*invXformL13 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*invXformL23 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*invXformL33) + invXformL12*((csetemp61 + csetemp62 + 
      csetemp63 + csetemp64)*invXformL03 + (csetemp73 + csetemp74 + csetemp75 
      + csetemp76)*invXformL13 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*invXformL23 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*invXformL33) + invXformL22*((csetemp65 + csetemp66 + 
      csetemp67 + csetemp68)*invXformL03 + (csetemp77 + csetemp78 + csetemp79 
      + csetemp80)*invXformL13 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*invXformL23 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL33) + invXformL32*((csetemp69 + csetemp70 + 
      csetemp71 + csetemp72)*invXformL03 + (csetemp81 + csetemp82 + csetemp83 
      + csetemp84)*invXformL13 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL23 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*invXformL33);
    
    CCTK_REAL dg4231 = invXformL02*((csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL03 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*invXformL13 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*invXformL23 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL33) + invXformL12*((csetemp101 + csetemp102 + 
      csetemp103 + csetemp104)*invXformL03 + (csetemp113 + csetemp114 + 
      csetemp115 + csetemp116)*invXformL13 + (csetemp117 + csetemp118 + 
      csetemp119 + csetemp120)*invXformL23 + (csetemp121 + csetemp122 + 
      csetemp123 + csetemp124)*invXformL33) + invXformL22*((csetemp105 + 
      csetemp106 + csetemp107 + csetemp108)*invXformL03 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*invXformL13 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*invXformL23 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*invXformL33) + 
      invXformL32*((csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL03 + (csetemp121 + csetemp122 + csetemp123 + 
      csetemp124)*invXformL13 + (csetemp129 + csetemp130 + csetemp131 + 
      csetemp132)*invXformL23 + (csetemp133 + csetemp134 + csetemp135 + 
      csetemp136)*invXformL33);
    
    CCTK_REAL dg4232 = invXformL02*((csetemp137 + csetemp138 + csetemp139 
      + csetemp140)*invXformL03 + (csetemp141 + csetemp142 + csetemp143 + 
      csetemp144)*invXformL13 + (csetemp145 + csetemp146 + csetemp147 + 
      csetemp148)*invXformL23 + (csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL33) + invXformL12*((csetemp141 + csetemp142 + 
      csetemp143 + csetemp144)*invXformL03 + (csetemp153 + csetemp154 + 
      csetemp155 + csetemp156)*invXformL13 + (csetemp157 + csetemp158 + 
      csetemp159 + csetemp160)*invXformL23 + (csetemp161 + csetemp162 + 
      csetemp163 + csetemp164)*invXformL33) + invXformL22*((csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*invXformL03 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*invXformL13 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*invXformL23 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*invXformL33) + 
      invXformL32*((csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL03 + (csetemp161 + csetemp162 + csetemp163 + 
      csetemp164)*invXformL13 + (csetemp169 + csetemp170 + csetemp171 + 
      csetemp172)*invXformL23 + (csetemp173 + csetemp174 + csetemp175 + 
      csetemp176)*invXformL33);
    
    CCTK_REAL dg4233 = invXformL02*((csetemp177 + csetemp178 + csetemp179 
      + csetemp180)*invXformL03 + (csetemp181 + csetemp182 + csetemp183 + 
      csetemp184)*invXformL13 + (csetemp185 + csetemp186 + csetemp187 + 
      csetemp188)*invXformL23 + (csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL33) + invXformL12*((csetemp181 + csetemp182 + 
      csetemp183 + csetemp184)*invXformL03 + (csetemp193 + csetemp194 + 
      csetemp195 + csetemp196)*invXformL13 + (csetemp197 + csetemp198 + 
      csetemp199 + csetemp200)*invXformL23 + (csetemp201 + csetemp202 + 
      csetemp203 + csetemp204)*invXformL33) + invXformL22*((csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*invXformL03 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*invXformL13 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*invXformL23 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*invXformL33) + 
      invXformL32*((csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL03 + (csetemp201 + csetemp202 + csetemp203 + 
      csetemp204)*invXformL13 + (csetemp209 + csetemp210 + csetemp211 + 
      csetemp212)*invXformL23 + (csetemp213 + csetemp214 + csetemp215 + 
      csetemp216)*invXformL33);
    
    CCTK_REAL dg4330 = invXformL03*((csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*invXformL03 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*invXformL13 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*invXformL23 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*invXformL33) + invXformL13*((csetemp61 + csetemp62 + 
      csetemp63 + csetemp64)*invXformL03 + (csetemp73 + csetemp74 + csetemp75 
      + csetemp76)*invXformL13 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*invXformL23 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*invXformL33) + invXformL23*((csetemp65 + csetemp66 + 
      csetemp67 + csetemp68)*invXformL03 + (csetemp77 + csetemp78 + csetemp79 
      + csetemp80)*invXformL13 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*invXformL23 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL33) + invXformL33*((csetemp69 + csetemp70 + 
      csetemp71 + csetemp72)*invXformL03 + (csetemp81 + csetemp82 + csetemp83 
      + csetemp84)*invXformL13 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*invXformL23 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*invXformL33);
    
    CCTK_REAL dg4331 = invXformL03*((csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL03 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*invXformL13 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*invXformL23 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL33) + invXformL13*((csetemp101 + csetemp102 + 
      csetemp103 + csetemp104)*invXformL03 + (csetemp113 + csetemp114 + 
      csetemp115 + csetemp116)*invXformL13 + (csetemp117 + csetemp118 + 
      csetemp119 + csetemp120)*invXformL23 + (csetemp121 + csetemp122 + 
      csetemp123 + csetemp124)*invXformL33) + invXformL23*((csetemp105 + 
      csetemp106 + csetemp107 + csetemp108)*invXformL03 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*invXformL13 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*invXformL23 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*invXformL33) + 
      invXformL33*((csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*invXformL03 + (csetemp121 + csetemp122 + csetemp123 + 
      csetemp124)*invXformL13 + (csetemp129 + csetemp130 + csetemp131 + 
      csetemp132)*invXformL23 + (csetemp133 + csetemp134 + csetemp135 + 
      csetemp136)*invXformL33);
    
    CCTK_REAL dg4332 = invXformL03*((csetemp137 + csetemp138 + csetemp139 
      + csetemp140)*invXformL03 + (csetemp141 + csetemp142 + csetemp143 + 
      csetemp144)*invXformL13 + (csetemp145 + csetemp146 + csetemp147 + 
      csetemp148)*invXformL23 + (csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL33) + invXformL13*((csetemp141 + csetemp142 + 
      csetemp143 + csetemp144)*invXformL03 + (csetemp153 + csetemp154 + 
      csetemp155 + csetemp156)*invXformL13 + (csetemp157 + csetemp158 + 
      csetemp159 + csetemp160)*invXformL23 + (csetemp161 + csetemp162 + 
      csetemp163 + csetemp164)*invXformL33) + invXformL23*((csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*invXformL03 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*invXformL13 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*invXformL23 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*invXformL33) + 
      invXformL33*((csetemp149 + csetemp150 + csetemp151 + 
      csetemp152)*invXformL03 + (csetemp161 + csetemp162 + csetemp163 + 
      csetemp164)*invXformL13 + (csetemp169 + csetemp170 + csetemp171 + 
      csetemp172)*invXformL23 + (csetemp173 + csetemp174 + csetemp175 + 
      csetemp176)*invXformL33);
    
    CCTK_REAL dg4333 = invXformL03*((csetemp177 + csetemp178 + csetemp179 
      + csetemp180)*invXformL03 + (csetemp181 + csetemp182 + csetemp183 + 
      csetemp184)*invXformL13 + (csetemp185 + csetemp186 + csetemp187 + 
      csetemp188)*invXformL23 + (csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL33) + invXformL13*((csetemp181 + csetemp182 + 
      csetemp183 + csetemp184)*invXformL03 + (csetemp193 + csetemp194 + 
      csetemp195 + csetemp196)*invXformL13 + (csetemp197 + csetemp198 + 
      csetemp199 + csetemp200)*invXformL23 + (csetemp201 + csetemp202 + 
      csetemp203 + csetemp204)*invXformL33) + invXformL23*((csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*invXformL03 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*invXformL13 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*invXformL23 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*invXformL33) + 
      invXformL33*((csetemp189 + csetemp190 + csetemp191 + 
      csetemp192)*invXformL03 + (csetemp201 + csetemp202 + csetemp203 + 
      csetemp204)*invXformL13 + (csetemp209 + csetemp210 + csetemp211 + 
      csetemp212)*invXformL23 + (csetemp213 + csetemp214 + csetemp215 + 
      csetemp216)*invXformL33);
    
    CCTK_REAL betal1 = g401;
    
    CCTK_REAL betal2 = g402;
    
    CCTK_REAL betal3 = g403;
    
    gxxL = g411;
    
    gxyL = g412;
    
    gxzL = g413;
    
    gyyL = g422;
    
    gyzL = g423;
    
    gzzL = g433;
    
    CCTK_REAL csetemp217 = SQR(gxzL);
    
    CCTK_REAL csetemp218 = SQR(gyzL);
    
    CCTK_REAL csetemp219 = SQR(gxyL);
    
    CCTK_REAL detg = 2*gxyL*gxzL*gyzL + gyyL*(gxxL*gzzL - 
      csetemp217) - gxxL*csetemp218 - gzzL*csetemp219;
    
    CCTK_REAL csetemp220 = INV(detg);
    
    CCTK_REAL gu11 = (gyyL*gzzL - csetemp218)*csetemp220;
    
    CCTK_REAL gu12 = (gxzL*gyzL - gxyL*gzzL)*csetemp220;
    
    CCTK_REAL gu13 = (-(gxzL*gyyL) + gxyL*gyzL)*csetemp220;
    
    CCTK_REAL gu22 = (gxxL*gzzL - csetemp217)*csetemp220;
    
    CCTK_REAL gu23 = (gxyL*gxzL - gxxL*gyzL)*csetemp220;
    
    CCTK_REAL gu33 = (gxxL*gyyL - csetemp219)*csetemp220;
    
    betaxL = betal1*gu11 + betal2*gu12 + betal3*gu13;
    
    betayL = betal1*gu12 + betal2*gu22 + betal3*gu23;
    
    betazL = betal1*gu13 + betal2*gu23 + betal3*gu33;
    
    CCTK_REAL betasq = betaxL*betal1 + betayL*betal2 + 
      betazL*betal3;
    
    alpL = sqrt(betasq - g400);
    
    CCTK_REAL dtg11 = dg4110;
    
    CCTK_REAL dtg12 = dg4120;
    
    CCTK_REAL dtg13 = dg4130;
    
    CCTK_REAL dtg22 = dg4220;
    
    CCTK_REAL dtg23 = dg4230;
    
    CCTK_REAL dtg33 = dg4330;
    
    CCTK_REAL dg111 = dg4111;
    
    CCTK_REAL dg112 = dg4112;
    
    CCTK_REAL dg113 = dg4113;
    
    CCTK_REAL dg121 = dg4121;
    
    CCTK_REAL dg122 = dg4122;
    
    CCTK_REAL dg123 = dg4123;
    
    CCTK_REAL dg131 = dg4131;
    
    CCTK_REAL dg132 = dg4132;
    
    CCTK_REAL dg133 = dg4133;
    
    CCTK_REAL dg221 = dg4221;
    
    CCTK_REAL dg222 = dg4222;
    
    CCTK_REAL dg223 = dg4223;
    
    CCTK_REAL dg231 = dg4231;
    
    CCTK_REAL dg232 = dg4232;
    
    CCTK_REAL dg233 = dg4233;
    
    CCTK_REAL dg331 = dg4331;
    
    CCTK_REAL dg332 = dg4332;
    
    CCTK_REAL dg333 = dg4333;
    
    CCTK_REAL csetemp221 = SQR(gu11);
    
    CCTK_REAL csetemp222 = SQR(gu12);
    
    CCTK_REAL csetemp223 = SQR(gu13);
    
    CCTK_REAL dtgu11 = -(csetemp221*dtg11) - csetemp222*dtg22 - 
      csetemp223*dtg33 - 2*dtg12*gu11*gu12 - 2*dtg13*gu11*gu13 - 
      2*dtg23*gu12*gu13;
    
    CCTK_REAL dtgu12 = gu12*(-(dtg11*gu11) - dtg13*gu13 - dtg22*gu22) + 
      dtg12*(-csetemp222 - gu11*gu22) + (-(dtg13*gu11) - dtg33*gu13)*gu23 + 
      dtg23*(-(gu13*gu22) - gu12*gu23);
    
    CCTK_REAL dtgu13 = (-(dtg12*gu11) - dtg22*gu12)*gu23 - dtg23*gu12*gu33 
      + gu13*(-(dtg11*gu11) - dtg12*gu12 - dtg23*gu23 - dtg33*gu33) + 
      dtg13*(-csetemp223 - gu11*gu33);
    
    CCTK_REAL csetemp224 = SQR(gu22);
    
    CCTK_REAL csetemp225 = SQR(gu23);
    
    CCTK_REAL dtgu22 = -(csetemp222*dtg11) - csetemp224*dtg22 - 
      csetemp225*dtg33 - 2*dtg12*gu12*gu22 - 2*dtg13*gu12*gu23 - 
      2*dtg23*gu22*gu23;
    
    CCTK_REAL dtgu23 = gu13*(-(dtg11*gu12) - dtg12*gu22 - dtg13*gu23) - 
      dtg13*gu12*gu33 + gu23*(-(dtg12*gu12) - dtg22*gu22 - dtg33*gu33) + 
      dtg23*(-csetemp225 - gu22*gu33);
    
    CCTK_REAL csetemp226 = SQR(gu33);
    
    CCTK_REAL dtgu33 = -(csetemp223*dtg11) - csetemp225*dtg22 - 
      csetemp226*dtg33 - 2*dtg12*gu13*gu23 - 2*dtg13*gu13*gu33 - 
      2*dtg23*gu23*gu33;
    
    CCTK_REAL dgu111 = -(csetemp221*dg111) - csetemp222*dg221 - 
      csetemp223*dg331 - 2*dg121*gu11*gu12 - 2*dg131*gu11*gu13 - 
      2*dg231*gu12*gu13;
    
    CCTK_REAL dgu121 = gu12*(-(dg111*gu11) - dg131*gu13 - dg221*gu22) + 
      dg121*(-csetemp222 - gu11*gu22) + (-(dg131*gu11) - dg331*gu13)*gu23 + 
      dg231*(-(gu13*gu22) - gu12*gu23);
    
    CCTK_REAL dgu131 = (-(dg121*gu11) - dg221*gu12)*gu23 - dg231*gu12*gu33 
      + gu13*(-(dg111*gu11) - dg121*gu12 - dg231*gu23 - dg331*gu33) + 
      dg131*(-csetemp223 - gu11*gu33);
    
    CCTK_REAL dgu221 = -(csetemp222*dg111) - csetemp224*dg221 - 
      csetemp225*dg331 - 2*dg121*gu12*gu22 - 2*dg131*gu12*gu23 - 
      2*dg231*gu22*gu23;
    
    CCTK_REAL dgu231 = gu13*(-(dg111*gu12) - dg121*gu22 - dg131*gu23) - 
      dg131*gu12*gu33 + gu23*(-(dg121*gu12) - dg221*gu22 - dg331*gu33) + 
      dg231*(-csetemp225 - gu22*gu33);
    
    CCTK_REAL dgu331 = -(csetemp223*dg111) - csetemp225*dg221 - 
      csetemp226*dg331 - 2*dg121*gu13*gu23 - 2*dg131*gu13*gu33 - 
      2*dg231*gu23*gu33;
    
    CCTK_REAL dgu112 = -(csetemp221*dg112) - csetemp222*dg222 - 
      csetemp223*dg332 - 2*dg122*gu11*gu12 - 2*dg132*gu11*gu13 - 
      2*dg232*gu12*gu13;
    
    CCTK_REAL dgu122 = gu12*(-(dg112*gu11) - dg132*gu13 - dg222*gu22) + 
      dg122*(-csetemp222 - gu11*gu22) + (-(dg132*gu11) - dg332*gu13)*gu23 + 
      dg232*(-(gu13*gu22) - gu12*gu23);
    
    CCTK_REAL dgu132 = (-(dg122*gu11) - dg222*gu12)*gu23 - dg232*gu12*gu33 
      + gu13*(-(dg112*gu11) - dg122*gu12 - dg232*gu23 - dg332*gu33) + 
      dg132*(-csetemp223 - gu11*gu33);
    
    CCTK_REAL dgu222 = -(csetemp222*dg112) - csetemp224*dg222 - 
      csetemp225*dg332 - 2*dg122*gu12*gu22 - 2*dg132*gu12*gu23 - 
      2*dg232*gu22*gu23;
    
    CCTK_REAL dgu232 = gu13*(-(dg112*gu12) - dg122*gu22 - dg132*gu23) - 
      dg132*gu12*gu33 + gu23*(-(dg122*gu12) - dg222*gu22 - dg332*gu33) + 
      dg232*(-csetemp225 - gu22*gu33);
    
    CCTK_REAL dgu332 = -(csetemp223*dg112) - csetemp225*dg222 - 
      csetemp226*dg332 - 2*dg122*gu13*gu23 - 2*dg132*gu13*gu33 - 
      2*dg232*gu23*gu33;
    
    CCTK_REAL dgu113 = -(csetemp221*dg113) - csetemp222*dg223 - 
      csetemp223*dg333 - 2*dg123*gu11*gu12 - 2*dg133*gu11*gu13 - 
      2*dg233*gu12*gu13;
    
    CCTK_REAL dgu123 = gu12*(-(dg113*gu11) - dg133*gu13 - dg223*gu22) + 
      dg123*(-csetemp222 - gu11*gu22) + (-(dg133*gu11) - dg333*gu13)*gu23 + 
      dg233*(-(gu13*gu22) - gu12*gu23);
    
    CCTK_REAL dgu133 = (-(dg123*gu11) - dg223*gu12)*gu23 - dg233*gu12*gu33 
      + gu13*(-(dg113*gu11) - dg123*gu12 - dg233*gu23 - dg333*gu33) + 
      dg133*(-csetemp223 - gu11*gu33);
    
    CCTK_REAL dgu223 = -(csetemp222*dg113) - csetemp224*dg223 - 
      csetemp225*dg333 - 2*dg123*gu12*gu22 - 2*dg133*gu12*gu23 - 
      2*dg233*gu22*gu23;
    
    CCTK_REAL dgu233 = gu13*(-(dg113*gu12) - dg123*gu22 - dg133*gu23) - 
      dg133*gu12*gu33 + gu23*(-(dg123*gu12) - dg223*gu22 - dg333*gu33) + 
      dg233*(-csetemp225 - gu22*gu33);
    
    CCTK_REAL dgu333 = -(csetemp223*dg113) - csetemp225*dg223 - 
      csetemp226*dg333 - 2*dg123*gu13*gu23 - 2*dg133*gu13*gu33 - 
      2*dg233*gu23*gu33;
    
    CCTK_REAL dtbetal1 = dg4010;
    
    CCTK_REAL dtbetal2 = dg4020;
    
    CCTK_REAL dtbetal3 = dg4030;
    
    CCTK_REAL dbetal11 = dg4011;
    
    CCTK_REAL dbetal12 = dg4012;
    
    CCTK_REAL dbetal13 = dg4013;
    
    CCTK_REAL dbetal21 = dg4021;
    
    CCTK_REAL dbetal22 = dg4022;
    
    CCTK_REAL dbetal23 = dg4023;
    
    CCTK_REAL dbetal31 = dg4031;
    
    CCTK_REAL dbetal32 = dg4032;
    
    CCTK_REAL dbetal33 = dg4033;
    
    dtbetaxL = betal1*dtgu11 + betal2*dtgu12 + betal3*dtgu13 + 
      dtbetal1*gu11 + dtbetal2*gu12 + dtbetal3*gu13;
    
    dtbetayL = betal1*dtgu12 + betal2*dtgu22 + betal3*dtgu23 + 
      dtbetal1*gu12 + dtbetal2*gu22 + dtbetal3*gu23;
    
    dtbetazL = betal1*dtgu13 + betal2*dtgu23 + betal3*dtgu33 + 
      dtbetal1*gu13 + dtbetal2*gu23 + dtbetal3*gu33;
    
    CCTK_REAL dbeta11 = betal1*dgu111 + betal2*dgu121 + betal3*dgu131 + 
      dbetal11*gu11 + dbetal21*gu12 + dbetal31*gu13;
    
    CCTK_REAL dbeta21 = betal1*dgu121 + betal2*dgu221 + betal3*dgu231 + 
      dbetal11*gu12 + dbetal21*gu22 + dbetal31*gu23;
    
    CCTK_REAL dbeta31 = betal1*dgu131 + betal2*dgu231 + betal3*dgu331 + 
      dbetal11*gu13 + dbetal21*gu23 + dbetal31*gu33;
    
    CCTK_REAL dbeta12 = betal1*dgu112 + betal2*dgu122 + betal3*dgu132 + 
      dbetal12*gu11 + dbetal22*gu12 + dbetal32*gu13;
    
    CCTK_REAL dbeta22 = betal1*dgu122 + betal2*dgu222 + betal3*dgu232 + 
      dbetal12*gu12 + dbetal22*gu22 + dbetal32*gu23;
    
    CCTK_REAL dbeta32 = betal1*dgu132 + betal2*dgu232 + betal3*dgu332 + 
      dbetal12*gu13 + dbetal22*gu23 + dbetal32*gu33;
    
    CCTK_REAL dbeta13 = betal1*dgu113 + betal2*dgu123 + betal3*dgu133 + 
      dbetal13*gu11 + dbetal23*gu12 + dbetal33*gu13;
    
    CCTK_REAL dbeta23 = betal1*dgu123 + betal2*dgu223 + betal3*dgu233 + 
      dbetal13*gu12 + dbetal23*gu22 + dbetal33*gu23;
    
    CCTK_REAL dbeta33 = betal1*dgu133 + betal2*dgu233 + betal3*dgu333 + 
      dbetal13*gu13 + dbetal23*gu23 + dbetal33*gu33;
    
    CCTK_REAL dtbetasq = dtbetaxL*betal1 + dtbetayL*betal2 + 
      dtbetazL*betal3 + betaxL*dtbetal1 + betayL*dtbetal2 + 
      betazL*dtbetal3;
    
    CCTK_REAL csetemp227 = INV(alpL);
    
    CCTK_REAL dtalpL = 0.5*csetemp227*(-dg4000 + dtbetasq);
    
    CCTK_REAL kxxL = 0.5*csetemp227*(2*(gxxL*dbeta11 + gxyL*dbeta21 + 
      gxzL*dbeta31) + betaxL*dg111 + betayL*dg112 + betazL*dg113 - 
      dtg11);
    
    CCTK_REAL kxyL = 0.5*csetemp227*(gxxL*dbeta12 + gyyL*dbeta21 + 
      gxyL*(dbeta11 + dbeta22) + gyzL*dbeta31 + gxzL*dbeta32 + 
      betaxL*dg121 + betayL*dg122 + betazL*dg123 - dtg12);
    
    CCTK_REAL kxzL = 0.5*csetemp227*(gxxL*dbeta13 + gyzL*dbeta21 + 
      gxyL*dbeta23 + gzzL*dbeta31 + gxzL*(dbeta11 + dbeta33) + 
      betaxL*dg131 + betayL*dg132 + betazL*dg133 - dtg13);
    
    CCTK_REAL kyyL = 0.5*csetemp227*(2*(gxyL*dbeta12 + gyyL*dbeta22 + 
      gyzL*dbeta32) + betaxL*dg221 + betayL*dg222 + betazL*dg223 - 
      dtg22);
    
    CCTK_REAL kyzL = 0.5*csetemp227*(gxzL*dbeta12 + gxyL*dbeta13 + 
      gyyL*dbeta23 + gzzL*dbeta32 + gyzL*(dbeta22 + dbeta33) + 
      betaxL*dg231 + betayL*dg232 + betazL*dg233 - dtg23);
    
    CCTK_REAL kzzL = 0.5*csetemp227*(2*(gxzL*dbeta13 + gyzL*dbeta23 + 
      gzzL*dbeta33) + betaxL*dg331 + betayL*dg332 + betazL*dg333 - 
      dtg33);
    
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
  CCTK_ENDLOOP3(Minkowski_initial);
}

extern "C" void Minkowski_initial(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering Minkowski_initial_Body");
  }
  
  if (cctk_iteration % Minkowski_initial_calc_every != Minkowski_initial_calc_offset)
  {
    return;
  }
  
  const char *const groups[] = {
    "admbase::curv",
    "admbase::dtlapse",
    "admbase::dtshift",
    "admbase::lapse",
    "admbase::metric",
    "admbase::shift"};
  GenericFD_AssertGroupStorage(cctkGH, "Minkowski_initial", 6, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, Minkowski_initial_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving Minkowski_initial_Body");
  }
}
