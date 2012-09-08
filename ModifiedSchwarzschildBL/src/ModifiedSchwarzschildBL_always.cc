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

static void ModifiedSchwarzschildBL_always_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_LOOP3(ModifiedSchwarzschildBL_always,
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
    CCTK_REAL xL = x[index];
    CCTK_REAL yL = y[index];
    CCTK_REAL zL = z[index];
    
    
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
    
    CCTK_REAL xx0 = t - ToReal(timeoffset);
    
    CCTK_REAL xx1 = xL - ToReal(positionx);
    
    CCTK_REAL xx2 = yL - ToReal(positiony);
    
    CCTK_REAL xx3 = zL - ToReal(positionz);
    
    CCTK_REAL txx1 = invXformL10*xx0 + invXformL11*xx1 + invXformL12*xx2 + 
      invXformL13*xx3;
    
    CCTK_REAL txx2 = invXformL20*xx0 + invXformL21*xx1 + invXformL22*xx2 + 
      invXformL23*xx3;
    
    CCTK_REAL txx3 = invXformL30*xx0 + invXformL31*xx1 + invXformL32*xx2 + 
      invXformL33*xx3;
    
    CCTK_REAL X = txx1;
    
    CCTK_REAL Y = txx2;
    
    CCTK_REAL Z = txx3;
    
    CCTK_REAL rXYZ = sqrt(SQR(X) + SQR(Y) + SQR(Z));
    
    CCTK_REAL tg400 = -1;
    
    CCTK_REAL tg401 = 0;
    
    CCTK_REAL tg402 = 0;
    
    CCTK_REAL tg403 = 0;
    
    CCTK_REAL csetemp9 = 0.5*ToReal(M);
    
    CCTK_REAL csetemp10 = INV(rXYZ);
    
    CCTK_REAL csetemp11 = INV(ToReal(M));
    
    CCTK_REAL csetemp12 = SQR(ToReal(M));
    
    CCTK_REAL csetemp13 = INV(SQR(Pi));
    
    CCTK_REAL tg411 = 0.00390625*QAD(4 + csetemp10*(8*fmin(csetemp9,rXYZ) 
      + csetemp11*(csetemp12*(csetemp13 - 
      csetemp13*cos(4*csetemp11*Pi*fmin(csetemp9,rXYZ))) - 
      8*SQR(fmin(csetemp9,rXYZ)))));
    
    CCTK_REAL tg412 = 0;
    
    CCTK_REAL tg413 = 0;
    
    CCTK_REAL tg422 = 0.00390625*QAD(4 + csetemp10*(8*fmin(csetemp9,rXYZ) 
      + csetemp11*(csetemp12*(csetemp13 - 
      csetemp13*cos(4*csetemp11*Pi*fmin(csetemp9,rXYZ))) - 
      8*SQR(fmin(csetemp9,rXYZ)))));
    
    CCTK_REAL tg423 = 0;
    
    CCTK_REAL tg433 = 0.00390625*QAD(4 + csetemp10*(8*fmin(csetemp9,rXYZ) 
      + csetemp11*(csetemp12*(csetemp13 - 
      csetemp13*cos(4*csetemp11*Pi*fmin(csetemp9,rXYZ))) - 
      8*SQR(fmin(csetemp9,rXYZ)))));
    
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
    
    CCTK_REAL csetemp14 = pow(rXYZ,-6);
    
    CCTK_REAL csetemp15 = 2*rXYZ;
    
    CCTK_REAL csetemp16 = INV(QAD(ToReal(M)));
    
    CCTK_REAL csetemp17 = pow(Pi,-8);
    
    CCTK_REAL csetemp18 = SQR(Pi);
    
    CCTK_REAL csetemp19 = SQR(rXYZ);
    
    CCTK_REAL tdg4111 = IfThen(ToReal(M) <= 
      csetemp15,-0.25*csetemp14*X*CUB(csetemp15 + 
      ToReal(M))*ToReal(M),-0.25*csetemp14*csetemp16*csetemp17*X*CUB(csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) 
      + csetemp18*(-4*csetemp19 + 6*rXYZ*ToReal(M)))*(4*csetemp18*csetemp19 + 
      csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp11*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4112 = IfThen(ToReal(M) <= 
      csetemp15,-0.25*csetemp14*Y*CUB(csetemp15 + 
      ToReal(M))*ToReal(M),-0.25*csetemp14*csetemp16*csetemp17*Y*CUB(csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) 
      + csetemp18*(-4*csetemp19 + 6*rXYZ*ToReal(M)))*(4*csetemp18*csetemp19 + 
      csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp11*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4113 = IfThen(ToReal(M) <= 
      csetemp15,-0.25*csetemp14*Z*CUB(csetemp15 + 
      ToReal(M))*ToReal(M),-0.25*csetemp14*csetemp16*csetemp17*Z*CUB(csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) 
      + csetemp18*(-4*csetemp19 + 6*rXYZ*ToReal(M)))*(4*csetemp18*csetemp19 + 
      csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp11*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4120 = 0;
    
    CCTK_REAL tdg4121 = 0;
    
    CCTK_REAL tdg4122 = 0;
    
    CCTK_REAL tdg4123 = 0;
    
    CCTK_REAL tdg4130 = 0;
    
    CCTK_REAL tdg4131 = 0;
    
    CCTK_REAL tdg4132 = 0;
    
    CCTK_REAL tdg4133 = 0;
    
    CCTK_REAL tdg4220 = 0;
    
    CCTK_REAL tdg4221 = IfThen(ToReal(M) <= 
      csetemp15,-0.25*csetemp14*X*CUB(csetemp15 + 
      ToReal(M))*ToReal(M),-0.25*csetemp14*csetemp16*csetemp17*X*CUB(csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) 
      + csetemp18*(-4*csetemp19 + 6*rXYZ*ToReal(M)))*(4*csetemp18*csetemp19 + 
      csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp11*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4222 = IfThen(ToReal(M) <= 
      csetemp15,-0.25*csetemp14*Y*CUB(csetemp15 + 
      ToReal(M))*ToReal(M),-0.25*csetemp14*csetemp16*csetemp17*Y*CUB(csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) 
      + csetemp18*(-4*csetemp19 + 6*rXYZ*ToReal(M)))*(4*csetemp18*csetemp19 + 
      csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp11*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4223 = IfThen(ToReal(M) <= 
      csetemp15,-0.25*csetemp14*Z*CUB(csetemp15 + 
      ToReal(M))*ToReal(M),-0.25*csetemp14*csetemp16*csetemp17*Z*CUB(csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) 
      + csetemp18*(-4*csetemp19 + 6*rXYZ*ToReal(M)))*(4*csetemp18*csetemp19 + 
      csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp11*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4230 = 0;
    
    CCTK_REAL tdg4231 = 0;
    
    CCTK_REAL tdg4232 = 0;
    
    CCTK_REAL tdg4233 = 0;
    
    CCTK_REAL tdg4330 = 0;
    
    CCTK_REAL tdg4331 = IfThen(ToReal(M) <= 
      csetemp15,-0.25*csetemp14*X*CUB(csetemp15 + 
      ToReal(M))*ToReal(M),-0.25*csetemp14*csetemp16*csetemp17*X*CUB(csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) 
      + csetemp18*(-4*csetemp19 + 6*rXYZ*ToReal(M)))*(4*csetemp18*csetemp19 + 
      csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp11*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4332 = IfThen(ToReal(M) <= 
      csetemp15,-0.25*csetemp14*Y*CUB(csetemp15 + 
      ToReal(M))*ToReal(M),-0.25*csetemp14*csetemp16*csetemp17*Y*CUB(csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) 
      + csetemp18*(-4*csetemp19 + 6*rXYZ*ToReal(M)))*(4*csetemp18*csetemp19 + 
      csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp11*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4333 = IfThen(ToReal(M) <= 
      csetemp15,-0.25*csetemp14*Z*CUB(csetemp15 + 
      ToReal(M))*ToReal(M),-0.25*csetemp14*csetemp16*csetemp17*Z*CUB(csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) 
      + csetemp18*(-4*csetemp19 + 6*rXYZ*ToReal(M)))*(4*csetemp18*csetemp19 + 
      csetemp12*SQR(sin(2*csetemp11*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp11*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL g400 = 2*(invXformL00*(invXformL10*tg401 + invXformL20*tg402 
      + invXformL30*tg403) + invXformL10*(invXformL20*tg412 + 
      invXformL30*tg413) + invXformL20*invXformL30*tg423) + 
      tg400*SQR(invXformL00) + tg411*SQR(invXformL10) + 
      tg422*SQR(invXformL20) + tg433*SQR(invXformL30);
    
    CCTK_REAL csetemp20 = invXformL01*tg400;
    
    CCTK_REAL csetemp21 = invXformL11*tg401;
    
    CCTK_REAL csetemp22 = invXformL21*tg402;
    
    CCTK_REAL csetemp23 = invXformL31*tg403;
    
    CCTK_REAL csetemp24 = invXformL01*tg401;
    
    CCTK_REAL csetemp25 = invXformL11*tg411;
    
    CCTK_REAL csetemp26 = invXformL21*tg412;
    
    CCTK_REAL csetemp27 = invXformL31*tg413;
    
    CCTK_REAL csetemp28 = invXformL01*tg402;
    
    CCTK_REAL csetemp29 = invXformL11*tg412;
    
    CCTK_REAL csetemp30 = invXformL21*tg422;
    
    CCTK_REAL csetemp31 = invXformL31*tg423;
    
    CCTK_REAL csetemp32 = invXformL01*tg403;
    
    CCTK_REAL csetemp33 = invXformL11*tg413;
    
    CCTK_REAL csetemp34 = invXformL21*tg423;
    
    CCTK_REAL csetemp35 = invXformL31*tg433;
    
    CCTK_REAL g401 = (csetemp20 + csetemp21 + csetemp22 + 
      csetemp23)*invXformL00 + (csetemp24 + csetemp25 + csetemp26 + 
      csetemp27)*invXformL10 + (csetemp28 + csetemp29 + csetemp30 + 
      csetemp31)*invXformL20 + (csetemp32 + csetemp33 + csetemp34 + 
      csetemp35)*invXformL30;
    
    CCTK_REAL csetemp36 = invXformL02*tg400;
    
    CCTK_REAL csetemp37 = invXformL12*tg401;
    
    CCTK_REAL csetemp38 = invXformL22*tg402;
    
    CCTK_REAL csetemp39 = invXformL32*tg403;
    
    CCTK_REAL csetemp40 = invXformL02*tg401;
    
    CCTK_REAL csetemp41 = invXformL12*tg411;
    
    CCTK_REAL csetemp42 = invXformL22*tg412;
    
    CCTK_REAL csetemp43 = invXformL32*tg413;
    
    CCTK_REAL csetemp44 = invXformL02*tg402;
    
    CCTK_REAL csetemp45 = invXformL12*tg412;
    
    CCTK_REAL csetemp46 = invXformL22*tg422;
    
    CCTK_REAL csetemp47 = invXformL32*tg423;
    
    CCTK_REAL csetemp48 = invXformL02*tg403;
    
    CCTK_REAL csetemp49 = invXformL12*tg413;
    
    CCTK_REAL csetemp50 = invXformL22*tg423;
    
    CCTK_REAL csetemp51 = invXformL32*tg433;
    
    CCTK_REAL g402 = (csetemp36 + csetemp37 + csetemp38 + 
      csetemp39)*invXformL00 + (csetemp40 + csetemp41 + csetemp42 + 
      csetemp43)*invXformL10 + (csetemp44 + csetemp45 + csetemp46 + 
      csetemp47)*invXformL20 + (csetemp48 + csetemp49 + csetemp50 + 
      csetemp51)*invXformL30;
    
    CCTK_REAL csetemp52 = invXformL03*tg400;
    
    CCTK_REAL csetemp53 = invXformL13*tg401;
    
    CCTK_REAL csetemp54 = invXformL23*tg402;
    
    CCTK_REAL csetemp55 = invXformL33*tg403;
    
    CCTK_REAL csetemp56 = invXformL03*tg401;
    
    CCTK_REAL csetemp57 = invXformL13*tg411;
    
    CCTK_REAL csetemp58 = invXformL23*tg412;
    
    CCTK_REAL csetemp59 = invXformL33*tg413;
    
    CCTK_REAL csetemp60 = invXformL03*tg402;
    
    CCTK_REAL csetemp61 = invXformL13*tg412;
    
    CCTK_REAL csetemp62 = invXformL23*tg422;
    
    CCTK_REAL csetemp63 = invXformL33*tg423;
    
    CCTK_REAL csetemp64 = invXformL03*tg403;
    
    CCTK_REAL csetemp65 = invXformL13*tg413;
    
    CCTK_REAL csetemp66 = invXformL23*tg423;
    
    CCTK_REAL csetemp67 = invXformL33*tg433;
    
    CCTK_REAL g403 = (csetemp52 + csetemp53 + csetemp54 + 
      csetemp55)*invXformL00 + (csetemp56 + csetemp57 + csetemp58 + 
      csetemp59)*invXformL10 + (csetemp60 + csetemp61 + csetemp62 + 
      csetemp63)*invXformL20 + (csetemp64 + csetemp65 + csetemp66 + 
      csetemp67)*invXformL30;
    
    CCTK_REAL g411 = (csetemp20 + csetemp21 + csetemp22 + 
      csetemp23)*invXformL01 + (csetemp24 + csetemp25 + csetemp26 + 
      csetemp27)*invXformL11 + (csetemp28 + csetemp29 + csetemp30 + 
      csetemp31)*invXformL21 + (csetemp32 + csetemp33 + csetemp34 + 
      csetemp35)*invXformL31;
    
    CCTK_REAL g412 = (csetemp36 + csetemp37 + csetemp38 + 
      csetemp39)*invXformL01 + (csetemp40 + csetemp41 + csetemp42 + 
      csetemp43)*invXformL11 + (csetemp44 + csetemp45 + csetemp46 + 
      csetemp47)*invXformL21 + (csetemp48 + csetemp49 + csetemp50 + 
      csetemp51)*invXformL31;
    
    CCTK_REAL g413 = (csetemp52 + csetemp53 + csetemp54 + 
      csetemp55)*invXformL01 + (csetemp56 + csetemp57 + csetemp58 + 
      csetemp59)*invXformL11 + (csetemp60 + csetemp61 + csetemp62 + 
      csetemp63)*invXformL21 + (csetemp64 + csetemp65 + csetemp66 + 
      csetemp67)*invXformL31;
    
    CCTK_REAL g422 = (csetemp36 + csetemp37 + csetemp38 + 
      csetemp39)*invXformL02 + (csetemp40 + csetemp41 + csetemp42 + 
      csetemp43)*invXformL12 + (csetemp44 + csetemp45 + csetemp46 + 
      csetemp47)*invXformL22 + (csetemp48 + csetemp49 + csetemp50 + 
      csetemp51)*invXformL32;
    
    CCTK_REAL g423 = (csetemp52 + csetemp53 + csetemp54 + 
      csetemp55)*invXformL02 + (csetemp56 + csetemp57 + csetemp58 + 
      csetemp59)*invXformL12 + (csetemp60 + csetemp61 + csetemp62 + 
      csetemp63)*invXformL22 + (csetemp64 + csetemp65 + csetemp66 + 
      csetemp67)*invXformL32;
    
    CCTK_REAL g433 = (csetemp52 + csetemp53 + csetemp54 + 
      csetemp55)*invXformL03 + (csetemp56 + csetemp57 + csetemp58 + 
      csetemp59)*invXformL13 + (csetemp60 + csetemp61 + csetemp62 + 
      csetemp63)*invXformL23 + (csetemp64 + csetemp65 + csetemp66 + 
      csetemp67)*invXformL33;
    
    CCTK_REAL csetemp68 = invXformL00*tdg4000;
    
    CCTK_REAL csetemp69 = invXformL10*tdg4001;
    
    CCTK_REAL csetemp70 = invXformL20*tdg4002;
    
    CCTK_REAL csetemp71 = invXformL30*tdg4003;
    
    CCTK_REAL csetemp72 = invXformL00*tdg4010;
    
    CCTK_REAL csetemp73 = invXformL10*tdg4011;
    
    CCTK_REAL csetemp74 = invXformL20*tdg4012;
    
    CCTK_REAL csetemp75 = invXformL30*tdg4013;
    
    CCTK_REAL csetemp76 = invXformL00*tdg4020;
    
    CCTK_REAL csetemp77 = invXformL10*tdg4021;
    
    CCTK_REAL csetemp78 = invXformL20*tdg4022;
    
    CCTK_REAL csetemp79 = invXformL30*tdg4023;
    
    CCTK_REAL csetemp80 = invXformL00*tdg4030;
    
    CCTK_REAL csetemp81 = invXformL10*tdg4031;
    
    CCTK_REAL csetemp82 = invXformL20*tdg4032;
    
    CCTK_REAL csetemp83 = invXformL30*tdg4033;
    
    CCTK_REAL csetemp84 = invXformL00*tdg4110;
    
    CCTK_REAL csetemp85 = invXformL10*tdg4111;
    
    CCTK_REAL csetemp86 = invXformL20*tdg4112;
    
    CCTK_REAL csetemp87 = invXformL30*tdg4113;
    
    CCTK_REAL csetemp88 = invXformL00*tdg4120;
    
    CCTK_REAL csetemp89 = invXformL10*tdg4121;
    
    CCTK_REAL csetemp90 = invXformL20*tdg4122;
    
    CCTK_REAL csetemp91 = invXformL30*tdg4123;
    
    CCTK_REAL csetemp92 = invXformL00*tdg4130;
    
    CCTK_REAL csetemp93 = invXformL10*tdg4131;
    
    CCTK_REAL csetemp94 = invXformL20*tdg4132;
    
    CCTK_REAL csetemp95 = invXformL30*tdg4133;
    
    CCTK_REAL csetemp96 = invXformL00*tdg4220;
    
    CCTK_REAL csetemp97 = invXformL10*tdg4221;
    
    CCTK_REAL csetemp98 = invXformL20*tdg4222;
    
    CCTK_REAL csetemp99 = invXformL30*tdg4223;
    
    CCTK_REAL csetemp100 = invXformL00*tdg4230;
    
    CCTK_REAL csetemp101 = invXformL10*tdg4231;
    
    CCTK_REAL csetemp102 = invXformL20*tdg4232;
    
    CCTK_REAL csetemp103 = invXformL30*tdg4233;
    
    CCTK_REAL csetemp104 = invXformL00*tdg4330;
    
    CCTK_REAL csetemp105 = invXformL10*tdg4331;
    
    CCTK_REAL csetemp106 = invXformL20*tdg4332;
    
    CCTK_REAL csetemp107 = invXformL30*tdg4333;
    
    CCTK_REAL dg4000 = invXformL20*((csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL00 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL10 + (csetemp96 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL20 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL30) + invXformL30*((csetemp80 + csetemp81 + 
      csetemp82 + csetemp83)*invXformL00 + (csetemp92 + csetemp93 + csetemp94 
      + csetemp95)*invXformL10 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL20 + (csetemp104 + csetemp105 + csetemp106 + 
      csetemp107)*invXformL30) + invXformL00*((csetemp68 + csetemp69 + 
      csetemp70 + csetemp71)*invXformL00 + (csetemp72 + csetemp73 + csetemp74 
      + csetemp75)*invXformL10 + (csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL20 + (csetemp80 + csetemp81 + csetemp82 + 
      csetemp83)*invXformL30) + invXformL10*((csetemp72 + csetemp73 + 
      csetemp74 + csetemp75)*invXformL00 + (csetemp84 + csetemp85 + csetemp86 
      + csetemp87)*invXformL10 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL20 + (csetemp92 + csetemp93 + csetemp94 + 
      csetemp95)*invXformL30);
    
    CCTK_REAL dg4010 = invXformL20*((csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL01 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL11 + (csetemp96 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL21 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL31) + invXformL30*((csetemp80 + csetemp81 + 
      csetemp82 + csetemp83)*invXformL01 + (csetemp92 + csetemp93 + csetemp94 
      + csetemp95)*invXformL11 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL21 + (csetemp104 + csetemp105 + csetemp106 + 
      csetemp107)*invXformL31) + invXformL00*((csetemp68 + csetemp69 + 
      csetemp70 + csetemp71)*invXformL01 + (csetemp72 + csetemp73 + csetemp74 
      + csetemp75)*invXformL11 + (csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL21 + (csetemp80 + csetemp81 + csetemp82 + 
      csetemp83)*invXformL31) + invXformL10*((csetemp72 + csetemp73 + 
      csetemp74 + csetemp75)*invXformL01 + (csetemp84 + csetemp85 + csetemp86 
      + csetemp87)*invXformL11 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL21 + (csetemp92 + csetemp93 + csetemp94 + 
      csetemp95)*invXformL31);
    
    CCTK_REAL csetemp108 = invXformL01*tdg4000;
    
    CCTK_REAL csetemp109 = invXformL11*tdg4001;
    
    CCTK_REAL csetemp110 = invXformL21*tdg4002;
    
    CCTK_REAL csetemp111 = invXformL31*tdg4003;
    
    CCTK_REAL csetemp112 = invXformL01*tdg4010;
    
    CCTK_REAL csetemp113 = invXformL11*tdg4011;
    
    CCTK_REAL csetemp114 = invXformL21*tdg4012;
    
    CCTK_REAL csetemp115 = invXformL31*tdg4013;
    
    CCTK_REAL csetemp116 = invXformL01*tdg4020;
    
    CCTK_REAL csetemp117 = invXformL11*tdg4021;
    
    CCTK_REAL csetemp118 = invXformL21*tdg4022;
    
    CCTK_REAL csetemp119 = invXformL31*tdg4023;
    
    CCTK_REAL csetemp120 = invXformL01*tdg4030;
    
    CCTK_REAL csetemp121 = invXformL11*tdg4031;
    
    CCTK_REAL csetemp122 = invXformL21*tdg4032;
    
    CCTK_REAL csetemp123 = invXformL31*tdg4033;
    
    CCTK_REAL csetemp124 = invXformL01*tdg4110;
    
    CCTK_REAL csetemp125 = invXformL11*tdg4111;
    
    CCTK_REAL csetemp126 = invXformL21*tdg4112;
    
    CCTK_REAL csetemp127 = invXformL31*tdg4113;
    
    CCTK_REAL csetemp128 = invXformL01*tdg4120;
    
    CCTK_REAL csetemp129 = invXformL11*tdg4121;
    
    CCTK_REAL csetemp130 = invXformL21*tdg4122;
    
    CCTK_REAL csetemp131 = invXformL31*tdg4123;
    
    CCTK_REAL csetemp132 = invXformL01*tdg4130;
    
    CCTK_REAL csetemp133 = invXformL11*tdg4131;
    
    CCTK_REAL csetemp134 = invXformL21*tdg4132;
    
    CCTK_REAL csetemp135 = invXformL31*tdg4133;
    
    CCTK_REAL csetemp136 = invXformL01*tdg4220;
    
    CCTK_REAL csetemp137 = invXformL11*tdg4221;
    
    CCTK_REAL csetemp138 = invXformL21*tdg4222;
    
    CCTK_REAL csetemp139 = invXformL31*tdg4223;
    
    CCTK_REAL csetemp140 = invXformL01*tdg4230;
    
    CCTK_REAL csetemp141 = invXformL11*tdg4231;
    
    CCTK_REAL csetemp142 = invXformL21*tdg4232;
    
    CCTK_REAL csetemp143 = invXformL31*tdg4233;
    
    CCTK_REAL csetemp144 = invXformL01*tdg4330;
    
    CCTK_REAL csetemp145 = invXformL11*tdg4331;
    
    CCTK_REAL csetemp146 = invXformL21*tdg4332;
    
    CCTK_REAL csetemp147 = invXformL31*tdg4333;
    
    CCTK_REAL dg4011 = invXformL00*((csetemp108 + csetemp109 + csetemp110 
      + csetemp111)*invXformL01 + (csetemp112 + csetemp113 + csetemp114 + 
      csetemp115)*invXformL11 + (csetemp116 + csetemp117 + csetemp118 + 
      csetemp119)*invXformL21 + (csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL31) + invXformL10*((csetemp112 + csetemp113 + 
      csetemp114 + csetemp115)*invXformL01 + (csetemp124 + csetemp125 + 
      csetemp126 + csetemp127)*invXformL11 + (csetemp128 + csetemp129 + 
      csetemp130 + csetemp131)*invXformL21 + (csetemp132 + csetemp133 + 
      csetemp134 + csetemp135)*invXformL31) + invXformL20*((csetemp116 + 
      csetemp117 + csetemp118 + csetemp119)*invXformL01 + (csetemp128 + 
      csetemp129 + csetemp130 + csetemp131)*invXformL11 + (csetemp136 + 
      csetemp137 + csetemp138 + csetemp139)*invXformL21 + (csetemp140 + 
      csetemp141 + csetemp142 + csetemp143)*invXformL31) + 
      invXformL30*((csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL01 + (csetemp132 + csetemp133 + csetemp134 + 
      csetemp135)*invXformL11 + (csetemp140 + csetemp141 + csetemp142 + 
      csetemp143)*invXformL21 + (csetemp144 + csetemp145 + csetemp146 + 
      csetemp147)*invXformL31);
    
    CCTK_REAL csetemp148 = invXformL02*tdg4000;
    
    CCTK_REAL csetemp149 = invXformL12*tdg4001;
    
    CCTK_REAL csetemp150 = invXformL22*tdg4002;
    
    CCTK_REAL csetemp151 = invXformL32*tdg4003;
    
    CCTK_REAL csetemp152 = invXformL02*tdg4010;
    
    CCTK_REAL csetemp153 = invXformL12*tdg4011;
    
    CCTK_REAL csetemp154 = invXformL22*tdg4012;
    
    CCTK_REAL csetemp155 = invXformL32*tdg4013;
    
    CCTK_REAL csetemp156 = invXformL02*tdg4020;
    
    CCTK_REAL csetemp157 = invXformL12*tdg4021;
    
    CCTK_REAL csetemp158 = invXformL22*tdg4022;
    
    CCTK_REAL csetemp159 = invXformL32*tdg4023;
    
    CCTK_REAL csetemp160 = invXformL02*tdg4030;
    
    CCTK_REAL csetemp161 = invXformL12*tdg4031;
    
    CCTK_REAL csetemp162 = invXformL22*tdg4032;
    
    CCTK_REAL csetemp163 = invXformL32*tdg4033;
    
    CCTK_REAL csetemp164 = invXformL02*tdg4110;
    
    CCTK_REAL csetemp165 = invXformL12*tdg4111;
    
    CCTK_REAL csetemp166 = invXformL22*tdg4112;
    
    CCTK_REAL csetemp167 = invXformL32*tdg4113;
    
    CCTK_REAL csetemp168 = invXformL02*tdg4120;
    
    CCTK_REAL csetemp169 = invXformL12*tdg4121;
    
    CCTK_REAL csetemp170 = invXformL22*tdg4122;
    
    CCTK_REAL csetemp171 = invXformL32*tdg4123;
    
    CCTK_REAL csetemp172 = invXformL02*tdg4130;
    
    CCTK_REAL csetemp173 = invXformL12*tdg4131;
    
    CCTK_REAL csetemp174 = invXformL22*tdg4132;
    
    CCTK_REAL csetemp175 = invXformL32*tdg4133;
    
    CCTK_REAL csetemp176 = invXformL02*tdg4220;
    
    CCTK_REAL csetemp177 = invXformL12*tdg4221;
    
    CCTK_REAL csetemp178 = invXformL22*tdg4222;
    
    CCTK_REAL csetemp179 = invXformL32*tdg4223;
    
    CCTK_REAL csetemp180 = invXformL02*tdg4230;
    
    CCTK_REAL csetemp181 = invXformL12*tdg4231;
    
    CCTK_REAL csetemp182 = invXformL22*tdg4232;
    
    CCTK_REAL csetemp183 = invXformL32*tdg4233;
    
    CCTK_REAL csetemp184 = invXformL02*tdg4330;
    
    CCTK_REAL csetemp185 = invXformL12*tdg4331;
    
    CCTK_REAL csetemp186 = invXformL22*tdg4332;
    
    CCTK_REAL csetemp187 = invXformL32*tdg4333;
    
    CCTK_REAL dg4012 = invXformL00*((csetemp148 + csetemp149 + csetemp150 
      + csetemp151)*invXformL01 + (csetemp152 + csetemp153 + csetemp154 + 
      csetemp155)*invXformL11 + (csetemp156 + csetemp157 + csetemp158 + 
      csetemp159)*invXformL21 + (csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL31) + invXformL10*((csetemp152 + csetemp153 + 
      csetemp154 + csetemp155)*invXformL01 + (csetemp164 + csetemp165 + 
      csetemp166 + csetemp167)*invXformL11 + (csetemp168 + csetemp169 + 
      csetemp170 + csetemp171)*invXformL21 + (csetemp172 + csetemp173 + 
      csetemp174 + csetemp175)*invXformL31) + invXformL20*((csetemp156 + 
      csetemp157 + csetemp158 + csetemp159)*invXformL01 + (csetemp168 + 
      csetemp169 + csetemp170 + csetemp171)*invXformL11 + (csetemp176 + 
      csetemp177 + csetemp178 + csetemp179)*invXformL21 + (csetemp180 + 
      csetemp181 + csetemp182 + csetemp183)*invXformL31) + 
      invXformL30*((csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL01 + (csetemp172 + csetemp173 + csetemp174 + 
      csetemp175)*invXformL11 + (csetemp180 + csetemp181 + csetemp182 + 
      csetemp183)*invXformL21 + (csetemp184 + csetemp185 + csetemp186 + 
      csetemp187)*invXformL31);
    
    CCTK_REAL csetemp188 = invXformL03*tdg4000;
    
    CCTK_REAL csetemp189 = invXformL13*tdg4001;
    
    CCTK_REAL csetemp190 = invXformL23*tdg4002;
    
    CCTK_REAL csetemp191 = invXformL33*tdg4003;
    
    CCTK_REAL csetemp192 = invXformL03*tdg4010;
    
    CCTK_REAL csetemp193 = invXformL13*tdg4011;
    
    CCTK_REAL csetemp194 = invXformL23*tdg4012;
    
    CCTK_REAL csetemp195 = invXformL33*tdg4013;
    
    CCTK_REAL csetemp196 = invXformL03*tdg4020;
    
    CCTK_REAL csetemp197 = invXformL13*tdg4021;
    
    CCTK_REAL csetemp198 = invXformL23*tdg4022;
    
    CCTK_REAL csetemp199 = invXformL33*tdg4023;
    
    CCTK_REAL csetemp200 = invXformL03*tdg4030;
    
    CCTK_REAL csetemp201 = invXformL13*tdg4031;
    
    CCTK_REAL csetemp202 = invXformL23*tdg4032;
    
    CCTK_REAL csetemp203 = invXformL33*tdg4033;
    
    CCTK_REAL csetemp204 = invXformL03*tdg4110;
    
    CCTK_REAL csetemp205 = invXformL13*tdg4111;
    
    CCTK_REAL csetemp206 = invXformL23*tdg4112;
    
    CCTK_REAL csetemp207 = invXformL33*tdg4113;
    
    CCTK_REAL csetemp208 = invXformL03*tdg4120;
    
    CCTK_REAL csetemp209 = invXformL13*tdg4121;
    
    CCTK_REAL csetemp210 = invXformL23*tdg4122;
    
    CCTK_REAL csetemp211 = invXformL33*tdg4123;
    
    CCTK_REAL csetemp212 = invXformL03*tdg4130;
    
    CCTK_REAL csetemp213 = invXformL13*tdg4131;
    
    CCTK_REAL csetemp214 = invXformL23*tdg4132;
    
    CCTK_REAL csetemp215 = invXformL33*tdg4133;
    
    CCTK_REAL csetemp216 = invXformL03*tdg4220;
    
    CCTK_REAL csetemp217 = invXformL13*tdg4221;
    
    CCTK_REAL csetemp218 = invXformL23*tdg4222;
    
    CCTK_REAL csetemp219 = invXformL33*tdg4223;
    
    CCTK_REAL csetemp220 = invXformL03*tdg4230;
    
    CCTK_REAL csetemp221 = invXformL13*tdg4231;
    
    CCTK_REAL csetemp222 = invXformL23*tdg4232;
    
    CCTK_REAL csetemp223 = invXformL33*tdg4233;
    
    CCTK_REAL csetemp224 = invXformL03*tdg4330;
    
    CCTK_REAL csetemp225 = invXformL13*tdg4331;
    
    CCTK_REAL csetemp226 = invXformL23*tdg4332;
    
    CCTK_REAL csetemp227 = invXformL33*tdg4333;
    
    CCTK_REAL dg4013 = invXformL00*((csetemp188 + csetemp189 + csetemp190 
      + csetemp191)*invXformL01 + (csetemp192 + csetemp193 + csetemp194 + 
      csetemp195)*invXformL11 + (csetemp196 + csetemp197 + csetemp198 + 
      csetemp199)*invXformL21 + (csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL31) + invXformL10*((csetemp192 + csetemp193 + 
      csetemp194 + csetemp195)*invXformL01 + (csetemp204 + csetemp205 + 
      csetemp206 + csetemp207)*invXformL11 + (csetemp208 + csetemp209 + 
      csetemp210 + csetemp211)*invXformL21 + (csetemp212 + csetemp213 + 
      csetemp214 + csetemp215)*invXformL31) + invXformL20*((csetemp196 + 
      csetemp197 + csetemp198 + csetemp199)*invXformL01 + (csetemp208 + 
      csetemp209 + csetemp210 + csetemp211)*invXformL11 + (csetemp216 + 
      csetemp217 + csetemp218 + csetemp219)*invXformL21 + (csetemp220 + 
      csetemp221 + csetemp222 + csetemp223)*invXformL31) + 
      invXformL30*((csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL01 + (csetemp212 + csetemp213 + csetemp214 + 
      csetemp215)*invXformL11 + (csetemp220 + csetemp221 + csetemp222 + 
      csetemp223)*invXformL21 + (csetemp224 + csetemp225 + csetemp226 + 
      csetemp227)*invXformL31);
    
    CCTK_REAL dg4020 = invXformL20*((csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL02 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL12 + (csetemp96 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL22 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL32) + invXformL30*((csetemp80 + csetemp81 + 
      csetemp82 + csetemp83)*invXformL02 + (csetemp92 + csetemp93 + csetemp94 
      + csetemp95)*invXformL12 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL22 + (csetemp104 + csetemp105 + csetemp106 + 
      csetemp107)*invXformL32) + invXformL00*((csetemp68 + csetemp69 + 
      csetemp70 + csetemp71)*invXformL02 + (csetemp72 + csetemp73 + csetemp74 
      + csetemp75)*invXformL12 + (csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL22 + (csetemp80 + csetemp81 + csetemp82 + 
      csetemp83)*invXformL32) + invXformL10*((csetemp72 + csetemp73 + 
      csetemp74 + csetemp75)*invXformL02 + (csetemp84 + csetemp85 + csetemp86 
      + csetemp87)*invXformL12 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL22 + (csetemp92 + csetemp93 + csetemp94 + 
      csetemp95)*invXformL32);
    
    CCTK_REAL dg4021 = invXformL00*((csetemp108 + csetemp109 + csetemp110 
      + csetemp111)*invXformL02 + (csetemp112 + csetemp113 + csetemp114 + 
      csetemp115)*invXformL12 + (csetemp116 + csetemp117 + csetemp118 + 
      csetemp119)*invXformL22 + (csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL32) + invXformL10*((csetemp112 + csetemp113 + 
      csetemp114 + csetemp115)*invXformL02 + (csetemp124 + csetemp125 + 
      csetemp126 + csetemp127)*invXformL12 + (csetemp128 + csetemp129 + 
      csetemp130 + csetemp131)*invXformL22 + (csetemp132 + csetemp133 + 
      csetemp134 + csetemp135)*invXformL32) + invXformL20*((csetemp116 + 
      csetemp117 + csetemp118 + csetemp119)*invXformL02 + (csetemp128 + 
      csetemp129 + csetemp130 + csetemp131)*invXformL12 + (csetemp136 + 
      csetemp137 + csetemp138 + csetemp139)*invXformL22 + (csetemp140 + 
      csetemp141 + csetemp142 + csetemp143)*invXformL32) + 
      invXformL30*((csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL02 + (csetemp132 + csetemp133 + csetemp134 + 
      csetemp135)*invXformL12 + (csetemp140 + csetemp141 + csetemp142 + 
      csetemp143)*invXformL22 + (csetemp144 + csetemp145 + csetemp146 + 
      csetemp147)*invXformL32);
    
    CCTK_REAL dg4022 = invXformL00*((csetemp148 + csetemp149 + csetemp150 
      + csetemp151)*invXformL02 + (csetemp152 + csetemp153 + csetemp154 + 
      csetemp155)*invXformL12 + (csetemp156 + csetemp157 + csetemp158 + 
      csetemp159)*invXformL22 + (csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL32) + invXformL10*((csetemp152 + csetemp153 + 
      csetemp154 + csetemp155)*invXformL02 + (csetemp164 + csetemp165 + 
      csetemp166 + csetemp167)*invXformL12 + (csetemp168 + csetemp169 + 
      csetemp170 + csetemp171)*invXformL22 + (csetemp172 + csetemp173 + 
      csetemp174 + csetemp175)*invXformL32) + invXformL20*((csetemp156 + 
      csetemp157 + csetemp158 + csetemp159)*invXformL02 + (csetemp168 + 
      csetemp169 + csetemp170 + csetemp171)*invXformL12 + (csetemp176 + 
      csetemp177 + csetemp178 + csetemp179)*invXformL22 + (csetemp180 + 
      csetemp181 + csetemp182 + csetemp183)*invXformL32) + 
      invXformL30*((csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL02 + (csetemp172 + csetemp173 + csetemp174 + 
      csetemp175)*invXformL12 + (csetemp180 + csetemp181 + csetemp182 + 
      csetemp183)*invXformL22 + (csetemp184 + csetemp185 + csetemp186 + 
      csetemp187)*invXformL32);
    
    CCTK_REAL dg4023 = invXformL00*((csetemp188 + csetemp189 + csetemp190 
      + csetemp191)*invXformL02 + (csetemp192 + csetemp193 + csetemp194 + 
      csetemp195)*invXformL12 + (csetemp196 + csetemp197 + csetemp198 + 
      csetemp199)*invXformL22 + (csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL32) + invXformL10*((csetemp192 + csetemp193 + 
      csetemp194 + csetemp195)*invXformL02 + (csetemp204 + csetemp205 + 
      csetemp206 + csetemp207)*invXformL12 + (csetemp208 + csetemp209 + 
      csetemp210 + csetemp211)*invXformL22 + (csetemp212 + csetemp213 + 
      csetemp214 + csetemp215)*invXformL32) + invXformL20*((csetemp196 + 
      csetemp197 + csetemp198 + csetemp199)*invXformL02 + (csetemp208 + 
      csetemp209 + csetemp210 + csetemp211)*invXformL12 + (csetemp216 + 
      csetemp217 + csetemp218 + csetemp219)*invXformL22 + (csetemp220 + 
      csetemp221 + csetemp222 + csetemp223)*invXformL32) + 
      invXformL30*((csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL02 + (csetemp212 + csetemp213 + csetemp214 + 
      csetemp215)*invXformL12 + (csetemp220 + csetemp221 + csetemp222 + 
      csetemp223)*invXformL22 + (csetemp224 + csetemp225 + csetemp226 + 
      csetemp227)*invXformL32);
    
    CCTK_REAL dg4030 = invXformL20*((csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL03 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL13 + (csetemp96 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL23 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL33) + invXformL30*((csetemp80 + csetemp81 + 
      csetemp82 + csetemp83)*invXformL03 + (csetemp92 + csetemp93 + csetemp94 
      + csetemp95)*invXformL13 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL23 + (csetemp104 + csetemp105 + csetemp106 + 
      csetemp107)*invXformL33) + invXformL00*((csetemp68 + csetemp69 + 
      csetemp70 + csetemp71)*invXformL03 + (csetemp72 + csetemp73 + csetemp74 
      + csetemp75)*invXformL13 + (csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL23 + (csetemp80 + csetemp81 + csetemp82 + 
      csetemp83)*invXformL33) + invXformL10*((csetemp72 + csetemp73 + 
      csetemp74 + csetemp75)*invXformL03 + (csetemp84 + csetemp85 + csetemp86 
      + csetemp87)*invXformL13 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL23 + (csetemp92 + csetemp93 + csetemp94 + 
      csetemp95)*invXformL33);
    
    CCTK_REAL dg4031 = invXformL00*((csetemp108 + csetemp109 + csetemp110 
      + csetemp111)*invXformL03 + (csetemp112 + csetemp113 + csetemp114 + 
      csetemp115)*invXformL13 + (csetemp116 + csetemp117 + csetemp118 + 
      csetemp119)*invXformL23 + (csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL33) + invXformL10*((csetemp112 + csetemp113 + 
      csetemp114 + csetemp115)*invXformL03 + (csetemp124 + csetemp125 + 
      csetemp126 + csetemp127)*invXformL13 + (csetemp128 + csetemp129 + 
      csetemp130 + csetemp131)*invXformL23 + (csetemp132 + csetemp133 + 
      csetemp134 + csetemp135)*invXformL33) + invXformL20*((csetemp116 + 
      csetemp117 + csetemp118 + csetemp119)*invXformL03 + (csetemp128 + 
      csetemp129 + csetemp130 + csetemp131)*invXformL13 + (csetemp136 + 
      csetemp137 + csetemp138 + csetemp139)*invXformL23 + (csetemp140 + 
      csetemp141 + csetemp142 + csetemp143)*invXformL33) + 
      invXformL30*((csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL03 + (csetemp132 + csetemp133 + csetemp134 + 
      csetemp135)*invXformL13 + (csetemp140 + csetemp141 + csetemp142 + 
      csetemp143)*invXformL23 + (csetemp144 + csetemp145 + csetemp146 + 
      csetemp147)*invXformL33);
    
    CCTK_REAL dg4032 = invXformL00*((csetemp148 + csetemp149 + csetemp150 
      + csetemp151)*invXformL03 + (csetemp152 + csetemp153 + csetemp154 + 
      csetemp155)*invXformL13 + (csetemp156 + csetemp157 + csetemp158 + 
      csetemp159)*invXformL23 + (csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL33) + invXformL10*((csetemp152 + csetemp153 + 
      csetemp154 + csetemp155)*invXformL03 + (csetemp164 + csetemp165 + 
      csetemp166 + csetemp167)*invXformL13 + (csetemp168 + csetemp169 + 
      csetemp170 + csetemp171)*invXformL23 + (csetemp172 + csetemp173 + 
      csetemp174 + csetemp175)*invXformL33) + invXformL20*((csetemp156 + 
      csetemp157 + csetemp158 + csetemp159)*invXformL03 + (csetemp168 + 
      csetemp169 + csetemp170 + csetemp171)*invXformL13 + (csetemp176 + 
      csetemp177 + csetemp178 + csetemp179)*invXformL23 + (csetemp180 + 
      csetemp181 + csetemp182 + csetemp183)*invXformL33) + 
      invXformL30*((csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL03 + (csetemp172 + csetemp173 + csetemp174 + 
      csetemp175)*invXformL13 + (csetemp180 + csetemp181 + csetemp182 + 
      csetemp183)*invXformL23 + (csetemp184 + csetemp185 + csetemp186 + 
      csetemp187)*invXformL33);
    
    CCTK_REAL dg4033 = invXformL00*((csetemp188 + csetemp189 + csetemp190 
      + csetemp191)*invXformL03 + (csetemp192 + csetemp193 + csetemp194 + 
      csetemp195)*invXformL13 + (csetemp196 + csetemp197 + csetemp198 + 
      csetemp199)*invXformL23 + (csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL33) + invXformL10*((csetemp192 + csetemp193 + 
      csetemp194 + csetemp195)*invXformL03 + (csetemp204 + csetemp205 + 
      csetemp206 + csetemp207)*invXformL13 + (csetemp208 + csetemp209 + 
      csetemp210 + csetemp211)*invXformL23 + (csetemp212 + csetemp213 + 
      csetemp214 + csetemp215)*invXformL33) + invXformL20*((csetemp196 + 
      csetemp197 + csetemp198 + csetemp199)*invXformL03 + (csetemp208 + 
      csetemp209 + csetemp210 + csetemp211)*invXformL13 + (csetemp216 + 
      csetemp217 + csetemp218 + csetemp219)*invXformL23 + (csetemp220 + 
      csetemp221 + csetemp222 + csetemp223)*invXformL33) + 
      invXformL30*((csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL03 + (csetemp212 + csetemp213 + csetemp214 + 
      csetemp215)*invXformL13 + (csetemp220 + csetemp221 + csetemp222 + 
      csetemp223)*invXformL23 + (csetemp224 + csetemp225 + csetemp226 + 
      csetemp227)*invXformL33);
    
    CCTK_REAL dg4110 = invXformL21*((csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL01 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL11 + (csetemp96 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL21 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL31) + invXformL31*((csetemp80 + csetemp81 + 
      csetemp82 + csetemp83)*invXformL01 + (csetemp92 + csetemp93 + csetemp94 
      + csetemp95)*invXformL11 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL21 + (csetemp104 + csetemp105 + csetemp106 + 
      csetemp107)*invXformL31) + invXformL01*((csetemp68 + csetemp69 + 
      csetemp70 + csetemp71)*invXformL01 + (csetemp72 + csetemp73 + csetemp74 
      + csetemp75)*invXformL11 + (csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL21 + (csetemp80 + csetemp81 + csetemp82 + 
      csetemp83)*invXformL31) + invXformL11*((csetemp72 + csetemp73 + 
      csetemp74 + csetemp75)*invXformL01 + (csetemp84 + csetemp85 + csetemp86 
      + csetemp87)*invXformL11 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL21 + (csetemp92 + csetemp93 + csetemp94 + 
      csetemp95)*invXformL31);
    
    CCTK_REAL dg4111 = invXformL01*((csetemp108 + csetemp109 + csetemp110 
      + csetemp111)*invXformL01 + (csetemp112 + csetemp113 + csetemp114 + 
      csetemp115)*invXformL11 + (csetemp116 + csetemp117 + csetemp118 + 
      csetemp119)*invXformL21 + (csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL31) + invXformL11*((csetemp112 + csetemp113 + 
      csetemp114 + csetemp115)*invXformL01 + (csetemp124 + csetemp125 + 
      csetemp126 + csetemp127)*invXformL11 + (csetemp128 + csetemp129 + 
      csetemp130 + csetemp131)*invXformL21 + (csetemp132 + csetemp133 + 
      csetemp134 + csetemp135)*invXformL31) + invXformL21*((csetemp116 + 
      csetemp117 + csetemp118 + csetemp119)*invXformL01 + (csetemp128 + 
      csetemp129 + csetemp130 + csetemp131)*invXformL11 + (csetemp136 + 
      csetemp137 + csetemp138 + csetemp139)*invXformL21 + (csetemp140 + 
      csetemp141 + csetemp142 + csetemp143)*invXformL31) + 
      invXformL31*((csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL01 + (csetemp132 + csetemp133 + csetemp134 + 
      csetemp135)*invXformL11 + (csetemp140 + csetemp141 + csetemp142 + 
      csetemp143)*invXformL21 + (csetemp144 + csetemp145 + csetemp146 + 
      csetemp147)*invXformL31);
    
    CCTK_REAL dg4112 = invXformL01*((csetemp148 + csetemp149 + csetemp150 
      + csetemp151)*invXformL01 + (csetemp152 + csetemp153 + csetemp154 + 
      csetemp155)*invXformL11 + (csetemp156 + csetemp157 + csetemp158 + 
      csetemp159)*invXformL21 + (csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL31) + invXformL11*((csetemp152 + csetemp153 + 
      csetemp154 + csetemp155)*invXformL01 + (csetemp164 + csetemp165 + 
      csetemp166 + csetemp167)*invXformL11 + (csetemp168 + csetemp169 + 
      csetemp170 + csetemp171)*invXformL21 + (csetemp172 + csetemp173 + 
      csetemp174 + csetemp175)*invXformL31) + invXformL21*((csetemp156 + 
      csetemp157 + csetemp158 + csetemp159)*invXformL01 + (csetemp168 + 
      csetemp169 + csetemp170 + csetemp171)*invXformL11 + (csetemp176 + 
      csetemp177 + csetemp178 + csetemp179)*invXformL21 + (csetemp180 + 
      csetemp181 + csetemp182 + csetemp183)*invXformL31) + 
      invXformL31*((csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL01 + (csetemp172 + csetemp173 + csetemp174 + 
      csetemp175)*invXformL11 + (csetemp180 + csetemp181 + csetemp182 + 
      csetemp183)*invXformL21 + (csetemp184 + csetemp185 + csetemp186 + 
      csetemp187)*invXformL31);
    
    CCTK_REAL dg4113 = invXformL01*((csetemp188 + csetemp189 + csetemp190 
      + csetemp191)*invXformL01 + (csetemp192 + csetemp193 + csetemp194 + 
      csetemp195)*invXformL11 + (csetemp196 + csetemp197 + csetemp198 + 
      csetemp199)*invXformL21 + (csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL31) + invXformL11*((csetemp192 + csetemp193 + 
      csetemp194 + csetemp195)*invXformL01 + (csetemp204 + csetemp205 + 
      csetemp206 + csetemp207)*invXformL11 + (csetemp208 + csetemp209 + 
      csetemp210 + csetemp211)*invXformL21 + (csetemp212 + csetemp213 + 
      csetemp214 + csetemp215)*invXformL31) + invXformL21*((csetemp196 + 
      csetemp197 + csetemp198 + csetemp199)*invXformL01 + (csetemp208 + 
      csetemp209 + csetemp210 + csetemp211)*invXformL11 + (csetemp216 + 
      csetemp217 + csetemp218 + csetemp219)*invXformL21 + (csetemp220 + 
      csetemp221 + csetemp222 + csetemp223)*invXformL31) + 
      invXformL31*((csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL01 + (csetemp212 + csetemp213 + csetemp214 + 
      csetemp215)*invXformL11 + (csetemp220 + csetemp221 + csetemp222 + 
      csetemp223)*invXformL21 + (csetemp224 + csetemp225 + csetemp226 + 
      csetemp227)*invXformL31);
    
    CCTK_REAL dg4120 = invXformL21*((csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL02 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL12 + (csetemp96 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL22 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL32) + invXformL31*((csetemp80 + csetemp81 + 
      csetemp82 + csetemp83)*invXformL02 + (csetemp92 + csetemp93 + csetemp94 
      + csetemp95)*invXformL12 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL22 + (csetemp104 + csetemp105 + csetemp106 + 
      csetemp107)*invXformL32) + invXformL01*((csetemp68 + csetemp69 + 
      csetemp70 + csetemp71)*invXformL02 + (csetemp72 + csetemp73 + csetemp74 
      + csetemp75)*invXformL12 + (csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL22 + (csetemp80 + csetemp81 + csetemp82 + 
      csetemp83)*invXformL32) + invXformL11*((csetemp72 + csetemp73 + 
      csetemp74 + csetemp75)*invXformL02 + (csetemp84 + csetemp85 + csetemp86 
      + csetemp87)*invXformL12 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL22 + (csetemp92 + csetemp93 + csetemp94 + 
      csetemp95)*invXformL32);
    
    CCTK_REAL dg4121 = invXformL01*((csetemp108 + csetemp109 + csetemp110 
      + csetemp111)*invXformL02 + (csetemp112 + csetemp113 + csetemp114 + 
      csetemp115)*invXformL12 + (csetemp116 + csetemp117 + csetemp118 + 
      csetemp119)*invXformL22 + (csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL32) + invXformL11*((csetemp112 + csetemp113 + 
      csetemp114 + csetemp115)*invXformL02 + (csetemp124 + csetemp125 + 
      csetemp126 + csetemp127)*invXformL12 + (csetemp128 + csetemp129 + 
      csetemp130 + csetemp131)*invXformL22 + (csetemp132 + csetemp133 + 
      csetemp134 + csetemp135)*invXformL32) + invXformL21*((csetemp116 + 
      csetemp117 + csetemp118 + csetemp119)*invXformL02 + (csetemp128 + 
      csetemp129 + csetemp130 + csetemp131)*invXformL12 + (csetemp136 + 
      csetemp137 + csetemp138 + csetemp139)*invXformL22 + (csetemp140 + 
      csetemp141 + csetemp142 + csetemp143)*invXformL32) + 
      invXformL31*((csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL02 + (csetemp132 + csetemp133 + csetemp134 + 
      csetemp135)*invXformL12 + (csetemp140 + csetemp141 + csetemp142 + 
      csetemp143)*invXformL22 + (csetemp144 + csetemp145 + csetemp146 + 
      csetemp147)*invXformL32);
    
    CCTK_REAL dg4122 = invXformL01*((csetemp148 + csetemp149 + csetemp150 
      + csetemp151)*invXformL02 + (csetemp152 + csetemp153 + csetemp154 + 
      csetemp155)*invXformL12 + (csetemp156 + csetemp157 + csetemp158 + 
      csetemp159)*invXformL22 + (csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL32) + invXformL11*((csetemp152 + csetemp153 + 
      csetemp154 + csetemp155)*invXformL02 + (csetemp164 + csetemp165 + 
      csetemp166 + csetemp167)*invXformL12 + (csetemp168 + csetemp169 + 
      csetemp170 + csetemp171)*invXformL22 + (csetemp172 + csetemp173 + 
      csetemp174 + csetemp175)*invXformL32) + invXformL21*((csetemp156 + 
      csetemp157 + csetemp158 + csetemp159)*invXformL02 + (csetemp168 + 
      csetemp169 + csetemp170 + csetemp171)*invXformL12 + (csetemp176 + 
      csetemp177 + csetemp178 + csetemp179)*invXformL22 + (csetemp180 + 
      csetemp181 + csetemp182 + csetemp183)*invXformL32) + 
      invXformL31*((csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL02 + (csetemp172 + csetemp173 + csetemp174 + 
      csetemp175)*invXformL12 + (csetemp180 + csetemp181 + csetemp182 + 
      csetemp183)*invXformL22 + (csetemp184 + csetemp185 + csetemp186 + 
      csetemp187)*invXformL32);
    
    CCTK_REAL dg4123 = invXformL01*((csetemp188 + csetemp189 + csetemp190 
      + csetemp191)*invXformL02 + (csetemp192 + csetemp193 + csetemp194 + 
      csetemp195)*invXformL12 + (csetemp196 + csetemp197 + csetemp198 + 
      csetemp199)*invXformL22 + (csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL32) + invXformL11*((csetemp192 + csetemp193 + 
      csetemp194 + csetemp195)*invXformL02 + (csetemp204 + csetemp205 + 
      csetemp206 + csetemp207)*invXformL12 + (csetemp208 + csetemp209 + 
      csetemp210 + csetemp211)*invXformL22 + (csetemp212 + csetemp213 + 
      csetemp214 + csetemp215)*invXformL32) + invXformL21*((csetemp196 + 
      csetemp197 + csetemp198 + csetemp199)*invXformL02 + (csetemp208 + 
      csetemp209 + csetemp210 + csetemp211)*invXformL12 + (csetemp216 + 
      csetemp217 + csetemp218 + csetemp219)*invXformL22 + (csetemp220 + 
      csetemp221 + csetemp222 + csetemp223)*invXformL32) + 
      invXformL31*((csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL02 + (csetemp212 + csetemp213 + csetemp214 + 
      csetemp215)*invXformL12 + (csetemp220 + csetemp221 + csetemp222 + 
      csetemp223)*invXformL22 + (csetemp224 + csetemp225 + csetemp226 + 
      csetemp227)*invXformL32);
    
    CCTK_REAL dg4130 = invXformL21*((csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL03 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL13 + (csetemp96 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL23 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL33) + invXformL31*((csetemp80 + csetemp81 + 
      csetemp82 + csetemp83)*invXformL03 + (csetemp92 + csetemp93 + csetemp94 
      + csetemp95)*invXformL13 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL23 + (csetemp104 + csetemp105 + csetemp106 + 
      csetemp107)*invXformL33) + invXformL01*((csetemp68 + csetemp69 + 
      csetemp70 + csetemp71)*invXformL03 + (csetemp72 + csetemp73 + csetemp74 
      + csetemp75)*invXformL13 + (csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL23 + (csetemp80 + csetemp81 + csetemp82 + 
      csetemp83)*invXformL33) + invXformL11*((csetemp72 + csetemp73 + 
      csetemp74 + csetemp75)*invXformL03 + (csetemp84 + csetemp85 + csetemp86 
      + csetemp87)*invXformL13 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL23 + (csetemp92 + csetemp93 + csetemp94 + 
      csetemp95)*invXformL33);
    
    CCTK_REAL dg4131 = invXformL01*((csetemp108 + csetemp109 + csetemp110 
      + csetemp111)*invXformL03 + (csetemp112 + csetemp113 + csetemp114 + 
      csetemp115)*invXformL13 + (csetemp116 + csetemp117 + csetemp118 + 
      csetemp119)*invXformL23 + (csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL33) + invXformL11*((csetemp112 + csetemp113 + 
      csetemp114 + csetemp115)*invXformL03 + (csetemp124 + csetemp125 + 
      csetemp126 + csetemp127)*invXformL13 + (csetemp128 + csetemp129 + 
      csetemp130 + csetemp131)*invXformL23 + (csetemp132 + csetemp133 + 
      csetemp134 + csetemp135)*invXformL33) + invXformL21*((csetemp116 + 
      csetemp117 + csetemp118 + csetemp119)*invXformL03 + (csetemp128 + 
      csetemp129 + csetemp130 + csetemp131)*invXformL13 + (csetemp136 + 
      csetemp137 + csetemp138 + csetemp139)*invXformL23 + (csetemp140 + 
      csetemp141 + csetemp142 + csetemp143)*invXformL33) + 
      invXformL31*((csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL03 + (csetemp132 + csetemp133 + csetemp134 + 
      csetemp135)*invXformL13 + (csetemp140 + csetemp141 + csetemp142 + 
      csetemp143)*invXformL23 + (csetemp144 + csetemp145 + csetemp146 + 
      csetemp147)*invXformL33);
    
    CCTK_REAL dg4132 = invXformL01*((csetemp148 + csetemp149 + csetemp150 
      + csetemp151)*invXformL03 + (csetemp152 + csetemp153 + csetemp154 + 
      csetemp155)*invXformL13 + (csetemp156 + csetemp157 + csetemp158 + 
      csetemp159)*invXformL23 + (csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL33) + invXformL11*((csetemp152 + csetemp153 + 
      csetemp154 + csetemp155)*invXformL03 + (csetemp164 + csetemp165 + 
      csetemp166 + csetemp167)*invXformL13 + (csetemp168 + csetemp169 + 
      csetemp170 + csetemp171)*invXformL23 + (csetemp172 + csetemp173 + 
      csetemp174 + csetemp175)*invXformL33) + invXformL21*((csetemp156 + 
      csetemp157 + csetemp158 + csetemp159)*invXformL03 + (csetemp168 + 
      csetemp169 + csetemp170 + csetemp171)*invXformL13 + (csetemp176 + 
      csetemp177 + csetemp178 + csetemp179)*invXformL23 + (csetemp180 + 
      csetemp181 + csetemp182 + csetemp183)*invXformL33) + 
      invXformL31*((csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL03 + (csetemp172 + csetemp173 + csetemp174 + 
      csetemp175)*invXformL13 + (csetemp180 + csetemp181 + csetemp182 + 
      csetemp183)*invXformL23 + (csetemp184 + csetemp185 + csetemp186 + 
      csetemp187)*invXformL33);
    
    CCTK_REAL dg4133 = invXformL01*((csetemp188 + csetemp189 + csetemp190 
      + csetemp191)*invXformL03 + (csetemp192 + csetemp193 + csetemp194 + 
      csetemp195)*invXformL13 + (csetemp196 + csetemp197 + csetemp198 + 
      csetemp199)*invXformL23 + (csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL33) + invXformL11*((csetemp192 + csetemp193 + 
      csetemp194 + csetemp195)*invXformL03 + (csetemp204 + csetemp205 + 
      csetemp206 + csetemp207)*invXformL13 + (csetemp208 + csetemp209 + 
      csetemp210 + csetemp211)*invXformL23 + (csetemp212 + csetemp213 + 
      csetemp214 + csetemp215)*invXformL33) + invXformL21*((csetemp196 + 
      csetemp197 + csetemp198 + csetemp199)*invXformL03 + (csetemp208 + 
      csetemp209 + csetemp210 + csetemp211)*invXformL13 + (csetemp216 + 
      csetemp217 + csetemp218 + csetemp219)*invXformL23 + (csetemp220 + 
      csetemp221 + csetemp222 + csetemp223)*invXformL33) + 
      invXformL31*((csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL03 + (csetemp212 + csetemp213 + csetemp214 + 
      csetemp215)*invXformL13 + (csetemp220 + csetemp221 + csetemp222 + 
      csetemp223)*invXformL23 + (csetemp224 + csetemp225 + csetemp226 + 
      csetemp227)*invXformL33);
    
    CCTK_REAL dg4220 = invXformL22*((csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL02 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL12 + (csetemp96 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL22 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL32) + invXformL32*((csetemp80 + csetemp81 + 
      csetemp82 + csetemp83)*invXformL02 + (csetemp92 + csetemp93 + csetemp94 
      + csetemp95)*invXformL12 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL22 + (csetemp104 + csetemp105 + csetemp106 + 
      csetemp107)*invXformL32) + invXformL02*((csetemp68 + csetemp69 + 
      csetemp70 + csetemp71)*invXformL02 + (csetemp72 + csetemp73 + csetemp74 
      + csetemp75)*invXformL12 + (csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL22 + (csetemp80 + csetemp81 + csetemp82 + 
      csetemp83)*invXformL32) + invXformL12*((csetemp72 + csetemp73 + 
      csetemp74 + csetemp75)*invXformL02 + (csetemp84 + csetemp85 + csetemp86 
      + csetemp87)*invXformL12 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL22 + (csetemp92 + csetemp93 + csetemp94 + 
      csetemp95)*invXformL32);
    
    CCTK_REAL dg4221 = invXformL02*((csetemp108 + csetemp109 + csetemp110 
      + csetemp111)*invXformL02 + (csetemp112 + csetemp113 + csetemp114 + 
      csetemp115)*invXformL12 + (csetemp116 + csetemp117 + csetemp118 + 
      csetemp119)*invXformL22 + (csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL32) + invXformL12*((csetemp112 + csetemp113 + 
      csetemp114 + csetemp115)*invXformL02 + (csetemp124 + csetemp125 + 
      csetemp126 + csetemp127)*invXformL12 + (csetemp128 + csetemp129 + 
      csetemp130 + csetemp131)*invXformL22 + (csetemp132 + csetemp133 + 
      csetemp134 + csetemp135)*invXformL32) + invXformL22*((csetemp116 + 
      csetemp117 + csetemp118 + csetemp119)*invXformL02 + (csetemp128 + 
      csetemp129 + csetemp130 + csetemp131)*invXformL12 + (csetemp136 + 
      csetemp137 + csetemp138 + csetemp139)*invXformL22 + (csetemp140 + 
      csetemp141 + csetemp142 + csetemp143)*invXformL32) + 
      invXformL32*((csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL02 + (csetemp132 + csetemp133 + csetemp134 + 
      csetemp135)*invXformL12 + (csetemp140 + csetemp141 + csetemp142 + 
      csetemp143)*invXformL22 + (csetemp144 + csetemp145 + csetemp146 + 
      csetemp147)*invXformL32);
    
    CCTK_REAL dg4222 = invXformL02*((csetemp148 + csetemp149 + csetemp150 
      + csetemp151)*invXformL02 + (csetemp152 + csetemp153 + csetemp154 + 
      csetemp155)*invXformL12 + (csetemp156 + csetemp157 + csetemp158 + 
      csetemp159)*invXformL22 + (csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL32) + invXformL12*((csetemp152 + csetemp153 + 
      csetemp154 + csetemp155)*invXformL02 + (csetemp164 + csetemp165 + 
      csetemp166 + csetemp167)*invXformL12 + (csetemp168 + csetemp169 + 
      csetemp170 + csetemp171)*invXformL22 + (csetemp172 + csetemp173 + 
      csetemp174 + csetemp175)*invXformL32) + invXformL22*((csetemp156 + 
      csetemp157 + csetemp158 + csetemp159)*invXformL02 + (csetemp168 + 
      csetemp169 + csetemp170 + csetemp171)*invXformL12 + (csetemp176 + 
      csetemp177 + csetemp178 + csetemp179)*invXformL22 + (csetemp180 + 
      csetemp181 + csetemp182 + csetemp183)*invXformL32) + 
      invXformL32*((csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL02 + (csetemp172 + csetemp173 + csetemp174 + 
      csetemp175)*invXformL12 + (csetemp180 + csetemp181 + csetemp182 + 
      csetemp183)*invXformL22 + (csetemp184 + csetemp185 + csetemp186 + 
      csetemp187)*invXformL32);
    
    CCTK_REAL dg4223 = invXformL02*((csetemp188 + csetemp189 + csetemp190 
      + csetemp191)*invXformL02 + (csetemp192 + csetemp193 + csetemp194 + 
      csetemp195)*invXformL12 + (csetemp196 + csetemp197 + csetemp198 + 
      csetemp199)*invXformL22 + (csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL32) + invXformL12*((csetemp192 + csetemp193 + 
      csetemp194 + csetemp195)*invXformL02 + (csetemp204 + csetemp205 + 
      csetemp206 + csetemp207)*invXformL12 + (csetemp208 + csetemp209 + 
      csetemp210 + csetemp211)*invXformL22 + (csetemp212 + csetemp213 + 
      csetemp214 + csetemp215)*invXformL32) + invXformL22*((csetemp196 + 
      csetemp197 + csetemp198 + csetemp199)*invXformL02 + (csetemp208 + 
      csetemp209 + csetemp210 + csetemp211)*invXformL12 + (csetemp216 + 
      csetemp217 + csetemp218 + csetemp219)*invXformL22 + (csetemp220 + 
      csetemp221 + csetemp222 + csetemp223)*invXformL32) + 
      invXformL32*((csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL02 + (csetemp212 + csetemp213 + csetemp214 + 
      csetemp215)*invXformL12 + (csetemp220 + csetemp221 + csetemp222 + 
      csetemp223)*invXformL22 + (csetemp224 + csetemp225 + csetemp226 + 
      csetemp227)*invXformL32);
    
    CCTK_REAL dg4230 = invXformL22*((csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL03 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL13 + (csetemp96 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL23 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL33) + invXformL32*((csetemp80 + csetemp81 + 
      csetemp82 + csetemp83)*invXformL03 + (csetemp92 + csetemp93 + csetemp94 
      + csetemp95)*invXformL13 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL23 + (csetemp104 + csetemp105 + csetemp106 + 
      csetemp107)*invXformL33) + invXformL02*((csetemp68 + csetemp69 + 
      csetemp70 + csetemp71)*invXformL03 + (csetemp72 + csetemp73 + csetemp74 
      + csetemp75)*invXformL13 + (csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL23 + (csetemp80 + csetemp81 + csetemp82 + 
      csetemp83)*invXformL33) + invXformL12*((csetemp72 + csetemp73 + 
      csetemp74 + csetemp75)*invXformL03 + (csetemp84 + csetemp85 + csetemp86 
      + csetemp87)*invXformL13 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL23 + (csetemp92 + csetemp93 + csetemp94 + 
      csetemp95)*invXformL33);
    
    CCTK_REAL dg4231 = invXformL02*((csetemp108 + csetemp109 + csetemp110 
      + csetemp111)*invXformL03 + (csetemp112 + csetemp113 + csetemp114 + 
      csetemp115)*invXformL13 + (csetemp116 + csetemp117 + csetemp118 + 
      csetemp119)*invXformL23 + (csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL33) + invXformL12*((csetemp112 + csetemp113 + 
      csetemp114 + csetemp115)*invXformL03 + (csetemp124 + csetemp125 + 
      csetemp126 + csetemp127)*invXformL13 + (csetemp128 + csetemp129 + 
      csetemp130 + csetemp131)*invXformL23 + (csetemp132 + csetemp133 + 
      csetemp134 + csetemp135)*invXformL33) + invXformL22*((csetemp116 + 
      csetemp117 + csetemp118 + csetemp119)*invXformL03 + (csetemp128 + 
      csetemp129 + csetemp130 + csetemp131)*invXformL13 + (csetemp136 + 
      csetemp137 + csetemp138 + csetemp139)*invXformL23 + (csetemp140 + 
      csetemp141 + csetemp142 + csetemp143)*invXformL33) + 
      invXformL32*((csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL03 + (csetemp132 + csetemp133 + csetemp134 + 
      csetemp135)*invXformL13 + (csetemp140 + csetemp141 + csetemp142 + 
      csetemp143)*invXformL23 + (csetemp144 + csetemp145 + csetemp146 + 
      csetemp147)*invXformL33);
    
    CCTK_REAL dg4232 = invXformL02*((csetemp148 + csetemp149 + csetemp150 
      + csetemp151)*invXformL03 + (csetemp152 + csetemp153 + csetemp154 + 
      csetemp155)*invXformL13 + (csetemp156 + csetemp157 + csetemp158 + 
      csetemp159)*invXformL23 + (csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL33) + invXformL12*((csetemp152 + csetemp153 + 
      csetemp154 + csetemp155)*invXformL03 + (csetemp164 + csetemp165 + 
      csetemp166 + csetemp167)*invXformL13 + (csetemp168 + csetemp169 + 
      csetemp170 + csetemp171)*invXformL23 + (csetemp172 + csetemp173 + 
      csetemp174 + csetemp175)*invXformL33) + invXformL22*((csetemp156 + 
      csetemp157 + csetemp158 + csetemp159)*invXformL03 + (csetemp168 + 
      csetemp169 + csetemp170 + csetemp171)*invXformL13 + (csetemp176 + 
      csetemp177 + csetemp178 + csetemp179)*invXformL23 + (csetemp180 + 
      csetemp181 + csetemp182 + csetemp183)*invXformL33) + 
      invXformL32*((csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL03 + (csetemp172 + csetemp173 + csetemp174 + 
      csetemp175)*invXformL13 + (csetemp180 + csetemp181 + csetemp182 + 
      csetemp183)*invXformL23 + (csetemp184 + csetemp185 + csetemp186 + 
      csetemp187)*invXformL33);
    
    CCTK_REAL dg4233 = invXformL02*((csetemp188 + csetemp189 + csetemp190 
      + csetemp191)*invXformL03 + (csetemp192 + csetemp193 + csetemp194 + 
      csetemp195)*invXformL13 + (csetemp196 + csetemp197 + csetemp198 + 
      csetemp199)*invXformL23 + (csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL33) + invXformL12*((csetemp192 + csetemp193 + 
      csetemp194 + csetemp195)*invXformL03 + (csetemp204 + csetemp205 + 
      csetemp206 + csetemp207)*invXformL13 + (csetemp208 + csetemp209 + 
      csetemp210 + csetemp211)*invXformL23 + (csetemp212 + csetemp213 + 
      csetemp214 + csetemp215)*invXformL33) + invXformL22*((csetemp196 + 
      csetemp197 + csetemp198 + csetemp199)*invXformL03 + (csetemp208 + 
      csetemp209 + csetemp210 + csetemp211)*invXformL13 + (csetemp216 + 
      csetemp217 + csetemp218 + csetemp219)*invXformL23 + (csetemp220 + 
      csetemp221 + csetemp222 + csetemp223)*invXformL33) + 
      invXformL32*((csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL03 + (csetemp212 + csetemp213 + csetemp214 + 
      csetemp215)*invXformL13 + (csetemp220 + csetemp221 + csetemp222 + 
      csetemp223)*invXformL23 + (csetemp224 + csetemp225 + csetemp226 + 
      csetemp227)*invXformL33);
    
    CCTK_REAL dg4330 = invXformL23*((csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL03 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL13 + (csetemp96 + csetemp97 + csetemp98 + 
      csetemp99)*invXformL23 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL33) + invXformL33*((csetemp80 + csetemp81 + 
      csetemp82 + csetemp83)*invXformL03 + (csetemp92 + csetemp93 + csetemp94 
      + csetemp95)*invXformL13 + (csetemp100 + csetemp101 + csetemp102 + 
      csetemp103)*invXformL23 + (csetemp104 + csetemp105 + csetemp106 + 
      csetemp107)*invXformL33) + invXformL03*((csetemp68 + csetemp69 + 
      csetemp70 + csetemp71)*invXformL03 + (csetemp72 + csetemp73 + csetemp74 
      + csetemp75)*invXformL13 + (csetemp76 + csetemp77 + csetemp78 + 
      csetemp79)*invXformL23 + (csetemp80 + csetemp81 + csetemp82 + 
      csetemp83)*invXformL33) + invXformL13*((csetemp72 + csetemp73 + 
      csetemp74 + csetemp75)*invXformL03 + (csetemp84 + csetemp85 + csetemp86 
      + csetemp87)*invXformL13 + (csetemp88 + csetemp89 + csetemp90 + 
      csetemp91)*invXformL23 + (csetemp92 + csetemp93 + csetemp94 + 
      csetemp95)*invXformL33);
    
    CCTK_REAL dg4331 = invXformL03*((csetemp108 + csetemp109 + csetemp110 
      + csetemp111)*invXformL03 + (csetemp112 + csetemp113 + csetemp114 + 
      csetemp115)*invXformL13 + (csetemp116 + csetemp117 + csetemp118 + 
      csetemp119)*invXformL23 + (csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL33) + invXformL13*((csetemp112 + csetemp113 + 
      csetemp114 + csetemp115)*invXformL03 + (csetemp124 + csetemp125 + 
      csetemp126 + csetemp127)*invXformL13 + (csetemp128 + csetemp129 + 
      csetemp130 + csetemp131)*invXformL23 + (csetemp132 + csetemp133 + 
      csetemp134 + csetemp135)*invXformL33) + invXformL23*((csetemp116 + 
      csetemp117 + csetemp118 + csetemp119)*invXformL03 + (csetemp128 + 
      csetemp129 + csetemp130 + csetemp131)*invXformL13 + (csetemp136 + 
      csetemp137 + csetemp138 + csetemp139)*invXformL23 + (csetemp140 + 
      csetemp141 + csetemp142 + csetemp143)*invXformL33) + 
      invXformL33*((csetemp120 + csetemp121 + csetemp122 + 
      csetemp123)*invXformL03 + (csetemp132 + csetemp133 + csetemp134 + 
      csetemp135)*invXformL13 + (csetemp140 + csetemp141 + csetemp142 + 
      csetemp143)*invXformL23 + (csetemp144 + csetemp145 + csetemp146 + 
      csetemp147)*invXformL33);
    
    CCTK_REAL dg4332 = invXformL03*((csetemp148 + csetemp149 + csetemp150 
      + csetemp151)*invXformL03 + (csetemp152 + csetemp153 + csetemp154 + 
      csetemp155)*invXformL13 + (csetemp156 + csetemp157 + csetemp158 + 
      csetemp159)*invXformL23 + (csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL33) + invXformL13*((csetemp152 + csetemp153 + 
      csetemp154 + csetemp155)*invXformL03 + (csetemp164 + csetemp165 + 
      csetemp166 + csetemp167)*invXformL13 + (csetemp168 + csetemp169 + 
      csetemp170 + csetemp171)*invXformL23 + (csetemp172 + csetemp173 + 
      csetemp174 + csetemp175)*invXformL33) + invXformL23*((csetemp156 + 
      csetemp157 + csetemp158 + csetemp159)*invXformL03 + (csetemp168 + 
      csetemp169 + csetemp170 + csetemp171)*invXformL13 + (csetemp176 + 
      csetemp177 + csetemp178 + csetemp179)*invXformL23 + (csetemp180 + 
      csetemp181 + csetemp182 + csetemp183)*invXformL33) + 
      invXformL33*((csetemp160 + csetemp161 + csetemp162 + 
      csetemp163)*invXformL03 + (csetemp172 + csetemp173 + csetemp174 + 
      csetemp175)*invXformL13 + (csetemp180 + csetemp181 + csetemp182 + 
      csetemp183)*invXformL23 + (csetemp184 + csetemp185 + csetemp186 + 
      csetemp187)*invXformL33);
    
    CCTK_REAL dg4333 = invXformL03*((csetemp188 + csetemp189 + csetemp190 
      + csetemp191)*invXformL03 + (csetemp192 + csetemp193 + csetemp194 + 
      csetemp195)*invXformL13 + (csetemp196 + csetemp197 + csetemp198 + 
      csetemp199)*invXformL23 + (csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL33) + invXformL13*((csetemp192 + csetemp193 + 
      csetemp194 + csetemp195)*invXformL03 + (csetemp204 + csetemp205 + 
      csetemp206 + csetemp207)*invXformL13 + (csetemp208 + csetemp209 + 
      csetemp210 + csetemp211)*invXformL23 + (csetemp212 + csetemp213 + 
      csetemp214 + csetemp215)*invXformL33) + invXformL23*((csetemp196 + 
      csetemp197 + csetemp198 + csetemp199)*invXformL03 + (csetemp208 + 
      csetemp209 + csetemp210 + csetemp211)*invXformL13 + (csetemp216 + 
      csetemp217 + csetemp218 + csetemp219)*invXformL23 + (csetemp220 + 
      csetemp221 + csetemp222 + csetemp223)*invXformL33) + 
      invXformL33*((csetemp200 + csetemp201 + csetemp202 + 
      csetemp203)*invXformL03 + (csetemp212 + csetemp213 + csetemp214 + 
      csetemp215)*invXformL13 + (csetemp220 + csetemp221 + csetemp222 + 
      csetemp223)*invXformL23 + (csetemp224 + csetemp225 + csetemp226 + 
      csetemp227)*invXformL33);
    
    CCTK_REAL betal1 = g401;
    
    CCTK_REAL betal2 = g402;
    
    CCTK_REAL betal3 = g403;
    
    gxxL = g411;
    
    gxyL = g412;
    
    gxzL = g413;
    
    gyyL = g422;
    
    gyzL = g423;
    
    gzzL = g433;
    
    CCTK_REAL csetemp228 = SQR(gxzL);
    
    CCTK_REAL csetemp229 = SQR(gyzL);
    
    CCTK_REAL csetemp230 = SQR(gxyL);
    
    CCTK_REAL detg = 2*gxyL*gxzL*gyzL + gyyL*(gxxL*gzzL - 
      csetemp228) - gxxL*csetemp229 - gzzL*csetemp230;
    
    CCTK_REAL csetemp231 = INV(detg);
    
    CCTK_REAL gu11 = (gyyL*gzzL - csetemp229)*csetemp231;
    
    CCTK_REAL gu12 = (gxzL*gyzL - gxyL*gzzL)*csetemp231;
    
    CCTK_REAL gu13 = (-(gxzL*gyyL) + gxyL*gyzL)*csetemp231;
    
    CCTK_REAL gu22 = (gxxL*gzzL - csetemp228)*csetemp231;
    
    CCTK_REAL gu23 = (gxyL*gxzL - gxxL*gyzL)*csetemp231;
    
    CCTK_REAL gu33 = (gxxL*gyyL - csetemp230)*csetemp231;
    
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
    
    CCTK_REAL csetemp232 = SQR(gu11);
    
    CCTK_REAL csetemp233 = SQR(gu12);
    
    CCTK_REAL csetemp234 = SQR(gu13);
    
    CCTK_REAL dtgu11 = -(csetemp232*dtg11) - csetemp233*dtg22 - 
      csetemp234*dtg33 - 2*dtg12*gu11*gu12 - 2*dtg13*gu11*gu13 - 
      2*dtg23*gu12*gu13;
    
    CCTK_REAL dtgu12 = gu12*(-(dtg11*gu11) - dtg13*gu13 - dtg22*gu22) + 
      dtg12*(-csetemp233 - gu11*gu22) + (-(dtg13*gu11) - dtg33*gu13)*gu23 + 
      dtg23*(-(gu13*gu22) - gu12*gu23);
    
    CCTK_REAL dtgu13 = (-(dtg12*gu11) - dtg22*gu12)*gu23 - dtg23*gu12*gu33 
      + gu13*(-(dtg11*gu11) - dtg12*gu12 - dtg23*gu23 - dtg33*gu33) + 
      dtg13*(-csetemp234 - gu11*gu33);
    
    CCTK_REAL csetemp235 = SQR(gu22);
    
    CCTK_REAL csetemp236 = SQR(gu23);
    
    CCTK_REAL dtgu22 = -(csetemp233*dtg11) - csetemp235*dtg22 - 
      csetemp236*dtg33 - 2*dtg12*gu12*gu22 - 2*dtg13*gu12*gu23 - 
      2*dtg23*gu22*gu23;
    
    CCTK_REAL dtgu23 = gu13*(-(dtg11*gu12) - dtg12*gu22 - dtg13*gu23) - 
      dtg13*gu12*gu33 + gu23*(-(dtg12*gu12) - dtg22*gu22 - dtg33*gu33) + 
      dtg23*(-csetemp236 - gu22*gu33);
    
    CCTK_REAL csetemp237 = SQR(gu33);
    
    CCTK_REAL dtgu33 = -(csetemp234*dtg11) - csetemp236*dtg22 - 
      csetemp237*dtg33 - 2*dtg12*gu13*gu23 - 2*dtg13*gu13*gu33 - 
      2*dtg23*gu23*gu33;
    
    CCTK_REAL dgu111 = -(csetemp232*dg111) - csetemp233*dg221 - 
      csetemp234*dg331 - 2*dg121*gu11*gu12 - 2*dg131*gu11*gu13 - 
      2*dg231*gu12*gu13;
    
    CCTK_REAL dgu121 = gu12*(-(dg111*gu11) - dg131*gu13 - dg221*gu22) + 
      dg121*(-csetemp233 - gu11*gu22) + (-(dg131*gu11) - dg331*gu13)*gu23 + 
      dg231*(-(gu13*gu22) - gu12*gu23);
    
    CCTK_REAL dgu131 = (-(dg121*gu11) - dg221*gu12)*gu23 - dg231*gu12*gu33 
      + gu13*(-(dg111*gu11) - dg121*gu12 - dg231*gu23 - dg331*gu33) + 
      dg131*(-csetemp234 - gu11*gu33);
    
    CCTK_REAL dgu221 = -(csetemp233*dg111) - csetemp235*dg221 - 
      csetemp236*dg331 - 2*dg121*gu12*gu22 - 2*dg131*gu12*gu23 - 
      2*dg231*gu22*gu23;
    
    CCTK_REAL dgu231 = gu13*(-(dg111*gu12) - dg121*gu22 - dg131*gu23) - 
      dg131*gu12*gu33 + gu23*(-(dg121*gu12) - dg221*gu22 - dg331*gu33) + 
      dg231*(-csetemp236 - gu22*gu33);
    
    CCTK_REAL dgu331 = -(csetemp234*dg111) - csetemp236*dg221 - 
      csetemp237*dg331 - 2*dg121*gu13*gu23 - 2*dg131*gu13*gu33 - 
      2*dg231*gu23*gu33;
    
    CCTK_REAL dgu112 = -(csetemp232*dg112) - csetemp233*dg222 - 
      csetemp234*dg332 - 2*dg122*gu11*gu12 - 2*dg132*gu11*gu13 - 
      2*dg232*gu12*gu13;
    
    CCTK_REAL dgu122 = gu12*(-(dg112*gu11) - dg132*gu13 - dg222*gu22) + 
      dg122*(-csetemp233 - gu11*gu22) + (-(dg132*gu11) - dg332*gu13)*gu23 + 
      dg232*(-(gu13*gu22) - gu12*gu23);
    
    CCTK_REAL dgu132 = (-(dg122*gu11) - dg222*gu12)*gu23 - dg232*gu12*gu33 
      + gu13*(-(dg112*gu11) - dg122*gu12 - dg232*gu23 - dg332*gu33) + 
      dg132*(-csetemp234 - gu11*gu33);
    
    CCTK_REAL dgu222 = -(csetemp233*dg112) - csetemp235*dg222 - 
      csetemp236*dg332 - 2*dg122*gu12*gu22 - 2*dg132*gu12*gu23 - 
      2*dg232*gu22*gu23;
    
    CCTK_REAL dgu232 = gu13*(-(dg112*gu12) - dg122*gu22 - dg132*gu23) - 
      dg132*gu12*gu33 + gu23*(-(dg122*gu12) - dg222*gu22 - dg332*gu33) + 
      dg232*(-csetemp236 - gu22*gu33);
    
    CCTK_REAL dgu332 = -(csetemp234*dg112) - csetemp236*dg222 - 
      csetemp237*dg332 - 2*dg122*gu13*gu23 - 2*dg132*gu13*gu33 - 
      2*dg232*gu23*gu33;
    
    CCTK_REAL dgu113 = -(csetemp232*dg113) - csetemp233*dg223 - 
      csetemp234*dg333 - 2*dg123*gu11*gu12 - 2*dg133*gu11*gu13 - 
      2*dg233*gu12*gu13;
    
    CCTK_REAL dgu123 = gu12*(-(dg113*gu11) - dg133*gu13 - dg223*gu22) + 
      dg123*(-csetemp233 - gu11*gu22) + (-(dg133*gu11) - dg333*gu13)*gu23 + 
      dg233*(-(gu13*gu22) - gu12*gu23);
    
    CCTK_REAL dgu133 = (-(dg123*gu11) - dg223*gu12)*gu23 - dg233*gu12*gu33 
      + gu13*(-(dg113*gu11) - dg123*gu12 - dg233*gu23 - dg333*gu33) + 
      dg133*(-csetemp234 - gu11*gu33);
    
    CCTK_REAL dgu223 = -(csetemp233*dg113) - csetemp235*dg223 - 
      csetemp236*dg333 - 2*dg123*gu12*gu22 - 2*dg133*gu12*gu23 - 
      2*dg233*gu22*gu23;
    
    CCTK_REAL dgu233 = gu13*(-(dg113*gu12) - dg123*gu22 - dg133*gu23) - 
      dg133*gu12*gu33 + gu23*(-(dg123*gu12) - dg223*gu22 - dg333*gu33) + 
      dg233*(-csetemp236 - gu22*gu33);
    
    CCTK_REAL dgu333 = -(csetemp234*dg113) - csetemp236*dg223 - 
      csetemp237*dg333 - 2*dg123*gu13*gu23 - 2*dg133*gu13*gu33 - 
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
    
    CCTK_REAL csetemp238 = INV(alpL);
    
    CCTK_REAL dtalpL = 0.5*csetemp238*(-dg4000 + dtbetasq);
    
    CCTK_REAL kxxL = 0.5*csetemp238*(2*(gxxL*dbeta11 + gxyL*dbeta21 + 
      gxzL*dbeta31) + betaxL*dg111 + betayL*dg112 + betazL*dg113 - 
      dtg11);
    
    CCTK_REAL kxyL = 0.5*csetemp238*(gxxL*dbeta12 + gyyL*dbeta21 + 
      gxyL*(dbeta11 + dbeta22) + gyzL*dbeta31 + gxzL*dbeta32 + 
      betaxL*dg121 + betayL*dg122 + betazL*dg123 - dtg12);
    
    CCTK_REAL kxzL = 0.5*csetemp238*(gxxL*dbeta13 + gyzL*dbeta21 + 
      gxyL*dbeta23 + gzzL*dbeta31 + gxzL*(dbeta11 + dbeta33) + 
      betaxL*dg131 + betayL*dg132 + betazL*dg133 - dtg13);
    
    CCTK_REAL kyyL = 0.5*csetemp238*(2*(gxyL*dbeta12 + gyyL*dbeta22 + 
      gyzL*dbeta32) + betaxL*dg221 + betayL*dg222 + betazL*dg223 - 
      dtg22);
    
    CCTK_REAL kyzL = 0.5*csetemp238*(gxzL*dbeta12 + gxyL*dbeta13 + 
      gyyL*dbeta23 + gzzL*dbeta32 + gyzL*(dbeta22 + dbeta33) + 
      betaxL*dg231 + betayL*dg232 + betazL*dg233 - dtg23);
    
    CCTK_REAL kzzL = 0.5*csetemp238*(2*(gxzL*dbeta13 + gyzL*dbeta23 + 
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
  CCTK_ENDLOOP3(ModifiedSchwarzschildBL_always);
}

extern "C" void ModifiedSchwarzschildBL_always(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ModifiedSchwarzschildBL_always_Body");
  }
  
  if (cctk_iteration % ModifiedSchwarzschildBL_always_calc_every != ModifiedSchwarzschildBL_always_calc_offset)
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
  GenericFD_AssertGroupStorage(cctkGH, "ModifiedSchwarzschildBL_always", 7, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, ModifiedSchwarzschildBL_always_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ModifiedSchwarzschildBL_always_Body");
  }
}
