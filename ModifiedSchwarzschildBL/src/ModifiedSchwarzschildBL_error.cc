/*  File produced by Kranc */

#define KRANC_C

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "Differencing.h"
#include "loopcontrol.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define INV(x) ((CCTK_REAL)1.0 / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * SQR(x))
#define QAD(x) (SQR(SQR(x)))

namespace ModifiedSchwarzschildBL {


static void ModifiedSchwarzschildBL_error_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);
  const CCTK_REAL cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_ORIGIN_SPACE(0));
  const CCTK_REAL cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_ORIGIN_SPACE(1));
  const CCTK_REAL cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_ORIGIN_SPACE(2));
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_TIME);
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(0));
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(1));
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(2));
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = INV(dx);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = INV(dy);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = INV(dz);
  const CCTK_REAL khalf CCTK_ATTRIBUTE_UNUSED = 0.5;
  const CCTK_REAL kthird CCTK_ATTRIBUTE_UNUSED = 
    0.333333333333333333333333333333;
  const CCTK_REAL ktwothird CCTK_ATTRIBUTE_UNUSED = 
    0.666666666666666666666666666667;
  const CCTK_REAL kfourthird CCTK_ATTRIBUTE_UNUSED = 
    1.33333333333333333333333333333;
  const CCTK_REAL hdxi CCTK_ATTRIBUTE_UNUSED = 0.5*dxi;
  const CCTK_REAL hdyi CCTK_ATTRIBUTE_UNUSED = 0.5*dyi;
  const CCTK_REAL hdzi CCTK_ATTRIBUTE_UNUSED = 0.5*dzi;
  /* Initialize predefined quantities */
  /* Assign local copies of arrays functions */
  
  
  /* Calculate temporaries and arrays functions */
  /* Copy local copies back to grid functions */
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  #pragma omp parallel
  CCTK_LOOP3(ModifiedSchwarzschildBL_error,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL CCTK_ATTRIBUTE_UNUSED = alp[index];
    CCTK_REAL betaxL CCTK_ATTRIBUTE_UNUSED = betax[index];
    CCTK_REAL betayL CCTK_ATTRIBUTE_UNUSED = betay[index];
    CCTK_REAL betazL CCTK_ATTRIBUTE_UNUSED = betaz[index];
    CCTK_REAL dtalpL CCTK_ATTRIBUTE_UNUSED = dtalp[index];
    CCTK_REAL dtbetaxL CCTK_ATTRIBUTE_UNUSED = dtbetax[index];
    CCTK_REAL dtbetayL CCTK_ATTRIBUTE_UNUSED = dtbetay[index];
    CCTK_REAL dtbetazL CCTK_ATTRIBUTE_UNUSED = dtbetaz[index];
    CCTK_REAL gxxL CCTK_ATTRIBUTE_UNUSED = gxx[index];
    CCTK_REAL gxyL CCTK_ATTRIBUTE_UNUSED = gxy[index];
    CCTK_REAL gxzL CCTK_ATTRIBUTE_UNUSED = gxz[index];
    CCTK_REAL gyyL CCTK_ATTRIBUTE_UNUSED = gyy[index];
    CCTK_REAL gyzL CCTK_ATTRIBUTE_UNUSED = gyz[index];
    CCTK_REAL gzzL CCTK_ATTRIBUTE_UNUSED = gzz[index];
    CCTK_REAL kxxL CCTK_ATTRIBUTE_UNUSED = kxx[index];
    CCTK_REAL kxyL CCTK_ATTRIBUTE_UNUSED = kxy[index];
    CCTK_REAL kxzL CCTK_ATTRIBUTE_UNUSED = kxz[index];
    CCTK_REAL kyyL CCTK_ATTRIBUTE_UNUSED = kyy[index];
    CCTK_REAL kyzL CCTK_ATTRIBUTE_UNUSED = kyz[index];
    CCTK_REAL kzzL CCTK_ATTRIBUTE_UNUSED = kzz[index];
    CCTK_REAL xL CCTK_ATTRIBUTE_UNUSED = x[index];
    CCTK_REAL yL CCTK_ATTRIBUTE_UNUSED = y[index];
    CCTK_REAL zL CCTK_ATTRIBUTE_UNUSED = z[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    /* Calculate temporaries and grid functions */
    CCTK_REAL xform1L00 CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL xform1L01 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL xform1L02 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL xform1L03 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL xform1L10 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL csetemp0 CCTK_ATTRIBUTE_UNUSED = cos(ToReal(phi));
    
    CCTK_REAL csetemp1 CCTK_ATTRIBUTE_UNUSED = cos(ToReal(psi));
    
    CCTK_REAL csetemp2 CCTK_ATTRIBUTE_UNUSED = cos(ToReal(theta));
    
    CCTK_REAL csetemp3 CCTK_ATTRIBUTE_UNUSED = sin(ToReal(phi));
    
    CCTK_REAL csetemp4 CCTK_ATTRIBUTE_UNUSED = sin(ToReal(psi));
    
    CCTK_REAL xform1L11 CCTK_ATTRIBUTE_UNUSED = csetemp0*csetemp1 - 
      csetemp2*csetemp3*csetemp4;
    
    CCTK_REAL xform1L12 CCTK_ATTRIBUTE_UNUSED = csetemp1*csetemp3 + 
      csetemp0*csetemp2*csetemp4;
    
    CCTK_REAL csetemp5 CCTK_ATTRIBUTE_UNUSED = sin(ToReal(theta));
    
    CCTK_REAL xform1L13 CCTK_ATTRIBUTE_UNUSED = csetemp4*csetemp5;
    
    CCTK_REAL xform1L20 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL xform1L21 CCTK_ATTRIBUTE_UNUSED = 
      -(csetemp1*csetemp2*csetemp3) - csetemp0*csetemp4;
    
    CCTK_REAL xform1L22 CCTK_ATTRIBUTE_UNUSED = csetemp0*csetemp1*csetemp2 
      - csetemp3*csetemp4;
    
    CCTK_REAL xform1L23 CCTK_ATTRIBUTE_UNUSED = csetemp1*csetemp5;
    
    CCTK_REAL xform1L30 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL xform1L31 CCTK_ATTRIBUTE_UNUSED = csetemp3*csetemp5;
    
    CCTK_REAL xform1L32 CCTK_ATTRIBUTE_UNUSED = -(csetemp0*csetemp5);
    
    CCTK_REAL xform1L33 CCTK_ATTRIBUTE_UNUSED = csetemp2;
    
    CCTK_REAL csetemp6 CCTK_ATTRIBUTE_UNUSED = SQR(ToReal(boostx));
    
    CCTK_REAL csetemp7 CCTK_ATTRIBUTE_UNUSED = SQR(ToReal(boosty));
    
    CCTK_REAL csetemp8 CCTK_ATTRIBUTE_UNUSED = SQR(ToReal(boostz));
    
    CCTK_REAL csetemp9 CCTK_ATTRIBUTE_UNUSED = INV(ToReal(lapsefactor));
    
    CCTK_REAL xform2L00 CCTK_ATTRIBUTE_UNUSED = csetemp9*INV((-1 + 
      csetemp6 + csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(1 - csetemp6 - csetemp7 - csetemp8 + sqrt(1 - csetemp6 - 
      csetemp7 - csetemp8))*(-1 + ToReal(boostx)*ToReal(shiftaddx) + 
      ToReal(boosty)*ToReal(shiftaddy) + ToReal(boostz)*ToReal(shiftaddz));
    
    CCTK_REAL xform2L01 CCTK_ATTRIBUTE_UNUSED = csetemp9*INV((-1 + 
      csetemp6 + csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(1 - csetemp6 - csetemp7 - csetemp8 + sqrt(1 - csetemp6 - 
      csetemp7 - csetemp8))*ToReal(boostx);
    
    CCTK_REAL xform2L02 CCTK_ATTRIBUTE_UNUSED = csetemp9*INV((-1 + 
      csetemp6 + csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(1 - csetemp6 - csetemp7 - csetemp8 + sqrt(1 - csetemp6 - 
      csetemp7 - csetemp8))*ToReal(boosty);
    
    CCTK_REAL xform2L03 CCTK_ATTRIBUTE_UNUSED = csetemp9*INV((-1 + 
      csetemp6 + csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(1 - csetemp6 - csetemp7 - csetemp8 + sqrt(1 - csetemp6 - 
      csetemp7 - csetemp8))*ToReal(boostz);
    
    CCTK_REAL xform2L10 CCTK_ATTRIBUTE_UNUSED = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*((csetemp6 + (-1 + csetemp7 + csetemp8)*(1 + sqrt(1 - 
      csetemp6 - csetemp7 - csetemp8)))*ToReal(shiftaddx) - 
      ToReal(boostx)*(-1 + csetemp6 + csetemp7 + csetemp8 + sqrt(1 - csetemp6 
      - csetemp7 - csetemp8)*(-1 + ToReal(boosty)*ToReal(shiftaddy) + 
      ToReal(boostz)*ToReal(shiftaddz))));
    
    CCTK_REAL xform2L11 CCTK_ATTRIBUTE_UNUSED = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(csetemp6 + (-1 + csetemp7 + csetemp8)*(1 + sqrt(1 - 
      csetemp6 - csetemp7 - csetemp8)));
    
    CCTK_REAL xform2L12 CCTK_ATTRIBUTE_UNUSED = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boostx)*ToReal(boosty));
    
    CCTK_REAL xform2L13 CCTK_ATTRIBUTE_UNUSED = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boostx)*ToReal(boostz));
    
    CCTK_REAL xform2L20 CCTK_ATTRIBUTE_UNUSED = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*((csetemp7 + (-1 + csetemp6 + csetemp8)*(1 + sqrt(1 - 
      csetemp6 - csetemp7 - csetemp8)))*ToReal(shiftaddy) - 
      ToReal(boosty)*(-1 + csetemp6 + csetemp7 + csetemp8 + sqrt(1 - csetemp6 
      - csetemp7 - csetemp8)*(-1 + ToReal(boostx)*ToReal(shiftaddx) + 
      ToReal(boostz)*ToReal(shiftaddz))));
    
    CCTK_REAL xform2L21 CCTK_ATTRIBUTE_UNUSED = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boostx)*ToReal(boosty));
    
    CCTK_REAL xform2L22 CCTK_ATTRIBUTE_UNUSED = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(csetemp7 + (-1 + csetemp6 + csetemp8)*(1 + sqrt(1 - 
      csetemp6 - csetemp7 - csetemp8)));
    
    CCTK_REAL xform2L23 CCTK_ATTRIBUTE_UNUSED = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boosty)*ToReal(boostz));
    
    CCTK_REAL xform2L30 CCTK_ATTRIBUTE_UNUSED = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(-(ToReal(boostz)*(-1 + csetemp6 + csetemp7 + csetemp8 + 
      sqrt(1 - csetemp6 - csetemp7 - csetemp8)*(-1 + 
      ToReal(boostx)*ToReal(shiftaddx) + ToReal(boosty)*ToReal(shiftaddy)))) 
      + (-1 + csetemp6 + csetemp7 + csetemp8 + (-1 + csetemp6 + 
      csetemp7)*sqrt(1 - csetemp6 - csetemp7 - csetemp8))*ToReal(shiftaddz));
    
    CCTK_REAL xform2L31 CCTK_ATTRIBUTE_UNUSED = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boostx)*ToReal(boostz));
    
    CCTK_REAL xform2L32 CCTK_ATTRIBUTE_UNUSED = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boosty)*ToReal(boostz));
    
    CCTK_REAL xform2L33 CCTK_ATTRIBUTE_UNUSED = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(-1 + csetemp6 + csetemp7 + csetemp8 + (-1 + csetemp6 + 
      csetemp7)*sqrt(1 - csetemp6 - csetemp7 - csetemp8));
    
    CCTK_REAL xformL00 CCTK_ATTRIBUTE_UNUSED = xform1L00*xform2L00 + 
      xform1L01*xform2L10 + xform1L02*xform2L20 + xform1L03*xform2L30;
    
    CCTK_REAL xformL01 CCTK_ATTRIBUTE_UNUSED = xform1L00*xform2L01 + 
      xform1L01*xform2L11 + xform1L02*xform2L21 + xform1L03*xform2L31;
    
    CCTK_REAL xformL02 CCTK_ATTRIBUTE_UNUSED = xform1L00*xform2L02 + 
      xform1L01*xform2L12 + xform1L02*xform2L22 + xform1L03*xform2L32;
    
    CCTK_REAL xformL03 CCTK_ATTRIBUTE_UNUSED = xform1L00*xform2L03 + 
      xform1L01*xform2L13 + xform1L02*xform2L23 + xform1L03*xform2L33;
    
    CCTK_REAL xformL10 CCTK_ATTRIBUTE_UNUSED = xform1L10*xform2L00 + 
      xform1L11*xform2L10 + xform1L12*xform2L20 + xform1L13*xform2L30;
    
    CCTK_REAL xformL11 CCTK_ATTRIBUTE_UNUSED = xform1L10*xform2L01 + 
      xform1L11*xform2L11 + xform1L12*xform2L21 + xform1L13*xform2L31;
    
    CCTK_REAL xformL12 CCTK_ATTRIBUTE_UNUSED = xform1L10*xform2L02 + 
      xform1L11*xform2L12 + xform1L12*xform2L22 + xform1L13*xform2L32;
    
    CCTK_REAL xformL13 CCTK_ATTRIBUTE_UNUSED = xform1L10*xform2L03 + 
      xform1L11*xform2L13 + xform1L12*xform2L23 + xform1L13*xform2L33;
    
    CCTK_REAL xformL20 CCTK_ATTRIBUTE_UNUSED = xform1L20*xform2L00 + 
      xform1L21*xform2L10 + xform1L22*xform2L20 + xform1L23*xform2L30;
    
    CCTK_REAL xformL21 CCTK_ATTRIBUTE_UNUSED = xform1L20*xform2L01 + 
      xform1L21*xform2L11 + xform1L22*xform2L21 + xform1L23*xform2L31;
    
    CCTK_REAL xformL22 CCTK_ATTRIBUTE_UNUSED = xform1L20*xform2L02 + 
      xform1L21*xform2L12 + xform1L22*xform2L22 + xform1L23*xform2L32;
    
    CCTK_REAL xformL23 CCTK_ATTRIBUTE_UNUSED = xform1L20*xform2L03 + 
      xform1L21*xform2L13 + xform1L22*xform2L23 + xform1L23*xform2L33;
    
    CCTK_REAL xformL30 CCTK_ATTRIBUTE_UNUSED = xform1L30*xform2L00 + 
      xform1L31*xform2L10 + xform1L32*xform2L20 + xform1L33*xform2L30;
    
    CCTK_REAL xformL31 CCTK_ATTRIBUTE_UNUSED = xform1L30*xform2L01 + 
      xform1L31*xform2L11 + xform1L32*xform2L21 + xform1L33*xform2L31;
    
    CCTK_REAL xformL32 CCTK_ATTRIBUTE_UNUSED = xform1L30*xform2L02 + 
      xform1L31*xform2L12 + xform1L32*xform2L22 + xform1L33*xform2L32;
    
    CCTK_REAL xformL33 CCTK_ATTRIBUTE_UNUSED = xform1L30*xform2L03 + 
      xform1L31*xform2L13 + xform1L32*xform2L23 + xform1L33*xform2L33;
    
    CCTK_REAL xx0 CCTK_ATTRIBUTE_UNUSED = t - ToReal(timeoffset);
    
    CCTK_REAL xx1 CCTK_ATTRIBUTE_UNUSED = xL - ToReal(positionx);
    
    CCTK_REAL xx2 CCTK_ATTRIBUTE_UNUSED = yL - ToReal(positiony);
    
    CCTK_REAL xx3 CCTK_ATTRIBUTE_UNUSED = zL - ToReal(positionz);
    
    CCTK_REAL txx1 CCTK_ATTRIBUTE_UNUSED = xformL10*xx0 + xformL11*xx1 + 
      xformL12*xx2 + xformL13*xx3;
    
    CCTK_REAL txx2 CCTK_ATTRIBUTE_UNUSED = xformL20*xx0 + xformL21*xx1 + 
      xformL22*xx2 + xformL23*xx3;
    
    CCTK_REAL txx3 CCTK_ATTRIBUTE_UNUSED = xformL30*xx0 + xformL31*xx1 + 
      xformL32*xx2 + xformL33*xx3;
    
    CCTK_REAL X CCTK_ATTRIBUTE_UNUSED = txx1;
    
    CCTK_REAL Y CCTK_ATTRIBUTE_UNUSED = txx2;
    
    CCTK_REAL Z CCTK_ATTRIBUTE_UNUSED = txx3;
    
    CCTK_REAL rXYZ CCTK_ATTRIBUTE_UNUSED = sqrt(SQR(X) + SQR(Y) + SQR(Z));
    
    CCTK_REAL tg400 CCTK_ATTRIBUTE_UNUSED = -1;
    
    CCTK_REAL tg401 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tg402 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tg403 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL csetemp10 CCTK_ATTRIBUTE_UNUSED = 0.5*ToReal(M);
    
    CCTK_REAL csetemp11 CCTK_ATTRIBUTE_UNUSED = INV(rXYZ);
    
    CCTK_REAL csetemp12 CCTK_ATTRIBUTE_UNUSED = INV(ToReal(M));
    
    CCTK_REAL csetemp13 CCTK_ATTRIBUTE_UNUSED = SQR(ToReal(M));
    
    CCTK_REAL csetemp14 CCTK_ATTRIBUTE_UNUSED = INV(SQR(Pi));
    
    CCTK_REAL tg411 CCTK_ATTRIBUTE_UNUSED = 0.00390625*QAD(4 + 
      csetemp11*(8*fmin(csetemp10,rXYZ) + csetemp12*(csetemp13*(csetemp14 - 
      csetemp14*cos(4*csetemp12*Pi*fmin(csetemp10,rXYZ))) - 
      8*SQR(fmin(csetemp10,rXYZ)))));
    
    CCTK_REAL tg412 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tg413 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tg422 CCTK_ATTRIBUTE_UNUSED = 0.00390625*QAD(4 + 
      csetemp11*(8*fmin(csetemp10,rXYZ) + csetemp12*(csetemp13*(csetemp14 - 
      csetemp14*cos(4*csetemp12*Pi*fmin(csetemp10,rXYZ))) - 
      8*SQR(fmin(csetemp10,rXYZ)))));
    
    CCTK_REAL tg423 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tg433 CCTK_ATTRIBUTE_UNUSED = 0.00390625*QAD(4 + 
      csetemp11*(8*fmin(csetemp10,rXYZ) + csetemp12*(csetemp13*(csetemp14 - 
      csetemp14*cos(4*csetemp12*Pi*fmin(csetemp10,rXYZ))) - 
      8*SQR(fmin(csetemp10,rXYZ)))));
    
    CCTK_REAL tdg4000 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4001 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4002 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4003 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4010 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4011 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4012 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4013 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4020 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4021 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4022 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4023 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4030 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4031 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4032 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4033 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4110 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL csetemp15 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,-6);
    
    CCTK_REAL csetemp16 CCTK_ATTRIBUTE_UNUSED = 2*rXYZ;
    
    CCTK_REAL csetemp17 CCTK_ATTRIBUTE_UNUSED = INV(QAD(ToReal(M)));
    
    CCTK_REAL csetemp18 CCTK_ATTRIBUTE_UNUSED = pow(Pi,-8);
    
    CCTK_REAL csetemp19 CCTK_ATTRIBUTE_UNUSED = SQR(Pi);
    
    CCTK_REAL csetemp20 CCTK_ATTRIBUTE_UNUSED = SQR(rXYZ);
    
    CCTK_REAL tdg4111 CCTK_ATTRIBUTE_UNUSED = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*X*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*X*CUB(csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) 
      + csetemp19*(-4*csetemp20 + 6*rXYZ*ToReal(M)))*(4*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp12*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4112 CCTK_ATTRIBUTE_UNUSED = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Y*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Y*CUB(csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) 
      + csetemp19*(-4*csetemp20 + 6*rXYZ*ToReal(M)))*(4*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp12*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4113 CCTK_ATTRIBUTE_UNUSED = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Z*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Z*CUB(csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) 
      + csetemp19*(-4*csetemp20 + 6*rXYZ*ToReal(M)))*(4*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp12*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4120 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4121 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4122 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4123 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4130 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4131 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4132 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4133 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4220 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4221 CCTK_ATTRIBUTE_UNUSED = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*X*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*X*CUB(csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) 
      + csetemp19*(-4*csetemp20 + 6*rXYZ*ToReal(M)))*(4*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp12*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4222 CCTK_ATTRIBUTE_UNUSED = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Y*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Y*CUB(csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) 
      + csetemp19*(-4*csetemp20 + 6*rXYZ*ToReal(M)))*(4*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp12*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4223 CCTK_ATTRIBUTE_UNUSED = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Z*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Z*CUB(csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) 
      + csetemp19*(-4*csetemp20 + 6*rXYZ*ToReal(M)))*(4*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp12*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4230 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4231 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4232 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4233 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4330 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4331 CCTK_ATTRIBUTE_UNUSED = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*X*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*X*CUB(csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) 
      + csetemp19*(-4*csetemp20 + 6*rXYZ*ToReal(M)))*(4*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp12*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4332 CCTK_ATTRIBUTE_UNUSED = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Y*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Y*CUB(csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) 
      + csetemp19*(-4*csetemp20 + 6*rXYZ*ToReal(M)))*(4*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp12*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL tdg4333 CCTK_ATTRIBUTE_UNUSED = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Z*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Z*CUB(csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) 
      + csetemp19*(-4*csetemp20 + 6*rXYZ*ToReal(M)))*(4*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(2*csetemp12*Pi*rXYZ)) - 
      2*Pi*rXYZ*sin(4*csetemp12*Pi*rXYZ)*ToReal(M)));
    
    CCTK_REAL g400 CCTK_ATTRIBUTE_UNUSED = 2*(tg423*xformL20*xformL30 + 
      xformL00*(tg401*xformL10 + tg402*xformL20 + tg403*xformL30) + 
      xformL10*(tg412*xformL20 + tg413*xformL30)) + tg400*SQR(xformL00) + 
      tg411*SQR(xformL10) + tg422*SQR(xformL20) + tg433*SQR(xformL30);
    
    CCTK_REAL csetemp21 CCTK_ATTRIBUTE_UNUSED = tg400*xformL01;
    
    CCTK_REAL csetemp22 CCTK_ATTRIBUTE_UNUSED = tg401*xformL11;
    
    CCTK_REAL csetemp23 CCTK_ATTRIBUTE_UNUSED = tg402*xformL21;
    
    CCTK_REAL csetemp24 CCTK_ATTRIBUTE_UNUSED = tg403*xformL31;
    
    CCTK_REAL csetemp25 CCTK_ATTRIBUTE_UNUSED = tg401*xformL01;
    
    CCTK_REAL csetemp26 CCTK_ATTRIBUTE_UNUSED = tg411*xformL11;
    
    CCTK_REAL csetemp27 CCTK_ATTRIBUTE_UNUSED = tg412*xformL21;
    
    CCTK_REAL csetemp28 CCTK_ATTRIBUTE_UNUSED = tg413*xformL31;
    
    CCTK_REAL csetemp29 CCTK_ATTRIBUTE_UNUSED = tg402*xformL01;
    
    CCTK_REAL csetemp30 CCTK_ATTRIBUTE_UNUSED = tg412*xformL11;
    
    CCTK_REAL csetemp31 CCTK_ATTRIBUTE_UNUSED = tg422*xformL21;
    
    CCTK_REAL csetemp32 CCTK_ATTRIBUTE_UNUSED = tg423*xformL31;
    
    CCTK_REAL csetemp33 CCTK_ATTRIBUTE_UNUSED = tg403*xformL01;
    
    CCTK_REAL csetemp34 CCTK_ATTRIBUTE_UNUSED = tg413*xformL11;
    
    CCTK_REAL csetemp35 CCTK_ATTRIBUTE_UNUSED = tg423*xformL21;
    
    CCTK_REAL csetemp36 CCTK_ATTRIBUTE_UNUSED = tg433*xformL31;
    
    CCTK_REAL g401 CCTK_ATTRIBUTE_UNUSED = (csetemp21 + csetemp22 + 
      csetemp23 + csetemp24)*xformL00 + (csetemp25 + csetemp26 + csetemp27 + 
      csetemp28)*xformL10 + (csetemp29 + csetemp30 + csetemp31 + 
      csetemp32)*xformL20 + (csetemp33 + csetemp34 + csetemp35 + 
      csetemp36)*xformL30;
    
    CCTK_REAL csetemp37 CCTK_ATTRIBUTE_UNUSED = tg400*xformL02;
    
    CCTK_REAL csetemp38 CCTK_ATTRIBUTE_UNUSED = tg401*xformL12;
    
    CCTK_REAL csetemp39 CCTK_ATTRIBUTE_UNUSED = tg402*xformL22;
    
    CCTK_REAL csetemp40 CCTK_ATTRIBUTE_UNUSED = tg403*xformL32;
    
    CCTK_REAL csetemp41 CCTK_ATTRIBUTE_UNUSED = tg401*xformL02;
    
    CCTK_REAL csetemp42 CCTK_ATTRIBUTE_UNUSED = tg411*xformL12;
    
    CCTK_REAL csetemp43 CCTK_ATTRIBUTE_UNUSED = tg412*xformL22;
    
    CCTK_REAL csetemp44 CCTK_ATTRIBUTE_UNUSED = tg413*xformL32;
    
    CCTK_REAL csetemp45 CCTK_ATTRIBUTE_UNUSED = tg402*xformL02;
    
    CCTK_REAL csetemp46 CCTK_ATTRIBUTE_UNUSED = tg412*xformL12;
    
    CCTK_REAL csetemp47 CCTK_ATTRIBUTE_UNUSED = tg422*xformL22;
    
    CCTK_REAL csetemp48 CCTK_ATTRIBUTE_UNUSED = tg423*xformL32;
    
    CCTK_REAL csetemp49 CCTK_ATTRIBUTE_UNUSED = tg403*xformL02;
    
    CCTK_REAL csetemp50 CCTK_ATTRIBUTE_UNUSED = tg413*xformL12;
    
    CCTK_REAL csetemp51 CCTK_ATTRIBUTE_UNUSED = tg423*xformL22;
    
    CCTK_REAL csetemp52 CCTK_ATTRIBUTE_UNUSED = tg433*xformL32;
    
    CCTK_REAL g402 CCTK_ATTRIBUTE_UNUSED = (csetemp37 + csetemp38 + 
      csetemp39 + csetemp40)*xformL00 + (csetemp41 + csetemp42 + csetemp43 + 
      csetemp44)*xformL10 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL20 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL30;
    
    CCTK_REAL csetemp53 CCTK_ATTRIBUTE_UNUSED = tg400*xformL03;
    
    CCTK_REAL csetemp54 CCTK_ATTRIBUTE_UNUSED = tg401*xformL13;
    
    CCTK_REAL csetemp55 CCTK_ATTRIBUTE_UNUSED = tg402*xformL23;
    
    CCTK_REAL csetemp56 CCTK_ATTRIBUTE_UNUSED = tg403*xformL33;
    
    CCTK_REAL csetemp57 CCTK_ATTRIBUTE_UNUSED = tg401*xformL03;
    
    CCTK_REAL csetemp58 CCTK_ATTRIBUTE_UNUSED = tg411*xformL13;
    
    CCTK_REAL csetemp59 CCTK_ATTRIBUTE_UNUSED = tg412*xformL23;
    
    CCTK_REAL csetemp60 CCTK_ATTRIBUTE_UNUSED = tg413*xformL33;
    
    CCTK_REAL csetemp61 CCTK_ATTRIBUTE_UNUSED = tg402*xformL03;
    
    CCTK_REAL csetemp62 CCTK_ATTRIBUTE_UNUSED = tg412*xformL13;
    
    CCTK_REAL csetemp63 CCTK_ATTRIBUTE_UNUSED = tg422*xformL23;
    
    CCTK_REAL csetemp64 CCTK_ATTRIBUTE_UNUSED = tg423*xformL33;
    
    CCTK_REAL csetemp65 CCTK_ATTRIBUTE_UNUSED = tg403*xformL03;
    
    CCTK_REAL csetemp66 CCTK_ATTRIBUTE_UNUSED = tg413*xformL13;
    
    CCTK_REAL csetemp67 CCTK_ATTRIBUTE_UNUSED = tg423*xformL23;
    
    CCTK_REAL csetemp68 CCTK_ATTRIBUTE_UNUSED = tg433*xformL33;
    
    CCTK_REAL g403 CCTK_ATTRIBUTE_UNUSED = (csetemp53 + csetemp54 + 
      csetemp55 + csetemp56)*xformL00 + (csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*xformL10 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL20 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL30;
    
    CCTK_REAL g411 CCTK_ATTRIBUTE_UNUSED = (csetemp21 + csetemp22 + 
      csetemp23 + csetemp24)*xformL01 + (csetemp25 + csetemp26 + csetemp27 + 
      csetemp28)*xformL11 + (csetemp29 + csetemp30 + csetemp31 + 
      csetemp32)*xformL21 + (csetemp33 + csetemp34 + csetemp35 + 
      csetemp36)*xformL31;
    
    CCTK_REAL g412 CCTK_ATTRIBUTE_UNUSED = (csetemp37 + csetemp38 + 
      csetemp39 + csetemp40)*xformL01 + (csetemp41 + csetemp42 + csetemp43 + 
      csetemp44)*xformL11 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL21 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL31;
    
    CCTK_REAL g413 CCTK_ATTRIBUTE_UNUSED = (csetemp53 + csetemp54 + 
      csetemp55 + csetemp56)*xformL01 + (csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*xformL11 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL21 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL31;
    
    CCTK_REAL g422 CCTK_ATTRIBUTE_UNUSED = (csetemp37 + csetemp38 + 
      csetemp39 + csetemp40)*xformL02 + (csetemp41 + csetemp42 + csetemp43 + 
      csetemp44)*xformL12 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL22 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL32;
    
    CCTK_REAL g423 CCTK_ATTRIBUTE_UNUSED = (csetemp53 + csetemp54 + 
      csetemp55 + csetemp56)*xformL02 + (csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*xformL12 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL22 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL32;
    
    CCTK_REAL g433 CCTK_ATTRIBUTE_UNUSED = (csetemp53 + csetemp54 + 
      csetemp55 + csetemp56)*xformL03 + (csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*xformL13 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL23 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL33;
    
    CCTK_REAL csetemp69 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL00;
    
    CCTK_REAL csetemp70 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL10;
    
    CCTK_REAL csetemp71 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL20;
    
    CCTK_REAL csetemp72 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL30;
    
    CCTK_REAL csetemp73 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL00;
    
    CCTK_REAL csetemp74 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL10;
    
    CCTK_REAL csetemp75 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL20;
    
    CCTK_REAL csetemp76 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL30;
    
    CCTK_REAL csetemp77 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL00;
    
    CCTK_REAL csetemp78 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL10;
    
    CCTK_REAL csetemp79 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL20;
    
    CCTK_REAL csetemp80 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL30;
    
    CCTK_REAL csetemp81 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL00;
    
    CCTK_REAL csetemp82 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL10;
    
    CCTK_REAL csetemp83 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL20;
    
    CCTK_REAL csetemp84 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL30;
    
    CCTK_REAL csetemp85 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL00;
    
    CCTK_REAL csetemp86 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL10;
    
    CCTK_REAL csetemp87 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL20;
    
    CCTK_REAL csetemp88 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL30;
    
    CCTK_REAL csetemp89 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL00;
    
    CCTK_REAL csetemp90 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL10;
    
    CCTK_REAL csetemp91 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL20;
    
    CCTK_REAL csetemp92 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL30;
    
    CCTK_REAL csetemp93 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL00;
    
    CCTK_REAL csetemp94 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL10;
    
    CCTK_REAL csetemp95 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL20;
    
    CCTK_REAL csetemp96 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL30;
    
    CCTK_REAL csetemp97 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL00;
    
    CCTK_REAL csetemp98 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL10;
    
    CCTK_REAL csetemp99 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL20;
    
    CCTK_REAL csetemp100 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL30;
    
    CCTK_REAL csetemp101 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL00;
    
    CCTK_REAL csetemp102 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL10;
    
    CCTK_REAL csetemp103 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL20;
    
    CCTK_REAL csetemp104 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL30;
    
    CCTK_REAL csetemp105 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL00;
    
    CCTK_REAL csetemp106 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL10;
    
    CCTK_REAL csetemp107 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL20;
    
    CCTK_REAL csetemp108 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL30;
    
    CCTK_REAL dg4000 CCTK_ATTRIBUTE_UNUSED = xformL20*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL00 + (csetemp89 + csetemp90 + 
      csetemp91 + csetemp92)*xformL10 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL20 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL30) + xformL30*((csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL00 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL10 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL20 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL30) + xformL00*((csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL00 + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL10 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL20 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL30) + xformL10*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL00 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL10 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL20 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL30);
    
    CCTK_REAL dg4010 CCTK_ATTRIBUTE_UNUSED = xformL20*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL01 + (csetemp89 + csetemp90 + 
      csetemp91 + csetemp92)*xformL11 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL21 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL31) + xformL30*((csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL01 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL11 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL21 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL31) + xformL00*((csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL01 + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL11 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL21 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL31) + xformL10*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL01 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL11 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL21 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL31);
    
    CCTK_REAL csetemp109 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL01;
    
    CCTK_REAL csetemp110 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL11;
    
    CCTK_REAL csetemp111 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL21;
    
    CCTK_REAL csetemp112 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL31;
    
    CCTK_REAL csetemp113 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL01;
    
    CCTK_REAL csetemp114 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL11;
    
    CCTK_REAL csetemp115 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL21;
    
    CCTK_REAL csetemp116 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL31;
    
    CCTK_REAL csetemp117 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL01;
    
    CCTK_REAL csetemp118 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL11;
    
    CCTK_REAL csetemp119 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL21;
    
    CCTK_REAL csetemp120 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL31;
    
    CCTK_REAL csetemp121 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL01;
    
    CCTK_REAL csetemp122 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL11;
    
    CCTK_REAL csetemp123 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL21;
    
    CCTK_REAL csetemp124 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL31;
    
    CCTK_REAL csetemp125 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL01;
    
    CCTK_REAL csetemp126 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL11;
    
    CCTK_REAL csetemp127 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL21;
    
    CCTK_REAL csetemp128 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL31;
    
    CCTK_REAL csetemp129 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL01;
    
    CCTK_REAL csetemp130 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL11;
    
    CCTK_REAL csetemp131 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL21;
    
    CCTK_REAL csetemp132 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL31;
    
    CCTK_REAL csetemp133 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL01;
    
    CCTK_REAL csetemp134 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL11;
    
    CCTK_REAL csetemp135 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL21;
    
    CCTK_REAL csetemp136 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL31;
    
    CCTK_REAL csetemp137 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL01;
    
    CCTK_REAL csetemp138 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL11;
    
    CCTK_REAL csetemp139 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL21;
    
    CCTK_REAL csetemp140 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL31;
    
    CCTK_REAL csetemp141 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL01;
    
    CCTK_REAL csetemp142 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL11;
    
    CCTK_REAL csetemp143 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL21;
    
    CCTK_REAL csetemp144 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL31;
    
    CCTK_REAL csetemp145 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL01;
    
    CCTK_REAL csetemp146 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL11;
    
    CCTK_REAL csetemp147 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL21;
    
    CCTK_REAL csetemp148 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL31;
    
    CCTK_REAL dg4011 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp109 + 
      csetemp110 + csetemp111 + csetemp112)*xformL01 + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL11 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL21 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL31) + xformL10*((csetemp113 
      + csetemp114 + csetemp115 + csetemp116)*xformL01 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL11 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL21 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL31) + xformL20*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL01 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL11 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL21 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL31) + xformL30*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL01 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL11 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL21 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL31);
    
    CCTK_REAL csetemp149 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL02;
    
    CCTK_REAL csetemp150 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL12;
    
    CCTK_REAL csetemp151 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL22;
    
    CCTK_REAL csetemp152 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL32;
    
    CCTK_REAL csetemp153 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL02;
    
    CCTK_REAL csetemp154 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL12;
    
    CCTK_REAL csetemp155 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL22;
    
    CCTK_REAL csetemp156 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL32;
    
    CCTK_REAL csetemp157 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL02;
    
    CCTK_REAL csetemp158 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL12;
    
    CCTK_REAL csetemp159 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL22;
    
    CCTK_REAL csetemp160 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL32;
    
    CCTK_REAL csetemp161 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL02;
    
    CCTK_REAL csetemp162 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL12;
    
    CCTK_REAL csetemp163 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL22;
    
    CCTK_REAL csetemp164 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL32;
    
    CCTK_REAL csetemp165 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL02;
    
    CCTK_REAL csetemp166 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL12;
    
    CCTK_REAL csetemp167 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL22;
    
    CCTK_REAL csetemp168 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL32;
    
    CCTK_REAL csetemp169 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL02;
    
    CCTK_REAL csetemp170 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL12;
    
    CCTK_REAL csetemp171 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL22;
    
    CCTK_REAL csetemp172 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL32;
    
    CCTK_REAL csetemp173 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL02;
    
    CCTK_REAL csetemp174 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL12;
    
    CCTK_REAL csetemp175 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL22;
    
    CCTK_REAL csetemp176 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL32;
    
    CCTK_REAL csetemp177 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL02;
    
    CCTK_REAL csetemp178 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL12;
    
    CCTK_REAL csetemp179 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL22;
    
    CCTK_REAL csetemp180 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL32;
    
    CCTK_REAL csetemp181 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL02;
    
    CCTK_REAL csetemp182 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL12;
    
    CCTK_REAL csetemp183 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL22;
    
    CCTK_REAL csetemp184 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL32;
    
    CCTK_REAL csetemp185 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL02;
    
    CCTK_REAL csetemp186 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL12;
    
    CCTK_REAL csetemp187 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL22;
    
    CCTK_REAL csetemp188 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL32;
    
    CCTK_REAL dg4012 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL01 + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL11 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL21 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL31) + xformL10*((csetemp153 
      + csetemp154 + csetemp155 + csetemp156)*xformL01 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL11 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL21 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL31) + xformL20*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL01 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL11 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL21 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL31) + xformL30*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL01 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL11 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL21 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL31);
    
    CCTK_REAL csetemp189 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL03;
    
    CCTK_REAL csetemp190 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL13;
    
    CCTK_REAL csetemp191 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL23;
    
    CCTK_REAL csetemp192 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL33;
    
    CCTK_REAL csetemp193 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL03;
    
    CCTK_REAL csetemp194 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL13;
    
    CCTK_REAL csetemp195 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL23;
    
    CCTK_REAL csetemp196 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL33;
    
    CCTK_REAL csetemp197 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL03;
    
    CCTK_REAL csetemp198 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL13;
    
    CCTK_REAL csetemp199 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL23;
    
    CCTK_REAL csetemp200 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL33;
    
    CCTK_REAL csetemp201 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL03;
    
    CCTK_REAL csetemp202 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL13;
    
    CCTK_REAL csetemp203 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL23;
    
    CCTK_REAL csetemp204 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL33;
    
    CCTK_REAL csetemp205 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL03;
    
    CCTK_REAL csetemp206 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL13;
    
    CCTK_REAL csetemp207 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL23;
    
    CCTK_REAL csetemp208 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL33;
    
    CCTK_REAL csetemp209 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL03;
    
    CCTK_REAL csetemp210 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL13;
    
    CCTK_REAL csetemp211 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL23;
    
    CCTK_REAL csetemp212 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL33;
    
    CCTK_REAL csetemp213 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL03;
    
    CCTK_REAL csetemp214 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL13;
    
    CCTK_REAL csetemp215 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL23;
    
    CCTK_REAL csetemp216 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL33;
    
    CCTK_REAL csetemp217 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL03;
    
    CCTK_REAL csetemp218 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL13;
    
    CCTK_REAL csetemp219 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL23;
    
    CCTK_REAL csetemp220 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL33;
    
    CCTK_REAL csetemp221 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL03;
    
    CCTK_REAL csetemp222 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL13;
    
    CCTK_REAL csetemp223 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL23;
    
    CCTK_REAL csetemp224 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL33;
    
    CCTK_REAL csetemp225 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL03;
    
    CCTK_REAL csetemp226 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL13;
    
    CCTK_REAL csetemp227 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL23;
    
    CCTK_REAL csetemp228 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL33;
    
    CCTK_REAL dg4013 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL01 + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL11 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL21 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL31) + xformL10*((csetemp193 
      + csetemp194 + csetemp195 + csetemp196)*xformL01 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL11 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL21 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL31) + xformL20*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL01 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL11 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL21 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL31) + xformL30*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL01 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL11 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL21 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL31);
    
    CCTK_REAL dg4020 CCTK_ATTRIBUTE_UNUSED = xformL20*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL02 + (csetemp89 + csetemp90 + 
      csetemp91 + csetemp92)*xformL12 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL22 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL32) + xformL30*((csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL02 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL12 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL22 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL32) + xformL00*((csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL02 + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL12 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL22 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL32) + xformL10*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL02 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL12 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL22 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL32);
    
    CCTK_REAL dg4021 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp109 + 
      csetemp110 + csetemp111 + csetemp112)*xformL02 + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL12 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL22 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL32) + xformL10*((csetemp113 
      + csetemp114 + csetemp115 + csetemp116)*xformL02 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL12 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL22 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL32) + xformL20*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL02 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL12 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL22 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL32) + xformL30*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL02 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL12 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL22 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL32);
    
    CCTK_REAL dg4022 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL02 + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL12 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL22 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL32) + xformL10*((csetemp153 
      + csetemp154 + csetemp155 + csetemp156)*xformL02 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL12 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL22 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL32) + xformL20*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL02 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL12 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL22 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL32) + xformL30*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL02 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL12 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL22 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL32);
    
    CCTK_REAL dg4023 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL02 + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL12 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL22 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL32) + xformL10*((csetemp193 
      + csetemp194 + csetemp195 + csetemp196)*xformL02 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL12 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL22 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL32) + xformL20*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL02 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL12 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL22 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL32) + xformL30*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL02 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL12 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL22 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL32);
    
    CCTK_REAL dg4030 CCTK_ATTRIBUTE_UNUSED = xformL20*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL03 + (csetemp89 + csetemp90 + 
      csetemp91 + csetemp92)*xformL13 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL23 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL33) + xformL30*((csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL03 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL13 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL23 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL33) + xformL00*((csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL03 + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL13 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL23 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL33) + xformL10*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL03 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL13 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL23 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL33);
    
    CCTK_REAL dg4031 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp109 + 
      csetemp110 + csetemp111 + csetemp112)*xformL03 + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL13 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL23 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL33) + xformL10*((csetemp113 
      + csetemp114 + csetemp115 + csetemp116)*xformL03 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL13 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL23 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL33) + xformL20*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL03 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL13 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL23 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL33) + xformL30*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL03 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL13 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL23 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL33);
    
    CCTK_REAL dg4032 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL03 + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL13 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL23 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL33) + xformL10*((csetemp153 
      + csetemp154 + csetemp155 + csetemp156)*xformL03 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL13 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL23 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL33) + xformL20*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL03 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL13 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL23 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL33) + xformL30*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL03 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL13 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL23 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL33);
    
    CCTK_REAL dg4033 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL03 + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL13 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL23 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL33) + xformL10*((csetemp193 
      + csetemp194 + csetemp195 + csetemp196)*xformL03 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL13 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL23 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL33) + xformL20*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL03 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL13 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL23 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL33) + xformL30*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL03 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL13 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL23 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL33);
    
    CCTK_REAL dg4110 CCTK_ATTRIBUTE_UNUSED = xformL21*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL01 + (csetemp89 + csetemp90 + 
      csetemp91 + csetemp92)*xformL11 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL21 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL31) + xformL31*((csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL01 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL11 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL21 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL31) + xformL01*((csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL01 + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL11 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL21 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL31) + xformL11*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL01 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL11 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL21 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL31);
    
    CCTK_REAL dg4111 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp109 + 
      csetemp110 + csetemp111 + csetemp112)*xformL01 + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL11 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL21 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL31) + xformL11*((csetemp113 
      + csetemp114 + csetemp115 + csetemp116)*xformL01 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL11 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL21 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL31) + xformL21*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL01 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL11 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL21 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL31) + xformL31*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL01 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL11 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL21 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL31);
    
    CCTK_REAL dg4112 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL01 + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL11 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL21 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL31) + xformL11*((csetemp153 
      + csetemp154 + csetemp155 + csetemp156)*xformL01 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL11 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL21 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL31) + xformL21*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL01 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL11 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL21 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL31) + xformL31*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL01 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL11 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL21 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL31);
    
    CCTK_REAL dg4113 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL01 + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL11 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL21 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL31) + xformL11*((csetemp193 
      + csetemp194 + csetemp195 + csetemp196)*xformL01 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL11 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL21 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL31) + xformL21*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL01 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL11 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL21 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL31) + xformL31*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL01 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL11 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL21 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL31);
    
    CCTK_REAL dg4120 CCTK_ATTRIBUTE_UNUSED = xformL21*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL02 + (csetemp89 + csetemp90 + 
      csetemp91 + csetemp92)*xformL12 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL22 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL32) + xformL31*((csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL02 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL12 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL22 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL32) + xformL01*((csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL02 + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL12 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL22 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL32) + xformL11*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL02 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL12 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL22 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL32);
    
    CCTK_REAL dg4121 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp109 + 
      csetemp110 + csetemp111 + csetemp112)*xformL02 + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL12 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL22 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL32) + xformL11*((csetemp113 
      + csetemp114 + csetemp115 + csetemp116)*xformL02 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL12 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL22 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL32) + xformL21*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL02 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL12 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL22 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL32) + xformL31*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL02 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL12 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL22 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL32);
    
    CCTK_REAL dg4122 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL02 + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL12 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL22 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL32) + xformL11*((csetemp153 
      + csetemp154 + csetemp155 + csetemp156)*xformL02 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL12 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL22 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL32) + xformL21*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL02 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL12 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL22 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL32) + xformL31*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL02 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL12 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL22 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL32);
    
    CCTK_REAL dg4123 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL02 + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL12 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL22 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL32) + xformL11*((csetemp193 
      + csetemp194 + csetemp195 + csetemp196)*xformL02 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL12 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL22 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL32) + xformL21*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL02 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL12 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL22 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL32) + xformL31*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL02 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL12 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL22 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL32);
    
    CCTK_REAL dg4130 CCTK_ATTRIBUTE_UNUSED = xformL21*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL03 + (csetemp89 + csetemp90 + 
      csetemp91 + csetemp92)*xformL13 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL23 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL33) + xformL31*((csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL03 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL13 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL23 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL33) + xformL01*((csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL03 + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL13 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL23 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL33) + xformL11*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL03 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL13 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL23 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL33);
    
    CCTK_REAL dg4131 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp109 + 
      csetemp110 + csetemp111 + csetemp112)*xformL03 + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL13 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL23 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL33) + xformL11*((csetemp113 
      + csetemp114 + csetemp115 + csetemp116)*xformL03 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL13 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL23 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL33) + xformL21*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL03 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL13 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL23 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL33) + xformL31*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL03 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL13 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL23 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL33);
    
    CCTK_REAL dg4132 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL03 + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL13 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL23 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL33) + xformL11*((csetemp153 
      + csetemp154 + csetemp155 + csetemp156)*xformL03 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL13 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL23 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL33) + xformL21*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL03 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL13 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL23 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL33) + xformL31*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL03 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL13 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL23 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL33);
    
    CCTK_REAL dg4133 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL03 + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL13 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL23 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL33) + xformL11*((csetemp193 
      + csetemp194 + csetemp195 + csetemp196)*xformL03 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL13 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL23 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL33) + xformL21*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL03 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL13 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL23 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL33) + xformL31*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL03 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL13 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL23 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL33);
    
    CCTK_REAL dg4220 CCTK_ATTRIBUTE_UNUSED = xformL22*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL02 + (csetemp89 + csetemp90 + 
      csetemp91 + csetemp92)*xformL12 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL22 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL32) + xformL32*((csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL02 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL12 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL22 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL32) + xformL02*((csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL02 + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL12 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL22 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL32) + xformL12*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL02 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL12 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL22 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL32);
    
    CCTK_REAL dg4221 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp109 + 
      csetemp110 + csetemp111 + csetemp112)*xformL02 + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL12 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL22 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL32) + xformL12*((csetemp113 
      + csetemp114 + csetemp115 + csetemp116)*xformL02 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL12 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL22 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL32) + xformL22*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL02 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL12 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL22 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL32) + xformL32*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL02 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL12 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL22 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL32);
    
    CCTK_REAL dg4222 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL02 + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL12 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL22 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL32) + xformL12*((csetemp153 
      + csetemp154 + csetemp155 + csetemp156)*xformL02 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL12 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL22 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL32) + xformL22*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL02 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL12 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL22 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL32) + xformL32*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL02 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL12 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL22 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL32);
    
    CCTK_REAL dg4223 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL02 + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL12 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL22 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL32) + xformL12*((csetemp193 
      + csetemp194 + csetemp195 + csetemp196)*xformL02 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL12 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL22 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL32) + xformL22*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL02 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL12 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL22 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL32) + xformL32*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL02 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL12 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL22 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL32);
    
    CCTK_REAL dg4230 CCTK_ATTRIBUTE_UNUSED = xformL22*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL03 + (csetemp89 + csetemp90 + 
      csetemp91 + csetemp92)*xformL13 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL23 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL33) + xformL32*((csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL03 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL13 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL23 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL33) + xformL02*((csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL03 + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL13 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL23 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL33) + xformL12*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL03 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL13 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL23 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL33);
    
    CCTK_REAL dg4231 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp109 + 
      csetemp110 + csetemp111 + csetemp112)*xformL03 + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL13 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL23 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL33) + xformL12*((csetemp113 
      + csetemp114 + csetemp115 + csetemp116)*xformL03 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL13 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL23 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL33) + xformL22*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL03 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL13 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL23 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL33) + xformL32*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL03 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL13 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL23 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL33);
    
    CCTK_REAL dg4232 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL03 + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL13 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL23 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL33) + xformL12*((csetemp153 
      + csetemp154 + csetemp155 + csetemp156)*xformL03 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL13 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL23 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL33) + xformL22*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL03 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL13 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL23 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL33) + xformL32*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL03 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL13 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL23 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL33);
    
    CCTK_REAL dg4233 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL03 + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL13 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL23 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL33) + xformL12*((csetemp193 
      + csetemp194 + csetemp195 + csetemp196)*xformL03 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL13 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL23 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL33) + xformL22*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL03 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL13 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL23 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL33) + xformL32*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL03 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL13 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL23 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL33);
    
    CCTK_REAL dg4330 CCTK_ATTRIBUTE_UNUSED = xformL23*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL03 + (csetemp89 + csetemp90 + 
      csetemp91 + csetemp92)*xformL13 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL23 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL33) + xformL33*((csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL03 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL13 + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*xformL23 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL33) + xformL03*((csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL03 + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL13 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL23 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL33) + xformL13*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL03 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL13 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL23 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL33);
    
    CCTK_REAL dg4331 CCTK_ATTRIBUTE_UNUSED = xformL03*((csetemp109 + 
      csetemp110 + csetemp111 + csetemp112)*xformL03 + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL13 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL23 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL33) + xformL13*((csetemp113 
      + csetemp114 + csetemp115 + csetemp116)*xformL03 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL13 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL23 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL33) + xformL23*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL03 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL13 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL23 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL33) + xformL33*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL03 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL13 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL23 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL33);
    
    CCTK_REAL dg4332 CCTK_ATTRIBUTE_UNUSED = xformL03*((csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL03 + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL13 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL23 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL33) + xformL13*((csetemp153 
      + csetemp154 + csetemp155 + csetemp156)*xformL03 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL13 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL23 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL33) + xformL23*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL03 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL13 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL23 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL33) + xformL33*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL03 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL13 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL23 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL33);
    
    CCTK_REAL dg4333 CCTK_ATTRIBUTE_UNUSED = xformL03*((csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL03 + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL13 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL23 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL33) + xformL13*((csetemp193 
      + csetemp194 + csetemp195 + csetemp196)*xformL03 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL13 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL23 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL33) + xformL23*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL03 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL13 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL23 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL33) + xformL33*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL03 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL13 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL23 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL33);
    
    CCTK_REAL betal1 CCTK_ATTRIBUTE_UNUSED = g401;
    
    CCTK_REAL betal2 CCTK_ATTRIBUTE_UNUSED = g402;
    
    CCTK_REAL betal3 CCTK_ATTRIBUTE_UNUSED = g403;
    
    CCTK_REAL gerr11L CCTK_ATTRIBUTE_UNUSED = gxxL - g411;
    
    CCTK_REAL gerr12L CCTK_ATTRIBUTE_UNUSED = gxyL - g412;
    
    CCTK_REAL gerr13L CCTK_ATTRIBUTE_UNUSED = gxzL - g413;
    
    CCTK_REAL gerr22L CCTK_ATTRIBUTE_UNUSED = gyyL - g422;
    
    CCTK_REAL gerr23L CCTK_ATTRIBUTE_UNUSED = gyzL - g423;
    
    CCTK_REAL gerr33L CCTK_ATTRIBUTE_UNUSED = gzzL - g433;
    
    CCTK_REAL csetemp229 CCTK_ATTRIBUTE_UNUSED = SQR(gxzL);
    
    CCTK_REAL csetemp230 CCTK_ATTRIBUTE_UNUSED = SQR(gyzL);
    
    CCTK_REAL csetemp231 CCTK_ATTRIBUTE_UNUSED = SQR(gxyL);
    
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = 2*gxyL*gxzL*gyzL + 
      gyyL*(gxxL*gzzL - csetemp229) - gxxL*csetemp230 - gzzL*csetemp231;
    
    CCTK_REAL csetemp232 CCTK_ATTRIBUTE_UNUSED = INV(detg);
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = (gyyL*gzzL - 
      csetemp230)*csetemp232;
    
    CCTK_REAL gu12 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyzL - 
      gxyL*gzzL)*csetemp232;
    
    CCTK_REAL gu13 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyyL) + 
      gxyL*gyzL)*csetemp232;
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = (gxxL*gzzL - 
      csetemp229)*csetemp232;
    
    CCTK_REAL gu23 CCTK_ATTRIBUTE_UNUSED = (gxyL*gxzL - 
      gxxL*gyzL)*csetemp232;
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = (gxxL*gyyL - 
      csetemp231)*csetemp232;
    
    CCTK_REAL betaerr1L CCTK_ATTRIBUTE_UNUSED = betaxL - betal1*gu11 - 
      betal2*gu12 - betal3*gu13;
    
    CCTK_REAL betaerr2L CCTK_ATTRIBUTE_UNUSED = betayL - betal1*gu12 - 
      betal2*gu22 - betal3*gu23;
    
    CCTK_REAL betaerr3L CCTK_ATTRIBUTE_UNUSED = betazL - betal1*gu13 - 
      betal2*gu23 - betal3*gu33;
    
    CCTK_REAL betasq CCTK_ATTRIBUTE_UNUSED = betaxL*betal1 + betayL*betal2 
      + betazL*betal3;
    
    CCTK_REAL alperrL CCTK_ATTRIBUTE_UNUSED = alpL - sqrt(betasq - g400);
    
    CCTK_REAL dtg11 CCTK_ATTRIBUTE_UNUSED = dg4110;
    
    CCTK_REAL dtg12 CCTK_ATTRIBUTE_UNUSED = dg4120;
    
    CCTK_REAL dtg13 CCTK_ATTRIBUTE_UNUSED = dg4130;
    
    CCTK_REAL dtg22 CCTK_ATTRIBUTE_UNUSED = dg4220;
    
    CCTK_REAL dtg23 CCTK_ATTRIBUTE_UNUSED = dg4230;
    
    CCTK_REAL dtg33 CCTK_ATTRIBUTE_UNUSED = dg4330;
    
    CCTK_REAL dg111 CCTK_ATTRIBUTE_UNUSED = dg4111;
    
    CCTK_REAL dg112 CCTK_ATTRIBUTE_UNUSED = dg4112;
    
    CCTK_REAL dg113 CCTK_ATTRIBUTE_UNUSED = dg4113;
    
    CCTK_REAL dg121 CCTK_ATTRIBUTE_UNUSED = dg4121;
    
    CCTK_REAL dg122 CCTK_ATTRIBUTE_UNUSED = dg4122;
    
    CCTK_REAL dg123 CCTK_ATTRIBUTE_UNUSED = dg4123;
    
    CCTK_REAL dg131 CCTK_ATTRIBUTE_UNUSED = dg4131;
    
    CCTK_REAL dg132 CCTK_ATTRIBUTE_UNUSED = dg4132;
    
    CCTK_REAL dg133 CCTK_ATTRIBUTE_UNUSED = dg4133;
    
    CCTK_REAL dg221 CCTK_ATTRIBUTE_UNUSED = dg4221;
    
    CCTK_REAL dg222 CCTK_ATTRIBUTE_UNUSED = dg4222;
    
    CCTK_REAL dg223 CCTK_ATTRIBUTE_UNUSED = dg4223;
    
    CCTK_REAL dg231 CCTK_ATTRIBUTE_UNUSED = dg4231;
    
    CCTK_REAL dg232 CCTK_ATTRIBUTE_UNUSED = dg4232;
    
    CCTK_REAL dg233 CCTK_ATTRIBUTE_UNUSED = dg4233;
    
    CCTK_REAL dg331 CCTK_ATTRIBUTE_UNUSED = dg4331;
    
    CCTK_REAL dg332 CCTK_ATTRIBUTE_UNUSED = dg4332;
    
    CCTK_REAL dg333 CCTK_ATTRIBUTE_UNUSED = dg4333;
    
    CCTK_REAL csetemp233 CCTK_ATTRIBUTE_UNUSED = dtg11*gu11;
    
    CCTK_REAL csetemp234 CCTK_ATTRIBUTE_UNUSED = dtg12*gu12;
    
    CCTK_REAL csetemp235 CCTK_ATTRIBUTE_UNUSED = dtg13*gu13;
    
    CCTK_REAL csetemp236 CCTK_ATTRIBUTE_UNUSED = dtg12*gu11;
    
    CCTK_REAL csetemp237 CCTK_ATTRIBUTE_UNUSED = dtg22*gu12;
    
    CCTK_REAL csetemp238 CCTK_ATTRIBUTE_UNUSED = dtg23*gu13;
    
    CCTK_REAL csetemp239 CCTK_ATTRIBUTE_UNUSED = dtg13*gu11;
    
    CCTK_REAL csetemp240 CCTK_ATTRIBUTE_UNUSED = dtg23*gu12;
    
    CCTK_REAL csetemp241 CCTK_ATTRIBUTE_UNUSED = dtg33*gu13;
    
    CCTK_REAL dtgu11 CCTK_ATTRIBUTE_UNUSED = -((csetemp233 + csetemp234 + 
      csetemp235)*gu11) - (csetemp236 + csetemp237 + csetemp238)*gu12 - 
      (csetemp239 + csetemp240 + csetemp241)*gu13;
    
    CCTK_REAL dtgu12 CCTK_ATTRIBUTE_UNUSED = -((csetemp233 + csetemp234 + 
      csetemp235)*gu12) - (csetemp236 + csetemp237 + csetemp238)*gu22 - 
      (csetemp239 + csetemp240 + csetemp241)*gu23;
    
    CCTK_REAL dtgu13 CCTK_ATTRIBUTE_UNUSED = -((csetemp233 + csetemp234 + 
      csetemp235)*gu13) - (csetemp236 + csetemp237 + csetemp238)*gu23 - 
      (csetemp239 + csetemp240 + csetemp241)*gu33;
    
    CCTK_REAL csetemp242 CCTK_ATTRIBUTE_UNUSED = dtg11*gu12;
    
    CCTK_REAL csetemp243 CCTK_ATTRIBUTE_UNUSED = dtg12*gu22;
    
    CCTK_REAL csetemp244 CCTK_ATTRIBUTE_UNUSED = dtg13*gu23;
    
    CCTK_REAL csetemp245 CCTK_ATTRIBUTE_UNUSED = dtg22*gu22;
    
    CCTK_REAL csetemp246 CCTK_ATTRIBUTE_UNUSED = dtg23*gu23;
    
    CCTK_REAL csetemp247 CCTK_ATTRIBUTE_UNUSED = dtg13*gu12;
    
    CCTK_REAL csetemp248 CCTK_ATTRIBUTE_UNUSED = dtg23*gu22;
    
    CCTK_REAL csetemp249 CCTK_ATTRIBUTE_UNUSED = dtg33*gu23;
    
    CCTK_REAL dtgu22 CCTK_ATTRIBUTE_UNUSED = -((csetemp242 + csetemp243 + 
      csetemp244)*gu12) - (csetemp234 + csetemp245 + csetemp246)*gu22 - 
      (csetemp247 + csetemp248 + csetemp249)*gu23;
    
    CCTK_REAL dtgu23 CCTK_ATTRIBUTE_UNUSED = -((csetemp242 + csetemp243 + 
      csetemp244)*gu13) - (csetemp234 + csetemp245 + csetemp246)*gu23 - 
      (csetemp247 + csetemp248 + csetemp249)*gu33;
    
    CCTK_REAL dtgu33 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dtg11*gu13 + 
      dtg12*gu23 + dtg13*gu33)) - gu23*(dtg12*gu13 + dtg22*gu23 + dtg23*gu33) 
      - gu33*(csetemp235 + csetemp246 + dtg33*gu33);
    
    CCTK_REAL csetemp250 CCTK_ATTRIBUTE_UNUSED = dg111*gu11;
    
    CCTK_REAL csetemp251 CCTK_ATTRIBUTE_UNUSED = dg121*gu12;
    
    CCTK_REAL csetemp252 CCTK_ATTRIBUTE_UNUSED = dg131*gu13;
    
    CCTK_REAL csetemp253 CCTK_ATTRIBUTE_UNUSED = dg121*gu11;
    
    CCTK_REAL csetemp254 CCTK_ATTRIBUTE_UNUSED = dg221*gu12;
    
    CCTK_REAL csetemp255 CCTK_ATTRIBUTE_UNUSED = dg231*gu13;
    
    CCTK_REAL csetemp256 CCTK_ATTRIBUTE_UNUSED = dg131*gu11;
    
    CCTK_REAL csetemp257 CCTK_ATTRIBUTE_UNUSED = dg231*gu12;
    
    CCTK_REAL csetemp258 CCTK_ATTRIBUTE_UNUSED = dg331*gu13;
    
    CCTK_REAL dgu111 CCTK_ATTRIBUTE_UNUSED = -((csetemp250 + csetemp251 + 
      csetemp252)*gu11) - (csetemp253 + csetemp254 + csetemp255)*gu12 - 
      (csetemp256 + csetemp257 + csetemp258)*gu13;
    
    CCTK_REAL dgu121 CCTK_ATTRIBUTE_UNUSED = -((csetemp250 + csetemp251 + 
      csetemp252)*gu12) - (csetemp253 + csetemp254 + csetemp255)*gu22 - 
      (csetemp256 + csetemp257 + csetemp258)*gu23;
    
    CCTK_REAL dgu131 CCTK_ATTRIBUTE_UNUSED = -((csetemp250 + csetemp251 + 
      csetemp252)*gu13) - (csetemp253 + csetemp254 + csetemp255)*gu23 - 
      (csetemp256 + csetemp257 + csetemp258)*gu33;
    
    CCTK_REAL csetemp259 CCTK_ATTRIBUTE_UNUSED = dg111*gu12;
    
    CCTK_REAL csetemp260 CCTK_ATTRIBUTE_UNUSED = dg121*gu22;
    
    CCTK_REAL csetemp261 CCTK_ATTRIBUTE_UNUSED = dg131*gu23;
    
    CCTK_REAL csetemp262 CCTK_ATTRIBUTE_UNUSED = dg221*gu22;
    
    CCTK_REAL csetemp263 CCTK_ATTRIBUTE_UNUSED = dg231*gu23;
    
    CCTK_REAL csetemp264 CCTK_ATTRIBUTE_UNUSED = dg131*gu12;
    
    CCTK_REAL csetemp265 CCTK_ATTRIBUTE_UNUSED = dg231*gu22;
    
    CCTK_REAL csetemp266 CCTK_ATTRIBUTE_UNUSED = dg331*gu23;
    
    CCTK_REAL dgu221 CCTK_ATTRIBUTE_UNUSED = -((csetemp259 + csetemp260 + 
      csetemp261)*gu12) - (csetemp251 + csetemp262 + csetemp263)*gu22 - 
      (csetemp264 + csetemp265 + csetemp266)*gu23;
    
    CCTK_REAL dgu231 CCTK_ATTRIBUTE_UNUSED = -((csetemp259 + csetemp260 + 
      csetemp261)*gu13) - (csetemp251 + csetemp262 + csetemp263)*gu23 - 
      (csetemp264 + csetemp265 + csetemp266)*gu33;
    
    CCTK_REAL dgu331 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dg111*gu13 + 
      dg121*gu23 + dg131*gu33)) - gu23*(dg121*gu13 + dg221*gu23 + dg231*gu33) 
      - gu33*(csetemp252 + csetemp263 + dg331*gu33);
    
    CCTK_REAL csetemp267 CCTK_ATTRIBUTE_UNUSED = dg112*gu11;
    
    CCTK_REAL csetemp268 CCTK_ATTRIBUTE_UNUSED = dg122*gu12;
    
    CCTK_REAL csetemp269 CCTK_ATTRIBUTE_UNUSED = dg132*gu13;
    
    CCTK_REAL csetemp270 CCTK_ATTRIBUTE_UNUSED = dg122*gu11;
    
    CCTK_REAL csetemp271 CCTK_ATTRIBUTE_UNUSED = dg222*gu12;
    
    CCTK_REAL csetemp272 CCTK_ATTRIBUTE_UNUSED = dg232*gu13;
    
    CCTK_REAL csetemp273 CCTK_ATTRIBUTE_UNUSED = dg132*gu11;
    
    CCTK_REAL csetemp274 CCTK_ATTRIBUTE_UNUSED = dg232*gu12;
    
    CCTK_REAL csetemp275 CCTK_ATTRIBUTE_UNUSED = dg332*gu13;
    
    CCTK_REAL dgu112 CCTK_ATTRIBUTE_UNUSED = -((csetemp267 + csetemp268 + 
      csetemp269)*gu11) - (csetemp270 + csetemp271 + csetemp272)*gu12 - 
      (csetemp273 + csetemp274 + csetemp275)*gu13;
    
    CCTK_REAL dgu122 CCTK_ATTRIBUTE_UNUSED = -((csetemp267 + csetemp268 + 
      csetemp269)*gu12) - (csetemp270 + csetemp271 + csetemp272)*gu22 - 
      (csetemp273 + csetemp274 + csetemp275)*gu23;
    
    CCTK_REAL dgu132 CCTK_ATTRIBUTE_UNUSED = -((csetemp267 + csetemp268 + 
      csetemp269)*gu13) - (csetemp270 + csetemp271 + csetemp272)*gu23 - 
      (csetemp273 + csetemp274 + csetemp275)*gu33;
    
    CCTK_REAL csetemp276 CCTK_ATTRIBUTE_UNUSED = dg112*gu12;
    
    CCTK_REAL csetemp277 CCTK_ATTRIBUTE_UNUSED = dg122*gu22;
    
    CCTK_REAL csetemp278 CCTK_ATTRIBUTE_UNUSED = dg132*gu23;
    
    CCTK_REAL csetemp279 CCTK_ATTRIBUTE_UNUSED = dg222*gu22;
    
    CCTK_REAL csetemp280 CCTK_ATTRIBUTE_UNUSED = dg232*gu23;
    
    CCTK_REAL csetemp281 CCTK_ATTRIBUTE_UNUSED = dg132*gu12;
    
    CCTK_REAL csetemp282 CCTK_ATTRIBUTE_UNUSED = dg232*gu22;
    
    CCTK_REAL csetemp283 CCTK_ATTRIBUTE_UNUSED = dg332*gu23;
    
    CCTK_REAL dgu222 CCTK_ATTRIBUTE_UNUSED = -((csetemp276 + csetemp277 + 
      csetemp278)*gu12) - (csetemp268 + csetemp279 + csetemp280)*gu22 - 
      (csetemp281 + csetemp282 + csetemp283)*gu23;
    
    CCTK_REAL dgu232 CCTK_ATTRIBUTE_UNUSED = -((csetemp276 + csetemp277 + 
      csetemp278)*gu13) - (csetemp268 + csetemp279 + csetemp280)*gu23 - 
      (csetemp281 + csetemp282 + csetemp283)*gu33;
    
    CCTK_REAL dgu332 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dg112*gu13 + 
      dg122*gu23 + dg132*gu33)) - gu23*(dg122*gu13 + dg222*gu23 + dg232*gu33) 
      - gu33*(csetemp269 + csetemp280 + dg332*gu33);
    
    CCTK_REAL csetemp284 CCTK_ATTRIBUTE_UNUSED = dg113*gu11;
    
    CCTK_REAL csetemp285 CCTK_ATTRIBUTE_UNUSED = dg123*gu12;
    
    CCTK_REAL csetemp286 CCTK_ATTRIBUTE_UNUSED = dg133*gu13;
    
    CCTK_REAL csetemp287 CCTK_ATTRIBUTE_UNUSED = dg123*gu11;
    
    CCTK_REAL csetemp288 CCTK_ATTRIBUTE_UNUSED = dg223*gu12;
    
    CCTK_REAL csetemp289 CCTK_ATTRIBUTE_UNUSED = dg233*gu13;
    
    CCTK_REAL csetemp290 CCTK_ATTRIBUTE_UNUSED = dg133*gu11;
    
    CCTK_REAL csetemp291 CCTK_ATTRIBUTE_UNUSED = dg233*gu12;
    
    CCTK_REAL csetemp292 CCTK_ATTRIBUTE_UNUSED = dg333*gu13;
    
    CCTK_REAL dgu113 CCTK_ATTRIBUTE_UNUSED = -((csetemp284 + csetemp285 + 
      csetemp286)*gu11) - (csetemp287 + csetemp288 + csetemp289)*gu12 - 
      (csetemp290 + csetemp291 + csetemp292)*gu13;
    
    CCTK_REAL dgu123 CCTK_ATTRIBUTE_UNUSED = -((csetemp284 + csetemp285 + 
      csetemp286)*gu12) - (csetemp287 + csetemp288 + csetemp289)*gu22 - 
      (csetemp290 + csetemp291 + csetemp292)*gu23;
    
    CCTK_REAL dgu133 CCTK_ATTRIBUTE_UNUSED = -((csetemp284 + csetemp285 + 
      csetemp286)*gu13) - (csetemp287 + csetemp288 + csetemp289)*gu23 - 
      (csetemp290 + csetemp291 + csetemp292)*gu33;
    
    CCTK_REAL csetemp293 CCTK_ATTRIBUTE_UNUSED = dg113*gu12;
    
    CCTK_REAL csetemp294 CCTK_ATTRIBUTE_UNUSED = dg123*gu22;
    
    CCTK_REAL csetemp295 CCTK_ATTRIBUTE_UNUSED = dg133*gu23;
    
    CCTK_REAL csetemp296 CCTK_ATTRIBUTE_UNUSED = dg223*gu22;
    
    CCTK_REAL csetemp297 CCTK_ATTRIBUTE_UNUSED = dg233*gu23;
    
    CCTK_REAL csetemp298 CCTK_ATTRIBUTE_UNUSED = dg133*gu12;
    
    CCTK_REAL csetemp299 CCTK_ATTRIBUTE_UNUSED = dg233*gu22;
    
    CCTK_REAL csetemp300 CCTK_ATTRIBUTE_UNUSED = dg333*gu23;
    
    CCTK_REAL dgu223 CCTK_ATTRIBUTE_UNUSED = -((csetemp293 + csetemp294 + 
      csetemp295)*gu12) - (csetemp285 + csetemp296 + csetemp297)*gu22 - 
      (csetemp298 + csetemp299 + csetemp300)*gu23;
    
    CCTK_REAL dgu233 CCTK_ATTRIBUTE_UNUSED = -((csetemp293 + csetemp294 + 
      csetemp295)*gu13) - (csetemp285 + csetemp296 + csetemp297)*gu23 - 
      (csetemp298 + csetemp299 + csetemp300)*gu33;
    
    CCTK_REAL dgu333 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dg113*gu13 + 
      dg123*gu23 + dg133*gu33)) - gu23*(dg123*gu13 + dg223*gu23 + dg233*gu33) 
      - gu33*(csetemp286 + csetemp297 + dg333*gu33);
    
    CCTK_REAL dtbetal1 CCTK_ATTRIBUTE_UNUSED = dg4010;
    
    CCTK_REAL dtbetal2 CCTK_ATTRIBUTE_UNUSED = dg4020;
    
    CCTK_REAL dtbetal3 CCTK_ATTRIBUTE_UNUSED = dg4030;
    
    CCTK_REAL dbetal11 CCTK_ATTRIBUTE_UNUSED = dg4011;
    
    CCTK_REAL dbetal12 CCTK_ATTRIBUTE_UNUSED = dg4012;
    
    CCTK_REAL dbetal13 CCTK_ATTRIBUTE_UNUSED = dg4013;
    
    CCTK_REAL dbetal21 CCTK_ATTRIBUTE_UNUSED = dg4021;
    
    CCTK_REAL dbetal22 CCTK_ATTRIBUTE_UNUSED = dg4022;
    
    CCTK_REAL dbetal23 CCTK_ATTRIBUTE_UNUSED = dg4023;
    
    CCTK_REAL dbetal31 CCTK_ATTRIBUTE_UNUSED = dg4031;
    
    CCTK_REAL dbetal32 CCTK_ATTRIBUTE_UNUSED = dg4032;
    
    CCTK_REAL dbetal33 CCTK_ATTRIBUTE_UNUSED = dg4033;
    
    CCTK_REAL dtbetaerr1L CCTK_ATTRIBUTE_UNUSED = dtbetaxL - betal1*dtgu11 
      - betal2*dtgu12 - betal3*dtgu13 - dtbetal1*gu11 - dtbetal2*gu12 - 
      dtbetal3*gu13;
    
    CCTK_REAL dtbetaerr2L CCTK_ATTRIBUTE_UNUSED = dtbetayL - betal1*dtgu12 
      - betal2*dtgu22 - betal3*dtgu23 - dtbetal1*gu12 - dtbetal2*gu22 - 
      dtbetal3*gu23;
    
    CCTK_REAL dtbetaerr3L CCTK_ATTRIBUTE_UNUSED = dtbetazL - betal1*dtgu13 
      - betal2*dtgu23 - betal3*dtgu33 - dtbetal1*gu13 - dtbetal2*gu23 - 
      dtbetal3*gu33;
    
    CCTK_REAL dbeta11 CCTK_ATTRIBUTE_UNUSED = betal1*dgu111 + 
      betal2*dgu121 + betal3*dgu131 + dbetal11*gu11 + dbetal21*gu12 + 
      dbetal31*gu13;
    
    CCTK_REAL dbeta21 CCTK_ATTRIBUTE_UNUSED = betal1*dgu121 + 
      betal2*dgu221 + betal3*dgu231 + dbetal11*gu12 + dbetal21*gu22 + 
      dbetal31*gu23;
    
    CCTK_REAL dbeta31 CCTK_ATTRIBUTE_UNUSED = betal1*dgu131 + 
      betal2*dgu231 + betal3*dgu331 + dbetal11*gu13 + dbetal21*gu23 + 
      dbetal31*gu33;
    
    CCTK_REAL dbeta12 CCTK_ATTRIBUTE_UNUSED = betal1*dgu112 + 
      betal2*dgu122 + betal3*dgu132 + dbetal12*gu11 + dbetal22*gu12 + 
      dbetal32*gu13;
    
    CCTK_REAL dbeta22 CCTK_ATTRIBUTE_UNUSED = betal1*dgu122 + 
      betal2*dgu222 + betal3*dgu232 + dbetal12*gu12 + dbetal22*gu22 + 
      dbetal32*gu23;
    
    CCTK_REAL dbeta32 CCTK_ATTRIBUTE_UNUSED = betal1*dgu132 + 
      betal2*dgu232 + betal3*dgu332 + dbetal12*gu13 + dbetal22*gu23 + 
      dbetal32*gu33;
    
    CCTK_REAL dbeta13 CCTK_ATTRIBUTE_UNUSED = betal1*dgu113 + 
      betal2*dgu123 + betal3*dgu133 + dbetal13*gu11 + dbetal23*gu12 + 
      dbetal33*gu13;
    
    CCTK_REAL dbeta23 CCTK_ATTRIBUTE_UNUSED = betal1*dgu123 + 
      betal2*dgu223 + betal3*dgu233 + dbetal13*gu12 + dbetal23*gu22 + 
      dbetal33*gu23;
    
    CCTK_REAL dbeta33 CCTK_ATTRIBUTE_UNUSED = betal1*dgu133 + 
      betal2*dgu233 + betal3*dgu333 + dbetal13*gu13 + dbetal23*gu23 + 
      dbetal33*gu33;
    
    CCTK_REAL dtbetasq CCTK_ATTRIBUTE_UNUSED = dtbetaxL*betal1 + 
      dtbetayL*betal2 + dtbetazL*betal3 + betaxL*dtbetal1 + betayL*dtbetal2 + 
      betazL*dtbetal3;
    
    CCTK_REAL csetemp301 CCTK_ATTRIBUTE_UNUSED = INV(alpL);
    
    CCTK_REAL dtalperrL CCTK_ATTRIBUTE_UNUSED = dtalpL + 
      0.5*csetemp301*(dg4000 - dtbetasq);
    
    CCTK_REAL kerr11L CCTK_ATTRIBUTE_UNUSED = kxxL - 
      0.5*csetemp301*(2*(gxxL*dbeta11 + gxyL*dbeta21 + gxzL*dbeta31) + 
      betaxL*dg111 + betayL*dg112 + betazL*dg113 - dtg11);
    
    CCTK_REAL kerr12L CCTK_ATTRIBUTE_UNUSED = kxyL - 
      0.5*csetemp301*(gxxL*dbeta12 + gyyL*dbeta21 + gxyL*(dbeta11 + dbeta22) 
      + gyzL*dbeta31 + gxzL*dbeta32 + betaxL*dg121 + betayL*dg122 + 
      betazL*dg123 - dtg12);
    
    CCTK_REAL kerr13L CCTK_ATTRIBUTE_UNUSED = kxzL - 
      0.5*csetemp301*(gxxL*dbeta13 + gyzL*dbeta21 + gxyL*dbeta23 + 
      gzzL*dbeta31 + gxzL*(dbeta11 + dbeta33) + betaxL*dg131 + betayL*dg132 + 
      betazL*dg133 - dtg13);
    
    CCTK_REAL kerr22L CCTK_ATTRIBUTE_UNUSED = kyyL - 
      0.5*csetemp301*(2*(gxyL*dbeta12 + gyyL*dbeta22 + gyzL*dbeta32) + 
      betaxL*dg221 + betayL*dg222 + betazL*dg223 - dtg22);
    
    CCTK_REAL kerr23L CCTK_ATTRIBUTE_UNUSED = kyzL - 
      0.5*csetemp301*(gxzL*dbeta12 + gxyL*dbeta13 + gyyL*dbeta23 + 
      gzzL*dbeta32 + gyzL*(dbeta22 + dbeta33) + betaxL*dg231 + betayL*dg232 + 
      betazL*dg233 - dtg23);
    
    CCTK_REAL kerr33L CCTK_ATTRIBUTE_UNUSED = kzzL - 
      0.5*csetemp301*(2*(gxzL*dbeta13 + gyzL*dbeta23 + gzzL*dbeta33) + 
      betaxL*dg331 + betayL*dg332 + betazL*dg333 - dtg33);
    /* Copy local copies back to grid functions */
    alperr[index] = alperrL;
    betaerr1[index] = betaerr1L;
    betaerr2[index] = betaerr2L;
    betaerr3[index] = betaerr3L;
    dtalperr[index] = dtalperrL;
    dtbetaerr1[index] = dtbetaerr1L;
    dtbetaerr2[index] = dtbetaerr2L;
    dtbetaerr3[index] = dtbetaerr3L;
    gerr11[index] = gerr11L;
    gerr12[index] = gerr12L;
    gerr13[index] = gerr13L;
    gerr22[index] = gerr22L;
    gerr23[index] = gerr23L;
    gerr33[index] = gerr33L;
    kerr11[index] = kerr11L;
    kerr12[index] = kerr12L;
    kerr13[index] = kerr13L;
    kerr22[index] = kerr22L;
    kerr23[index] = kerr23L;
    kerr33[index] = kerr33L;
  }
  CCTK_ENDLOOP3(ModifiedSchwarzschildBL_error);
}
extern "C" void ModifiedSchwarzschildBL_error(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ModifiedSchwarzschildBL_error_Body");
  }
  if (cctk_iteration % ModifiedSchwarzschildBL_error_calc_every != ModifiedSchwarzschildBL_error_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "admbase::curv",
    "admbase::dtlapse",
    "admbase::dtshift",
    "admbase::lapse",
    "admbase::metric",
    "admbase::shift",
    "ModifiedSchwarzschildBL::curv_error",
    "ModifiedSchwarzschildBL::dtlapse_error",
    "ModifiedSchwarzschildBL::dtshift_error",
    "grid::coordinates",
    "ModifiedSchwarzschildBL::lapse_error",
    "ModifiedSchwarzschildBL::metric_error",
    "ModifiedSchwarzschildBL::shift_error"};
  AssertGroupStorage(cctkGH, "ModifiedSchwarzschildBL_error", 13, groups);
  
  
  LoopOverEverything(cctkGH, ModifiedSchwarzschildBL_error_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ModifiedSchwarzschildBL_error_Body");
  }
}

} // namespace ModifiedSchwarzschildBL
