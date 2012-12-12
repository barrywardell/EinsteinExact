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

static void KerrSchild_initial_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Include user-supplied include files */
  
  /* Initialise finite differencing variables */
  ptrdiff_t const di CCTK_ATTRIBUTE_UNUSED  = 1;
  ptrdiff_t const dj CCTK_ATTRIBUTE_UNUSED  = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk CCTK_ATTRIBUTE_UNUSED  = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const cdi CCTK_ATTRIBUTE_UNUSED  = sizeof(CCTK_REAL) * di;
  ptrdiff_t const cdj CCTK_ATTRIBUTE_UNUSED  = sizeof(CCTK_REAL) * dj;
  ptrdiff_t const cdk CCTK_ATTRIBUTE_UNUSED  = sizeof(CCTK_REAL) * dk;
  CCTK_REAL const dx CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL const dy CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL const dz CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL const dt CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL const t CCTK_ATTRIBUTE_UNUSED  = ToReal(cctk_time);
  CCTK_REAL const dxi CCTK_ATTRIBUTE_UNUSED  = INV(dx);
  CCTK_REAL const dyi CCTK_ATTRIBUTE_UNUSED  = INV(dy);
  CCTK_REAL const dzi CCTK_ATTRIBUTE_UNUSED  = INV(dz);
  CCTK_REAL const khalf CCTK_ATTRIBUTE_UNUSED  = 0.5;
  CCTK_REAL const kthird CCTK_ATTRIBUTE_UNUSED  = 1/3.0;
  CCTK_REAL const ktwothird CCTK_ATTRIBUTE_UNUSED  = 2.0/3.0;
  CCTK_REAL const kfourthird CCTK_ATTRIBUTE_UNUSED  = 4.0/3.0;
  CCTK_REAL const keightthird CCTK_ATTRIBUTE_UNUSED  = 8.0/3.0;
  CCTK_REAL const hdxi CCTK_ATTRIBUTE_UNUSED  = 0.5 * dxi;
  CCTK_REAL const hdyi CCTK_ATTRIBUTE_UNUSED  = 0.5 * dyi;
  CCTK_REAL const hdzi CCTK_ATTRIBUTE_UNUSED  = 0.5 * dzi;
  
  /* Initialize predefined quantities */
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  CCTK_LOOP3(KerrSchild_initial,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    ptrdiff_t const index CCTK_ATTRIBUTE_UNUSED  = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL CCTK_ATTRIBUTE_UNUSED = alp[index];
    CCTK_REAL betaxL CCTK_ATTRIBUTE_UNUSED = betax[index];
    CCTK_REAL betayL CCTK_ATTRIBUTE_UNUSED = betay[index];
    CCTK_REAL betazL CCTK_ATTRIBUTE_UNUSED = betaz[index];
    CCTK_REAL dtbetaxL CCTK_ATTRIBUTE_UNUSED = dtbetax[index];
    CCTK_REAL dtbetayL CCTK_ATTRIBUTE_UNUSED = dtbetay[index];
    CCTK_REAL dtbetazL CCTK_ATTRIBUTE_UNUSED = dtbetaz[index];
    CCTK_REAL gxxL CCTK_ATTRIBUTE_UNUSED = gxx[index];
    CCTK_REAL gxyL CCTK_ATTRIBUTE_UNUSED = gxy[index];
    CCTK_REAL gxzL CCTK_ATTRIBUTE_UNUSED = gxz[index];
    CCTK_REAL gyyL CCTK_ATTRIBUTE_UNUSED = gyy[index];
    CCTK_REAL gyzL CCTK_ATTRIBUTE_UNUSED = gyz[index];
    CCTK_REAL gzzL CCTK_ATTRIBUTE_UNUSED = gzz[index];
    CCTK_REAL xL CCTK_ATTRIBUTE_UNUSED = x[index];
    CCTK_REAL yL CCTK_ATTRIBUTE_UNUSED = y[index];
    CCTK_REAL zL CCTK_ATTRIBUTE_UNUSED = z[index];
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L00 = 1;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L01 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L02 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L03 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L10 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp0 = cos(ToReal(phi));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp1 = cos(ToReal(psi));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp2 = cos(ToReal(theta));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp3 = sin(ToReal(phi));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp4 = sin(ToReal(psi));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L11 = csetemp0*csetemp1 - 
      csetemp2*csetemp3*csetemp4;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L12 = csetemp1*csetemp3 + 
      csetemp0*csetemp2*csetemp4;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp5 = sin(ToReal(theta));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L13 = csetemp4*csetemp5;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L20 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L21 = 
      -(csetemp1*csetemp2*csetemp3) - csetemp0*csetemp4;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L22 = csetemp0*csetemp1*csetemp2 
      - csetemp3*csetemp4;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L23 = csetemp1*csetemp5;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L30 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L31 = csetemp3*csetemp5;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L32 = -(csetemp0*csetemp5);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L33 = csetemp2;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp6 = SQR(ToReal(boostx));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp7 = SQR(ToReal(boosty));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp8 = SQR(ToReal(boostz));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp9 = INV(ToReal(lapsefactor));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L00 = csetemp9*INV((-1 + 
      csetemp6 + csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(1 - csetemp6 - csetemp7 - csetemp8 + sqrt(1 - csetemp6 - 
      csetemp7 - csetemp8))*(-1 + ToReal(boostx)*ToReal(shiftaddx) + 
      ToReal(boosty)*ToReal(shiftaddy) + ToReal(boostz)*ToReal(shiftaddz));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L01 = csetemp9*INV((-1 + 
      csetemp6 + csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(1 - csetemp6 - csetemp7 - csetemp8 + sqrt(1 - csetemp6 - 
      csetemp7 - csetemp8))*ToReal(boostx);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L02 = csetemp9*INV((-1 + 
      csetemp6 + csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(1 - csetemp6 - csetemp7 - csetemp8 + sqrt(1 - csetemp6 - 
      csetemp7 - csetemp8))*ToReal(boosty);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L03 = csetemp9*INV((-1 + 
      csetemp6 + csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(1 - csetemp6 - csetemp7 - csetemp8 + sqrt(1 - csetemp6 - 
      csetemp7 - csetemp8))*ToReal(boostz);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L10 = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*((csetemp6 + (-1 + csetemp7 + csetemp8)*(1 + sqrt(1 - 
      csetemp6 - csetemp7 - csetemp8)))*ToReal(shiftaddx) - 
      ToReal(boostx)*(-1 + csetemp6 + csetemp7 + csetemp8 + sqrt(1 - csetemp6 
      - csetemp7 - csetemp8)*(-1 + ToReal(boosty)*ToReal(shiftaddy) + 
      ToReal(boostz)*ToReal(shiftaddz))));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L11 = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(csetemp6 + (-1 + csetemp7 + csetemp8)*(1 + sqrt(1 - 
      csetemp6 - csetemp7 - csetemp8)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L12 = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boostx)*ToReal(boosty));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L13 = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boostx)*ToReal(boostz));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L20 = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*((csetemp7 + (-1 + csetemp6 + csetemp8)*(1 + sqrt(1 - 
      csetemp6 - csetemp7 - csetemp8)))*ToReal(shiftaddy) - 
      ToReal(boosty)*(-1 + csetemp6 + csetemp7 + csetemp8 + sqrt(1 - csetemp6 
      - csetemp7 - csetemp8)*(-1 + ToReal(boostx)*ToReal(shiftaddx) + 
      ToReal(boostz)*ToReal(shiftaddz))));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L21 = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boostx)*ToReal(boosty));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L22 = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(csetemp7 + (-1 + csetemp6 + csetemp8)*(1 + sqrt(1 - 
      csetemp6 - csetemp7 - csetemp8)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L23 = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boosty)*ToReal(boostz));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L30 = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(-(ToReal(boostz)*(-1 + csetemp6 + csetemp7 + csetemp8 + 
      sqrt(1 - csetemp6 - csetemp7 - csetemp8)*(-1 + 
      ToReal(boostx)*ToReal(shiftaddx) + ToReal(boosty)*ToReal(shiftaddy)))) 
      + (-1 + csetemp6 + csetemp7 + csetemp8 + (-1 + csetemp6 + 
      csetemp7)*sqrt(1 - csetemp6 - csetemp7 - csetemp8))*ToReal(shiftaddz));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L31 = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boostx)*ToReal(boostz));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L32 = -(INV(-1 + csetemp6 + 
      csetemp7 + csetemp8 - sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8))*ToReal(boosty)*ToReal(boostz));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L33 = INV((-1 + csetemp6 + 
      csetemp7 + csetemp8)*(1 + sqrt(1 - csetemp6 - csetemp7 - 
      csetemp8)))*(-1 + csetemp6 + csetemp7 + csetemp8 + (-1 + csetemp6 + 
      csetemp7)*sqrt(1 - csetemp6 - csetemp7 - csetemp8));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL00 = xform1L00*xform2L00 + 
      xform1L01*xform2L10 + xform1L02*xform2L20 + xform1L03*xform2L30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL01 = xform1L00*xform2L01 + 
      xform1L01*xform2L11 + xform1L02*xform2L21 + xform1L03*xform2L31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL02 = xform1L00*xform2L02 + 
      xform1L01*xform2L12 + xform1L02*xform2L22 + xform1L03*xform2L32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL03 = xform1L00*xform2L03 + 
      xform1L01*xform2L13 + xform1L02*xform2L23 + xform1L03*xform2L33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL10 = xform1L10*xform2L00 + 
      xform1L11*xform2L10 + xform1L12*xform2L20 + xform1L13*xform2L30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL11 = xform1L10*xform2L01 + 
      xform1L11*xform2L11 + xform1L12*xform2L21 + xform1L13*xform2L31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL12 = xform1L10*xform2L02 + 
      xform1L11*xform2L12 + xform1L12*xform2L22 + xform1L13*xform2L32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL13 = xform1L10*xform2L03 + 
      xform1L11*xform2L13 + xform1L12*xform2L23 + xform1L13*xform2L33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL20 = xform1L20*xform2L00 + 
      xform1L21*xform2L10 + xform1L22*xform2L20 + xform1L23*xform2L30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL21 = xform1L20*xform2L01 + 
      xform1L21*xform2L11 + xform1L22*xform2L21 + xform1L23*xform2L31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL22 = xform1L20*xform2L02 + 
      xform1L21*xform2L12 + xform1L22*xform2L22 + xform1L23*xform2L32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL23 = xform1L20*xform2L03 + 
      xform1L21*xform2L13 + xform1L22*xform2L23 + xform1L23*xform2L33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL30 = xform1L30*xform2L00 + 
      xform1L31*xform2L10 + xform1L32*xform2L20 + xform1L33*xform2L30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL31 = xform1L30*xform2L01 + 
      xform1L31*xform2L11 + xform1L32*xform2L21 + xform1L33*xform2L31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL32 = xform1L30*xform2L02 + 
      xform1L31*xform2L12 + xform1L32*xform2L22 + xform1L33*xform2L32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xformL33 = xform1L30*xform2L03 + 
      xform1L31*xform2L13 + xform1L32*xform2L23 + xform1L33*xform2L33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xx0 = t - ToReal(timeoffset);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xx1 = xL - ToReal(positionx);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xx2 = yL - ToReal(positiony);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xx3 = zL - ToReal(positionz);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED txx1 = xformL10*xx0 + xformL11*xx1 + 
      xformL12*xx2 + xformL13*xx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED txx2 = xformL20*xx0 + xformL21*xx1 + 
      xformL22*xx2 + xformL23*xx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED txx3 = xformL30*xx0 + xformL31*xx1 + 
      xformL32*xx2 + xformL33*xx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED X = txx1;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED Y = txx2;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED Z = txx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp10 = SQR(ToReal(a));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp11 = SQR(X);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp12 = SQR(Y);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp13 = SQR(Z);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED rXYZ = sqrt(INV(2)*(-csetemp10 + 
      csetemp11 + csetemp12 + csetemp13 + sqrt(4*csetemp10*csetemp13 + 
      SQR(-csetemp10 + csetemp11 + csetemp12 + csetemp13))));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp14 = CUB(rXYZ);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp15 = QAD(rXYZ);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg400 = -1 + 
      2*csetemp14*INV(csetemp10*csetemp13 + csetemp15)*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp16 = SQR(rXYZ);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp17 = rXYZ*X;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp18 = Y*ToReal(a);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg401 = 2*csetemp14*(csetemp17 + 
      csetemp18)*INV((csetemp10*csetemp13 + csetemp15)*(csetemp10 + 
      csetemp16))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp19 = X*ToReal(a);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp20 = -csetemp19;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp21 = rXYZ*Y;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg402 = 2*csetemp14*(csetemp20 + 
      csetemp21)*INV((csetemp10*csetemp13 + csetemp15)*(csetemp10 + 
      csetemp16))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg403 = 
      2*csetemp16*Z*INV(csetemp10*csetemp13 + csetemp15)*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg411 = 1 + 
      2*csetemp14*INV((csetemp10*csetemp13 + csetemp15)*SQR(csetemp10 + 
      csetemp16))*SQR(csetemp17 + csetemp18)*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg412 = 2*csetemp14*(csetemp17 + 
      csetemp18)*(csetemp20 + csetemp21)*INV((csetemp10*csetemp13 + 
      csetemp15)*SQR(csetemp10 + csetemp16))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg413 = 2*csetemp16*(csetemp17 + 
      csetemp18)*Z*INV((csetemp10*csetemp13 + csetemp15)*(csetemp10 + 
      csetemp16))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg422 = 1 + 
      2*csetemp14*INV((csetemp10*csetemp13 + csetemp15)*SQR(csetemp10 + 
      csetemp16))*SQR(csetemp20 + csetemp21)*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg423 = 2*csetemp16*(csetemp20 + 
      csetemp21)*Z*INV((csetemp10*csetemp13 + csetemp15)*(csetemp10 + 
      csetemp16))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg433 = 1 + 
      2*csetemp13*rXYZ*INV(csetemp10*csetemp13 + csetemp15)*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4000 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp22 = pow(rXYZ,7);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4001 = 
      2*(3*csetemp10*csetemp13*csetemp14 - 
      csetemp22)*X*INV(SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4002 = 
      2*(3*csetemp10*csetemp13*csetemp14 - 
      csetemp22)*Y*INV(SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp23 = pow(rXYZ,6);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp24 = QAD(ToReal(a));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4003 = 2*(csetemp10*(-5*csetemp15 + 
      (-2*csetemp10 + 2*(csetemp11 + csetemp12) + 5*csetemp13)*csetemp16) - 
      csetemp23 + 3*csetemp13*csetemp24)*rXYZ*Z*INV(SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4010 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp25 = pow(rXYZ,9);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp26 = CUB(ToReal(a));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp27 = pow(ToReal(a),5);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp28 = pow(rXYZ,5);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4011 = 
      2*csetemp14*INV(SQR(csetemp10*csetemp13 + csetemp15)*SQR(csetemp10 + 
      csetemp16)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((3*csetemp10 - 3*csetemp11 - csetemp12 - 
      csetemp13)*csetemp22 + 2*csetemp25 + csetemp10*((3*csetemp10 + 
      csetemp11 - csetemp12 - csetemp13)*csetemp13*csetemp14 + (csetemp10 - 
      csetemp11 - csetemp12 + csetemp13)*csetemp28) + csetemp13*((csetemp10 + 
      3*csetemp11 - csetemp12 - csetemp13)*csetemp24*rXYZ + 
      csetemp16*csetemp26*X*Y) + X*Y*(-(csetemp15*csetemp26) + 
      3*csetemp13*csetemp27 - 3*csetemp23*ToReal(a)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp29 = pow(rXYZ,8);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4012 = 
      2*csetemp14*INV(SQR(csetemp10*csetemp13 + csetemp15)*SQR(csetemp10 + 
      csetemp16)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*(((csetemp10 - csetemp11 - 2*csetemp12 + 
      csetemp13)*csetemp15 + (3*csetemp10 - csetemp11 - 
      csetemp13)*csetemp13*csetemp16)*csetemp26 + (-2*csetemp22 + 
      4*csetemp13*csetemp24*rXYZ)*X*Y + csetemp13*((csetemp10 - csetemp11 + 
      2*csetemp12 - csetemp13)*csetemp27 + 2*csetemp10*csetemp14*X*Y) + 
      ((3*csetemp10 - csetemp11 - 4*csetemp12 - csetemp13)*csetemp23 + 
      2*csetemp29)*ToReal(a))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4013 = 2*rXYZ*Z*INV((csetemp10 + 
      csetemp16)*SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((-2*csetemp22 + csetemp10*(2*(-csetemp10 + 
      csetemp11 + csetemp12 + 2*csetemp13)*csetemp14 - 4*csetemp28) + 
      4*csetemp13*csetemp24*rXYZ)*X + Y*((-5*csetemp15 + (-2*csetemp10 + 
      2*(csetemp11 + csetemp12) + 3*csetemp13)*csetemp16)*csetemp26 + 
      3*csetemp13*csetemp27 - 3*csetemp23*ToReal(a)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4020 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4021 = 
      -2*csetemp14*INV(SQR(csetemp10*csetemp13 + csetemp15)*SQR(csetemp10 + 
      csetemp16)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*(((csetemp10 - 2*csetemp11 - csetemp12 + 
      csetemp13)*csetemp15 + (3*csetemp10 - csetemp12 - 
      csetemp13)*csetemp13*csetemp16)*csetemp26 + (2*csetemp22 - 
      4*csetemp13*csetemp24*rXYZ)*X*Y + csetemp13*((csetemp10 + 2*csetemp11 - 
      csetemp12 - csetemp13)*csetemp27 - 2*csetemp10*csetemp14*X*Y) + 
      ((3*csetemp10 - 4*csetemp11 - csetemp12 - csetemp13)*csetemp23 + 
      2*csetemp29)*ToReal(a))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4022 = 
      2*csetemp14*INV(SQR(csetemp10*csetemp13 + csetemp15)*SQR(csetemp10 + 
      csetemp16)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((3*csetemp10 - csetemp11 - 3*csetemp12 - 
      csetemp13)*csetemp22 + 2*csetemp25 + csetemp10*((3*csetemp10 - 
      csetemp11 + csetemp12 - csetemp13)*csetemp13*csetemp14 + (csetemp10 - 
      csetemp11 - csetemp12 + csetemp13)*csetemp28) + csetemp13*((csetemp10 - 
      csetemp11 + 3*csetemp12 - csetemp13)*csetemp24*rXYZ - 
      csetemp16*csetemp26*X*Y) + X*Y*(csetemp15*csetemp26 - 
      3*csetemp13*csetemp27 + 3*csetemp23*ToReal(a)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4023 = 2*rXYZ*Z*INV((csetemp10 + 
      csetemp16)*SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((-2*csetemp22 + csetemp10*(2*(-csetemp10 + 
      csetemp11 + csetemp12 + 2*csetemp13)*csetemp14 - 4*csetemp28))*Y + 
      csetemp13*(-3*csetemp27*X + 4*csetemp24*rXYZ*Y) + X*((5*csetemp15 + 
      (2*csetemp10 - 2*(csetemp11 + csetemp12) - 
      3*csetemp13)*csetemp16)*csetemp26 + 3*csetemp23*ToReal(a)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4030 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4031 = 
      4*(csetemp10*csetemp13*csetemp16 - 
      csetemp23)*X*Z*INV(SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4032 = 
      4*(csetemp10*csetemp13*csetemp16 - 
      csetemp23)*Y*Z*INV(SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4033 = 2*(csetemp10*csetemp13 - 
      csetemp15)*(2*csetemp10*csetemp13 - 2*csetemp15 + (-csetemp10 + 
      csetemp11 + csetemp12 + 
      3*csetemp13)*csetemp16)*INV(SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4110 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4111 = 2*csetemp14*(csetemp17 + 
      csetemp18)*INV(CUB(csetemp10 + csetemp16)*SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((6*csetemp10 - 5*csetemp11 - 2*(csetemp12 + 
      csetemp13))*csetemp22 + 4*csetemp25 + csetemp10*(csetemp13*(6*csetemp10 
      - csetemp11 - 2*(csetemp12 + csetemp13))*csetemp14 + (-csetemp11 - 
      2*csetemp12 + 2*(csetemp10 + csetemp13))*csetemp28) + 
      csetemp13*((2*csetemp10 + 3*csetemp11 - 2*(csetemp12 + 
      csetemp13))*csetemp24*rXYZ - csetemp16*csetemp26*X*Y) + 
      X*Y*(-(csetemp15*csetemp26) + 3*csetemp13*csetemp27 - 
      5*csetemp23*ToReal(a)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4112 = 2*csetemp14*(csetemp17 + 
      csetemp18)*INV(CUB(csetemp10 + csetemp16)*SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*(((-2*csetemp11 - 3*csetemp12 + 2*(csetemp10 + 
      csetemp13))*csetemp15 + csetemp13*(6*csetemp10 - 3*csetemp12 - 
      2*(csetemp11 + csetemp13))*csetemp16)*csetemp26 + (-3*csetemp22 + 
      csetemp10*csetemp28 + 5*csetemp13*csetemp24*rXYZ)*X*Y + 
      csetemp13*((2*csetemp10 + csetemp12 - 2*(csetemp11 + 
      csetemp13))*csetemp27 + csetemp10*csetemp14*X*Y) + ((6*csetemp10 - 
      7*csetemp12 - 2*(csetemp11 + csetemp13))*csetemp23 + 
      4*csetemp29)*ToReal(a))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4113 = 2*(csetemp17 + 
      csetemp18)*rXYZ*Z*INV(SQR(csetemp10*csetemp13 + 
      csetemp15)*SQR(csetemp10 + csetemp16)*sqrt(4*csetemp10*csetemp13 + 
      SQR(-csetemp10 + csetemp11 + csetemp12 + csetemp13)))*((-3*csetemp22 + 
      csetemp10*((-2*csetemp10 + 2*(csetemp11 + csetemp12) + 
      3*csetemp13)*csetemp14 - 3*csetemp28) + 5*csetemp13*csetemp24*rXYZ)*X + 
      Y*((-5*csetemp15 + (-2*csetemp10 + 2*(csetemp11 + csetemp12) + 
      csetemp13)*csetemp16)*csetemp26 + 3*csetemp13*csetemp27 - 
      5*csetemp23*ToReal(a)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4120 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp30 = pow(rXYZ,10);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp31 = pow(ToReal(a),6);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4121 = 
      2*csetemp14*INV(CUB(csetemp10 + csetemp16)*SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((csetemp10*(2*(-csetemp10 + 3*csetemp11 + 
      csetemp13)*csetemp23 - csetemp15*(csetemp13*(csetemp12 + csetemp13) - 
      csetemp10*(csetemp12 + 2*(csetemp11 + csetemp13)) + csetemp24)) - 
      (-csetemp10 + 4*csetemp11 + csetemp12 + csetemp13)*csetemp29 + 
      2*csetemp30 + csetemp13*(-csetemp10 - 2*csetemp11 + csetemp12 + 
      csetemp13)*csetemp31)*Y - 2*csetemp13*((csetemp10 + csetemp11 - 
      3*csetemp12 - csetemp13)*csetemp27*rXYZ*X + (csetemp10 - 
      3*csetemp11)*csetemp16*csetemp24*Y) + 
      X*(csetemp26*(2*csetemp13*(-3*csetemp10 + csetemp11 + csetemp12 + 
      csetemp13)*csetemp14 - 2*(csetemp10 - csetemp11 - csetemp12 + 
      csetemp13)*csetemp28) + (2*(-3*csetemp10 + 3*csetemp11 - csetemp12 + 
      csetemp13)*csetemp22 - 4*csetemp25)*ToReal(a)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4122 = 
      -2*csetemp14*INV(CUB(csetemp10 + csetemp16)*SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((csetemp10*csetemp15*(csetemp13*(csetemp11 + 
      csetemp13) - csetemp10*(csetemp11 + 2*(csetemp12 + csetemp13)) + 
      csetemp24) + (-csetemp10 + csetemp11 + 4*csetemp12 + 
      csetemp13)*csetemp29 - 2*csetemp30 + (csetemp10 - csetemp11 + 
      2*csetemp12 - csetemp13)*csetemp13*csetemp31)*X + Y*(-2*((csetemp10 - 
      csetemp11 - csetemp12 + csetemp13)*csetemp26*csetemp28 + (csetemp10 - 
      3*csetemp11 + csetemp12 - csetemp13)*csetemp13*csetemp27*rXYZ) - 
      4*csetemp25*ToReal(a)) + 2*((csetemp10*(csetemp10 - 3*csetemp12 - 
      csetemp13)*csetemp23 + (csetemp10 - 
      3*csetemp12)*csetemp13*csetemp16*csetemp24)*X + 
      Y*(csetemp13*(-3*csetemp10 + csetemp11 + csetemp12 + 
      csetemp13)*csetemp14*csetemp26 + (-3*csetemp10 - csetemp11 + 
      3*csetemp12 + csetemp13)*csetemp22*ToReal(a))))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4123 = 
      2*rXYZ*Z*INV(SQR(csetemp10*csetemp13 + csetemp15)*SQR(csetemp10 + 
      csetemp16)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((csetemp10*((2*(csetemp11 + csetemp12) + 
      3*(csetemp10 + csetemp13))*csetemp15 + 2*csetemp23) + 2*(csetemp10 - 
      csetemp11 - csetemp12 + 2*csetemp13)*csetemp16*csetemp24 - 
      3*csetemp29)*X*Y + csetemp13*(4*(-csetemp11 + csetemp12)*csetemp27*rXYZ 
      - 3*csetemp31*X*Y) + (csetemp11 - csetemp12)*(csetemp26*(-2*(-csetemp10 
      + csetemp11 + csetemp12 + csetemp13)*csetemp14 + 4*csetemp28) + 
      4*csetemp22*ToReal(a)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4130 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4131 = 
      -2*csetemp16*Z*INV(SQR(csetemp10*csetemp13 + csetemp15)*SQR(csetemp10 + 
      csetemp16)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((-3*csetemp10 + 4*csetemp11 + csetemp12 + 
      csetemp13)*csetemp22 + csetemp10*(csetemp13*(-3*csetemp10 + csetemp12 + 
      csetemp13)*csetemp14 - (csetemp10 - 2*csetemp11 - csetemp12 + 
      csetemp13)*csetemp28) + csetemp13*(-csetemp10 - 2*csetemp11 + csetemp12 
      + csetemp13)*csetemp24*rXYZ - 2*(csetemp25 + csetemp13*csetemp27*X*Y) + 
      X*Y*(2*csetemp15*csetemp26 + 4*csetemp23*ToReal(a)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4132 = 
      -2*csetemp16*Z*INV(SQR(csetemp10*csetemp13 + csetemp15)*SQR(csetemp10 + 
      csetemp16)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((-((csetemp10 - csetemp11 - 3*csetemp12 + 
      csetemp13)*csetemp15) + csetemp13*(-3*csetemp10 + csetemp11 + csetemp12 
      + csetemp13)*csetemp16)*csetemp26 + (3*csetemp22 + csetemp10*csetemp28 
      - 3*csetemp13*csetemp24*rXYZ)*X*Y + csetemp13*((-csetemp10 + csetemp11 
      - csetemp12 + csetemp13)*csetemp27 - csetemp10*csetemp14*X*Y) + 
      ((-3*csetemp10 + csetemp11 + 5*csetemp12 + csetemp13)*csetemp23 - 
      2*csetemp29)*ToReal(a))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp32 = QAD(Z);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4133 = 2*INV((csetemp10 + 
      csetemp16)*SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((-((-csetemp10 + csetemp11 + csetemp12 + 
      4*csetemp13)*csetemp22) + 2*csetemp25 + 
      csetemp10*csetemp13*((-csetemp10 + csetemp11 + csetemp12 + 
      2*csetemp13)*csetemp14 - 3*csetemp28) + 3*csetemp24*csetemp32*rXYZ)*X + 
      Y*(csetemp13*(-4*csetemp15 + (-csetemp10 + csetemp11 + csetemp12 + 
      csetemp13)*csetemp16)*csetemp26 + (csetemp10 - csetemp11 - csetemp12 - 
      5*csetemp13)*csetemp23*ToReal(a) + 2*(csetemp27*csetemp32 + 
      csetemp29*ToReal(a))))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4220 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp33 = -csetemp21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4221 = 2*csetemp14*(csetemp19 + 
      csetemp33)*INV(CUB(csetemp10 + csetemp16)*SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*(((-3*csetemp11 - 2*csetemp12 + 2*(csetemp10 + 
      csetemp13))*csetemp15 + csetemp13*(6*csetemp10 - 3*csetemp11 - 
      2*(csetemp12 + csetemp13))*csetemp16)*csetemp26 + (3*csetemp22 - 
      csetemp10*csetemp28 - 5*csetemp13*csetemp24*rXYZ)*X*Y + 
      csetemp13*((2*csetemp10 + csetemp11 - 2*(csetemp12 + 
      csetemp13))*csetemp27 - csetemp10*csetemp14*X*Y) + ((6*csetemp10 - 
      7*csetemp11 - 2*(csetemp12 + csetemp13))*csetemp23 + 
      4*csetemp29)*ToReal(a))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4222 = 2*csetemp14*(csetemp20 + 
      csetemp21)*INV(CUB(csetemp10 + csetemp16)*SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((6*csetemp10 - 5*csetemp12 - 2*(csetemp11 + 
      csetemp13))*csetemp22 + 4*csetemp25 + csetemp10*(csetemp13*(6*csetemp10 
      - csetemp12 - 2*(csetemp11 + csetemp13))*csetemp14 + (-2*csetemp11 - 
      csetemp12 + 2*(csetemp10 + csetemp13))*csetemp28) + 
      csetemp13*((2*csetemp10 + 3*csetemp12 - 2*(csetemp11 + 
      csetemp13))*csetemp24*rXYZ + csetemp16*csetemp26*X*Y) + 
      X*Y*(csetemp15*csetemp26 - 3*csetemp13*csetemp27 + 
      5*csetemp23*ToReal(a)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4223 = 2*(csetemp19 + 
      csetemp33)*rXYZ*Z*INV(SQR(csetemp10*csetemp13 + 
      csetemp15)*SQR(csetemp10 + csetemp16)*sqrt(4*csetemp10*csetemp13 + 
      SQR(-csetemp10 + csetemp11 + csetemp12 + csetemp13)))*((-5*csetemp15 + 
      (-2*csetemp10 + 2*(csetemp11 + csetemp12) + 
      csetemp13)*csetemp16)*csetemp26*X + csetemp10*((2*csetemp10 - 
      2*(csetemp11 + csetemp12) - 3*csetemp13)*csetemp14 + 3*csetemp28)*Y + 
      3*(csetemp13*csetemp27*X + csetemp22*Y) - 5*(csetemp13*csetemp24*rXYZ*Y 
      + csetemp23*X*ToReal(a)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4230 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4231 = 
      2*csetemp16*Z*INV(SQR(csetemp10*csetemp13 + csetemp15)*SQR(csetemp10 + 
      csetemp16)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((-((csetemp10 - 3*csetemp11 - csetemp12 + 
      csetemp13)*csetemp15) + csetemp13*(-3*csetemp10 + csetemp11 + csetemp12 
      + csetemp13)*csetemp16)*csetemp26 + (-3*csetemp22 - csetemp10*csetemp28 
      + 3*csetemp13*csetemp24*rXYZ)*X*Y + csetemp13*((-csetemp10 - csetemp11 
      + csetemp12 + csetemp13)*csetemp27 + csetemp10*csetemp14*X*Y) + 
      ((-3*csetemp10 + 5*csetemp11 + csetemp12 + csetemp13)*csetemp23 - 
      2*csetemp29)*ToReal(a))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4232 = 
      -2*csetemp16*Z*INV(SQR(csetemp10*csetemp13 + csetemp15)*SQR(csetemp10 + 
      csetemp16)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*((-3*csetemp10 + csetemp11 + 4*csetemp12 + 
      csetemp13)*csetemp22 + csetemp10*(csetemp13*(-3*csetemp10 + csetemp11 + 
      csetemp13)*csetemp14 - (csetemp10 - csetemp11 - 2*csetemp12 + 
      csetemp13)*csetemp28) - 2*(csetemp25 + csetemp15*csetemp26*X*Y) + 
      csetemp13*((-csetemp10 + csetemp11 - 2*csetemp12 + 
      csetemp13)*csetemp24*rXYZ + 2*csetemp27*X*Y) - 
      4*csetemp23*X*Y*ToReal(a))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4233 = -2*INV((csetemp10 + 
      csetemp16)*SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*(((-csetemp10 + csetemp11 + csetemp12 + 
      4*csetemp13)*csetemp22 - 2*csetemp25 + 
      3*csetemp10*csetemp13*csetemp28)*Y + csetemp13*((-4*csetemp15 + 
      (-csetemp10 + csetemp11 + csetemp12 + csetemp13)*csetemp16)*csetemp26*X 
      + csetemp10*(csetemp10 - csetemp11 - csetemp12 - 
      2*csetemp13)*csetemp14*Y) + csetemp32*(2*csetemp27*X - 
      3*csetemp24*rXYZ*Y) + ((csetemp10 - csetemp11 - csetemp12 - 
      5*csetemp13)*csetemp23 + 2*csetemp29)*X*ToReal(a))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4330 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4331 = 2*csetemp13*(-3*csetemp28 + 
      csetemp10*csetemp13*rXYZ)*X*INV(SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4332 = 2*csetemp13*(-3*csetemp28 + 
      csetemp10*csetemp13*rXYZ)*Y*INV(SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4333 = 2*((2*csetemp10 - 
      2*(csetemp11 + csetemp12) - 5*csetemp13)*csetemp23 + 4*csetemp29 + 
      csetemp24*csetemp32 + csetemp10*(-3*csetemp13*csetemp15 + 
      csetemp16*csetemp32))*Z*INV(rXYZ*SQR(csetemp10*csetemp13 + 
      csetemp15)*sqrt(4*csetemp10*csetemp13 + SQR(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g400 = 2*(tg423*xformL20*xformL30 + 
      xformL00*(tg401*xformL10 + tg402*xformL20 + tg403*xformL30) + 
      xformL10*(tg412*xformL20 + tg413*xformL30)) + tg400*SQR(xformL00) + 
      tg411*SQR(xformL10) + tg422*SQR(xformL20) + tg433*SQR(xformL30);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp34 = tg400*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp35 = tg401*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp36 = tg402*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp37 = tg403*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp38 = tg401*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp39 = tg411*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp40 = tg412*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp41 = tg413*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp42 = tg402*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp43 = tg412*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp44 = tg422*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp45 = tg423*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp46 = tg403*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp47 = tg413*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp48 = tg423*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp49 = tg433*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g401 = (csetemp34 + csetemp35 + 
      csetemp36 + csetemp37)*xformL00 + (csetemp38 + csetemp39 + csetemp40 + 
      csetemp41)*xformL10 + (csetemp42 + csetemp43 + csetemp44 + 
      csetemp45)*xformL20 + (csetemp46 + csetemp47 + csetemp48 + 
      csetemp49)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp50 = tg400*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp51 = tg401*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp52 = tg402*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp53 = tg403*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp54 = tg401*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp55 = tg411*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp56 = tg412*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp57 = tg413*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp58 = tg402*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp59 = tg412*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp60 = tg422*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp61 = tg423*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp62 = tg403*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp63 = tg413*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp64 = tg423*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp65 = tg433*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g402 = (csetemp50 + csetemp51 + 
      csetemp52 + csetemp53)*xformL00 + (csetemp54 + csetemp55 + csetemp56 + 
      csetemp57)*xformL10 + (csetemp58 + csetemp59 + csetemp60 + 
      csetemp61)*xformL20 + (csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp66 = tg400*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp67 = tg401*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp68 = tg402*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp69 = tg403*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp70 = tg401*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp71 = tg411*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp72 = tg412*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp73 = tg413*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp74 = tg402*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp75 = tg412*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp76 = tg422*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp77 = tg423*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp78 = tg403*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp79 = tg413*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp80 = tg423*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp81 = tg433*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g403 = (csetemp66 + csetemp67 + 
      csetemp68 + csetemp69)*xformL00 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL10 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL20 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g411 = (csetemp34 + csetemp35 + 
      csetemp36 + csetemp37)*xformL01 + (csetemp38 + csetemp39 + csetemp40 + 
      csetemp41)*xformL11 + (csetemp42 + csetemp43 + csetemp44 + 
      csetemp45)*xformL21 + (csetemp46 + csetemp47 + csetemp48 + 
      csetemp49)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g412 = (csetemp50 + csetemp51 + 
      csetemp52 + csetemp53)*xformL01 + (csetemp54 + csetemp55 + csetemp56 + 
      csetemp57)*xformL11 + (csetemp58 + csetemp59 + csetemp60 + 
      csetemp61)*xformL21 + (csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g413 = (csetemp66 + csetemp67 + 
      csetemp68 + csetemp69)*xformL01 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL11 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL21 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g422 = (csetemp50 + csetemp51 + 
      csetemp52 + csetemp53)*xformL02 + (csetemp54 + csetemp55 + csetemp56 + 
      csetemp57)*xformL12 + (csetemp58 + csetemp59 + csetemp60 + 
      csetemp61)*xformL22 + (csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g423 = (csetemp66 + csetemp67 + 
      csetemp68 + csetemp69)*xformL02 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL12 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL22 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g433 = (csetemp66 + csetemp67 + 
      csetemp68 + csetemp69)*xformL03 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL13 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL23 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp82 = tdg4000*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp83 = tdg4001*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp84 = tdg4002*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp85 = tdg4003*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp86 = tdg4010*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp87 = tdg4011*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp88 = tdg4012*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp89 = tdg4013*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp90 = tdg4020*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp91 = tdg4021*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp92 = tdg4022*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp93 = tdg4023*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp94 = tdg4030*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp95 = tdg4031*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp96 = tdg4032*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp97 = tdg4033*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp98 = tdg4110*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp99 = tdg4111*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp100 = tdg4112*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp101 = tdg4113*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp102 = tdg4120*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp103 = tdg4121*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp104 = tdg4122*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp105 = tdg4123*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp106 = tdg4130*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp107 = tdg4131*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp108 = tdg4132*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp109 = tdg4133*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp110 = tdg4220*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp111 = tdg4221*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp112 = tdg4222*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp113 = tdg4223*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp114 = tdg4230*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp115 = tdg4231*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp116 = tdg4232*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp117 = tdg4233*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp118 = tdg4330*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp119 = tdg4331*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp120 = tdg4332*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp121 = tdg4333*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4000 = xformL10*((csetemp86 + 
      csetemp87 + csetemp88 + csetemp89)*xformL00 + (csetemp100 + csetemp101 
      + csetemp98 + csetemp99)*xformL10 + (csetemp102 + csetemp103 + 
      csetemp104 + csetemp105)*xformL20 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL30) + xformL20*((csetemp90 + csetemp91 + 
      csetemp92 + csetemp93)*xformL00 + (csetemp102 + csetemp103 + csetemp104 
      + csetemp105)*xformL10 + (csetemp110 + csetemp111 + csetemp112 + 
      csetemp113)*xformL20 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL30) + xformL30*((csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL00 + (csetemp106 + csetemp107 + csetemp108 + 
      csetemp109)*xformL10 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL20 + (csetemp118 + csetemp119 + csetemp120 + 
      csetemp121)*xformL30) + xformL00*((csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL00 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL10 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL20 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL30);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4010 = xformL10*((csetemp86 + 
      csetemp87 + csetemp88 + csetemp89)*xformL01 + (csetemp100 + csetemp101 
      + csetemp98 + csetemp99)*xformL11 + (csetemp102 + csetemp103 + 
      csetemp104 + csetemp105)*xformL21 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL31) + xformL20*((csetemp90 + csetemp91 + 
      csetemp92 + csetemp93)*xformL01 + (csetemp102 + csetemp103 + csetemp104 
      + csetemp105)*xformL11 + (csetemp110 + csetemp111 + csetemp112 + 
      csetemp113)*xformL21 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL31) + xformL30*((csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL01 + (csetemp106 + csetemp107 + csetemp108 + 
      csetemp109)*xformL11 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL21 + (csetemp118 + csetemp119 + csetemp120 + 
      csetemp121)*xformL31) + xformL00*((csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL01 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL11 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL21 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp122 = tdg4000*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp123 = tdg4001*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp124 = tdg4002*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp125 = tdg4003*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp126 = tdg4010*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp127 = tdg4011*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp128 = tdg4012*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp129 = tdg4013*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp130 = tdg4020*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp131 = tdg4021*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp132 = tdg4022*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp133 = tdg4023*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp134 = tdg4030*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp135 = tdg4031*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp136 = tdg4032*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp137 = tdg4033*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp138 = tdg4110*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp139 = tdg4111*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp140 = tdg4112*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp141 = tdg4113*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp142 = tdg4120*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp143 = tdg4121*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp144 = tdg4122*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp145 = tdg4123*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp146 = tdg4130*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp147 = tdg4131*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp148 = tdg4132*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp149 = tdg4133*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp150 = tdg4220*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp151 = tdg4221*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp152 = tdg4222*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp153 = tdg4223*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp154 = tdg4230*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp155 = tdg4231*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp156 = tdg4232*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp157 = tdg4233*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp158 = tdg4330*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp159 = tdg4331*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp160 = tdg4332*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp161 = tdg4333*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4011 = xformL00*((csetemp122 + 
      csetemp123 + csetemp124 + csetemp125)*xformL01 + (csetemp126 + 
      csetemp127 + csetemp128 + csetemp129)*xformL11 + (csetemp130 + 
      csetemp131 + csetemp132 + csetemp133)*xformL21 + (csetemp134 + 
      csetemp135 + csetemp136 + csetemp137)*xformL31) + xformL10*((csetemp126 
      + csetemp127 + csetemp128 + csetemp129)*xformL01 + (csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL11 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL21 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL31) + xformL20*((csetemp130 
      + csetemp131 + csetemp132 + csetemp133)*xformL01 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL11 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL21 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL31) + xformL30*((csetemp134 
      + csetemp135 + csetemp136 + csetemp137)*xformL01 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL11 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL21 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp162 = tdg4000*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp163 = tdg4001*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp164 = tdg4002*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp165 = tdg4003*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp166 = tdg4010*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp167 = tdg4011*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp168 = tdg4012*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp169 = tdg4013*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp170 = tdg4020*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp171 = tdg4021*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp172 = tdg4022*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp173 = tdg4023*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp174 = tdg4030*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp175 = tdg4031*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp176 = tdg4032*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp177 = tdg4033*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp178 = tdg4110*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp179 = tdg4111*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp180 = tdg4112*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp181 = tdg4113*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp182 = tdg4120*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp183 = tdg4121*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp184 = tdg4122*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp185 = tdg4123*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp186 = tdg4130*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp187 = tdg4131*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp188 = tdg4132*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp189 = tdg4133*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp190 = tdg4220*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp191 = tdg4221*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp192 = tdg4222*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp193 = tdg4223*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp194 = tdg4230*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp195 = tdg4231*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp196 = tdg4232*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp197 = tdg4233*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp198 = tdg4330*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp199 = tdg4331*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp200 = tdg4332*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp201 = tdg4333*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4012 = xformL00*((csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL01 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL11 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL21 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL31) + xformL10*((csetemp166 
      + csetemp167 + csetemp168 + csetemp169)*xformL01 + (csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL11 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL21 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL31) + xformL20*((csetemp170 
      + csetemp171 + csetemp172 + csetemp173)*xformL01 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL11 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL21 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL31) + xformL30*((csetemp174 
      + csetemp175 + csetemp176 + csetemp177)*xformL01 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL11 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL21 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp202 = tdg4000*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp203 = tdg4001*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp204 = tdg4002*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp205 = tdg4003*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp206 = tdg4010*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp207 = tdg4011*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp208 = tdg4012*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp209 = tdg4013*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp210 = tdg4020*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp211 = tdg4021*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp212 = tdg4022*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp213 = tdg4023*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp214 = tdg4030*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp215 = tdg4031*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp216 = tdg4032*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp217 = tdg4033*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp218 = tdg4110*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp219 = tdg4111*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp220 = tdg4112*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp221 = tdg4113*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp222 = tdg4120*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp223 = tdg4121*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp224 = tdg4122*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp225 = tdg4123*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp226 = tdg4130*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp227 = tdg4131*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp228 = tdg4132*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp229 = tdg4133*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp230 = tdg4220*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp231 = tdg4221*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp232 = tdg4222*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp233 = tdg4223*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp234 = tdg4230*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp235 = tdg4231*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp236 = tdg4232*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp237 = tdg4233*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp238 = tdg4330*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp239 = tdg4331*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp240 = tdg4332*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp241 = tdg4333*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4013 = xformL00*((csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL01 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL11 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL21 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL31) + xformL10*((csetemp206 
      + csetemp207 + csetemp208 + csetemp209)*xformL01 + (csetemp218 + 
      csetemp219 + csetemp220 + csetemp221)*xformL11 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL21 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL31) + xformL20*((csetemp210 
      + csetemp211 + csetemp212 + csetemp213)*xformL01 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL11 + (csetemp230 + 
      csetemp231 + csetemp232 + csetemp233)*xformL21 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL31) + xformL30*((csetemp214 
      + csetemp215 + csetemp216 + csetemp217)*xformL01 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL11 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL21 + (csetemp238 + 
      csetemp239 + csetemp240 + csetemp241)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4020 = xformL10*((csetemp86 + 
      csetemp87 + csetemp88 + csetemp89)*xformL02 + (csetemp100 + csetemp101 
      + csetemp98 + csetemp99)*xformL12 + (csetemp102 + csetemp103 + 
      csetemp104 + csetemp105)*xformL22 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL32) + xformL20*((csetemp90 + csetemp91 + 
      csetemp92 + csetemp93)*xformL02 + (csetemp102 + csetemp103 + csetemp104 
      + csetemp105)*xformL12 + (csetemp110 + csetemp111 + csetemp112 + 
      csetemp113)*xformL22 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL32) + xformL30*((csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL02 + (csetemp106 + csetemp107 + csetemp108 + 
      csetemp109)*xformL12 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL22 + (csetemp118 + csetemp119 + csetemp120 + 
      csetemp121)*xformL32) + xformL00*((csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL02 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL12 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL22 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4021 = xformL00*((csetemp122 + 
      csetemp123 + csetemp124 + csetemp125)*xformL02 + (csetemp126 + 
      csetemp127 + csetemp128 + csetemp129)*xformL12 + (csetemp130 + 
      csetemp131 + csetemp132 + csetemp133)*xformL22 + (csetemp134 + 
      csetemp135 + csetemp136 + csetemp137)*xformL32) + xformL10*((csetemp126 
      + csetemp127 + csetemp128 + csetemp129)*xformL02 + (csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL12 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL22 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL32) + xformL20*((csetemp130 
      + csetemp131 + csetemp132 + csetemp133)*xformL02 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL12 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL22 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL32) + xformL30*((csetemp134 
      + csetemp135 + csetemp136 + csetemp137)*xformL02 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL12 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL22 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4022 = xformL00*((csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL02 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL12 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL22 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL32) + xformL10*((csetemp166 
      + csetemp167 + csetemp168 + csetemp169)*xformL02 + (csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL12 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL22 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL32) + xformL20*((csetemp170 
      + csetemp171 + csetemp172 + csetemp173)*xformL02 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL12 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL22 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL32) + xformL30*((csetemp174 
      + csetemp175 + csetemp176 + csetemp177)*xformL02 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL12 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL22 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4023 = xformL00*((csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL02 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL12 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL22 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL32) + xformL10*((csetemp206 
      + csetemp207 + csetemp208 + csetemp209)*xformL02 + (csetemp218 + 
      csetemp219 + csetemp220 + csetemp221)*xformL12 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL22 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL32) + xformL20*((csetemp210 
      + csetemp211 + csetemp212 + csetemp213)*xformL02 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL12 + (csetemp230 + 
      csetemp231 + csetemp232 + csetemp233)*xformL22 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL32) + xformL30*((csetemp214 
      + csetemp215 + csetemp216 + csetemp217)*xformL02 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL12 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL22 + (csetemp238 + 
      csetemp239 + csetemp240 + csetemp241)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4030 = xformL10*((csetemp86 + 
      csetemp87 + csetemp88 + csetemp89)*xformL03 + (csetemp100 + csetemp101 
      + csetemp98 + csetemp99)*xformL13 + (csetemp102 + csetemp103 + 
      csetemp104 + csetemp105)*xformL23 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL33) + xformL20*((csetemp90 + csetemp91 + 
      csetemp92 + csetemp93)*xformL03 + (csetemp102 + csetemp103 + csetemp104 
      + csetemp105)*xformL13 + (csetemp110 + csetemp111 + csetemp112 + 
      csetemp113)*xformL23 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL33) + xformL30*((csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL03 + (csetemp106 + csetemp107 + csetemp108 + 
      csetemp109)*xformL13 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL23 + (csetemp118 + csetemp119 + csetemp120 + 
      csetemp121)*xformL33) + xformL00*((csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL03 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL13 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL23 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4031 = xformL00*((csetemp122 + 
      csetemp123 + csetemp124 + csetemp125)*xformL03 + (csetemp126 + 
      csetemp127 + csetemp128 + csetemp129)*xformL13 + (csetemp130 + 
      csetemp131 + csetemp132 + csetemp133)*xformL23 + (csetemp134 + 
      csetemp135 + csetemp136 + csetemp137)*xformL33) + xformL10*((csetemp126 
      + csetemp127 + csetemp128 + csetemp129)*xformL03 + (csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL13 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL23 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL33) + xformL20*((csetemp130 
      + csetemp131 + csetemp132 + csetemp133)*xformL03 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL13 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL23 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL33) + xformL30*((csetemp134 
      + csetemp135 + csetemp136 + csetemp137)*xformL03 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL13 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL23 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4032 = xformL00*((csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL03 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL13 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL23 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL33) + xformL10*((csetemp166 
      + csetemp167 + csetemp168 + csetemp169)*xformL03 + (csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL13 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL23 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL33) + xformL20*((csetemp170 
      + csetemp171 + csetemp172 + csetemp173)*xformL03 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL13 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL23 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL33) + xformL30*((csetemp174 
      + csetemp175 + csetemp176 + csetemp177)*xformL03 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL13 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL23 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4033 = xformL00*((csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL03 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL13 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL23 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL33) + xformL10*((csetemp206 
      + csetemp207 + csetemp208 + csetemp209)*xformL03 + (csetemp218 + 
      csetemp219 + csetemp220 + csetemp221)*xformL13 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL23 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL33) + xformL20*((csetemp210 
      + csetemp211 + csetemp212 + csetemp213)*xformL03 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL13 + (csetemp230 + 
      csetemp231 + csetemp232 + csetemp233)*xformL23 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL33) + xformL30*((csetemp214 
      + csetemp215 + csetemp216 + csetemp217)*xformL03 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL13 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL23 + (csetemp238 + 
      csetemp239 + csetemp240 + csetemp241)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4110 = xformL11*((csetemp86 + 
      csetemp87 + csetemp88 + csetemp89)*xformL01 + (csetemp100 + csetemp101 
      + csetemp98 + csetemp99)*xformL11 + (csetemp102 + csetemp103 + 
      csetemp104 + csetemp105)*xformL21 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL31) + xformL21*((csetemp90 + csetemp91 + 
      csetemp92 + csetemp93)*xformL01 + (csetemp102 + csetemp103 + csetemp104 
      + csetemp105)*xformL11 + (csetemp110 + csetemp111 + csetemp112 + 
      csetemp113)*xformL21 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL31) + xformL31*((csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL01 + (csetemp106 + csetemp107 + csetemp108 + 
      csetemp109)*xformL11 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL21 + (csetemp118 + csetemp119 + csetemp120 + 
      csetemp121)*xformL31) + xformL01*((csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL01 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL11 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL21 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4111 = xformL01*((csetemp122 + 
      csetemp123 + csetemp124 + csetemp125)*xformL01 + (csetemp126 + 
      csetemp127 + csetemp128 + csetemp129)*xformL11 + (csetemp130 + 
      csetemp131 + csetemp132 + csetemp133)*xformL21 + (csetemp134 + 
      csetemp135 + csetemp136 + csetemp137)*xformL31) + xformL11*((csetemp126 
      + csetemp127 + csetemp128 + csetemp129)*xformL01 + (csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL11 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL21 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL31) + xformL21*((csetemp130 
      + csetemp131 + csetemp132 + csetemp133)*xformL01 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL11 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL21 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL31) + xformL31*((csetemp134 
      + csetemp135 + csetemp136 + csetemp137)*xformL01 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL11 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL21 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4112 = xformL01*((csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL01 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL11 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL21 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL31) + xformL11*((csetemp166 
      + csetemp167 + csetemp168 + csetemp169)*xformL01 + (csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL11 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL21 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL31) + xformL21*((csetemp170 
      + csetemp171 + csetemp172 + csetemp173)*xformL01 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL11 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL21 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL31) + xformL31*((csetemp174 
      + csetemp175 + csetemp176 + csetemp177)*xformL01 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL11 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL21 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4113 = xformL01*((csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL01 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL11 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL21 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL31) + xformL11*((csetemp206 
      + csetemp207 + csetemp208 + csetemp209)*xformL01 + (csetemp218 + 
      csetemp219 + csetemp220 + csetemp221)*xformL11 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL21 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL31) + xformL21*((csetemp210 
      + csetemp211 + csetemp212 + csetemp213)*xformL01 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL11 + (csetemp230 + 
      csetemp231 + csetemp232 + csetemp233)*xformL21 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL31) + xformL31*((csetemp214 
      + csetemp215 + csetemp216 + csetemp217)*xformL01 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL11 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL21 + (csetemp238 + 
      csetemp239 + csetemp240 + csetemp241)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4120 = xformL11*((csetemp86 + 
      csetemp87 + csetemp88 + csetemp89)*xformL02 + (csetemp100 + csetemp101 
      + csetemp98 + csetemp99)*xformL12 + (csetemp102 + csetemp103 + 
      csetemp104 + csetemp105)*xformL22 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL32) + xformL21*((csetemp90 + csetemp91 + 
      csetemp92 + csetemp93)*xformL02 + (csetemp102 + csetemp103 + csetemp104 
      + csetemp105)*xformL12 + (csetemp110 + csetemp111 + csetemp112 + 
      csetemp113)*xformL22 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL32) + xformL31*((csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL02 + (csetemp106 + csetemp107 + csetemp108 + 
      csetemp109)*xformL12 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL22 + (csetemp118 + csetemp119 + csetemp120 + 
      csetemp121)*xformL32) + xformL01*((csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL02 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL12 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL22 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4121 = xformL01*((csetemp122 + 
      csetemp123 + csetemp124 + csetemp125)*xformL02 + (csetemp126 + 
      csetemp127 + csetemp128 + csetemp129)*xformL12 + (csetemp130 + 
      csetemp131 + csetemp132 + csetemp133)*xformL22 + (csetemp134 + 
      csetemp135 + csetemp136 + csetemp137)*xformL32) + xformL11*((csetemp126 
      + csetemp127 + csetemp128 + csetemp129)*xformL02 + (csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL12 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL22 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL32) + xformL21*((csetemp130 
      + csetemp131 + csetemp132 + csetemp133)*xformL02 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL12 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL22 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL32) + xformL31*((csetemp134 
      + csetemp135 + csetemp136 + csetemp137)*xformL02 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL12 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL22 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4122 = xformL01*((csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL02 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL12 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL22 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL32) + xformL11*((csetemp166 
      + csetemp167 + csetemp168 + csetemp169)*xformL02 + (csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL12 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL22 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL32) + xformL21*((csetemp170 
      + csetemp171 + csetemp172 + csetemp173)*xformL02 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL12 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL22 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL32) + xformL31*((csetemp174 
      + csetemp175 + csetemp176 + csetemp177)*xformL02 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL12 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL22 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4123 = xformL01*((csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL02 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL12 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL22 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL32) + xformL11*((csetemp206 
      + csetemp207 + csetemp208 + csetemp209)*xformL02 + (csetemp218 + 
      csetemp219 + csetemp220 + csetemp221)*xformL12 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL22 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL32) + xformL21*((csetemp210 
      + csetemp211 + csetemp212 + csetemp213)*xformL02 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL12 + (csetemp230 + 
      csetemp231 + csetemp232 + csetemp233)*xformL22 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL32) + xformL31*((csetemp214 
      + csetemp215 + csetemp216 + csetemp217)*xformL02 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL12 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL22 + (csetemp238 + 
      csetemp239 + csetemp240 + csetemp241)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4130 = xformL11*((csetemp86 + 
      csetemp87 + csetemp88 + csetemp89)*xformL03 + (csetemp100 + csetemp101 
      + csetemp98 + csetemp99)*xformL13 + (csetemp102 + csetemp103 + 
      csetemp104 + csetemp105)*xformL23 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL33) + xformL21*((csetemp90 + csetemp91 + 
      csetemp92 + csetemp93)*xformL03 + (csetemp102 + csetemp103 + csetemp104 
      + csetemp105)*xformL13 + (csetemp110 + csetemp111 + csetemp112 + 
      csetemp113)*xformL23 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL33) + xformL31*((csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL03 + (csetemp106 + csetemp107 + csetemp108 + 
      csetemp109)*xformL13 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL23 + (csetemp118 + csetemp119 + csetemp120 + 
      csetemp121)*xformL33) + xformL01*((csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL03 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL13 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL23 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4131 = xformL01*((csetemp122 + 
      csetemp123 + csetemp124 + csetemp125)*xformL03 + (csetemp126 + 
      csetemp127 + csetemp128 + csetemp129)*xformL13 + (csetemp130 + 
      csetemp131 + csetemp132 + csetemp133)*xformL23 + (csetemp134 + 
      csetemp135 + csetemp136 + csetemp137)*xformL33) + xformL11*((csetemp126 
      + csetemp127 + csetemp128 + csetemp129)*xformL03 + (csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL13 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL23 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL33) + xformL21*((csetemp130 
      + csetemp131 + csetemp132 + csetemp133)*xformL03 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL13 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL23 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL33) + xformL31*((csetemp134 
      + csetemp135 + csetemp136 + csetemp137)*xformL03 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL13 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL23 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4132 = xformL01*((csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL03 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL13 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL23 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL33) + xformL11*((csetemp166 
      + csetemp167 + csetemp168 + csetemp169)*xformL03 + (csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL13 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL23 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL33) + xformL21*((csetemp170 
      + csetemp171 + csetemp172 + csetemp173)*xformL03 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL13 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL23 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL33) + xformL31*((csetemp174 
      + csetemp175 + csetemp176 + csetemp177)*xformL03 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL13 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL23 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4133 = xformL01*((csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL03 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL13 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL23 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL33) + xformL11*((csetemp206 
      + csetemp207 + csetemp208 + csetemp209)*xformL03 + (csetemp218 + 
      csetemp219 + csetemp220 + csetemp221)*xformL13 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL23 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL33) + xformL21*((csetemp210 
      + csetemp211 + csetemp212 + csetemp213)*xformL03 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL13 + (csetemp230 + 
      csetemp231 + csetemp232 + csetemp233)*xformL23 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL33) + xformL31*((csetemp214 
      + csetemp215 + csetemp216 + csetemp217)*xformL03 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL13 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL23 + (csetemp238 + 
      csetemp239 + csetemp240 + csetemp241)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4220 = xformL12*((csetemp86 + 
      csetemp87 + csetemp88 + csetemp89)*xformL02 + (csetemp100 + csetemp101 
      + csetemp98 + csetemp99)*xformL12 + (csetemp102 + csetemp103 + 
      csetemp104 + csetemp105)*xformL22 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL32) + xformL22*((csetemp90 + csetemp91 + 
      csetemp92 + csetemp93)*xformL02 + (csetemp102 + csetemp103 + csetemp104 
      + csetemp105)*xformL12 + (csetemp110 + csetemp111 + csetemp112 + 
      csetemp113)*xformL22 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL32) + xformL32*((csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL02 + (csetemp106 + csetemp107 + csetemp108 + 
      csetemp109)*xformL12 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL22 + (csetemp118 + csetemp119 + csetemp120 + 
      csetemp121)*xformL32) + xformL02*((csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL02 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL12 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL22 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4221 = xformL02*((csetemp122 + 
      csetemp123 + csetemp124 + csetemp125)*xformL02 + (csetemp126 + 
      csetemp127 + csetemp128 + csetemp129)*xformL12 + (csetemp130 + 
      csetemp131 + csetemp132 + csetemp133)*xformL22 + (csetemp134 + 
      csetemp135 + csetemp136 + csetemp137)*xformL32) + xformL12*((csetemp126 
      + csetemp127 + csetemp128 + csetemp129)*xformL02 + (csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL12 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL22 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL32) + xformL22*((csetemp130 
      + csetemp131 + csetemp132 + csetemp133)*xformL02 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL12 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL22 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL32) + xformL32*((csetemp134 
      + csetemp135 + csetemp136 + csetemp137)*xformL02 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL12 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL22 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4222 = xformL02*((csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL02 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL12 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL22 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL32) + xformL12*((csetemp166 
      + csetemp167 + csetemp168 + csetemp169)*xformL02 + (csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL12 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL22 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL32) + xformL22*((csetemp170 
      + csetemp171 + csetemp172 + csetemp173)*xformL02 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL12 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL22 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL32) + xformL32*((csetemp174 
      + csetemp175 + csetemp176 + csetemp177)*xformL02 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL12 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL22 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4223 = xformL02*((csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL02 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL12 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL22 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL32) + xformL12*((csetemp206 
      + csetemp207 + csetemp208 + csetemp209)*xformL02 + (csetemp218 + 
      csetemp219 + csetemp220 + csetemp221)*xformL12 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL22 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL32) + xformL22*((csetemp210 
      + csetemp211 + csetemp212 + csetemp213)*xformL02 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL12 + (csetemp230 + 
      csetemp231 + csetemp232 + csetemp233)*xformL22 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL32) + xformL32*((csetemp214 
      + csetemp215 + csetemp216 + csetemp217)*xformL02 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL12 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL22 + (csetemp238 + 
      csetemp239 + csetemp240 + csetemp241)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4230 = xformL12*((csetemp86 + 
      csetemp87 + csetemp88 + csetemp89)*xformL03 + (csetemp100 + csetemp101 
      + csetemp98 + csetemp99)*xformL13 + (csetemp102 + csetemp103 + 
      csetemp104 + csetemp105)*xformL23 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL33) + xformL22*((csetemp90 + csetemp91 + 
      csetemp92 + csetemp93)*xformL03 + (csetemp102 + csetemp103 + csetemp104 
      + csetemp105)*xformL13 + (csetemp110 + csetemp111 + csetemp112 + 
      csetemp113)*xformL23 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL33) + xformL32*((csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL03 + (csetemp106 + csetemp107 + csetemp108 + 
      csetemp109)*xformL13 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL23 + (csetemp118 + csetemp119 + csetemp120 + 
      csetemp121)*xformL33) + xformL02*((csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL03 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL13 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL23 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4231 = xformL02*((csetemp122 + 
      csetemp123 + csetemp124 + csetemp125)*xformL03 + (csetemp126 + 
      csetemp127 + csetemp128 + csetemp129)*xformL13 + (csetemp130 + 
      csetemp131 + csetemp132 + csetemp133)*xformL23 + (csetemp134 + 
      csetemp135 + csetemp136 + csetemp137)*xformL33) + xformL12*((csetemp126 
      + csetemp127 + csetemp128 + csetemp129)*xformL03 + (csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL13 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL23 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL33) + xformL22*((csetemp130 
      + csetemp131 + csetemp132 + csetemp133)*xformL03 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL13 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL23 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL33) + xformL32*((csetemp134 
      + csetemp135 + csetemp136 + csetemp137)*xformL03 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL13 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL23 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4232 = xformL02*((csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL03 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL13 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL23 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL33) + xformL12*((csetemp166 
      + csetemp167 + csetemp168 + csetemp169)*xformL03 + (csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL13 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL23 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL33) + xformL22*((csetemp170 
      + csetemp171 + csetemp172 + csetemp173)*xformL03 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL13 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL23 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL33) + xformL32*((csetemp174 
      + csetemp175 + csetemp176 + csetemp177)*xformL03 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL13 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL23 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4233 = xformL02*((csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL03 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL13 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL23 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL33) + xformL12*((csetemp206 
      + csetemp207 + csetemp208 + csetemp209)*xformL03 + (csetemp218 + 
      csetemp219 + csetemp220 + csetemp221)*xformL13 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL23 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL33) + xformL22*((csetemp210 
      + csetemp211 + csetemp212 + csetemp213)*xformL03 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL13 + (csetemp230 + 
      csetemp231 + csetemp232 + csetemp233)*xformL23 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL33) + xformL32*((csetemp214 
      + csetemp215 + csetemp216 + csetemp217)*xformL03 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL13 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL23 + (csetemp238 + 
      csetemp239 + csetemp240 + csetemp241)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4330 = xformL13*((csetemp86 + 
      csetemp87 + csetemp88 + csetemp89)*xformL03 + (csetemp100 + csetemp101 
      + csetemp98 + csetemp99)*xformL13 + (csetemp102 + csetemp103 + 
      csetemp104 + csetemp105)*xformL23 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL33) + xformL23*((csetemp90 + csetemp91 + 
      csetemp92 + csetemp93)*xformL03 + (csetemp102 + csetemp103 + csetemp104 
      + csetemp105)*xformL13 + (csetemp110 + csetemp111 + csetemp112 + 
      csetemp113)*xformL23 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL33) + xformL33*((csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL03 + (csetemp106 + csetemp107 + csetemp108 + 
      csetemp109)*xformL13 + (csetemp114 + csetemp115 + csetemp116 + 
      csetemp117)*xformL23 + (csetemp118 + csetemp119 + csetemp120 + 
      csetemp121)*xformL33) + xformL03*((csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL03 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL13 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL23 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4331 = xformL03*((csetemp122 + 
      csetemp123 + csetemp124 + csetemp125)*xformL03 + (csetemp126 + 
      csetemp127 + csetemp128 + csetemp129)*xformL13 + (csetemp130 + 
      csetemp131 + csetemp132 + csetemp133)*xformL23 + (csetemp134 + 
      csetemp135 + csetemp136 + csetemp137)*xformL33) + xformL13*((csetemp126 
      + csetemp127 + csetemp128 + csetemp129)*xformL03 + (csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL13 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL23 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL33) + xformL23*((csetemp130 
      + csetemp131 + csetemp132 + csetemp133)*xformL03 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL13 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL23 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL33) + xformL33*((csetemp134 
      + csetemp135 + csetemp136 + csetemp137)*xformL03 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL13 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL23 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4332 = xformL03*((csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL03 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL13 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL23 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL33) + xformL13*((csetemp166 
      + csetemp167 + csetemp168 + csetemp169)*xformL03 + (csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL13 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL23 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL33) + xformL23*((csetemp170 
      + csetemp171 + csetemp172 + csetemp173)*xformL03 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL13 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL23 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL33) + xformL33*((csetemp174 
      + csetemp175 + csetemp176 + csetemp177)*xformL03 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL13 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL23 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4333 = xformL03*((csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL03 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL13 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL23 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL33) + xformL13*((csetemp206 
      + csetemp207 + csetemp208 + csetemp209)*xformL03 + (csetemp218 + 
      csetemp219 + csetemp220 + csetemp221)*xformL13 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL23 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL33) + xformL23*((csetemp210 
      + csetemp211 + csetemp212 + csetemp213)*xformL03 + (csetemp222 + 
      csetemp223 + csetemp224 + csetemp225)*xformL13 + (csetemp230 + 
      csetemp231 + csetemp232 + csetemp233)*xformL23 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL33) + xformL33*((csetemp214 
      + csetemp215 + csetemp216 + csetemp217)*xformL03 + (csetemp226 + 
      csetemp227 + csetemp228 + csetemp229)*xformL13 + (csetemp234 + 
      csetemp235 + csetemp236 + csetemp237)*xformL23 + (csetemp238 + 
      csetemp239 + csetemp240 + csetemp241)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal1 = g401;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal2 = g402;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal3 = g403;
    
    gxxL = g411;
    
    gxyL = g412;
    
    gxzL = g413;
    
    gyyL = g422;
    
    gyzL = g423;
    
    gzzL = g433;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp242 = SQR(gxzL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp243 = SQR(gyzL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp244 = SQR(gxyL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED detg = 2*gxyL*gxzL*gyzL + 
      gyyL*(gxxL*gzzL - csetemp242) - gxxL*csetemp243 - 
      gzzL*csetemp244;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp245 = INV(detg);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu11 = (gyyL*gzzL - 
      csetemp243)*csetemp245;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu12 = (gxzL*gyzL - 
      gxyL*gzzL)*csetemp245;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu13 = (-(gxzL*gyyL) + 
      gxyL*gyzL)*csetemp245;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu22 = (gxxL*gzzL - 
      csetemp242)*csetemp245;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu23 = (gxyL*gxzL - 
      gxxL*gyzL)*csetemp245;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu33 = (gxxL*gyyL - 
      csetemp244)*csetemp245;
    
    betaxL = betal1*gu11 + betal2*gu12 + betal3*gu13;
    
    betayL = betal1*gu12 + betal2*gu22 + betal3*gu23;
    
    betazL = betal1*gu13 + betal2*gu23 + betal3*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betasq = betaxL*betal1 + 
      betayL*betal2 + betazL*betal3;
    
    alpL = sqrt(betasq - g400);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtg11 = dg4110;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtg12 = dg4120;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtg13 = dg4130;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtg22 = dg4220;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtg23 = dg4230;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtg33 = dg4330;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg111 = dg4111;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg112 = dg4112;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg113 = dg4113;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg121 = dg4121;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg122 = dg4122;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg123 = dg4123;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg131 = dg4131;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg132 = dg4132;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg133 = dg4133;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg221 = dg4221;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg222 = dg4222;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg223 = dg4223;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg231 = dg4231;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg232 = dg4232;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg233 = dg4233;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg331 = dg4331;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg332 = dg4332;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg333 = dg4333;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp246 = dtg11*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp247 = dtg12*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp248 = dtg13*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp249 = dtg12*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp250 = dtg22*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp251 = dtg23*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp252 = dtg13*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp253 = dtg23*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp254 = dtg33*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu11 = -((csetemp246 + csetemp247 + 
      csetemp248)*gu11) - (csetemp249 + csetemp250 + csetemp251)*gu12 - 
      (csetemp252 + csetemp253 + csetemp254)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu12 = -((csetemp246 + csetemp247 + 
      csetemp248)*gu12) - (csetemp249 + csetemp250 + csetemp251)*gu22 - 
      (csetemp252 + csetemp253 + csetemp254)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu13 = -((csetemp246 + csetemp247 + 
      csetemp248)*gu13) - (csetemp249 + csetemp250 + csetemp251)*gu23 - 
      (csetemp252 + csetemp253 + csetemp254)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp255 = dtg11*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp256 = dtg12*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp257 = dtg13*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp258 = dtg22*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp259 = dtg23*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp260 = dtg13*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp261 = dtg23*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp262 = dtg33*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu22 = -((csetemp255 + csetemp256 + 
      csetemp257)*gu12) - (csetemp247 + csetemp258 + csetemp259)*gu22 - 
      (csetemp260 + csetemp261 + csetemp262)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu23 = -((csetemp255 + csetemp256 + 
      csetemp257)*gu13) - (csetemp247 + csetemp258 + csetemp259)*gu23 - 
      (csetemp260 + csetemp261 + csetemp262)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu33 = -(gu13*(dtg11*gu13 + 
      dtg12*gu23 + dtg13*gu33)) - gu23*(dtg12*gu13 + dtg22*gu23 + dtg23*gu33) 
      - gu33*(csetemp248 + csetemp259 + dtg33*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp263 = dg111*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp264 = dg121*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp265 = dg131*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp266 = dg121*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp267 = dg221*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp268 = dg231*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp269 = dg131*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp270 = dg231*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp271 = dg331*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu111 = -((csetemp263 + csetemp264 + 
      csetemp265)*gu11) - (csetemp266 + csetemp267 + csetemp268)*gu12 - 
      (csetemp269 + csetemp270 + csetemp271)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu121 = -((csetemp263 + csetemp264 + 
      csetemp265)*gu12) - (csetemp266 + csetemp267 + csetemp268)*gu22 - 
      (csetemp269 + csetemp270 + csetemp271)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu131 = -((csetemp263 + csetemp264 + 
      csetemp265)*gu13) - (csetemp266 + csetemp267 + csetemp268)*gu23 - 
      (csetemp269 + csetemp270 + csetemp271)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp272 = dg111*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp273 = dg121*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp274 = dg131*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp275 = dg221*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp276 = dg231*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp277 = dg131*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp278 = dg231*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp279 = dg331*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu221 = -((csetemp272 + csetemp273 + 
      csetemp274)*gu12) - (csetemp264 + csetemp275 + csetemp276)*gu22 - 
      (csetemp277 + csetemp278 + csetemp279)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu231 = -((csetemp272 + csetemp273 + 
      csetemp274)*gu13) - (csetemp264 + csetemp275 + csetemp276)*gu23 - 
      (csetemp277 + csetemp278 + csetemp279)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu331 = -(gu13*(dg111*gu13 + 
      dg121*gu23 + dg131*gu33)) - gu23*(dg121*gu13 + dg221*gu23 + dg231*gu33) 
      - gu33*(csetemp265 + csetemp276 + dg331*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp280 = dg112*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp281 = dg122*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp282 = dg132*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp283 = dg122*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp284 = dg222*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp285 = dg232*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp286 = dg132*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp287 = dg232*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp288 = dg332*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu112 = -((csetemp280 + csetemp281 + 
      csetemp282)*gu11) - (csetemp283 + csetemp284 + csetemp285)*gu12 - 
      (csetemp286 + csetemp287 + csetemp288)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu122 = -((csetemp280 + csetemp281 + 
      csetemp282)*gu12) - (csetemp283 + csetemp284 + csetemp285)*gu22 - 
      (csetemp286 + csetemp287 + csetemp288)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu132 = -((csetemp280 + csetemp281 + 
      csetemp282)*gu13) - (csetemp283 + csetemp284 + csetemp285)*gu23 - 
      (csetemp286 + csetemp287 + csetemp288)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp289 = dg112*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp290 = dg122*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp291 = dg132*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp292 = dg222*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp293 = dg232*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp294 = dg132*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp295 = dg232*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp296 = dg332*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu222 = -((csetemp289 + csetemp290 + 
      csetemp291)*gu12) - (csetemp281 + csetemp292 + csetemp293)*gu22 - 
      (csetemp294 + csetemp295 + csetemp296)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu232 = -((csetemp289 + csetemp290 + 
      csetemp291)*gu13) - (csetemp281 + csetemp292 + csetemp293)*gu23 - 
      (csetemp294 + csetemp295 + csetemp296)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu332 = -(gu13*(dg112*gu13 + 
      dg122*gu23 + dg132*gu33)) - gu23*(dg122*gu13 + dg222*gu23 + dg232*gu33) 
      - gu33*(csetemp282 + csetemp293 + dg332*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp297 = dg113*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp298 = dg123*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp299 = dg133*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp300 = dg123*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp301 = dg223*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp302 = dg233*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp303 = dg133*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp304 = dg233*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp305 = dg333*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu113 = -((csetemp297 + csetemp298 + 
      csetemp299)*gu11) - (csetemp300 + csetemp301 + csetemp302)*gu12 - 
      (csetemp303 + csetemp304 + csetemp305)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu123 = -((csetemp297 + csetemp298 + 
      csetemp299)*gu12) - (csetemp300 + csetemp301 + csetemp302)*gu22 - 
      (csetemp303 + csetemp304 + csetemp305)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu133 = -((csetemp297 + csetemp298 + 
      csetemp299)*gu13) - (csetemp300 + csetemp301 + csetemp302)*gu23 - 
      (csetemp303 + csetemp304 + csetemp305)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp306 = dg113*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp307 = dg123*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp308 = dg133*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp309 = dg223*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp310 = dg233*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp311 = dg133*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp312 = dg233*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp313 = dg333*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu223 = -((csetemp306 + csetemp307 + 
      csetemp308)*gu12) - (csetemp298 + csetemp309 + csetemp310)*gu22 - 
      (csetemp311 + csetemp312 + csetemp313)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu233 = -((csetemp306 + csetemp307 + 
      csetemp308)*gu13) - (csetemp298 + csetemp309 + csetemp310)*gu23 - 
      (csetemp311 + csetemp312 + csetemp313)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu333 = -(gu13*(dg113*gu13 + 
      dg123*gu23 + dg133*gu33)) - gu23*(dg123*gu13 + dg223*gu23 + dg233*gu33) 
      - gu33*(csetemp299 + csetemp310 + dg333*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtbetal1 = dg4010;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtbetal2 = dg4020;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtbetal3 = dg4030;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbetal11 = dg4011;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbetal12 = dg4012;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbetal13 = dg4013;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbetal21 = dg4021;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbetal22 = dg4022;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbetal23 = dg4023;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbetal31 = dg4031;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbetal32 = dg4032;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbetal33 = dg4033;
    
    dtbetaxL = betal1*dtgu11 + betal2*dtgu12 + betal3*dtgu13 + 
      dtbetal1*gu11 + dtbetal2*gu12 + dtbetal3*gu13;
    
    dtbetayL = betal1*dtgu12 + betal2*dtgu22 + betal3*dtgu23 + 
      dtbetal1*gu12 + dtbetal2*gu22 + dtbetal3*gu23;
    
    dtbetazL = betal1*dtgu13 + betal2*dtgu23 + betal3*dtgu33 + 
      dtbetal1*gu13 + dtbetal2*gu23 + dtbetal3*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbeta11 = betal1*dgu111 + 
      betal2*dgu121 + betal3*dgu131 + dbetal11*gu11 + dbetal21*gu12 + 
      dbetal31*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbeta21 = betal1*dgu121 + 
      betal2*dgu221 + betal3*dgu231 + dbetal11*gu12 + dbetal21*gu22 + 
      dbetal31*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbeta31 = betal1*dgu131 + 
      betal2*dgu231 + betal3*dgu331 + dbetal11*gu13 + dbetal21*gu23 + 
      dbetal31*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbeta12 = betal1*dgu112 + 
      betal2*dgu122 + betal3*dgu132 + dbetal12*gu11 + dbetal22*gu12 + 
      dbetal32*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbeta22 = betal1*dgu122 + 
      betal2*dgu222 + betal3*dgu232 + dbetal12*gu12 + dbetal22*gu22 + 
      dbetal32*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbeta32 = betal1*dgu132 + 
      betal2*dgu232 + betal3*dgu332 + dbetal12*gu13 + dbetal22*gu23 + 
      dbetal32*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbeta13 = betal1*dgu113 + 
      betal2*dgu123 + betal3*dgu133 + dbetal13*gu11 + dbetal23*gu12 + 
      dbetal33*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbeta23 = betal1*dgu123 + 
      betal2*dgu223 + betal3*dgu233 + dbetal13*gu12 + dbetal23*gu22 + 
      dbetal33*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dbeta33 = betal1*dgu133 + 
      betal2*dgu233 + betal3*dgu333 + dbetal13*gu13 + dbetal23*gu23 + 
      dbetal33*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtbetasq = dtbetaxL*betal1 + 
      dtbetayL*betal2 + dtbetazL*betal3 + betaxL*dtbetal1 + 
      betayL*dtbetal2 + betazL*dtbetal3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp314 = INV(alpL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtalpL = 0.5*csetemp314*(-dg4000 + 
      dtbetasq);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxxL = 
      0.5*csetemp314*(2*(gxxL*dbeta11 + gxyL*dbeta21 + gxzL*dbeta31) + 
      betaxL*dg111 + betayL*dg112 + betazL*dg113 - dtg11);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxyL = 0.5*csetemp314*(gxxL*dbeta12 
      + gyyL*dbeta21 + gxyL*(dbeta11 + dbeta22) + gyzL*dbeta31 + 
      gxzL*dbeta32 + betaxL*dg121 + betayL*dg122 + betazL*dg123 - 
      dtg12);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxzL = 0.5*csetemp314*(gxxL*dbeta13 
      + gyzL*dbeta21 + gxyL*dbeta23 + gzzL*dbeta31 + gxzL*(dbeta11 + 
      dbeta33) + betaxL*dg131 + betayL*dg132 + betazL*dg133 - dtg13);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kyyL = 
      0.5*csetemp314*(2*(gxyL*dbeta12 + gyyL*dbeta22 + gyzL*dbeta32) + 
      betaxL*dg221 + betayL*dg222 + betazL*dg223 - dtg22);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kyzL = 0.5*csetemp314*(gxzL*dbeta12 
      + gxyL*dbeta13 + gyyL*dbeta23 + gzzL*dbeta32 + gyzL*(dbeta22 + 
      dbeta33) + betaxL*dg231 + betayL*dg232 + betazL*dg233 - dtg23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kzzL = 
      0.5*csetemp314*(2*(gxzL*dbeta13 + gyzL*dbeta23 + gzzL*dbeta33) + 
      betaxL*dg331 + betayL*dg332 + betazL*dg333 - dtg33);
    
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
  CCTK_ENDLOOP3(KerrSchild_initial);
}

extern "C" void KerrSchild_initial(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering KerrSchild_initial_Body");
  }
  
  if (cctk_iteration % KerrSchild_initial_calc_every != KerrSchild_initial_calc_offset)
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
  GenericFD_AssertGroupStorage(cctkGH, "KerrSchild_initial", 7, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, KerrSchild_initial_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving KerrSchild_initial_Body");
  }
}
