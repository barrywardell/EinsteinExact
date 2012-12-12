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

static void Minkowski_always_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_LOOP3(Minkowski_always,
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg400 = -1;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg401 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg402 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg403 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg411 = 1;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg412 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg413 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg422 = 1;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg423 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg433 = 1;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4000 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4001 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4002 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4003 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4010 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4011 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4012 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4013 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4020 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4021 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4022 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4023 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4030 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4031 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4032 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4033 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4110 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4111 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4112 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4113 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4120 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4121 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4122 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4123 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4130 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4131 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4132 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4133 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4220 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4221 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4222 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4223 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4230 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4231 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4232 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4233 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4330 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4331 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4332 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4333 = 0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g400 = 2*(tg423*xformL20*xformL30 + 
      xformL00*(tg401*xformL10 + tg402*xformL20 + tg403*xformL30) + 
      xformL10*(tg412*xformL20 + tg413*xformL30)) + tg400*SQR(xformL00) + 
      tg411*SQR(xformL10) + tg422*SQR(xformL20) + tg433*SQR(xformL30);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp10 = tg400*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp11 = tg401*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp12 = tg402*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp13 = tg403*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp14 = tg401*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp15 = tg411*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp16 = tg412*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp17 = tg413*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp18 = tg402*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp19 = tg412*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp20 = tg422*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp21 = tg423*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp22 = tg403*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp23 = tg413*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp24 = tg423*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp25 = tg433*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g401 = (csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)*xformL00 + (csetemp14 + csetemp15 + csetemp16 + 
      csetemp17)*xformL10 + (csetemp18 + csetemp19 + csetemp20 + 
      csetemp21)*xformL20 + (csetemp22 + csetemp23 + csetemp24 + 
      csetemp25)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp26 = tg400*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp27 = tg401*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp28 = tg402*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp29 = tg403*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp30 = tg401*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp31 = tg411*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp32 = tg412*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp33 = tg413*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp34 = tg402*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp35 = tg412*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp36 = tg422*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp37 = tg423*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp38 = tg403*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp39 = tg413*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp40 = tg423*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp41 = tg433*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g402 = (csetemp26 + csetemp27 + 
      csetemp28 + csetemp29)*xformL00 + (csetemp30 + csetemp31 + csetemp32 + 
      csetemp33)*xformL10 + (csetemp34 + csetemp35 + csetemp36 + 
      csetemp37)*xformL20 + (csetemp38 + csetemp39 + csetemp40 + 
      csetemp41)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp42 = tg400*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp43 = tg401*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp44 = tg402*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp45 = tg403*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp46 = tg401*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp47 = tg411*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp48 = tg412*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp49 = tg413*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp50 = tg402*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp51 = tg412*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp52 = tg422*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp53 = tg423*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp54 = tg403*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp55 = tg413*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp56 = tg423*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp57 = tg433*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g403 = (csetemp42 + csetemp43 + 
      csetemp44 + csetemp45)*xformL00 + (csetemp46 + csetemp47 + csetemp48 + 
      csetemp49)*xformL10 + (csetemp50 + csetemp51 + csetemp52 + 
      csetemp53)*xformL20 + (csetemp54 + csetemp55 + csetemp56 + 
      csetemp57)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g411 = (csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)*xformL01 + (csetemp14 + csetemp15 + csetemp16 + 
      csetemp17)*xformL11 + (csetemp18 + csetemp19 + csetemp20 + 
      csetemp21)*xformL21 + (csetemp22 + csetemp23 + csetemp24 + 
      csetemp25)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g412 = (csetemp26 + csetemp27 + 
      csetemp28 + csetemp29)*xformL01 + (csetemp30 + csetemp31 + csetemp32 + 
      csetemp33)*xformL11 + (csetemp34 + csetemp35 + csetemp36 + 
      csetemp37)*xformL21 + (csetemp38 + csetemp39 + csetemp40 + 
      csetemp41)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g413 = (csetemp42 + csetemp43 + 
      csetemp44 + csetemp45)*xformL01 + (csetemp46 + csetemp47 + csetemp48 + 
      csetemp49)*xformL11 + (csetemp50 + csetemp51 + csetemp52 + 
      csetemp53)*xformL21 + (csetemp54 + csetemp55 + csetemp56 + 
      csetemp57)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g422 = (csetemp26 + csetemp27 + 
      csetemp28 + csetemp29)*xformL02 + (csetemp30 + csetemp31 + csetemp32 + 
      csetemp33)*xformL12 + (csetemp34 + csetemp35 + csetemp36 + 
      csetemp37)*xformL22 + (csetemp38 + csetemp39 + csetemp40 + 
      csetemp41)*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g423 = (csetemp42 + csetemp43 + 
      csetemp44 + csetemp45)*xformL02 + (csetemp46 + csetemp47 + csetemp48 + 
      csetemp49)*xformL12 + (csetemp50 + csetemp51 + csetemp52 + 
      csetemp53)*xformL22 + (csetemp54 + csetemp55 + csetemp56 + 
      csetemp57)*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g433 = (csetemp42 + csetemp43 + 
      csetemp44 + csetemp45)*xformL03 + (csetemp46 + csetemp47 + csetemp48 + 
      csetemp49)*xformL13 + (csetemp50 + csetemp51 + csetemp52 + 
      csetemp53)*xformL23 + (csetemp54 + csetemp55 + csetemp56 + 
      csetemp57)*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp58 = tdg4000*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp59 = tdg4001*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp60 = tdg4002*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp61 = tdg4003*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp62 = tdg4010*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp63 = tdg4011*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp64 = tdg4012*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp65 = tdg4013*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp66 = tdg4020*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp67 = tdg4021*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp68 = tdg4022*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp69 = tdg4023*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp70 = tdg4030*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp71 = tdg4031*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp72 = tdg4032*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp73 = tdg4033*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp74 = tdg4110*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp75 = tdg4111*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp76 = tdg4112*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp77 = tdg4113*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp78 = tdg4120*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp79 = tdg4121*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp80 = tdg4122*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp81 = tdg4123*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp82 = tdg4130*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp83 = tdg4131*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp84 = tdg4132*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp85 = tdg4133*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp86 = tdg4220*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp87 = tdg4221*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp88 = tdg4222*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp89 = tdg4223*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp90 = tdg4230*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp91 = tdg4231*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp92 = tdg4232*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp93 = tdg4233*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp94 = tdg4330*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp95 = tdg4331*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp96 = tdg4332*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp97 = tdg4333*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4000 = xformL00*((csetemp58 + 
      csetemp59 + csetemp60 + csetemp61)*xformL00 + (csetemp62 + csetemp63 + 
      csetemp64 + csetemp65)*xformL10 + (csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL20 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL30) + xformL10*((csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL00 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL10 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL20 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL30) + xformL20*((csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL00 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL10 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL20 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL30) + xformL30*((csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL00 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL10 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL20 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL30);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4010 = xformL00*((csetemp58 + 
      csetemp59 + csetemp60 + csetemp61)*xformL01 + (csetemp62 + csetemp63 + 
      csetemp64 + csetemp65)*xformL11 + (csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL21 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL31) + xformL10*((csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL01 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL11 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL21 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL31) + xformL20*((csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL01 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL11 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL21 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL31) + xformL30*((csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL01 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL11 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL21 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp98 = tdg4000*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp99 = tdg4001*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp100 = tdg4002*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp101 = tdg4003*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp102 = tdg4010*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp103 = tdg4011*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp104 = tdg4012*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp105 = tdg4013*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp106 = tdg4020*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp107 = tdg4021*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp108 = tdg4022*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp109 = tdg4023*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp110 = tdg4030*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp111 = tdg4031*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp112 = tdg4032*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp113 = tdg4033*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp114 = tdg4110*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp115 = tdg4111*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp116 = tdg4112*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp117 = tdg4113*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp118 = tdg4120*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp119 = tdg4121*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp120 = tdg4122*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp121 = tdg4123*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp122 = tdg4130*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp123 = tdg4131*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp124 = tdg4132*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp125 = tdg4133*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp126 = tdg4220*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp127 = tdg4221*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp128 = tdg4222*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp129 = tdg4223*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp130 = tdg4230*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp131 = tdg4231*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp132 = tdg4232*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp133 = tdg4233*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp134 = tdg4330*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp135 = tdg4331*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp136 = tdg4332*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp137 = tdg4333*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4011 = xformL00*((csetemp100 + 
      csetemp101 + csetemp98 + csetemp99)*xformL01 + (csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL11 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL21 + (csetemp110 + csetemp111 + 
      csetemp112 + csetemp113)*xformL31) + xformL10*((csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL01 + (csetemp114 + csetemp115 + 
      csetemp116 + csetemp117)*xformL11 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL21 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL31) + xformL20*((csetemp106 + csetemp107 
      + csetemp108 + csetemp109)*xformL01 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL11 + (csetemp126 + csetemp127 + 
      csetemp128 + csetemp129)*xformL21 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL31) + xformL30*((csetemp110 + csetemp111 
      + csetemp112 + csetemp113)*xformL01 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL11 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL21 + (csetemp134 + csetemp135 + 
      csetemp136 + csetemp137)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp138 = tdg4000*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp139 = tdg4001*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp140 = tdg4002*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp141 = tdg4003*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp142 = tdg4010*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp143 = tdg4011*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp144 = tdg4012*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp145 = tdg4013*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp146 = tdg4020*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp147 = tdg4021*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp148 = tdg4022*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp149 = tdg4023*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp150 = tdg4030*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp151 = tdg4031*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp152 = tdg4032*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp153 = tdg4033*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp154 = tdg4110*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp155 = tdg4111*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp156 = tdg4112*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp157 = tdg4113*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp158 = tdg4120*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp159 = tdg4121*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp160 = tdg4122*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp161 = tdg4123*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp162 = tdg4130*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp163 = tdg4131*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp164 = tdg4132*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp165 = tdg4133*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp166 = tdg4220*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp167 = tdg4221*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp168 = tdg4222*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp169 = tdg4223*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp170 = tdg4230*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp171 = tdg4231*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp172 = tdg4232*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp173 = tdg4233*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp174 = tdg4330*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp175 = tdg4331*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp176 = tdg4332*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp177 = tdg4333*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4012 = xformL00*((csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL01 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL11 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL21 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL31) + xformL10*((csetemp142 
      + csetemp143 + csetemp144 + csetemp145)*xformL01 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL11 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL21 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL31) + xformL20*((csetemp146 
      + csetemp147 + csetemp148 + csetemp149)*xformL01 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL11 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL21 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL31) + xformL30*((csetemp150 
      + csetemp151 + csetemp152 + csetemp153)*xformL01 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL11 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL21 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp178 = tdg4000*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp179 = tdg4001*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp180 = tdg4002*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp181 = tdg4003*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp182 = tdg4010*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp183 = tdg4011*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp184 = tdg4012*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp185 = tdg4013*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp186 = tdg4020*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp187 = tdg4021*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp188 = tdg4022*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp189 = tdg4023*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp190 = tdg4030*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp191 = tdg4031*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp192 = tdg4032*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp193 = tdg4033*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp194 = tdg4110*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp195 = tdg4111*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp196 = tdg4112*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp197 = tdg4113*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp198 = tdg4120*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp199 = tdg4121*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp200 = tdg4122*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp201 = tdg4123*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp202 = tdg4130*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp203 = tdg4131*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp204 = tdg4132*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp205 = tdg4133*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp206 = tdg4220*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp207 = tdg4221*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp208 = tdg4222*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp209 = tdg4223*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp210 = tdg4230*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp211 = tdg4231*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp212 = tdg4232*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp213 = tdg4233*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp214 = tdg4330*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp215 = tdg4331*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp216 = tdg4332*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp217 = tdg4333*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4013 = xformL00*((csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL01 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL11 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL21 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL31) + xformL10*((csetemp182 
      + csetemp183 + csetemp184 + csetemp185)*xformL01 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL11 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL21 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL31) + xformL20*((csetemp186 
      + csetemp187 + csetemp188 + csetemp189)*xformL01 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL11 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL21 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL31) + xformL30*((csetemp190 
      + csetemp191 + csetemp192 + csetemp193)*xformL01 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL11 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL21 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4020 = xformL00*((csetemp58 + 
      csetemp59 + csetemp60 + csetemp61)*xformL02 + (csetemp62 + csetemp63 + 
      csetemp64 + csetemp65)*xformL12 + (csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL22 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL32) + xformL10*((csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL02 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL12 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL22 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL32) + xformL20*((csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL02 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL12 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL22 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL32) + xformL30*((csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL02 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL12 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL22 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4021 = xformL00*((csetemp100 + 
      csetemp101 + csetemp98 + csetemp99)*xformL02 + (csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL12 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL22 + (csetemp110 + csetemp111 + 
      csetemp112 + csetemp113)*xformL32) + xformL10*((csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL02 + (csetemp114 + csetemp115 + 
      csetemp116 + csetemp117)*xformL12 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL22 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL32) + xformL20*((csetemp106 + csetemp107 
      + csetemp108 + csetemp109)*xformL02 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL12 + (csetemp126 + csetemp127 + 
      csetemp128 + csetemp129)*xformL22 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL32) + xformL30*((csetemp110 + csetemp111 
      + csetemp112 + csetemp113)*xformL02 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL12 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL22 + (csetemp134 + csetemp135 + 
      csetemp136 + csetemp137)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4022 = xformL00*((csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL02 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL12 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL22 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL32) + xformL10*((csetemp142 
      + csetemp143 + csetemp144 + csetemp145)*xformL02 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL12 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL22 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL32) + xformL20*((csetemp146 
      + csetemp147 + csetemp148 + csetemp149)*xformL02 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL12 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL22 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL32) + xformL30*((csetemp150 
      + csetemp151 + csetemp152 + csetemp153)*xformL02 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL12 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL22 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4023 = xformL00*((csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL02 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL12 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL22 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL32) + xformL10*((csetemp182 
      + csetemp183 + csetemp184 + csetemp185)*xformL02 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL12 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL22 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL32) + xformL20*((csetemp186 
      + csetemp187 + csetemp188 + csetemp189)*xformL02 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL12 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL22 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL32) + xformL30*((csetemp190 
      + csetemp191 + csetemp192 + csetemp193)*xformL02 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL12 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL22 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4030 = xformL00*((csetemp58 + 
      csetemp59 + csetemp60 + csetemp61)*xformL03 + (csetemp62 + csetemp63 + 
      csetemp64 + csetemp65)*xformL13 + (csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL23 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL33) + xformL10*((csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL03 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL13 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL23 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL33) + xformL20*((csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL03 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL13 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL23 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL33) + xformL30*((csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL03 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL13 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL23 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4031 = xformL00*((csetemp100 + 
      csetemp101 + csetemp98 + csetemp99)*xformL03 + (csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL13 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL23 + (csetemp110 + csetemp111 + 
      csetemp112 + csetemp113)*xformL33) + xformL10*((csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL03 + (csetemp114 + csetemp115 + 
      csetemp116 + csetemp117)*xformL13 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL23 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL33) + xformL20*((csetemp106 + csetemp107 
      + csetemp108 + csetemp109)*xformL03 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL13 + (csetemp126 + csetemp127 + 
      csetemp128 + csetemp129)*xformL23 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL33) + xformL30*((csetemp110 + csetemp111 
      + csetemp112 + csetemp113)*xformL03 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL13 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL23 + (csetemp134 + csetemp135 + 
      csetemp136 + csetemp137)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4032 = xformL00*((csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL03 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL13 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL23 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL33) + xformL10*((csetemp142 
      + csetemp143 + csetemp144 + csetemp145)*xformL03 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL13 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL23 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL33) + xformL20*((csetemp146 
      + csetemp147 + csetemp148 + csetemp149)*xformL03 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL13 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL23 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL33) + xformL30*((csetemp150 
      + csetemp151 + csetemp152 + csetemp153)*xformL03 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL13 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL23 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4033 = xformL00*((csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL03 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL13 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL23 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL33) + xformL10*((csetemp182 
      + csetemp183 + csetemp184 + csetemp185)*xformL03 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL13 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL23 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL33) + xformL20*((csetemp186 
      + csetemp187 + csetemp188 + csetemp189)*xformL03 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL13 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL23 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL33) + xformL30*((csetemp190 
      + csetemp191 + csetemp192 + csetemp193)*xformL03 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL13 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL23 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4110 = xformL01*((csetemp58 + 
      csetemp59 + csetemp60 + csetemp61)*xformL01 + (csetemp62 + csetemp63 + 
      csetemp64 + csetemp65)*xformL11 + (csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL21 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL31) + xformL11*((csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL01 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL11 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL21 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL31) + xformL21*((csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL01 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL11 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL21 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL31) + xformL31*((csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL01 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL11 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL21 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4111 = xformL01*((csetemp100 + 
      csetemp101 + csetemp98 + csetemp99)*xformL01 + (csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL11 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL21 + (csetemp110 + csetemp111 + 
      csetemp112 + csetemp113)*xformL31) + xformL11*((csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL01 + (csetemp114 + csetemp115 + 
      csetemp116 + csetemp117)*xformL11 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL21 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL31) + xformL21*((csetemp106 + csetemp107 
      + csetemp108 + csetemp109)*xformL01 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL11 + (csetemp126 + csetemp127 + 
      csetemp128 + csetemp129)*xformL21 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL31) + xformL31*((csetemp110 + csetemp111 
      + csetemp112 + csetemp113)*xformL01 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL11 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL21 + (csetemp134 + csetemp135 + 
      csetemp136 + csetemp137)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4112 = xformL01*((csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL01 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL11 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL21 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL31) + xformL11*((csetemp142 
      + csetemp143 + csetemp144 + csetemp145)*xformL01 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL11 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL21 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL31) + xformL21*((csetemp146 
      + csetemp147 + csetemp148 + csetemp149)*xformL01 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL11 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL21 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL31) + xformL31*((csetemp150 
      + csetemp151 + csetemp152 + csetemp153)*xformL01 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL11 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL21 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4113 = xformL01*((csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL01 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL11 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL21 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL31) + xformL11*((csetemp182 
      + csetemp183 + csetemp184 + csetemp185)*xformL01 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL11 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL21 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL31) + xformL21*((csetemp186 
      + csetemp187 + csetemp188 + csetemp189)*xformL01 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL11 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL21 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL31) + xformL31*((csetemp190 
      + csetemp191 + csetemp192 + csetemp193)*xformL01 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL11 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL21 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4120 = xformL01*((csetemp58 + 
      csetemp59 + csetemp60 + csetemp61)*xformL02 + (csetemp62 + csetemp63 + 
      csetemp64 + csetemp65)*xformL12 + (csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL22 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL32) + xformL11*((csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL02 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL12 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL22 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL32) + xformL21*((csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL02 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL12 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL22 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL32) + xformL31*((csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL02 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL12 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL22 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4121 = xformL01*((csetemp100 + 
      csetemp101 + csetemp98 + csetemp99)*xformL02 + (csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL12 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL22 + (csetemp110 + csetemp111 + 
      csetemp112 + csetemp113)*xformL32) + xformL11*((csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL02 + (csetemp114 + csetemp115 + 
      csetemp116 + csetemp117)*xformL12 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL22 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL32) + xformL21*((csetemp106 + csetemp107 
      + csetemp108 + csetemp109)*xformL02 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL12 + (csetemp126 + csetemp127 + 
      csetemp128 + csetemp129)*xformL22 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL32) + xformL31*((csetemp110 + csetemp111 
      + csetemp112 + csetemp113)*xformL02 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL12 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL22 + (csetemp134 + csetemp135 + 
      csetemp136 + csetemp137)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4122 = xformL01*((csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL02 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL12 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL22 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL32) + xformL11*((csetemp142 
      + csetemp143 + csetemp144 + csetemp145)*xformL02 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL12 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL22 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL32) + xformL21*((csetemp146 
      + csetemp147 + csetemp148 + csetemp149)*xformL02 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL12 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL22 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL32) + xformL31*((csetemp150 
      + csetemp151 + csetemp152 + csetemp153)*xformL02 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL12 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL22 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4123 = xformL01*((csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL02 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL12 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL22 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL32) + xformL11*((csetemp182 
      + csetemp183 + csetemp184 + csetemp185)*xformL02 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL12 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL22 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL32) + xformL21*((csetemp186 
      + csetemp187 + csetemp188 + csetemp189)*xformL02 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL12 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL22 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL32) + xformL31*((csetemp190 
      + csetemp191 + csetemp192 + csetemp193)*xformL02 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL12 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL22 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4130 = xformL01*((csetemp58 + 
      csetemp59 + csetemp60 + csetemp61)*xformL03 + (csetemp62 + csetemp63 + 
      csetemp64 + csetemp65)*xformL13 + (csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL23 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL33) + xformL11*((csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL03 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL13 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL23 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL33) + xformL21*((csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL03 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL13 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL23 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL33) + xformL31*((csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL03 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL13 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL23 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4131 = xformL01*((csetemp100 + 
      csetemp101 + csetemp98 + csetemp99)*xformL03 + (csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL13 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL23 + (csetemp110 + csetemp111 + 
      csetemp112 + csetemp113)*xformL33) + xformL11*((csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL03 + (csetemp114 + csetemp115 + 
      csetemp116 + csetemp117)*xformL13 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL23 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL33) + xformL21*((csetemp106 + csetemp107 
      + csetemp108 + csetemp109)*xformL03 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL13 + (csetemp126 + csetemp127 + 
      csetemp128 + csetemp129)*xformL23 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL33) + xformL31*((csetemp110 + csetemp111 
      + csetemp112 + csetemp113)*xformL03 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL13 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL23 + (csetemp134 + csetemp135 + 
      csetemp136 + csetemp137)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4132 = xformL01*((csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL03 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL13 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL23 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL33) + xformL11*((csetemp142 
      + csetemp143 + csetemp144 + csetemp145)*xformL03 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL13 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL23 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL33) + xformL21*((csetemp146 
      + csetemp147 + csetemp148 + csetemp149)*xformL03 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL13 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL23 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL33) + xformL31*((csetemp150 
      + csetemp151 + csetemp152 + csetemp153)*xformL03 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL13 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL23 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4133 = xformL01*((csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL03 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL13 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL23 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL33) + xformL11*((csetemp182 
      + csetemp183 + csetemp184 + csetemp185)*xformL03 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL13 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL23 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL33) + xformL21*((csetemp186 
      + csetemp187 + csetemp188 + csetemp189)*xformL03 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL13 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL23 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL33) + xformL31*((csetemp190 
      + csetemp191 + csetemp192 + csetemp193)*xformL03 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL13 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL23 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4220 = xformL02*((csetemp58 + 
      csetemp59 + csetemp60 + csetemp61)*xformL02 + (csetemp62 + csetemp63 + 
      csetemp64 + csetemp65)*xformL12 + (csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL22 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL32) + xformL12*((csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL02 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL12 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL22 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL32) + xformL22*((csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL02 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL12 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL22 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL32) + xformL32*((csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL02 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL12 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL22 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4221 = xformL02*((csetemp100 + 
      csetemp101 + csetemp98 + csetemp99)*xformL02 + (csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL12 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL22 + (csetemp110 + csetemp111 + 
      csetemp112 + csetemp113)*xformL32) + xformL12*((csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL02 + (csetemp114 + csetemp115 + 
      csetemp116 + csetemp117)*xformL12 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL22 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL32) + xformL22*((csetemp106 + csetemp107 
      + csetemp108 + csetemp109)*xformL02 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL12 + (csetemp126 + csetemp127 + 
      csetemp128 + csetemp129)*xformL22 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL32) + xformL32*((csetemp110 + csetemp111 
      + csetemp112 + csetemp113)*xformL02 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL12 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL22 + (csetemp134 + csetemp135 + 
      csetemp136 + csetemp137)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4222 = xformL02*((csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL02 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL12 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL22 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL32) + xformL12*((csetemp142 
      + csetemp143 + csetemp144 + csetemp145)*xformL02 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL12 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL22 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL32) + xformL22*((csetemp146 
      + csetemp147 + csetemp148 + csetemp149)*xformL02 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL12 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL22 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL32) + xformL32*((csetemp150 
      + csetemp151 + csetemp152 + csetemp153)*xformL02 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL12 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL22 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4223 = xformL02*((csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL02 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL12 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL22 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL32) + xformL12*((csetemp182 
      + csetemp183 + csetemp184 + csetemp185)*xformL02 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL12 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL22 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL32) + xformL22*((csetemp186 
      + csetemp187 + csetemp188 + csetemp189)*xformL02 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL12 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL22 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL32) + xformL32*((csetemp190 
      + csetemp191 + csetemp192 + csetemp193)*xformL02 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL12 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL22 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4230 = xformL02*((csetemp58 + 
      csetemp59 + csetemp60 + csetemp61)*xformL03 + (csetemp62 + csetemp63 + 
      csetemp64 + csetemp65)*xformL13 + (csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL23 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL33) + xformL12*((csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL03 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL13 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL23 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL33) + xformL22*((csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL03 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL13 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL23 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL33) + xformL32*((csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL03 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL13 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL23 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4231 = xformL02*((csetemp100 + 
      csetemp101 + csetemp98 + csetemp99)*xformL03 + (csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL13 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL23 + (csetemp110 + csetemp111 + 
      csetemp112 + csetemp113)*xformL33) + xformL12*((csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL03 + (csetemp114 + csetemp115 + 
      csetemp116 + csetemp117)*xformL13 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL23 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL33) + xformL22*((csetemp106 + csetemp107 
      + csetemp108 + csetemp109)*xformL03 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL13 + (csetemp126 + csetemp127 + 
      csetemp128 + csetemp129)*xformL23 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL33) + xformL32*((csetemp110 + csetemp111 
      + csetemp112 + csetemp113)*xformL03 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL13 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL23 + (csetemp134 + csetemp135 + 
      csetemp136 + csetemp137)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4232 = xformL02*((csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL03 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL13 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL23 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL33) + xformL12*((csetemp142 
      + csetemp143 + csetemp144 + csetemp145)*xformL03 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL13 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL23 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL33) + xformL22*((csetemp146 
      + csetemp147 + csetemp148 + csetemp149)*xformL03 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL13 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL23 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL33) + xformL32*((csetemp150 
      + csetemp151 + csetemp152 + csetemp153)*xformL03 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL13 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL23 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4233 = xformL02*((csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL03 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL13 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL23 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL33) + xformL12*((csetemp182 
      + csetemp183 + csetemp184 + csetemp185)*xformL03 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL13 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL23 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL33) + xformL22*((csetemp186 
      + csetemp187 + csetemp188 + csetemp189)*xformL03 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL13 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL23 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL33) + xformL32*((csetemp190 
      + csetemp191 + csetemp192 + csetemp193)*xformL03 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL13 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL23 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4330 = xformL03*((csetemp58 + 
      csetemp59 + csetemp60 + csetemp61)*xformL03 + (csetemp62 + csetemp63 + 
      csetemp64 + csetemp65)*xformL13 + (csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL23 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL33) + xformL13*((csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL03 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL13 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL23 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL33) + xformL23*((csetemp66 + csetemp67 + csetemp68 + 
      csetemp69)*xformL03 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL13 + (csetemp86 + csetemp87 + csetemp88 + 
      csetemp89)*xformL23 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL33) + xformL33*((csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL03 + (csetemp82 + csetemp83 + csetemp84 + 
      csetemp85)*xformL13 + (csetemp90 + csetemp91 + csetemp92 + 
      csetemp93)*xformL23 + (csetemp94 + csetemp95 + csetemp96 + 
      csetemp97)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4331 = xformL03*((csetemp100 + 
      csetemp101 + csetemp98 + csetemp99)*xformL03 + (csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL13 + (csetemp106 + csetemp107 + 
      csetemp108 + csetemp109)*xformL23 + (csetemp110 + csetemp111 + 
      csetemp112 + csetemp113)*xformL33) + xformL13*((csetemp102 + csetemp103 
      + csetemp104 + csetemp105)*xformL03 + (csetemp114 + csetemp115 + 
      csetemp116 + csetemp117)*xformL13 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL23 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL33) + xformL23*((csetemp106 + csetemp107 
      + csetemp108 + csetemp109)*xformL03 + (csetemp118 + csetemp119 + 
      csetemp120 + csetemp121)*xformL13 + (csetemp126 + csetemp127 + 
      csetemp128 + csetemp129)*xformL23 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL33) + xformL33*((csetemp110 + csetemp111 
      + csetemp112 + csetemp113)*xformL03 + (csetemp122 + csetemp123 + 
      csetemp124 + csetemp125)*xformL13 + (csetemp130 + csetemp131 + 
      csetemp132 + csetemp133)*xformL23 + (csetemp134 + csetemp135 + 
      csetemp136 + csetemp137)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4332 = xformL03*((csetemp138 + 
      csetemp139 + csetemp140 + csetemp141)*xformL03 + (csetemp142 + 
      csetemp143 + csetemp144 + csetemp145)*xformL13 + (csetemp146 + 
      csetemp147 + csetemp148 + csetemp149)*xformL23 + (csetemp150 + 
      csetemp151 + csetemp152 + csetemp153)*xformL33) + xformL13*((csetemp142 
      + csetemp143 + csetemp144 + csetemp145)*xformL03 + (csetemp154 + 
      csetemp155 + csetemp156 + csetemp157)*xformL13 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL23 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL33) + xformL23*((csetemp146 
      + csetemp147 + csetemp148 + csetemp149)*xformL03 + (csetemp158 + 
      csetemp159 + csetemp160 + csetemp161)*xformL13 + (csetemp166 + 
      csetemp167 + csetemp168 + csetemp169)*xformL23 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL33) + xformL33*((csetemp150 
      + csetemp151 + csetemp152 + csetemp153)*xformL03 + (csetemp162 + 
      csetemp163 + csetemp164 + csetemp165)*xformL13 + (csetemp170 + 
      csetemp171 + csetemp172 + csetemp173)*xformL23 + (csetemp174 + 
      csetemp175 + csetemp176 + csetemp177)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4333 = xformL03*((csetemp178 + 
      csetemp179 + csetemp180 + csetemp181)*xformL03 + (csetemp182 + 
      csetemp183 + csetemp184 + csetemp185)*xformL13 + (csetemp186 + 
      csetemp187 + csetemp188 + csetemp189)*xformL23 + (csetemp190 + 
      csetemp191 + csetemp192 + csetemp193)*xformL33) + xformL13*((csetemp182 
      + csetemp183 + csetemp184 + csetemp185)*xformL03 + (csetemp194 + 
      csetemp195 + csetemp196 + csetemp197)*xformL13 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL23 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL33) + xformL23*((csetemp186 
      + csetemp187 + csetemp188 + csetemp189)*xformL03 + (csetemp198 + 
      csetemp199 + csetemp200 + csetemp201)*xformL13 + (csetemp206 + 
      csetemp207 + csetemp208 + csetemp209)*xformL23 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL33) + xformL33*((csetemp190 
      + csetemp191 + csetemp192 + csetemp193)*xformL03 + (csetemp202 + 
      csetemp203 + csetemp204 + csetemp205)*xformL13 + (csetemp210 + 
      csetemp211 + csetemp212 + csetemp213)*xformL23 + (csetemp214 + 
      csetemp215 + csetemp216 + csetemp217)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal1 = g401;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal2 = g402;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal3 = g403;
    
    gxxL = g411;
    
    gxyL = g412;
    
    gxzL = g413;
    
    gyyL = g422;
    
    gyzL = g423;
    
    gzzL = g433;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp218 = SQR(gxzL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp219 = SQR(gyzL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp220 = SQR(gxyL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED detg = 2*gxyL*gxzL*gyzL + 
      gyyL*(gxxL*gzzL - csetemp218) - gxxL*csetemp219 - 
      gzzL*csetemp220;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp221 = INV(detg);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu11 = (gyyL*gzzL - 
      csetemp219)*csetemp221;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu12 = (gxzL*gyzL - 
      gxyL*gzzL)*csetemp221;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu13 = (-(gxzL*gyyL) + 
      gxyL*gyzL)*csetemp221;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu22 = (gxxL*gzzL - 
      csetemp218)*csetemp221;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu23 = (gxyL*gxzL - 
      gxxL*gyzL)*csetemp221;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu33 = (gxxL*gyyL - 
      csetemp220)*csetemp221;
    
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp222 = dtg11*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp223 = dtg12*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp224 = dtg13*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp225 = dtg12*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp226 = dtg22*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp227 = dtg23*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp228 = dtg13*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp229 = dtg23*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp230 = dtg33*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu11 = -((csetemp222 + csetemp223 + 
      csetemp224)*gu11) - (csetemp225 + csetemp226 + csetemp227)*gu12 - 
      (csetemp228 + csetemp229 + csetemp230)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu12 = -((csetemp222 + csetemp223 + 
      csetemp224)*gu12) - (csetemp225 + csetemp226 + csetemp227)*gu22 - 
      (csetemp228 + csetemp229 + csetemp230)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu13 = -((csetemp222 + csetemp223 + 
      csetemp224)*gu13) - (csetemp225 + csetemp226 + csetemp227)*gu23 - 
      (csetemp228 + csetemp229 + csetemp230)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp231 = dtg11*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp232 = dtg12*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp233 = dtg13*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp234 = dtg22*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp235 = dtg23*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp236 = dtg13*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp237 = dtg23*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp238 = dtg33*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu22 = -((csetemp231 + csetemp232 + 
      csetemp233)*gu12) - (csetemp223 + csetemp234 + csetemp235)*gu22 - 
      (csetemp236 + csetemp237 + csetemp238)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu23 = -((csetemp231 + csetemp232 + 
      csetemp233)*gu13) - (csetemp223 + csetemp234 + csetemp235)*gu23 - 
      (csetemp236 + csetemp237 + csetemp238)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu33 = -(gu13*(dtg11*gu13 + 
      dtg12*gu23 + dtg13*gu33)) - gu23*(dtg12*gu13 + dtg22*gu23 + dtg23*gu33) 
      - gu33*(csetemp224 + csetemp235 + dtg33*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp239 = dg111*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp240 = dg121*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp241 = dg131*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp242 = dg121*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp243 = dg221*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp244 = dg231*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp245 = dg131*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp246 = dg231*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp247 = dg331*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu111 = -((csetemp239 + csetemp240 + 
      csetemp241)*gu11) - (csetemp242 + csetemp243 + csetemp244)*gu12 - 
      (csetemp245 + csetemp246 + csetemp247)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu121 = -((csetemp239 + csetemp240 + 
      csetemp241)*gu12) - (csetemp242 + csetemp243 + csetemp244)*gu22 - 
      (csetemp245 + csetemp246 + csetemp247)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu131 = -((csetemp239 + csetemp240 + 
      csetemp241)*gu13) - (csetemp242 + csetemp243 + csetemp244)*gu23 - 
      (csetemp245 + csetemp246 + csetemp247)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp248 = dg111*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp249 = dg121*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp250 = dg131*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp251 = dg221*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp252 = dg231*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp253 = dg131*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp254 = dg231*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp255 = dg331*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu221 = -((csetemp248 + csetemp249 + 
      csetemp250)*gu12) - (csetemp240 + csetemp251 + csetemp252)*gu22 - 
      (csetemp253 + csetemp254 + csetemp255)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu231 = -((csetemp248 + csetemp249 + 
      csetemp250)*gu13) - (csetemp240 + csetemp251 + csetemp252)*gu23 - 
      (csetemp253 + csetemp254 + csetemp255)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu331 = -(gu13*(dg111*gu13 + 
      dg121*gu23 + dg131*gu33)) - gu23*(dg121*gu13 + dg221*gu23 + dg231*gu33) 
      - gu33*(csetemp241 + csetemp252 + dg331*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp256 = dg112*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp257 = dg122*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp258 = dg132*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp259 = dg122*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp260 = dg222*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp261 = dg232*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp262 = dg132*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp263 = dg232*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp264 = dg332*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu112 = -((csetemp256 + csetemp257 + 
      csetemp258)*gu11) - (csetemp259 + csetemp260 + csetemp261)*gu12 - 
      (csetemp262 + csetemp263 + csetemp264)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu122 = -((csetemp256 + csetemp257 + 
      csetemp258)*gu12) - (csetemp259 + csetemp260 + csetemp261)*gu22 - 
      (csetemp262 + csetemp263 + csetemp264)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu132 = -((csetemp256 + csetemp257 + 
      csetemp258)*gu13) - (csetemp259 + csetemp260 + csetemp261)*gu23 - 
      (csetemp262 + csetemp263 + csetemp264)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp265 = dg112*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp266 = dg122*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp267 = dg132*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp268 = dg222*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp269 = dg232*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp270 = dg132*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp271 = dg232*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp272 = dg332*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu222 = -((csetemp265 + csetemp266 + 
      csetemp267)*gu12) - (csetemp257 + csetemp268 + csetemp269)*gu22 - 
      (csetemp270 + csetemp271 + csetemp272)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu232 = -((csetemp265 + csetemp266 + 
      csetemp267)*gu13) - (csetemp257 + csetemp268 + csetemp269)*gu23 - 
      (csetemp270 + csetemp271 + csetemp272)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu332 = -(gu13*(dg112*gu13 + 
      dg122*gu23 + dg132*gu33)) - gu23*(dg122*gu13 + dg222*gu23 + dg232*gu33) 
      - gu33*(csetemp258 + csetemp269 + dg332*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp273 = dg113*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp274 = dg123*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp275 = dg133*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp276 = dg123*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp277 = dg223*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp278 = dg233*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp279 = dg133*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp280 = dg233*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp281 = dg333*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu113 = -((csetemp273 + csetemp274 + 
      csetemp275)*gu11) - (csetemp276 + csetemp277 + csetemp278)*gu12 - 
      (csetemp279 + csetemp280 + csetemp281)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu123 = -((csetemp273 + csetemp274 + 
      csetemp275)*gu12) - (csetemp276 + csetemp277 + csetemp278)*gu22 - 
      (csetemp279 + csetemp280 + csetemp281)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu133 = -((csetemp273 + csetemp274 + 
      csetemp275)*gu13) - (csetemp276 + csetemp277 + csetemp278)*gu23 - 
      (csetemp279 + csetemp280 + csetemp281)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp282 = dg113*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp283 = dg123*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp284 = dg133*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp285 = dg223*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp286 = dg233*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp287 = dg133*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp288 = dg233*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp289 = dg333*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu223 = -((csetemp282 + csetemp283 + 
      csetemp284)*gu12) - (csetemp274 + csetemp285 + csetemp286)*gu22 - 
      (csetemp287 + csetemp288 + csetemp289)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu233 = -((csetemp282 + csetemp283 + 
      csetemp284)*gu13) - (csetemp274 + csetemp285 + csetemp286)*gu23 - 
      (csetemp287 + csetemp288 + csetemp289)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu333 = -(gu13*(dg113*gu13 + 
      dg123*gu23 + dg133*gu33)) - gu23*(dg123*gu13 + dg223*gu23 + dg233*gu33) 
      - gu33*(csetemp275 + csetemp286 + dg333*gu33);
    
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp290 = INV(alpL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtalpL = 0.5*csetemp290*(-dg4000 + 
      dtbetasq);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxxL = 
      0.5*csetemp290*(2*(gxxL*dbeta11 + gxyL*dbeta21 + gxzL*dbeta31) + 
      betaxL*dg111 + betayL*dg112 + betazL*dg113 - dtg11);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxyL = 0.5*csetemp290*(gxxL*dbeta12 
      + gyyL*dbeta21 + gxyL*(dbeta11 + dbeta22) + gyzL*dbeta31 + 
      gxzL*dbeta32 + betaxL*dg121 + betayL*dg122 + betazL*dg123 - 
      dtg12);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxzL = 0.5*csetemp290*(gxxL*dbeta13 
      + gyzL*dbeta21 + gxyL*dbeta23 + gzzL*dbeta31 + gxzL*(dbeta11 + 
      dbeta33) + betaxL*dg131 + betayL*dg132 + betazL*dg133 - dtg13);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kyyL = 
      0.5*csetemp290*(2*(gxyL*dbeta12 + gyyL*dbeta22 + gyzL*dbeta32) + 
      betaxL*dg221 + betayL*dg222 + betazL*dg223 - dtg22);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kyzL = 0.5*csetemp290*(gxzL*dbeta12 
      + gxyL*dbeta13 + gyyL*dbeta23 + gzzL*dbeta32 + gyzL*(dbeta22 + 
      dbeta33) + betaxL*dg231 + betayL*dg232 + betazL*dg233 - dtg23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kzzL = 
      0.5*csetemp290*(2*(gxzL*dbeta13 + gyzL*dbeta23 + gzzL*dbeta33) + 
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
  CCTK_ENDLOOP3(Minkowski_always);
}

extern "C" void Minkowski_always(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering Minkowski_always_Body");
  }
  
  if (cctk_iteration % Minkowski_always_calc_every != Minkowski_always_calc_offset)
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
  GenericFD_AssertGroupStorage(cctkGH, "Minkowski_always", 6, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, Minkowski_always_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving Minkowski_always_Body");
  }
}
