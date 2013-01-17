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
  ptrdiff_t /*const*/ di CCTK_ATTRIBUTE_UNUSED  = 1;
  ptrdiff_t /*const*/ dj CCTK_ATTRIBUTE_UNUSED  = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t /*const*/ dk CCTK_ATTRIBUTE_UNUSED  = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t /*const*/ cdi CCTK_ATTRIBUTE_UNUSED  = sizeof(CCTK_REAL) * di;
  ptrdiff_t /*const*/ cdj CCTK_ATTRIBUTE_UNUSED  = sizeof(CCTK_REAL) * dj;
  ptrdiff_t /*const*/ cdk CCTK_ATTRIBUTE_UNUSED  = sizeof(CCTK_REAL) * dk;
  CCTK_REAL /*const*/ dx CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL /*const*/ dy CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL /*const*/ dz CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL /*const*/ dt CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL /*const*/ t CCTK_ATTRIBUTE_UNUSED  = ToReal(cctk_time);
  CCTK_REAL /*const*/ dxi CCTK_ATTRIBUTE_UNUSED  = INV(dx);
  CCTK_REAL /*const*/ dyi CCTK_ATTRIBUTE_UNUSED  = INV(dy);
  CCTK_REAL /*const*/ dzi CCTK_ATTRIBUTE_UNUSED  = INV(dz);
  CCTK_REAL /*const*/ khalf CCTK_ATTRIBUTE_UNUSED  = 0.5;
  CCTK_REAL /*const*/ kthird CCTK_ATTRIBUTE_UNUSED  = 1/3.0;
  CCTK_REAL /*const*/ ktwothird CCTK_ATTRIBUTE_UNUSED  = 2.0/3.0;
  CCTK_REAL /*const*/ kfourthird CCTK_ATTRIBUTE_UNUSED  = 4.0/3.0;
  CCTK_REAL /*const*/ keightthird CCTK_ATTRIBUTE_UNUSED  = 8.0/3.0;
  CCTK_REAL /*const*/ hdxi CCTK_ATTRIBUTE_UNUSED  = 0.5 * dxi;
  CCTK_REAL /*const*/ hdyi CCTK_ATTRIBUTE_UNUSED  = 0.5 * dyi;
  CCTK_REAL /*const*/ hdzi CCTK_ATTRIBUTE_UNUSED  = 0.5 * dzi;
  
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
    ptrdiff_t /*const*/ index CCTK_ATTRIBUTE_UNUSED  = di*i + dj*j + dk*k;
    
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
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L00 = 1.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L01 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L02 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L03 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L10 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp0 = cos(ToReal(phi));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp1 = cos(ToReal(psi));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp2 = cos(ToReal(theta));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp3 = sin(ToReal(phi));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp4 = sin(ToReal(psi));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L11 = csetemp0*csetemp1 - 
      1.*csetemp2*csetemp3*csetemp4;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L12 = csetemp1*csetemp3 + 
      csetemp0*csetemp2*csetemp4;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp5 = sin(ToReal(theta));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L13 = csetemp4*csetemp5;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L20 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L21 = 
      -1.*(csetemp1*csetemp2*csetemp3 + csetemp0*csetemp4);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L22 = csetemp0*csetemp1*csetemp2 
      - 1.*csetemp3*csetemp4;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L23 = csetemp1*csetemp5;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L30 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L31 = csetemp3*csetemp5;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L32 = -1.*csetemp0*csetemp5;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform1L33 = csetemp2;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp6 = SQR(ToReal(boostx));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp7 = SQR(ToReal(boosty));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp8 = SQR(ToReal(boostz));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp9 = INV(ToReal(lapsefactor));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L00 = csetemp9*INV((-1. + 
      csetemp6 + csetemp7 + csetemp8)*(1. + sqrt(1. - 1.*(csetemp6 + csetemp7 
      + csetemp8))))*(1. - 1.*(csetemp6 + csetemp7 + csetemp8) + sqrt(1. - 
      1.*(csetemp6 + csetemp7 + csetemp8)))*(-1. + 
      ToReal(boostx)*ToReal(shiftaddx) + ToReal(boosty)*ToReal(shiftaddy) + 
      ToReal(boostz)*ToReal(shiftaddz));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L01 = csetemp9*INV((-1. + 
      csetemp6 + csetemp7 + csetemp8)*(1. + sqrt(1. - 1.*(csetemp6 + csetemp7 
      + csetemp8))))*(1. - 1.*(csetemp6 + csetemp7 + csetemp8) + sqrt(1. - 
      1.*(csetemp6 + csetemp7 + csetemp8)))*ToReal(boostx);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L02 = csetemp9*INV((-1. + 
      csetemp6 + csetemp7 + csetemp8)*(1. + sqrt(1. - 1.*(csetemp6 + csetemp7 
      + csetemp8))))*(1. - 1.*(csetemp6 + csetemp7 + csetemp8) + sqrt(1. - 
      1.*(csetemp6 + csetemp7 + csetemp8)))*ToReal(boosty);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L03 = csetemp9*INV((-1. + 
      csetemp6 + csetemp7 + csetemp8)*(1. + sqrt(1. - 1.*(csetemp6 + csetemp7 
      + csetemp8))))*(1. - 1.*(csetemp6 + csetemp7 + csetemp8) + sqrt(1. - 
      1.*(csetemp6 + csetemp7 + csetemp8)))*ToReal(boostz);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L10 = INV((-1. + csetemp6 + 
      csetemp7 + csetemp8)*(1. + sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8))))*((csetemp6 + (-1. + csetemp7 + csetemp8)*(1. + sqrt(1. - 
      1.*(csetemp6 + csetemp7 + csetemp8))))*ToReal(shiftaddx) - 
      1.*ToReal(boostx)*(-1. + csetemp6 + csetemp7 + csetemp8 + sqrt(1. - 
      1.*(csetemp6 + csetemp7 + csetemp8))*(-1. + 
      ToReal(boosty)*ToReal(shiftaddy) + ToReal(boostz)*ToReal(shiftaddz))));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L11 = INV((-1. + csetemp6 + 
      csetemp7 + csetemp8)*(1. + sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8))))*(csetemp6 + (-1. + csetemp7 + csetemp8)*(1. + sqrt(1. - 
      1.*(csetemp6 + csetemp7 + csetemp8))));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L12 = -1.*INV(-1. + csetemp6 + 
      csetemp7 + csetemp8 - 1.*sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8)))*ToReal(boostx)*ToReal(boosty);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L13 = -1.*INV(-1. + csetemp6 + 
      csetemp7 + csetemp8 - 1.*sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8)))*ToReal(boostx)*ToReal(boostz);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L20 = INV((-1. + csetemp6 + 
      csetemp7 + csetemp8)*(1. + sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8))))*((csetemp7 + (-1. + csetemp6 + csetemp8)*(1. + sqrt(1. - 
      1.*(csetemp6 + csetemp7 + csetemp8))))*ToReal(shiftaddy) - 
      1.*ToReal(boosty)*(-1. + csetemp6 + csetemp7 + csetemp8 + sqrt(1. - 
      1.*(csetemp6 + csetemp7 + csetemp8))*(-1. + 
      ToReal(boostx)*ToReal(shiftaddx) + ToReal(boostz)*ToReal(shiftaddz))));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L21 = -1.*INV(-1. + csetemp6 + 
      csetemp7 + csetemp8 - 1.*sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8)))*ToReal(boostx)*ToReal(boosty);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L22 = INV((-1. + csetemp6 + 
      csetemp7 + csetemp8)*(1. + sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8))))*(csetemp7 + (-1. + csetemp6 + csetemp8)*(1. + sqrt(1. - 
      1.*(csetemp6 + csetemp7 + csetemp8))));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L23 = -1.*INV(-1. + csetemp6 + 
      csetemp7 + csetemp8 - 1.*sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8)))*ToReal(boosty)*ToReal(boostz);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L30 = INV((-1. + csetemp6 + 
      csetemp7 + csetemp8)*(1. + sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8))))*(-1.*ToReal(boostz)*(-1. + csetemp6 + csetemp7 + csetemp8 
      + sqrt(1. - 1.*(csetemp6 + csetemp7 + csetemp8))*(-1. + 
      ToReal(boostx)*ToReal(shiftaddx) + ToReal(boosty)*ToReal(shiftaddy))) + 
      (-1. + csetemp6 + csetemp7 + csetemp8 + (-1. + csetemp6 + 
      csetemp7)*sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8)))*ToReal(shiftaddz));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L31 = -1.*INV(-1. + csetemp6 + 
      csetemp7 + csetemp8 - 1.*sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8)))*ToReal(boostx)*ToReal(boostz);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L32 = -1.*INV(-1. + csetemp6 + 
      csetemp7 + csetemp8 - 1.*sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8)))*ToReal(boosty)*ToReal(boostz);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xform2L33 = INV((-1. + csetemp6 + 
      csetemp7 + csetemp8)*(1. + sqrt(1. - 1.*(csetemp6 + csetemp7 + 
      csetemp8))))*(-1. + csetemp6 + csetemp7 + csetemp8 + (-1. + csetemp6 + 
      csetemp7)*sqrt(1. - 1.*(csetemp6 + csetemp7 + csetemp8)));
    
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xx0 = t - 1.*ToReal(timeoffset);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xx1 = xL - 1.*ToReal(positionx);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xx2 = yL - 1.*ToReal(positiony);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED xx3 = zL - 1.*ToReal(positionz);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED txx1 = xformL10*xx0 + xformL11*xx1 + 
      xformL12*xx2 + xformL13*xx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED txx2 = xformL20*xx0 + xformL21*xx1 + 
      xformL22*xx2 + xformL23*xx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED txx3 = xformL30*xx0 + xformL31*xx1 + 
      xformL32*xx2 + xformL33*xx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED X = txx1;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED Y = txx2;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED Z = txx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED rXYZ = sqrt(SQR(X) + SQR(Y) + SQR(Z));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg400 = -1.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg401 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg402 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg403 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp10 = 0.5*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp11 = INV(rXYZ);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp12 = INV(ToReal(M));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp13 = SQR(ToReal(M));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp14 = INV(SQR(Pi));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg411 = 0.00390625*QAD(4. + 
      csetemp11*(8.*fmin(csetemp10,rXYZ) + csetemp12*(csetemp13*(csetemp14 - 
      1.*csetemp14*cos(12.566370614359172*csetemp12*fmin(csetemp10,rXYZ))) - 
      8.*SQR(fmin(csetemp10,rXYZ)))));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg412 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg413 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg422 = 0.00390625*QAD(4. + 
      csetemp11*(8.*fmin(csetemp10,rXYZ) + csetemp12*(csetemp13*(csetemp14 - 
      1.*csetemp14*cos(12.566370614359172*csetemp12*fmin(csetemp10,rXYZ))) - 
      8.*SQR(fmin(csetemp10,rXYZ)))));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg423 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg433 = 0.00390625*QAD(4. + 
      csetemp11*(8.*fmin(csetemp10,rXYZ) + csetemp12*(csetemp13*(csetemp14 - 
      1.*csetemp14*cos(12.566370614359172*csetemp12*fmin(csetemp10,rXYZ))) - 
      8.*SQR(fmin(csetemp10,rXYZ)))));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4000 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4001 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4002 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4003 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4010 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4011 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4012 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4013 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4020 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4021 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4022 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4023 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4030 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4031 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4032 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4033 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4110 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp15 = pow(rXYZ,-6.);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp16 = 2.*rXYZ;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp17 = INV(QAD(ToReal(M)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp18 = 0.0001053903916534937;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp19 = SQR(Pi);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp20 = SQR(rXYZ);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4111 = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*X*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*X*CUB(csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) 
      + csetemp19*(-4.*csetemp20 + 
      6.*rXYZ*ToReal(M)))*(4.*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) - 
      6.283185307179586*rXYZ*sin(12.566370614359172*csetemp12*rXYZ)*ToReal(M)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4112 = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Y*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Y*CUB(csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) 
      + csetemp19*(-4.*csetemp20 + 
      6.*rXYZ*ToReal(M)))*(4.*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) - 
      6.283185307179586*rXYZ*sin(12.566370614359172*csetemp12*rXYZ)*ToReal(M)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4113 = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Z*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Z*CUB(csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) 
      + csetemp19*(-4.*csetemp20 + 
      6.*rXYZ*ToReal(M)))*(4.*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) - 
      6.283185307179586*rXYZ*sin(12.566370614359172*csetemp12*rXYZ)*ToReal(M)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4120 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4121 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4122 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4123 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4130 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4131 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4132 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4133 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4220 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4221 = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*X*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*X*CUB(csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) 
      + csetemp19*(-4.*csetemp20 + 
      6.*rXYZ*ToReal(M)))*(4.*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) - 
      6.283185307179586*rXYZ*sin(12.566370614359172*csetemp12*rXYZ)*ToReal(M)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4222 = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Y*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Y*CUB(csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) 
      + csetemp19*(-4.*csetemp20 + 
      6.*rXYZ*ToReal(M)))*(4.*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) - 
      6.283185307179586*rXYZ*sin(12.566370614359172*csetemp12*rXYZ)*ToReal(M)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4223 = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Z*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Z*CUB(csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) 
      + csetemp19*(-4.*csetemp20 + 
      6.*rXYZ*ToReal(M)))*(4.*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) - 
      6.283185307179586*rXYZ*sin(12.566370614359172*csetemp12*rXYZ)*ToReal(M)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4230 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4231 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4232 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4233 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4330 = 0.;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4331 = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*X*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*X*CUB(csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) 
      + csetemp19*(-4.*csetemp20 + 
      6.*rXYZ*ToReal(M)))*(4.*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) - 
      6.283185307179586*rXYZ*sin(12.566370614359172*csetemp12*rXYZ)*ToReal(M)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4332 = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Y*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Y*CUB(csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) 
      + csetemp19*(-4.*csetemp20 + 
      6.*rXYZ*ToReal(M)))*(4.*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) - 
      6.283185307179586*rXYZ*sin(12.566370614359172*csetemp12*rXYZ)*ToReal(M)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4333 = IfThen(ToReal(M) <= 
      csetemp16,-0.25*csetemp15*Z*CUB(csetemp16 + 
      ToReal(M))*ToReal(M),-0.25*csetemp15*csetemp17*csetemp18*Z*CUB(csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) 
      + csetemp19*(-4.*csetemp20 + 
      6.*rXYZ*ToReal(M)))*(4.*csetemp19*csetemp20 + 
      csetemp13*SQR(sin(6.283185307179586*csetemp12*rXYZ)) - 
      6.283185307179586*rXYZ*sin(12.566370614359172*csetemp12*rXYZ)*ToReal(M)));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g400 = 2.*(tg423*xformL20*xformL30 + 
      xformL00*(tg401*xformL10 + tg402*xformL20 + tg403*xformL30) + 
      xformL10*(tg412*xformL20 + tg413*xformL30)) + tg400*SQR(xformL00) + 
      tg411*SQR(xformL10) + tg422*SQR(xformL20) + tg433*SQR(xformL30);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp21 = tg400*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp22 = tg401*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp23 = tg402*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp24 = tg403*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp25 = tg401*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp26 = tg411*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp27 = tg412*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp28 = tg413*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp29 = tg402*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp30 = tg412*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp31 = tg422*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp32 = tg423*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp33 = tg403*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp34 = tg413*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp35 = tg423*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp36 = tg433*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g401 = (csetemp21 + csetemp22 + 
      csetemp23 + csetemp24)*xformL00 + (csetemp25 + csetemp26 + csetemp27 + 
      csetemp28)*xformL10 + (csetemp29 + csetemp30 + csetemp31 + 
      csetemp32)*xformL20 + (csetemp33 + csetemp34 + csetemp35 + 
      csetemp36)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp37 = tg400*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp38 = tg401*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp39 = tg402*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp40 = tg403*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp41 = tg401*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp42 = tg411*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp43 = tg412*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp44 = tg413*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp45 = tg402*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp46 = tg412*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp47 = tg422*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp48 = tg423*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp49 = tg403*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp50 = tg413*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp51 = tg423*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp52 = tg433*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g402 = (csetemp37 + csetemp38 + 
      csetemp39 + csetemp40)*xformL00 + (csetemp41 + csetemp42 + csetemp43 + 
      csetemp44)*xformL10 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL20 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp53 = tg400*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp54 = tg401*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp55 = tg402*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp56 = tg403*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp57 = tg401*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp58 = tg411*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp59 = tg412*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp60 = tg413*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp61 = tg402*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp62 = tg412*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp63 = tg422*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp64 = tg423*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp65 = tg403*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp66 = tg413*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp67 = tg423*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp68 = tg433*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g403 = (csetemp53 + csetemp54 + 
      csetemp55 + csetemp56)*xformL00 + (csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*xformL10 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL20 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g411 = (csetemp21 + csetemp22 + 
      csetemp23 + csetemp24)*xformL01 + (csetemp25 + csetemp26 + csetemp27 + 
      csetemp28)*xformL11 + (csetemp29 + csetemp30 + csetemp31 + 
      csetemp32)*xformL21 + (csetemp33 + csetemp34 + csetemp35 + 
      csetemp36)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g412 = (csetemp37 + csetemp38 + 
      csetemp39 + csetemp40)*xformL01 + (csetemp41 + csetemp42 + csetemp43 + 
      csetemp44)*xformL11 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL21 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g413 = (csetemp53 + csetemp54 + 
      csetemp55 + csetemp56)*xformL01 + (csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*xformL11 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL21 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g422 = (csetemp37 + csetemp38 + 
      csetemp39 + csetemp40)*xformL02 + (csetemp41 + csetemp42 + csetemp43 + 
      csetemp44)*xformL12 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL22 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g423 = (csetemp53 + csetemp54 + 
      csetemp55 + csetemp56)*xformL02 + (csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*xformL12 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL22 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g433 = (csetemp53 + csetemp54 + 
      csetemp55 + csetemp56)*xformL03 + (csetemp57 + csetemp58 + csetemp59 + 
      csetemp60)*xformL13 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL23 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp69 = tdg4000*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp70 = tdg4001*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp71 = tdg4002*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp72 = tdg4003*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp73 = tdg4010*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp74 = tdg4011*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp75 = tdg4012*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp76 = tdg4013*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp77 = tdg4020*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp78 = tdg4021*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp79 = tdg4022*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp80 = tdg4023*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp81 = tdg4030*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp82 = tdg4031*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp83 = tdg4032*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp84 = tdg4033*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp85 = tdg4110*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp86 = tdg4111*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp87 = tdg4112*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp88 = tdg4113*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp89 = tdg4120*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp90 = tdg4121*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp91 = tdg4122*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp92 = tdg4123*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp93 = tdg4130*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp94 = tdg4131*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp95 = tdg4132*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp96 = tdg4133*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp97 = tdg4220*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp98 = tdg4221*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp99 = tdg4222*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp100 = tdg4223*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp101 = tdg4230*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp102 = tdg4231*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp103 = tdg4232*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp104 = tdg4233*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp105 = tdg4330*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp106 = tdg4331*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp107 = tdg4332*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp108 = tdg4333*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4000 = xformL20*((csetemp77 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4010 = xformL20*((csetemp77 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp109 = tdg4000*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp110 = tdg4001*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp111 = tdg4002*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp112 = tdg4003*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp113 = tdg4010*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp114 = tdg4011*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp115 = tdg4012*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp116 = tdg4013*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp117 = tdg4020*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp118 = tdg4021*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp119 = tdg4022*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp120 = tdg4023*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp121 = tdg4030*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp122 = tdg4031*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp123 = tdg4032*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp124 = tdg4033*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp125 = tdg4110*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp126 = tdg4111*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp127 = tdg4112*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp128 = tdg4113*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp129 = tdg4120*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp130 = tdg4121*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp131 = tdg4122*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp132 = tdg4123*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp133 = tdg4130*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp134 = tdg4131*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp135 = tdg4132*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp136 = tdg4133*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp137 = tdg4220*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp138 = tdg4221*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp139 = tdg4222*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp140 = tdg4223*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp141 = tdg4230*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp142 = tdg4231*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp143 = tdg4232*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp144 = tdg4233*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp145 = tdg4330*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp146 = tdg4331*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp147 = tdg4332*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp148 = tdg4333*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4011 = xformL00*((csetemp109 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp149 = tdg4000*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp150 = tdg4001*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp151 = tdg4002*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp152 = tdg4003*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp153 = tdg4010*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp154 = tdg4011*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp155 = tdg4012*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp156 = tdg4013*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp157 = tdg4020*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp158 = tdg4021*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp159 = tdg4022*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp160 = tdg4023*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp161 = tdg4030*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp162 = tdg4031*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp163 = tdg4032*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp164 = tdg4033*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp165 = tdg4110*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp166 = tdg4111*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp167 = tdg4112*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp168 = tdg4113*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp169 = tdg4120*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp170 = tdg4121*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp171 = tdg4122*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp172 = tdg4123*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp173 = tdg4130*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp174 = tdg4131*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp175 = tdg4132*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp176 = tdg4133*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp177 = tdg4220*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp178 = tdg4221*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp179 = tdg4222*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp180 = tdg4223*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp181 = tdg4230*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp182 = tdg4231*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp183 = tdg4232*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp184 = tdg4233*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp185 = tdg4330*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp186 = tdg4331*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp187 = tdg4332*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp188 = tdg4333*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4012 = xformL00*((csetemp149 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp189 = tdg4000*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp190 = tdg4001*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp191 = tdg4002*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp192 = tdg4003*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp193 = tdg4010*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp194 = tdg4011*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp195 = tdg4012*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp196 = tdg4013*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp197 = tdg4020*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp198 = tdg4021*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp199 = tdg4022*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp200 = tdg4023*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp201 = tdg4030*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp202 = tdg4031*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp203 = tdg4032*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp204 = tdg4033*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp205 = tdg4110*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp206 = tdg4111*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp207 = tdg4112*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp208 = tdg4113*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp209 = tdg4120*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp210 = tdg4121*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp211 = tdg4122*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp212 = tdg4123*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp213 = tdg4130*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp214 = tdg4131*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp215 = tdg4132*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp216 = tdg4133*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp217 = tdg4220*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp218 = tdg4221*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp219 = tdg4222*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp220 = tdg4223*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp221 = tdg4230*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp222 = tdg4231*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp223 = tdg4232*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp224 = tdg4233*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp225 = tdg4330*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp226 = tdg4331*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp227 = tdg4332*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp228 = tdg4333*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4013 = xformL00*((csetemp189 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4020 = xformL20*((csetemp77 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4021 = xformL00*((csetemp109 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4022 = xformL00*((csetemp149 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4023 = xformL00*((csetemp189 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4030 = xformL20*((csetemp77 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4031 = xformL00*((csetemp109 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4032 = xformL00*((csetemp149 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4033 = xformL00*((csetemp189 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4110 = xformL21*((csetemp77 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4111 = xformL01*((csetemp109 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4112 = xformL01*((csetemp149 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4113 = xformL01*((csetemp189 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4120 = xformL21*((csetemp77 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4121 = xformL01*((csetemp109 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4122 = xformL01*((csetemp149 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4123 = xformL01*((csetemp189 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4130 = xformL21*((csetemp77 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4131 = xformL01*((csetemp109 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4132 = xformL01*((csetemp149 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4133 = xformL01*((csetemp189 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4220 = xformL22*((csetemp77 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4221 = xformL02*((csetemp109 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4222 = xformL02*((csetemp149 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4223 = xformL02*((csetemp189 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4230 = xformL22*((csetemp77 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4231 = xformL02*((csetemp109 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4232 = xformL02*((csetemp149 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4233 = xformL02*((csetemp189 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4330 = xformL23*((csetemp77 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4331 = xformL03*((csetemp109 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4332 = xformL03*((csetemp149 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4333 = xformL03*((csetemp189 + 
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal1 = g401;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal2 = g402;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal3 = g403;
    
    gxxL = g411;
    
    gxyL = g412;
    
    gxzL = g413;
    
    gyyL = g422;
    
    gyzL = g423;
    
    gzzL = g433;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp229 = SQR(gxzL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp230 = SQR(gyzL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp231 = SQR(gxyL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED detg = 2.*gxyL*gxzL*gyzL + 
      gyyL*(gxxL*gzzL - 1.*csetemp229) - 1.*(gxxL*csetemp230 + 
      gzzL*csetemp231);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp232 = INV(detg);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu11 = (gyyL*gzzL - 
      1.*csetemp230)*csetemp232;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu12 = (gxzL*gyzL - 
      1.*gxyL*gzzL)*csetemp232;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu13 = (-1.*gxzL*gyyL + 
      gxyL*gyzL)*csetemp232;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu22 = (gxxL*gzzL - 
      1.*csetemp229)*csetemp232;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu23 = (gxyL*gxzL - 
      1.*gxxL*gyzL)*csetemp232;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu33 = (gxxL*gyyL - 
      1.*csetemp231)*csetemp232;
    
    betaxL = betal1*gu11 + betal2*gu12 + betal3*gu13;
    
    betayL = betal1*gu12 + betal2*gu22 + betal3*gu23;
    
    betazL = betal1*gu13 + betal2*gu23 + betal3*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betasq = betaxL*betal1 + 
      betayL*betal2 + betazL*betal3;
    
    alpL = sqrt(betasq - 1.*g400);
    
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp233 = dtg11*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp234 = dtg12*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp235 = dtg13*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp236 = dtg12*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp237 = dtg22*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp238 = dtg23*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp239 = dtg13*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp240 = dtg23*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp241 = dtg33*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu11 = -1.*((csetemp233 + csetemp234 
      + csetemp235)*gu11 + (csetemp236 + csetemp237 + csetemp238)*gu12 + 
      (csetemp239 + csetemp240 + csetemp241)*gu13);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu12 = -1.*((csetemp233 + csetemp234 
      + csetemp235)*gu12 + (csetemp236 + csetemp237 + csetemp238)*gu22 + 
      (csetemp239 + csetemp240 + csetemp241)*gu23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu13 = -1.*((csetemp233 + csetemp234 
      + csetemp235)*gu13 + (csetemp236 + csetemp237 + csetemp238)*gu23 + 
      (csetemp239 + csetemp240 + csetemp241)*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp242 = dtg11*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp243 = dtg12*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp244 = dtg13*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp245 = dtg22*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp246 = dtg23*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp247 = dtg13*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp248 = dtg23*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp249 = dtg33*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu22 = -1.*((csetemp242 + csetemp243 
      + csetemp244)*gu12 + (csetemp234 + csetemp245 + csetemp246)*gu22 + 
      (csetemp247 + csetemp248 + csetemp249)*gu23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu23 = -1.*((csetemp242 + csetemp243 
      + csetemp244)*gu13 + (csetemp234 + csetemp245 + csetemp246)*gu23 + 
      (csetemp247 + csetemp248 + csetemp249)*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu33 = -1.*(gu13*(dtg11*gu13 + 
      dtg12*gu23 + dtg13*gu33) + gu23*(dtg12*gu13 + dtg22*gu23 + dtg23*gu33) 
      + gu33*(csetemp235 + csetemp246 + dtg33*gu33));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp250 = dg111*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp251 = dg121*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp252 = dg131*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp253 = dg121*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp254 = dg221*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp255 = dg231*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp256 = dg131*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp257 = dg231*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp258 = dg331*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu111 = -1.*((csetemp250 + csetemp251 
      + csetemp252)*gu11 + (csetemp253 + csetemp254 + csetemp255)*gu12 + 
      (csetemp256 + csetemp257 + csetemp258)*gu13);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu121 = -1.*((csetemp250 + csetemp251 
      + csetemp252)*gu12 + (csetemp253 + csetemp254 + csetemp255)*gu22 + 
      (csetemp256 + csetemp257 + csetemp258)*gu23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu131 = -1.*((csetemp250 + csetemp251 
      + csetemp252)*gu13 + (csetemp253 + csetemp254 + csetemp255)*gu23 + 
      (csetemp256 + csetemp257 + csetemp258)*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp259 = dg111*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp260 = dg121*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp261 = dg131*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp262 = dg221*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp263 = dg231*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp264 = dg131*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp265 = dg231*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp266 = dg331*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu221 = -1.*((csetemp259 + csetemp260 
      + csetemp261)*gu12 + (csetemp251 + csetemp262 + csetemp263)*gu22 + 
      (csetemp264 + csetemp265 + csetemp266)*gu23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu231 = -1.*((csetemp259 + csetemp260 
      + csetemp261)*gu13 + (csetemp251 + csetemp262 + csetemp263)*gu23 + 
      (csetemp264 + csetemp265 + csetemp266)*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu331 = -1.*(gu13*(dg111*gu13 + 
      dg121*gu23 + dg131*gu33) + gu23*(dg121*gu13 + dg221*gu23 + dg231*gu33) 
      + gu33*(csetemp252 + csetemp263 + dg331*gu33));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp267 = dg112*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp268 = dg122*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp269 = dg132*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp270 = dg122*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp271 = dg222*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp272 = dg232*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp273 = dg132*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp274 = dg232*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp275 = dg332*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu112 = -1.*((csetemp267 + csetemp268 
      + csetemp269)*gu11 + (csetemp270 + csetemp271 + csetemp272)*gu12 + 
      (csetemp273 + csetemp274 + csetemp275)*gu13);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu122 = -1.*((csetemp267 + csetemp268 
      + csetemp269)*gu12 + (csetemp270 + csetemp271 + csetemp272)*gu22 + 
      (csetemp273 + csetemp274 + csetemp275)*gu23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu132 = -1.*((csetemp267 + csetemp268 
      + csetemp269)*gu13 + (csetemp270 + csetemp271 + csetemp272)*gu23 + 
      (csetemp273 + csetemp274 + csetemp275)*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp276 = dg112*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp277 = dg122*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp278 = dg132*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp279 = dg222*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp280 = dg232*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp281 = dg132*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp282 = dg232*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp283 = dg332*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu222 = -1.*((csetemp276 + csetemp277 
      + csetemp278)*gu12 + (csetemp268 + csetemp279 + csetemp280)*gu22 + 
      (csetemp281 + csetemp282 + csetemp283)*gu23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu232 = -1.*((csetemp276 + csetemp277 
      + csetemp278)*gu13 + (csetemp268 + csetemp279 + csetemp280)*gu23 + 
      (csetemp281 + csetemp282 + csetemp283)*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu332 = -1.*(gu13*(dg112*gu13 + 
      dg122*gu23 + dg132*gu33) + gu23*(dg122*gu13 + dg222*gu23 + dg232*gu33) 
      + gu33*(csetemp269 + csetemp280 + dg332*gu33));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp284 = dg113*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp285 = dg123*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp286 = dg133*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp287 = dg123*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp288 = dg223*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp289 = dg233*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp290 = dg133*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp291 = dg233*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp292 = dg333*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu113 = -1.*((csetemp284 + csetemp285 
      + csetemp286)*gu11 + (csetemp287 + csetemp288 + csetemp289)*gu12 + 
      (csetemp290 + csetemp291 + csetemp292)*gu13);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu123 = -1.*((csetemp284 + csetemp285 
      + csetemp286)*gu12 + (csetemp287 + csetemp288 + csetemp289)*gu22 + 
      (csetemp290 + csetemp291 + csetemp292)*gu23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu133 = -1.*((csetemp284 + csetemp285 
      + csetemp286)*gu13 + (csetemp287 + csetemp288 + csetemp289)*gu23 + 
      (csetemp290 + csetemp291 + csetemp292)*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp293 = dg113*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp294 = dg123*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp295 = dg133*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp296 = dg223*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp297 = dg233*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp298 = dg133*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp299 = dg233*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp300 = dg333*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu223 = -1.*((csetemp293 + csetemp294 
      + csetemp295)*gu12 + (csetemp285 + csetemp296 + csetemp297)*gu22 + 
      (csetemp298 + csetemp299 + csetemp300)*gu23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu233 = -1.*((csetemp293 + csetemp294 
      + csetemp295)*gu13 + (csetemp285 + csetemp296 + csetemp297)*gu23 + 
      (csetemp298 + csetemp299 + csetemp300)*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu333 = -1.*(gu13*(dg113*gu13 + 
      dg123*gu23 + dg133*gu33) + gu23*(dg123*gu13 + dg223*gu23 + dg233*gu33) 
      + gu33*(csetemp286 + csetemp297 + dg333*gu33));
    
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp301 = INV(alpL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtalpL = 0.5*csetemp301*(-1.*dg4000 + 
      dtbetasq);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxxL = 
      0.5*csetemp301*(2.*(gxxL*dbeta11 + gxyL*dbeta21 + gxzL*dbeta31) + 
      betaxL*dg111 + betayL*dg112 + betazL*dg113 - 1.*dtg11);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxyL = 0.5*csetemp301*(gxxL*dbeta12 
      + gyyL*dbeta21 + gxyL*(dbeta11 + dbeta22) + gyzL*dbeta31 + 
      gxzL*dbeta32 + betaxL*dg121 + betayL*dg122 + betazL*dg123 - 
      1.*dtg12);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxzL = 0.5*csetemp301*(gxxL*dbeta13 
      + gyzL*dbeta21 + gxyL*dbeta23 + gzzL*dbeta31 + gxzL*(dbeta11 + 
      dbeta33) + betaxL*dg131 + betayL*dg132 + betazL*dg133 - 
      1.*dtg13);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kyyL = 
      0.5*csetemp301*(2.*(gxyL*dbeta12 + gyyL*dbeta22 + gyzL*dbeta32) + 
      betaxL*dg221 + betayL*dg222 + betazL*dg223 - 1.*dtg22);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kyzL = 0.5*csetemp301*(gxzL*dbeta12 
      + gxyL*dbeta13 + gyyL*dbeta23 + gzzL*dbeta32 + gyzL*(dbeta22 + 
      dbeta33) + betaxL*dg231 + betayL*dg232 + betazL*dg233 - 
      1.*dtg23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kzzL = 
      0.5*csetemp301*(2.*(gxzL*dbeta13 + gyzL*dbeta23 + gzzL*dbeta33) + 
      betaxL*dg331 + betayL*dg332 + betazL*dg333 - 1.*dtg33);
    
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
