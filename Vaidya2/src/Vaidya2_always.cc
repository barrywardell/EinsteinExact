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

static void Vaidya2_always_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(0));
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(1));
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(2));
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_TIME);
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = INV(dx);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = INV(dy);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = INV(dz);
  const CCTK_REAL khalf CCTK_ATTRIBUTE_UNUSED = 0.5;
  const CCTK_REAL kthird CCTK_ATTRIBUTE_UNUSED = 1/3.0;
  const CCTK_REAL ktwothird CCTK_ATTRIBUTE_UNUSED = 2.0/3.0;
  const CCTK_REAL kfourthird CCTK_ATTRIBUTE_UNUSED = 4.0/3.0;
  const CCTK_REAL keightthird CCTK_ATTRIBUTE_UNUSED = 8.0/3.0;
  const CCTK_REAL hdxi CCTK_ATTRIBUTE_UNUSED = 0.5 * dxi;
  const CCTK_REAL hdyi CCTK_ATTRIBUTE_UNUSED = 0.5 * dyi;
  const CCTK_REAL hdzi CCTK_ATTRIBUTE_UNUSED = 0.5 * dzi;
  
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
  #pragma omp parallel // reduction(+: vec_iter_counter, vec_op_counter, vec_mem_counter)
  CCTK_LOOP3(Vaidya2_always,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    // ++vec_iter_counter;
    
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED txx0 = xformL00*xx0 + xformL01*xx1 + 
      xformL02*xx2 + xformL03*xx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED txx1 = xformL10*xx0 + xformL11*xx1 + 
      xformL12*xx2 + xformL13*xx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED txx2 = xformL20*xx0 + xformL21*xx1 + 
      xformL22*xx2 + xformL23*xx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED txx3 = xformL30*xx0 + xformL31*xx1 + 
      xformL32*xx2 + xformL33*xx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED T = txx0;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED X = txx1;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED Y = txx2;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED Z = txx3;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp10 = SQR(X);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp11 = SQR(Y);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp12 = SQR(Z);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED rXYZ = sqrt(csetemp10 + csetemp11 + 
      csetemp12);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp13 = INV(ToReal(M));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED mTXYZ = (1 + SQR(tanh(csetemp13*(T + 
      sqrt(csetemp10 + csetemp11 + csetemp12))*ToReal(dM))))*ToReal(M);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp14 = INV(rXYZ);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg400 = -1 + 2*csetemp14*mTXYZ;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp15 = INV(SQR(rXYZ));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg401 = 2*csetemp15*mTXYZ*X;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg402 = 2*csetemp15*mTXYZ*Y;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg403 = 2*csetemp15*mTXYZ*Z;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp16 = INV(CUB(rXYZ));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg411 = 1 + 
      2*csetemp10*csetemp16*mTXYZ;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg412 = 2*csetemp16*mTXYZ*X*Y;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg413 = 2*csetemp16*mTXYZ*X*Z;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg422 = 1 + 
      2*csetemp11*csetemp16*mTXYZ;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg423 = 2*csetemp16*mTXYZ*Y*Z;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tg433 = 1 + 
      2*csetemp12*csetemp16*mTXYZ;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp17 = rXYZ + T;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4000 = 
      4*csetemp14*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4001 = -2*csetemp16*X*(mTXYZ - 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4002 = -2*csetemp16*Y*(mTXYZ - 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4003 = -2*csetemp16*Z*(mTXYZ - 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4010 = 
      4*csetemp15*X*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp18 = INV(QAD(rXYZ));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp19 = SQR(rXYZ);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4011 = csetemp18*(2*(-2*csetemp10 + 
      csetemp19)*mTXYZ + 
      4*csetemp10*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4012 = -4*csetemp18*X*Y*(mTXYZ - 
      rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4013 = -4*csetemp18*X*Z*(mTXYZ - 
      rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4020 = 
      4*csetemp15*Y*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4021 = -4*csetemp18*X*Y*(mTXYZ - 
      rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4022 = csetemp18*(2*(-2*csetemp11 + 
      csetemp19)*mTXYZ + 
      4*csetemp11*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4023 = -4*csetemp18*Y*Z*(mTXYZ - 
      rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4030 = 
      4*csetemp15*Z*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4031 = -4*csetemp18*X*Z*(mTXYZ - 
      rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4032 = -4*csetemp18*Y*Z*(mTXYZ - 
      rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4033 = csetemp18*(2*(-2*csetemp12 + 
      csetemp19)*mTXYZ + 
      4*csetemp12*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4110 = 
      4*csetemp10*csetemp16*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp20 = pow(rXYZ,-5);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp21 = CUB(X);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4111 = 
      csetemp20*(mTXYZ*(-6*csetemp21 + 4*csetemp19*X) + 
      4*csetemp21*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp22 = -3*mTXYZ;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4112 = 
      2*csetemp10*csetemp20*Y*(csetemp22 + 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4113 = 
      2*csetemp10*csetemp20*Z*(csetemp22 + 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4120 = 
      4*csetemp16*X*Y*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4121 = 2*csetemp20*Y*((-2*csetemp10 
      + csetemp11 + csetemp12)*mTXYZ + 
      2*csetemp10*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4122 = 2*csetemp20*X*((csetemp10 - 
      2*csetemp11 + csetemp12)*mTXYZ + 
      2*csetemp11*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4123 = 2*csetemp20*X*Y*Z*(csetemp22 
      + 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4130 = 
      4*csetemp16*X*Z*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4131 = 2*csetemp20*Z*((-2*csetemp10 
      + csetemp11 + csetemp12)*mTXYZ + 
      2*csetemp10*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4132 = 2*csetemp20*X*Y*Z*(csetemp22 
      + 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4133 = 2*csetemp20*X*((csetemp10 + 
      csetemp11 - 2*csetemp12)*mTXYZ + 
      2*csetemp12*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4220 = 
      4*csetemp11*csetemp16*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4221 = 
      2*csetemp11*csetemp20*X*(csetemp22 + 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp23 = CUB(Y);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4222 = 
      csetemp20*(mTXYZ*(-6*csetemp23 + 4*csetemp19*Y) + 
      4*csetemp23*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4223 = 
      2*csetemp11*csetemp20*Z*(csetemp22 + 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4230 = 
      4*csetemp16*Y*Z*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4231 = 2*csetemp20*X*Y*Z*(csetemp22 
      + 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4232 = 2*csetemp20*Z*((csetemp10 - 
      2*csetemp11 + csetemp12)*mTXYZ + 
      2*csetemp11*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4233 = 2*csetemp20*Y*((csetemp10 + 
      csetemp11 - 2*csetemp12)*mTXYZ + 
      2*csetemp12*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4330 = 
      4*csetemp12*csetemp16*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4331 = 
      2*csetemp12*csetemp20*X*(csetemp22 + 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4332 = 
      2*csetemp12*csetemp20*Y*(csetemp22 + 
      2*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp24 = CUB(Z);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED tdg4333 = 
      csetemp20*(mTXYZ*(-6*csetemp24 + 4*csetemp19*Z) + 
      4*csetemp24*rXYZ*SQR(INV(cosh(csetemp13*csetemp17*ToReal(dM))))*tanh(csetemp13*csetemp17*ToReal(dM))*ToReal(dM));
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g400 = 2*(tg423*xformL20*xformL30 + 
      xformL00*(tg401*xformL10 + tg402*xformL20 + tg403*xformL30) + 
      xformL10*(tg412*xformL20 + tg413*xformL30)) + tg400*SQR(xformL00) + 
      tg411*SQR(xformL10) + tg422*SQR(xformL20) + tg433*SQR(xformL30);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp25 = tg400*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp26 = tg401*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp27 = tg402*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp28 = tg403*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp29 = tg401*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp30 = tg411*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp31 = tg412*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp32 = tg413*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp33 = tg402*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp34 = tg412*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp35 = tg422*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp36 = tg423*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp37 = tg403*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp38 = tg413*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp39 = tg423*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp40 = tg433*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g401 = (csetemp25 + csetemp26 + 
      csetemp27 + csetemp28)*xformL00 + (csetemp29 + csetemp30 + csetemp31 + 
      csetemp32)*xformL10 + (csetemp33 + csetemp34 + csetemp35 + 
      csetemp36)*xformL20 + (csetemp37 + csetemp38 + csetemp39 + 
      csetemp40)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp41 = tg400*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp42 = tg401*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp43 = tg402*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp44 = tg403*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp45 = tg401*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp46 = tg411*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp47 = tg412*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp48 = tg413*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp49 = tg402*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp50 = tg412*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp51 = tg422*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp52 = tg423*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp53 = tg403*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp54 = tg413*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp55 = tg423*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp56 = tg433*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g402 = (csetemp41 + csetemp42 + 
      csetemp43 + csetemp44)*xformL00 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL10 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL20 + (csetemp53 + csetemp54 + csetemp55 + 
      csetemp56)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp57 = tg400*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp58 = tg401*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp59 = tg402*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp60 = tg403*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp61 = tg401*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp62 = tg411*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp63 = tg412*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp64 = tg413*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp65 = tg402*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp66 = tg412*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp67 = tg422*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp68 = tg423*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp69 = tg403*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp70 = tg413*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp71 = tg423*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp72 = tg433*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g403 = (csetemp57 + csetemp58 + 
      csetemp59 + csetemp60)*xformL00 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL10 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL20 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g411 = (csetemp25 + csetemp26 + 
      csetemp27 + csetemp28)*xformL01 + (csetemp29 + csetemp30 + csetemp31 + 
      csetemp32)*xformL11 + (csetemp33 + csetemp34 + csetemp35 + 
      csetemp36)*xformL21 + (csetemp37 + csetemp38 + csetemp39 + 
      csetemp40)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g412 = (csetemp41 + csetemp42 + 
      csetemp43 + csetemp44)*xformL01 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL11 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL21 + (csetemp53 + csetemp54 + csetemp55 + 
      csetemp56)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g413 = (csetemp57 + csetemp58 + 
      csetemp59 + csetemp60)*xformL01 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL11 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL21 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g422 = (csetemp41 + csetemp42 + 
      csetemp43 + csetemp44)*xformL02 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL12 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL22 + (csetemp53 + csetemp54 + csetemp55 + 
      csetemp56)*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g423 = (csetemp57 + csetemp58 + 
      csetemp59 + csetemp60)*xformL02 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL12 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL22 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED g433 = (csetemp57 + csetemp58 + 
      csetemp59 + csetemp60)*xformL03 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL13 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL23 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp73 = tdg4000*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp74 = tdg4001*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp75 = tdg4002*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp76 = tdg4003*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp77 = tdg4010*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp78 = tdg4011*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp79 = tdg4012*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp80 = tdg4013*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp81 = tdg4020*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp82 = tdg4021*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp83 = tdg4022*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp84 = tdg4023*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp85 = tdg4030*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp86 = tdg4031*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp87 = tdg4032*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp88 = tdg4033*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp89 = tdg4110*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp90 = tdg4111*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp91 = tdg4112*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp92 = tdg4113*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp93 = tdg4120*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp94 = tdg4121*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp95 = tdg4122*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp96 = tdg4123*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp97 = tdg4130*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp98 = tdg4131*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp99 = tdg4132*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp100 = tdg4133*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp101 = tdg4220*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp102 = tdg4221*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp103 = tdg4222*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp104 = tdg4223*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp105 = tdg4230*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp106 = tdg4231*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp107 = tdg4232*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp108 = tdg4233*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp109 = tdg4330*xformL00;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp110 = tdg4331*xformL10;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp111 = tdg4332*xformL20;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp112 = tdg4333*xformL30;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4000 = xformL20*((csetemp81 + 
      csetemp82 + csetemp83 + csetemp84)*xformL00 + (csetemp93 + csetemp94 + 
      csetemp95 + csetemp96)*xformL10 + (csetemp101 + csetemp102 + csetemp103 
      + csetemp104)*xformL20 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL30) + xformL30*((csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL00 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL10 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL20 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*xformL30) + xformL00*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL00 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL10 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL20 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL30) + xformL10*((csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL00 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL10 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL20 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL30);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4010 = xformL20*((csetemp81 + 
      csetemp82 + csetemp83 + csetemp84)*xformL01 + (csetemp93 + csetemp94 + 
      csetemp95 + csetemp96)*xformL11 + (csetemp101 + csetemp102 + csetemp103 
      + csetemp104)*xformL21 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL31) + xformL30*((csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL01 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL11 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL21 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*xformL31) + xformL00*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL01 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL11 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL21 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL31) + xformL10*((csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL01 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL11 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL21 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp113 = tdg4000*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp114 = tdg4001*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp115 = tdg4002*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp116 = tdg4003*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp117 = tdg4010*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp118 = tdg4011*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp119 = tdg4012*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp120 = tdg4013*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp121 = tdg4020*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp122 = tdg4021*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp123 = tdg4022*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp124 = tdg4023*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp125 = tdg4030*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp126 = tdg4031*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp127 = tdg4032*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp128 = tdg4033*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp129 = tdg4110*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp130 = tdg4111*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp131 = tdg4112*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp132 = tdg4113*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp133 = tdg4120*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp134 = tdg4121*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp135 = tdg4122*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp136 = tdg4123*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp137 = tdg4130*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp138 = tdg4131*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp139 = tdg4132*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp140 = tdg4133*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp141 = tdg4220*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp142 = tdg4221*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp143 = tdg4222*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp144 = tdg4223*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp145 = tdg4230*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp146 = tdg4231*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp147 = tdg4232*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp148 = tdg4233*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp149 = tdg4330*xformL01;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp150 = tdg4331*xformL11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp151 = tdg4332*xformL21;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp152 = tdg4333*xformL31;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4011 = xformL00*((csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL01 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL11 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL21 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL31) + xformL10*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL01 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL11 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL21 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL31) + xformL20*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL01 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL11 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL21 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL31) + xformL30*((csetemp125 
      + csetemp126 + csetemp127 + csetemp128)*xformL01 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL11 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL21 + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp153 = tdg4000*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp154 = tdg4001*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp155 = tdg4002*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp156 = tdg4003*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp157 = tdg4010*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp158 = tdg4011*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp159 = tdg4012*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp160 = tdg4013*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp161 = tdg4020*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp162 = tdg4021*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp163 = tdg4022*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp164 = tdg4023*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp165 = tdg4030*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp166 = tdg4031*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp167 = tdg4032*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp168 = tdg4033*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp169 = tdg4110*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp170 = tdg4111*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp171 = tdg4112*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp172 = tdg4113*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp173 = tdg4120*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp174 = tdg4121*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp175 = tdg4122*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp176 = tdg4123*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp177 = tdg4130*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp178 = tdg4131*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp179 = tdg4132*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp180 = tdg4133*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp181 = tdg4220*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp182 = tdg4221*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp183 = tdg4222*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp184 = tdg4223*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp185 = tdg4230*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp186 = tdg4231*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp187 = tdg4232*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp188 = tdg4233*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp189 = tdg4330*xformL02;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp190 = tdg4331*xformL12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp191 = tdg4332*xformL22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp192 = tdg4333*xformL32;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4012 = xformL00*((csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL01 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL11 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL21 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL31) + xformL10*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL01 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL11 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL21 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL31) + xformL20*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL01 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL11 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL21 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL31) + xformL30*((csetemp165 
      + csetemp166 + csetemp167 + csetemp168)*xformL01 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL11 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL21 + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp193 = tdg4000*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp194 = tdg4001*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp195 = tdg4002*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp196 = tdg4003*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp197 = tdg4010*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp198 = tdg4011*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp199 = tdg4012*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp200 = tdg4013*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp201 = tdg4020*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp202 = tdg4021*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp203 = tdg4022*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp204 = tdg4023*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp205 = tdg4030*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp206 = tdg4031*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp207 = tdg4032*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp208 = tdg4033*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp209 = tdg4110*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp210 = tdg4111*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp211 = tdg4112*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp212 = tdg4113*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp213 = tdg4120*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp214 = tdg4121*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp215 = tdg4122*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp216 = tdg4123*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp217 = tdg4130*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp218 = tdg4131*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp219 = tdg4132*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp220 = tdg4133*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp221 = tdg4220*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp222 = tdg4221*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp223 = tdg4222*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp224 = tdg4223*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp225 = tdg4230*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp226 = tdg4231*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp227 = tdg4232*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp228 = tdg4233*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp229 = tdg4330*xformL03;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp230 = tdg4331*xformL13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp231 = tdg4332*xformL23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp232 = tdg4333*xformL33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4013 = xformL00*((csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL01 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL11 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL21 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL31) + xformL10*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL01 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL11 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL21 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL31) + xformL20*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL01 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL11 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL21 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL31) + xformL30*((csetemp205 
      + csetemp206 + csetemp207 + csetemp208)*xformL01 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL11 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL21 + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4020 = xformL20*((csetemp81 + 
      csetemp82 + csetemp83 + csetemp84)*xformL02 + (csetemp93 + csetemp94 + 
      csetemp95 + csetemp96)*xformL12 + (csetemp101 + csetemp102 + csetemp103 
      + csetemp104)*xformL22 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL32) + xformL30*((csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL02 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL12 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL22 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*xformL32) + xformL00*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL02 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL12 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL22 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL32) + xformL10*((csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL02 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL12 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL22 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4021 = xformL00*((csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL02 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL12 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL22 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL32) + xformL10*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL02 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL12 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL22 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL32) + xformL20*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL02 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL12 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL22 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL32) + xformL30*((csetemp125 
      + csetemp126 + csetemp127 + csetemp128)*xformL02 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL12 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL22 + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4022 = xformL00*((csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL02 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL12 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL22 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL32) + xformL10*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL02 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL12 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL22 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL32) + xformL20*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL02 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL12 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL22 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL32) + xformL30*((csetemp165 
      + csetemp166 + csetemp167 + csetemp168)*xformL02 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL12 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL22 + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4023 = xformL00*((csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL02 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL12 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL22 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL32) + xformL10*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL02 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL12 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL22 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL32) + xformL20*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL02 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL12 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL22 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL32) + xformL30*((csetemp205 
      + csetemp206 + csetemp207 + csetemp208)*xformL02 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL12 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL22 + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4030 = xformL20*((csetemp81 + 
      csetemp82 + csetemp83 + csetemp84)*xformL03 + (csetemp93 + csetemp94 + 
      csetemp95 + csetemp96)*xformL13 + (csetemp101 + csetemp102 + csetemp103 
      + csetemp104)*xformL23 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL33) + xformL30*((csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL03 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL13 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL23 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*xformL33) + xformL00*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL03 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL13 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL23 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL33) + xformL10*((csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL03 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL13 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL23 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4031 = xformL00*((csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL03 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL13 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL23 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL33) + xformL10*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL03 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL13 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL23 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL33) + xformL20*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL03 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL13 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL23 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL33) + xformL30*((csetemp125 
      + csetemp126 + csetemp127 + csetemp128)*xformL03 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL13 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL23 + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4032 = xformL00*((csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL03 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL13 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL23 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL33) + xformL10*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL03 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL13 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL23 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL33) + xformL20*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL03 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL13 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL23 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL33) + xformL30*((csetemp165 
      + csetemp166 + csetemp167 + csetemp168)*xformL03 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL13 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL23 + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4033 = xformL00*((csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL03 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL13 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL23 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL33) + xformL10*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL03 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL13 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL23 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL33) + xformL20*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL03 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL13 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL23 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL33) + xformL30*((csetemp205 
      + csetemp206 + csetemp207 + csetemp208)*xformL03 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL13 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL23 + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4110 = xformL21*((csetemp81 + 
      csetemp82 + csetemp83 + csetemp84)*xformL01 + (csetemp93 + csetemp94 + 
      csetemp95 + csetemp96)*xformL11 + (csetemp101 + csetemp102 + csetemp103 
      + csetemp104)*xformL21 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL31) + xformL31*((csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL01 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL11 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL21 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*xformL31) + xformL01*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL01 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL11 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL21 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL31) + xformL11*((csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL01 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL11 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL21 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4111 = xformL01*((csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL01 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL11 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL21 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL31) + xformL11*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL01 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL11 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL21 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL31) + xformL21*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL01 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL11 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL21 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL31) + xformL31*((csetemp125 
      + csetemp126 + csetemp127 + csetemp128)*xformL01 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL11 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL21 + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4112 = xformL01*((csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL01 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL11 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL21 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL31) + xformL11*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL01 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL11 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL21 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL31) + xformL21*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL01 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL11 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL21 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL31) + xformL31*((csetemp165 
      + csetemp166 + csetemp167 + csetemp168)*xformL01 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL11 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL21 + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4113 = xformL01*((csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL01 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL11 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL21 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL31) + xformL11*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL01 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL11 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL21 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL31) + xformL21*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL01 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL11 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL21 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL31) + xformL31*((csetemp205 
      + csetemp206 + csetemp207 + csetemp208)*xformL01 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL11 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL21 + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*xformL31);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4120 = xformL21*((csetemp81 + 
      csetemp82 + csetemp83 + csetemp84)*xformL02 + (csetemp93 + csetemp94 + 
      csetemp95 + csetemp96)*xformL12 + (csetemp101 + csetemp102 + csetemp103 
      + csetemp104)*xformL22 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL32) + xformL31*((csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL02 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL12 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL22 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*xformL32) + xformL01*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL02 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL12 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL22 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL32) + xformL11*((csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL02 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL12 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL22 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4121 = xformL01*((csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL02 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL12 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL22 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL32) + xformL11*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL02 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL12 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL22 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL32) + xformL21*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL02 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL12 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL22 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL32) + xformL31*((csetemp125 
      + csetemp126 + csetemp127 + csetemp128)*xformL02 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL12 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL22 + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4122 = xformL01*((csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL02 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL12 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL22 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL32) + xformL11*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL02 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL12 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL22 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL32) + xformL21*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL02 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL12 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL22 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL32) + xformL31*((csetemp165 
      + csetemp166 + csetemp167 + csetemp168)*xformL02 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL12 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL22 + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4123 = xformL01*((csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL02 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL12 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL22 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL32) + xformL11*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL02 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL12 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL22 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL32) + xformL21*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL02 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL12 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL22 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL32) + xformL31*((csetemp205 
      + csetemp206 + csetemp207 + csetemp208)*xformL02 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL12 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL22 + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4130 = xformL21*((csetemp81 + 
      csetemp82 + csetemp83 + csetemp84)*xformL03 + (csetemp93 + csetemp94 + 
      csetemp95 + csetemp96)*xformL13 + (csetemp101 + csetemp102 + csetemp103 
      + csetemp104)*xformL23 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL33) + xformL31*((csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL03 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL13 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL23 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*xformL33) + xformL01*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL03 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL13 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL23 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL33) + xformL11*((csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL03 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL13 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL23 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4131 = xformL01*((csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL03 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL13 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL23 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL33) + xformL11*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL03 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL13 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL23 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL33) + xformL21*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL03 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL13 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL23 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL33) + xformL31*((csetemp125 
      + csetemp126 + csetemp127 + csetemp128)*xformL03 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL13 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL23 + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4132 = xformL01*((csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL03 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL13 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL23 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL33) + xformL11*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL03 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL13 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL23 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL33) + xformL21*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL03 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL13 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL23 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL33) + xformL31*((csetemp165 
      + csetemp166 + csetemp167 + csetemp168)*xformL03 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL13 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL23 + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4133 = xformL01*((csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL03 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL13 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL23 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL33) + xformL11*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL03 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL13 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL23 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL33) + xformL21*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL03 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL13 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL23 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL33) + xformL31*((csetemp205 
      + csetemp206 + csetemp207 + csetemp208)*xformL03 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL13 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL23 + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4220 = xformL22*((csetemp81 + 
      csetemp82 + csetemp83 + csetemp84)*xformL02 + (csetemp93 + csetemp94 + 
      csetemp95 + csetemp96)*xformL12 + (csetemp101 + csetemp102 + csetemp103 
      + csetemp104)*xformL22 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL32) + xformL32*((csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL02 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL12 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL22 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*xformL32) + xformL02*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL02 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL12 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL22 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL32) + xformL12*((csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL02 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL12 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL22 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4221 = xformL02*((csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL02 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL12 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL22 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL32) + xformL12*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL02 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL12 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL22 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL32) + xformL22*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL02 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL12 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL22 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL32) + xformL32*((csetemp125 
      + csetemp126 + csetemp127 + csetemp128)*xformL02 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL12 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL22 + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4222 = xformL02*((csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL02 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL12 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL22 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL32) + xformL12*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL02 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL12 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL22 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL32) + xformL22*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL02 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL12 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL22 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL32) + xformL32*((csetemp165 
      + csetemp166 + csetemp167 + csetemp168)*xformL02 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL12 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL22 + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4223 = xformL02*((csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL02 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL12 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL22 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL32) + xformL12*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL02 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL12 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL22 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL32) + xformL22*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL02 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL12 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL22 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL32) + xformL32*((csetemp205 
      + csetemp206 + csetemp207 + csetemp208)*xformL02 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL12 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL22 + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*xformL32);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4230 = xformL22*((csetemp81 + 
      csetemp82 + csetemp83 + csetemp84)*xformL03 + (csetemp93 + csetemp94 + 
      csetemp95 + csetemp96)*xformL13 + (csetemp101 + csetemp102 + csetemp103 
      + csetemp104)*xformL23 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL33) + xformL32*((csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL03 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL13 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL23 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*xformL33) + xformL02*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL03 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL13 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL23 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL33) + xformL12*((csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL03 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL13 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL23 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4231 = xformL02*((csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL03 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL13 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL23 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL33) + xformL12*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL03 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL13 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL23 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL33) + xformL22*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL03 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL13 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL23 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL33) + xformL32*((csetemp125 
      + csetemp126 + csetemp127 + csetemp128)*xformL03 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL13 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL23 + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4232 = xformL02*((csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL03 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL13 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL23 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL33) + xformL12*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL03 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL13 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL23 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL33) + xformL22*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL03 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL13 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL23 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL33) + xformL32*((csetemp165 
      + csetemp166 + csetemp167 + csetemp168)*xformL03 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL13 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL23 + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4233 = xformL02*((csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL03 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL13 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL23 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL33) + xformL12*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL03 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL13 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL23 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL33) + xformL22*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL03 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL13 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL23 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL33) + xformL32*((csetemp205 
      + csetemp206 + csetemp207 + csetemp208)*xformL03 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL13 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL23 + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4330 = xformL23*((csetemp81 + 
      csetemp82 + csetemp83 + csetemp84)*xformL03 + (csetemp93 + csetemp94 + 
      csetemp95 + csetemp96)*xformL13 + (csetemp101 + csetemp102 + csetemp103 
      + csetemp104)*xformL23 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL33) + xformL33*((csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL03 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL13 + (csetemp105 + csetemp106 + csetemp107 + 
      csetemp108)*xformL23 + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*xformL33) + xformL03*((csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*xformL03 + (csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL13 + (csetemp81 + csetemp82 + csetemp83 + 
      csetemp84)*xformL23 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL33) + xformL13*((csetemp77 + csetemp78 + csetemp79 + 
      csetemp80)*xformL03 + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*xformL13 + (csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL23 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4331 = xformL03*((csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*xformL03 + (csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL13 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL23 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL33) + xformL13*((csetemp117 
      + csetemp118 + csetemp119 + csetemp120)*xformL03 + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*xformL13 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL23 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL33) + xformL23*((csetemp121 
      + csetemp122 + csetemp123 + csetemp124)*xformL03 + (csetemp133 + 
      csetemp134 + csetemp135 + csetemp136)*xformL13 + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*xformL23 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL33) + xformL33*((csetemp125 
      + csetemp126 + csetemp127 + csetemp128)*xformL03 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL13 + (csetemp145 + 
      csetemp146 + csetemp147 + csetemp148)*xformL23 + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4332 = xformL03*((csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*xformL03 + (csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL13 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL23 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL33) + xformL13*((csetemp157 
      + csetemp158 + csetemp159 + csetemp160)*xformL03 + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*xformL13 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL23 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL33) + xformL23*((csetemp161 
      + csetemp162 + csetemp163 + csetemp164)*xformL03 + (csetemp173 + 
      csetemp174 + csetemp175 + csetemp176)*xformL13 + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*xformL23 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL33) + xformL33*((csetemp165 
      + csetemp166 + csetemp167 + csetemp168)*xformL03 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL13 + (csetemp185 + 
      csetemp186 + csetemp187 + csetemp188)*xformL23 + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dg4333 = xformL03*((csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*xformL03 + (csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL13 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL23 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL33) + xformL13*((csetemp197 
      + csetemp198 + csetemp199 + csetemp200)*xformL03 + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*xformL13 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL23 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL33) + xformL23*((csetemp201 
      + csetemp202 + csetemp203 + csetemp204)*xformL03 + (csetemp213 + 
      csetemp214 + csetemp215 + csetemp216)*xformL13 + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*xformL23 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL33) + xformL33*((csetemp205 
      + csetemp206 + csetemp207 + csetemp208)*xformL03 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL13 + (csetemp225 + 
      csetemp226 + csetemp227 + csetemp228)*xformL23 + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*xformL33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal1 = g401;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal2 = g402;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED betal3 = g403;
    
    gxxL = g411;
    
    gxyL = g412;
    
    gxzL = g413;
    
    gyyL = g422;
    
    gyzL = g423;
    
    gzzL = g433;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp233 = SQR(gxzL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp234 = SQR(gyzL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp235 = SQR(gxyL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED detg = 2*gxyL*gxzL*gyzL + 
      gyyL*(gxxL*gzzL - csetemp233) - gxxL*csetemp234 - 
      gzzL*csetemp235;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp236 = INV(detg);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu11 = (gyyL*gzzL - 
      csetemp234)*csetemp236;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu12 = (gxzL*gyzL - 
      gxyL*gzzL)*csetemp236;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu13 = (-(gxzL*gyyL) + 
      gxyL*gyzL)*csetemp236;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu22 = (gxxL*gzzL - 
      csetemp233)*csetemp236;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu23 = (gxyL*gxzL - 
      gxxL*gyzL)*csetemp236;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED gu33 = (gxxL*gyyL - 
      csetemp235)*csetemp236;
    
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp237 = dtg11*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp238 = dtg12*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp239 = dtg13*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp240 = dtg12*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp241 = dtg22*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp242 = dtg23*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp243 = dtg13*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp244 = dtg23*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp245 = dtg33*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu11 = -((csetemp237 + csetemp238 + 
      csetemp239)*gu11) - (csetemp240 + csetemp241 + csetemp242)*gu12 - 
      (csetemp243 + csetemp244 + csetemp245)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu12 = -((csetemp237 + csetemp238 + 
      csetemp239)*gu12) - (csetemp240 + csetemp241 + csetemp242)*gu22 - 
      (csetemp243 + csetemp244 + csetemp245)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu13 = -((csetemp237 + csetemp238 + 
      csetemp239)*gu13) - (csetemp240 + csetemp241 + csetemp242)*gu23 - 
      (csetemp243 + csetemp244 + csetemp245)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp246 = dtg11*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp247 = dtg12*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp248 = dtg13*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp249 = dtg22*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp250 = dtg23*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp251 = dtg13*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp252 = dtg23*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp253 = dtg33*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu22 = -((csetemp246 + csetemp247 + 
      csetemp248)*gu12) - (csetemp238 + csetemp249 + csetemp250)*gu22 - 
      (csetemp251 + csetemp252 + csetemp253)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu23 = -((csetemp246 + csetemp247 + 
      csetemp248)*gu13) - (csetemp238 + csetemp249 + csetemp250)*gu23 - 
      (csetemp251 + csetemp252 + csetemp253)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtgu33 = -(gu13*(dtg11*gu13 + 
      dtg12*gu23 + dtg13*gu33)) - gu23*(dtg12*gu13 + dtg22*gu23 + dtg23*gu33) 
      - gu33*(csetemp239 + csetemp250 + dtg33*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp254 = dg111*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp255 = dg121*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp256 = dg131*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp257 = dg121*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp258 = dg221*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp259 = dg231*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp260 = dg131*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp261 = dg231*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp262 = dg331*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu111 = -((csetemp254 + csetemp255 + 
      csetemp256)*gu11) - (csetemp257 + csetemp258 + csetemp259)*gu12 - 
      (csetemp260 + csetemp261 + csetemp262)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu121 = -((csetemp254 + csetemp255 + 
      csetemp256)*gu12) - (csetemp257 + csetemp258 + csetemp259)*gu22 - 
      (csetemp260 + csetemp261 + csetemp262)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu131 = -((csetemp254 + csetemp255 + 
      csetemp256)*gu13) - (csetemp257 + csetemp258 + csetemp259)*gu23 - 
      (csetemp260 + csetemp261 + csetemp262)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp263 = dg111*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp264 = dg121*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp265 = dg131*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp266 = dg221*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp267 = dg231*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp268 = dg131*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp269 = dg231*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp270 = dg331*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu221 = -((csetemp263 + csetemp264 + 
      csetemp265)*gu12) - (csetemp255 + csetemp266 + csetemp267)*gu22 - 
      (csetemp268 + csetemp269 + csetemp270)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu231 = -((csetemp263 + csetemp264 + 
      csetemp265)*gu13) - (csetemp255 + csetemp266 + csetemp267)*gu23 - 
      (csetemp268 + csetemp269 + csetemp270)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu331 = -(gu13*(dg111*gu13 + 
      dg121*gu23 + dg131*gu33)) - gu23*(dg121*gu13 + dg221*gu23 + dg231*gu33) 
      - gu33*(csetemp256 + csetemp267 + dg331*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp271 = dg112*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp272 = dg122*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp273 = dg132*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp274 = dg122*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp275 = dg222*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp276 = dg232*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp277 = dg132*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp278 = dg232*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp279 = dg332*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu112 = -((csetemp271 + csetemp272 + 
      csetemp273)*gu11) - (csetemp274 + csetemp275 + csetemp276)*gu12 - 
      (csetemp277 + csetemp278 + csetemp279)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu122 = -((csetemp271 + csetemp272 + 
      csetemp273)*gu12) - (csetemp274 + csetemp275 + csetemp276)*gu22 - 
      (csetemp277 + csetemp278 + csetemp279)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu132 = -((csetemp271 + csetemp272 + 
      csetemp273)*gu13) - (csetemp274 + csetemp275 + csetemp276)*gu23 - 
      (csetemp277 + csetemp278 + csetemp279)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp280 = dg112*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp281 = dg122*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp282 = dg132*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp283 = dg222*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp284 = dg232*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp285 = dg132*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp286 = dg232*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp287 = dg332*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu222 = -((csetemp280 + csetemp281 + 
      csetemp282)*gu12) - (csetemp272 + csetemp283 + csetemp284)*gu22 - 
      (csetemp285 + csetemp286 + csetemp287)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu232 = -((csetemp280 + csetemp281 + 
      csetemp282)*gu13) - (csetemp272 + csetemp283 + csetemp284)*gu23 - 
      (csetemp285 + csetemp286 + csetemp287)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu332 = -(gu13*(dg112*gu13 + 
      dg122*gu23 + dg132*gu33)) - gu23*(dg122*gu13 + dg222*gu23 + dg232*gu33) 
      - gu33*(csetemp273 + csetemp284 + dg332*gu33);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp288 = dg113*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp289 = dg123*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp290 = dg133*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp291 = dg123*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp292 = dg223*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp293 = dg233*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp294 = dg133*gu11;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp295 = dg233*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp296 = dg333*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu113 = -((csetemp288 + csetemp289 + 
      csetemp290)*gu11) - (csetemp291 + csetemp292 + csetemp293)*gu12 - 
      (csetemp294 + csetemp295 + csetemp296)*gu13;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu123 = -((csetemp288 + csetemp289 + 
      csetemp290)*gu12) - (csetemp291 + csetemp292 + csetemp293)*gu22 - 
      (csetemp294 + csetemp295 + csetemp296)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu133 = -((csetemp288 + csetemp289 + 
      csetemp290)*gu13) - (csetemp291 + csetemp292 + csetemp293)*gu23 - 
      (csetemp294 + csetemp295 + csetemp296)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp297 = dg113*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp298 = dg123*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp299 = dg133*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp300 = dg223*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp301 = dg233*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp302 = dg133*gu12;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp303 = dg233*gu22;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp304 = dg333*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu223 = -((csetemp297 + csetemp298 + 
      csetemp299)*gu12) - (csetemp289 + csetemp300 + csetemp301)*gu22 - 
      (csetemp302 + csetemp303 + csetemp304)*gu23;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu233 = -((csetemp297 + csetemp298 + 
      csetemp299)*gu13) - (csetemp289 + csetemp300 + csetemp301)*gu23 - 
      (csetemp302 + csetemp303 + csetemp304)*gu33;
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dgu333 = -(gu13*(dg113*gu13 + 
      dg123*gu23 + dg133*gu33)) - gu23*(dg123*gu13 + dg223*gu23 + dg233*gu33) 
      - gu33*(csetemp290 + csetemp301 + dg333*gu33);
    
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
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED csetemp305 = INV(alpL);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED dtalpL = 0.5*csetemp305*(-dg4000 + 
      dtbetasq);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxxL = 
      0.5*csetemp305*(2*(gxxL*dbeta11 + gxyL*dbeta21 + gxzL*dbeta31) + 
      betaxL*dg111 + betayL*dg112 + betazL*dg113 - dtg11);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxyL = 0.5*csetemp305*(gxxL*dbeta12 
      + gyyL*dbeta21 + gxyL*(dbeta11 + dbeta22) + gyzL*dbeta31 + 
      gxzL*dbeta32 + betaxL*dg121 + betayL*dg122 + betazL*dg123 - 
      dtg12);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kxzL = 0.5*csetemp305*(gxxL*dbeta13 
      + gyzL*dbeta21 + gxyL*dbeta23 + gzzL*dbeta31 + gxzL*(dbeta11 + 
      dbeta33) + betaxL*dg131 + betayL*dg132 + betazL*dg133 - dtg13);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kyyL = 
      0.5*csetemp305*(2*(gxyL*dbeta12 + gyyL*dbeta22 + gyzL*dbeta32) + 
      betaxL*dg221 + betayL*dg222 + betazL*dg223 - dtg22);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kyzL = 0.5*csetemp305*(gxzL*dbeta12 
      + gxyL*dbeta13 + gyyL*dbeta23 + gzzL*dbeta32 + gyzL*(dbeta22 + 
      dbeta33) + betaxL*dg231 + betayL*dg232 + betazL*dg233 - dtg23);
    
    CCTK_REAL CCTK_ATTRIBUTE_UNUSED kzzL = 
      0.5*csetemp305*(2*(gxzL*dbeta13 + gyzL*dbeta23 + gzzL*dbeta33) + 
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
  
  const char* const groups[] = {
    "admbase::curv",
    "admbase::dtlapse",
    "admbase::dtshift",
    "admbase::lapse",
    "admbase::metric",
    "admbase::shift",
    "grid::coordinates"};
  GenericFD_AssertGroupStorage(cctkGH, "Vaidya2_always", 7, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, Vaidya2_always_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving Vaidya2_always_Body");
  }
}
