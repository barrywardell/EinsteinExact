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

namespace Vaidya2 {


static void Vaidya2_exact_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = cctk_time;
  const CCTK_REAL cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(0);
  const CCTK_REAL cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(1);
  const CCTK_REAL cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(2);
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_TIME;
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
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
  CCTK_LOOP3(Vaidya2_exact,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpExactL CCTK_ATTRIBUTE_UNUSED = alpExact[index];
    CCTK_REAL betaExact1L CCTK_ATTRIBUTE_UNUSED = betaExact1[index];
    CCTK_REAL betaExact2L CCTK_ATTRIBUTE_UNUSED = betaExact2[index];
    CCTK_REAL betaExact3L CCTK_ATTRIBUTE_UNUSED = betaExact3[index];
    CCTK_REAL dtbetaExact1L CCTK_ATTRIBUTE_UNUSED = dtbetaExact1[index];
    CCTK_REAL dtbetaExact2L CCTK_ATTRIBUTE_UNUSED = dtbetaExact2[index];
    CCTK_REAL dtbetaExact3L CCTK_ATTRIBUTE_UNUSED = dtbetaExact3[index];
    CCTK_REAL gExact11L CCTK_ATTRIBUTE_UNUSED = gExact11[index];
    CCTK_REAL gExact12L CCTK_ATTRIBUTE_UNUSED = gExact12[index];
    CCTK_REAL gExact13L CCTK_ATTRIBUTE_UNUSED = gExact13[index];
    CCTK_REAL gExact22L CCTK_ATTRIBUTE_UNUSED = gExact22[index];
    CCTK_REAL gExact23L CCTK_ATTRIBUTE_UNUSED = gExact23[index];
    CCTK_REAL gExact33L CCTK_ATTRIBUTE_UNUSED = gExact33[index];
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
    CCTK_REAL xform1L00 CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL xform1L01 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL xform1L02 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL xform1L03 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL xform1L10 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL csetemp0 CCTK_ATTRIBUTE_UNUSED = cos(phi);
    
    CCTK_REAL csetemp1 CCTK_ATTRIBUTE_UNUSED = cos(psi);
    
    CCTK_REAL csetemp2 CCTK_ATTRIBUTE_UNUSED = cos(theta);
    
    CCTK_REAL csetemp3 CCTK_ATTRIBUTE_UNUSED = sin(phi);
    
    CCTK_REAL csetemp4 CCTK_ATTRIBUTE_UNUSED = sin(psi);
    
    CCTK_REAL xform1L11 CCTK_ATTRIBUTE_UNUSED = csetemp0*csetemp1 - 
      csetemp2*csetemp3*csetemp4;
    
    CCTK_REAL xform1L12 CCTK_ATTRIBUTE_UNUSED = csetemp1*csetemp3 + 
      csetemp0*csetemp2*csetemp4;
    
    CCTK_REAL csetemp5 CCTK_ATTRIBUTE_UNUSED = sin(theta);
    
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
    
    CCTK_REAL csetemp6 CCTK_ATTRIBUTE_UNUSED = pow(boostx,2);
    
    CCTK_REAL csetemp7 CCTK_ATTRIBUTE_UNUSED = pow(boosty,2);
    
    CCTK_REAL csetemp8 CCTK_ATTRIBUTE_UNUSED = pow(boostz,2);
    
    CCTK_REAL csetemp9 CCTK_ATTRIBUTE_UNUSED = pow(lapsefactor,-1);
    
    CCTK_REAL xform2L00 CCTK_ATTRIBUTE_UNUSED = csetemp9*(-1 + 
      boostx*shiftaddx + boosty*shiftaddy + boostz*shiftaddz)*(1 - csetemp6 - 
      csetemp7 - csetemp8 + pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5))*pow(-1 + csetemp6 + csetemp7 + csetemp8,-1)*pow(1 + 
      pow(1 - csetemp6 - csetemp7 - csetemp8,0.5),-1);
    
    CCTK_REAL xform2L01 CCTK_ATTRIBUTE_UNUSED = boostx*csetemp9*(1 - 
      csetemp6 - csetemp7 - csetemp8 + pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5))*pow(-1 + csetemp6 + csetemp7 + csetemp8,-1)*pow(1 + 
      pow(1 - csetemp6 - csetemp7 - csetemp8,0.5),-1);
    
    CCTK_REAL xform2L02 CCTK_ATTRIBUTE_UNUSED = boosty*csetemp9*(1 - 
      csetemp6 - csetemp7 - csetemp8 + pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5))*pow(-1 + csetemp6 + csetemp7 + csetemp8,-1)*pow(1 + 
      pow(1 - csetemp6 - csetemp7 - csetemp8,0.5),-1);
    
    CCTK_REAL xform2L03 CCTK_ATTRIBUTE_UNUSED = boostz*csetemp9*(1 - 
      csetemp6 - csetemp7 - csetemp8 + pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5))*pow(-1 + csetemp6 + csetemp7 + csetemp8,-1)*pow(1 + 
      pow(1 - csetemp6 - csetemp7 - csetemp8,0.5),-1);
    
    CCTK_REAL xform2L10 CCTK_ATTRIBUTE_UNUSED = (-(boostx*(-1 + csetemp6 + 
      csetemp7 + csetemp8 + (-1 + boosty*shiftaddy + boostz*shiftaddz)*pow(1 
      - csetemp6 - csetemp7 - csetemp8,0.5))) + shiftaddx*(csetemp6 + (-1 + 
      csetemp7 + csetemp8)*(1 + pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5))))*pow(-1 + csetemp6 + csetemp7 + csetemp8,-1)*pow(1 + 
      pow(1 - csetemp6 - csetemp7 - csetemp8,0.5),-1);
    
    CCTK_REAL xform2L11 CCTK_ATTRIBUTE_UNUSED = pow(-1 + csetemp6 + 
      csetemp7 + csetemp8,-1)*(-1 + csetemp7 + csetemp8 + csetemp6*pow(1 + 
      pow(1 - csetemp6 - csetemp7 - csetemp8,0.5),-1));
    
    CCTK_REAL xform2L12 CCTK_ATTRIBUTE_UNUSED = -(boostx*boosty*pow(-1 + 
      csetemp6 + csetemp7 + csetemp8 - pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5),-1));
    
    CCTK_REAL xform2L13 CCTK_ATTRIBUTE_UNUSED = -(boostx*boostz*pow(-1 + 
      csetemp6 + csetemp7 + csetemp8 - pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5),-1));
    
    CCTK_REAL xform2L20 CCTK_ATTRIBUTE_UNUSED = (-(boosty*(-1 + csetemp6 + 
      csetemp7 + csetemp8 + (-1 + boostx*shiftaddx + boostz*shiftaddz)*pow(1 
      - csetemp6 - csetemp7 - csetemp8,0.5))) + shiftaddy*(csetemp7 + (-1 + 
      csetemp6 + csetemp8)*(1 + pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5))))*pow(-1 + csetemp6 + csetemp7 + csetemp8,-1)*pow(1 + 
      pow(1 - csetemp6 - csetemp7 - csetemp8,0.5),-1);
    
    CCTK_REAL xform2L21 CCTK_ATTRIBUTE_UNUSED = -(boostx*boosty*pow(-1 + 
      csetemp6 + csetemp7 + csetemp8 - pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5),-1));
    
    CCTK_REAL xform2L22 CCTK_ATTRIBUTE_UNUSED = pow(-1 + csetemp6 + 
      csetemp7 + csetemp8,-1)*(-1 + csetemp6 + csetemp8 + csetemp7*pow(1 + 
      pow(1 - csetemp6 - csetemp7 - csetemp8,0.5),-1));
    
    CCTK_REAL xform2L23 CCTK_ATTRIBUTE_UNUSED = -(boosty*boostz*pow(-1 + 
      csetemp6 + csetemp7 + csetemp8 - pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5),-1));
    
    CCTK_REAL xform2L30 CCTK_ATTRIBUTE_UNUSED = (shiftaddz*(-1 + csetemp6 
      + csetemp7 + csetemp8 + (-1 + csetemp6 + csetemp7)*pow(1 - csetemp6 - 
      csetemp7 - csetemp8,0.5)) - boostz*(-1 + csetemp6 + csetemp7 + csetemp8 
      + (-1 + boostx*shiftaddx + boosty*shiftaddy)*pow(1 - csetemp6 - 
      csetemp7 - csetemp8,0.5)))*pow(-1 + csetemp6 + csetemp7 + 
      csetemp8,-1)*pow(1 + pow(1 - csetemp6 - csetemp7 - csetemp8,0.5),-1);
    
    CCTK_REAL xform2L31 CCTK_ATTRIBUTE_UNUSED = -(boostx*boostz*pow(-1 + 
      csetemp6 + csetemp7 + csetemp8 - pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5),-1));
    
    CCTK_REAL xform2L32 CCTK_ATTRIBUTE_UNUSED = -(boosty*boostz*pow(-1 + 
      csetemp6 + csetemp7 + csetemp8 - pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5),-1));
    
    CCTK_REAL xform2L33 CCTK_ATTRIBUTE_UNUSED = (-1 + csetemp6 + csetemp7 
      + csetemp8 + (-1 + csetemp6 + csetemp7)*pow(1 - csetemp6 - csetemp7 - 
      csetemp8,0.5))*pow(-1 + csetemp6 + csetemp7 + csetemp8,-1)*pow(1 + 
      pow(1 - csetemp6 - csetemp7 - csetemp8,0.5),-1);
    
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
    
    CCTK_REAL xx0 CCTK_ATTRIBUTE_UNUSED = t - timeoffset;
    
    CCTK_REAL xx1 CCTK_ATTRIBUTE_UNUSED = xL - positionx;
    
    CCTK_REAL xx2 CCTK_ATTRIBUTE_UNUSED = yL - positiony;
    
    CCTK_REAL xx3 CCTK_ATTRIBUTE_UNUSED = zL - positionz;
    
    CCTK_REAL txx0 CCTK_ATTRIBUTE_UNUSED = xformL00*xx0 + xformL01*xx1 + 
      xformL02*xx2 + xformL03*xx3;
    
    CCTK_REAL txx1 CCTK_ATTRIBUTE_UNUSED = xformL10*xx0 + xformL11*xx1 + 
      xformL12*xx2 + xformL13*xx3;
    
    CCTK_REAL txx2 CCTK_ATTRIBUTE_UNUSED = xformL20*xx0 + xformL21*xx1 + 
      xformL22*xx2 + xformL23*xx3;
    
    CCTK_REAL txx3 CCTK_ATTRIBUTE_UNUSED = xformL30*xx0 + xformL31*xx1 + 
      xformL32*xx2 + xformL33*xx3;
    
    CCTK_REAL T CCTK_ATTRIBUTE_UNUSED = txx0;
    
    CCTK_REAL X CCTK_ATTRIBUTE_UNUSED = txx1;
    
    CCTK_REAL Y CCTK_ATTRIBUTE_UNUSED = txx2;
    
    CCTK_REAL Z CCTK_ATTRIBUTE_UNUSED = txx3;
    
    CCTK_REAL csetemp10 CCTK_ATTRIBUTE_UNUSED = pow(X,2);
    
    CCTK_REAL csetemp11 CCTK_ATTRIBUTE_UNUSED = pow(Y,2);
    
    CCTK_REAL csetemp12 CCTK_ATTRIBUTE_UNUSED = pow(Z,2);
    
    CCTK_REAL rXYZ CCTK_ATTRIBUTE_UNUSED = pow(csetemp10 + csetemp11 + 
      csetemp12,0.5);
    
    CCTK_REAL csetemp13 CCTK_ATTRIBUTE_UNUSED = pow(M,-1);
    
    CCTK_REAL mTXYZ CCTK_ATTRIBUTE_UNUSED = M*(1 + 
      pow(tanh(csetemp13*dM*(T + pow(csetemp10 + csetemp11 + 
      csetemp12,0.5))),2));
    
    CCTK_REAL csetemp14 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,-1);
    
    CCTK_REAL tg400 CCTK_ATTRIBUTE_UNUSED = -1 + 2*csetemp14*mTXYZ;
    
    CCTK_REAL csetemp15 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,-2);
    
    CCTK_REAL tg401 CCTK_ATTRIBUTE_UNUSED = 2*csetemp15*mTXYZ*X;
    
    CCTK_REAL tg402 CCTK_ATTRIBUTE_UNUSED = 2*csetemp15*mTXYZ*Y;
    
    CCTK_REAL tg403 CCTK_ATTRIBUTE_UNUSED = 2*csetemp15*mTXYZ*Z;
    
    CCTK_REAL csetemp16 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,-3);
    
    CCTK_REAL tg411 CCTK_ATTRIBUTE_UNUSED = 1 + 
      2*csetemp10*csetemp16*mTXYZ;
    
    CCTK_REAL tg412 CCTK_ATTRIBUTE_UNUSED = 2*csetemp16*mTXYZ*X*Y;
    
    CCTK_REAL tg413 CCTK_ATTRIBUTE_UNUSED = 2*csetemp16*mTXYZ*X*Z;
    
    CCTK_REAL tg422 CCTK_ATTRIBUTE_UNUSED = 1 + 
      2*csetemp11*csetemp16*mTXYZ;
    
    CCTK_REAL tg423 CCTK_ATTRIBUTE_UNUSED = 2*csetemp16*mTXYZ*Y*Z;
    
    CCTK_REAL tg433 CCTK_ATTRIBUTE_UNUSED = 1 + 
      2*csetemp12*csetemp16*mTXYZ;
    
    CCTK_REAL csetemp17 CCTK_ATTRIBUTE_UNUSED = rXYZ + T;
    
    CCTK_REAL tdg4000 CCTK_ATTRIBUTE_UNUSED = 
      4*csetemp14*dM*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM);
    
    CCTK_REAL tdg4001 CCTK_ATTRIBUTE_UNUSED = -2*csetemp16*X*(mTXYZ - 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4002 CCTK_ATTRIBUTE_UNUSED = -2*csetemp16*Y*(mTXYZ - 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4003 CCTK_ATTRIBUTE_UNUSED = -2*csetemp16*Z*(mTXYZ - 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4010 CCTK_ATTRIBUTE_UNUSED = 
      4*csetemp15*dM*X*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM);
    
    CCTK_REAL csetemp18 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,-4);
    
    CCTK_REAL csetemp19 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,2);
    
    CCTK_REAL tdg4011 CCTK_ATTRIBUTE_UNUSED = csetemp18*(2*(-2*csetemp10 + 
      csetemp19)*mTXYZ + 
      4*csetemp10*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4012 CCTK_ATTRIBUTE_UNUSED = -4*csetemp18*X*Y*(mTXYZ - 
      dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4013 CCTK_ATTRIBUTE_UNUSED = -4*csetemp18*X*Z*(mTXYZ - 
      dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4020 CCTK_ATTRIBUTE_UNUSED = 
      4*csetemp15*dM*Y*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM);
    
    CCTK_REAL tdg4021 CCTK_ATTRIBUTE_UNUSED = -4*csetemp18*X*Y*(mTXYZ - 
      dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4022 CCTK_ATTRIBUTE_UNUSED = csetemp18*(2*(-2*csetemp11 + 
      csetemp19)*mTXYZ + 
      4*csetemp11*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4023 CCTK_ATTRIBUTE_UNUSED = -4*csetemp18*Y*Z*(mTXYZ - 
      dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4030 CCTK_ATTRIBUTE_UNUSED = 
      4*csetemp15*dM*Z*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM);
    
    CCTK_REAL tdg4031 CCTK_ATTRIBUTE_UNUSED = -4*csetemp18*X*Z*(mTXYZ - 
      dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4032 CCTK_ATTRIBUTE_UNUSED = -4*csetemp18*Y*Z*(mTXYZ - 
      dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4033 CCTK_ATTRIBUTE_UNUSED = csetemp18*(2*(-2*csetemp12 + 
      csetemp19)*mTXYZ + 
      4*csetemp12*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4110 CCTK_ATTRIBUTE_UNUSED = 
      4*csetemp10*csetemp16*dM*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM);
    
    CCTK_REAL csetemp20 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,-5);
    
    CCTK_REAL csetemp21 CCTK_ATTRIBUTE_UNUSED = pow(X,3);
    
    CCTK_REAL tdg4111 CCTK_ATTRIBUTE_UNUSED = 
      csetemp20*(mTXYZ*(-6*csetemp21 + 4*csetemp19*X) + 
      4*csetemp21*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL csetemp22 CCTK_ATTRIBUTE_UNUSED = -3*mTXYZ;
    
    CCTK_REAL tdg4112 CCTK_ATTRIBUTE_UNUSED = 
      2*csetemp10*csetemp20*Y*(csetemp22 + 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4113 CCTK_ATTRIBUTE_UNUSED = 
      2*csetemp10*csetemp20*Z*(csetemp22 + 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4120 CCTK_ATTRIBUTE_UNUSED = 
      4*csetemp16*dM*X*Y*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM);
    
    CCTK_REAL tdg4121 CCTK_ATTRIBUTE_UNUSED = 2*csetemp20*Y*((-2*csetemp10 
      + csetemp11 + csetemp12)*mTXYZ + 
      2*csetemp10*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4122 CCTK_ATTRIBUTE_UNUSED = 2*csetemp20*X*((csetemp10 - 
      2*csetemp11 + csetemp12)*mTXYZ + 
      2*csetemp11*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4123 CCTK_ATTRIBUTE_UNUSED = 2*csetemp20*X*Y*Z*(csetemp22 
      + 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4130 CCTK_ATTRIBUTE_UNUSED = 
      4*csetemp16*dM*X*Z*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM);
    
    CCTK_REAL tdg4131 CCTK_ATTRIBUTE_UNUSED = 2*csetemp20*Z*((-2*csetemp10 
      + csetemp11 + csetemp12)*mTXYZ + 
      2*csetemp10*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4132 CCTK_ATTRIBUTE_UNUSED = 2*csetemp20*X*Y*Z*(csetemp22 
      + 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4133 CCTK_ATTRIBUTE_UNUSED = 2*csetemp20*X*((csetemp10 + 
      csetemp11 - 2*csetemp12)*mTXYZ + 
      2*csetemp12*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4220 CCTK_ATTRIBUTE_UNUSED = 
      4*csetemp11*csetemp16*dM*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM);
    
    CCTK_REAL tdg4221 CCTK_ATTRIBUTE_UNUSED = 
      2*csetemp11*csetemp20*X*(csetemp22 + 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL csetemp23 CCTK_ATTRIBUTE_UNUSED = pow(Y,3);
    
    CCTK_REAL tdg4222 CCTK_ATTRIBUTE_UNUSED = 
      csetemp20*(mTXYZ*(-6*csetemp23 + 4*csetemp19*Y) + 
      4*csetemp23*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4223 CCTK_ATTRIBUTE_UNUSED = 
      2*csetemp11*csetemp20*Z*(csetemp22 + 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4230 CCTK_ATTRIBUTE_UNUSED = 
      4*csetemp16*dM*Y*Z*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM);
    
    CCTK_REAL tdg4231 CCTK_ATTRIBUTE_UNUSED = 2*csetemp20*X*Y*Z*(csetemp22 
      + 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4232 CCTK_ATTRIBUTE_UNUSED = 2*csetemp20*Z*((csetemp10 - 
      2*csetemp11 + csetemp12)*mTXYZ + 
      2*csetemp11*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4233 CCTK_ATTRIBUTE_UNUSED = 2*csetemp20*Y*((csetemp10 + 
      csetemp11 - 2*csetemp12)*mTXYZ + 
      2*csetemp12*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4330 CCTK_ATTRIBUTE_UNUSED = 
      4*csetemp12*csetemp16*dM*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM);
    
    CCTK_REAL tdg4331 CCTK_ATTRIBUTE_UNUSED = 
      2*csetemp12*csetemp20*X*(csetemp22 + 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL tdg4332 CCTK_ATTRIBUTE_UNUSED = 
      2*csetemp12*csetemp20*Y*(csetemp22 + 
      2*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL csetemp24 CCTK_ATTRIBUTE_UNUSED = pow(Z,3);
    
    CCTK_REAL tdg4333 CCTK_ATTRIBUTE_UNUSED = 
      csetemp20*(mTXYZ*(-6*csetemp24 + 4*csetemp19*Z) + 
      4*csetemp24*dM*rXYZ*pow(pow(cosh(csetemp13*csetemp17*dM),-1),2)*tanh(csetemp13*csetemp17*dM));
    
    CCTK_REAL g400 CCTK_ATTRIBUTE_UNUSED = 2*(tg423*xformL20*xformL30 + 
      xformL00*(tg401*xformL10 + tg402*xformL20 + tg403*xformL30) + 
      xformL10*(tg412*xformL20 + tg413*xformL30)) + tg400*pow(xformL00,2) + 
      tg411*pow(xformL10,2) + tg422*pow(xformL20,2) + tg433*pow(xformL30,2);
    
    CCTK_REAL csetemp25 CCTK_ATTRIBUTE_UNUSED = tg400*xformL01;
    
    CCTK_REAL csetemp26 CCTK_ATTRIBUTE_UNUSED = tg401*xformL11;
    
    CCTK_REAL csetemp27 CCTK_ATTRIBUTE_UNUSED = tg402*xformL21;
    
    CCTK_REAL csetemp28 CCTK_ATTRIBUTE_UNUSED = tg403*xformL31;
    
    CCTK_REAL csetemp29 CCTK_ATTRIBUTE_UNUSED = tg401*xformL01;
    
    CCTK_REAL csetemp30 CCTK_ATTRIBUTE_UNUSED = tg411*xformL11;
    
    CCTK_REAL csetemp31 CCTK_ATTRIBUTE_UNUSED = tg412*xformL21;
    
    CCTK_REAL csetemp32 CCTK_ATTRIBUTE_UNUSED = tg413*xformL31;
    
    CCTK_REAL csetemp33 CCTK_ATTRIBUTE_UNUSED = tg402*xformL01;
    
    CCTK_REAL csetemp34 CCTK_ATTRIBUTE_UNUSED = tg412*xformL11;
    
    CCTK_REAL csetemp35 CCTK_ATTRIBUTE_UNUSED = tg422*xformL21;
    
    CCTK_REAL csetemp36 CCTK_ATTRIBUTE_UNUSED = tg423*xformL31;
    
    CCTK_REAL csetemp37 CCTK_ATTRIBUTE_UNUSED = tg403*xformL01;
    
    CCTK_REAL csetemp38 CCTK_ATTRIBUTE_UNUSED = tg413*xformL11;
    
    CCTK_REAL csetemp39 CCTK_ATTRIBUTE_UNUSED = tg423*xformL21;
    
    CCTK_REAL csetemp40 CCTK_ATTRIBUTE_UNUSED = tg433*xformL31;
    
    CCTK_REAL g401 CCTK_ATTRIBUTE_UNUSED = (csetemp25 + csetemp26 + 
      csetemp27 + csetemp28)*xformL00 + (csetemp29 + csetemp30 + csetemp31 + 
      csetemp32)*xformL10 + (csetemp33 + csetemp34 + csetemp35 + 
      csetemp36)*xformL20 + (csetemp37 + csetemp38 + csetemp39 + 
      csetemp40)*xformL30;
    
    CCTK_REAL csetemp41 CCTK_ATTRIBUTE_UNUSED = tg400*xformL02;
    
    CCTK_REAL csetemp42 CCTK_ATTRIBUTE_UNUSED = tg401*xformL12;
    
    CCTK_REAL csetemp43 CCTK_ATTRIBUTE_UNUSED = tg402*xformL22;
    
    CCTK_REAL csetemp44 CCTK_ATTRIBUTE_UNUSED = tg403*xformL32;
    
    CCTK_REAL csetemp45 CCTK_ATTRIBUTE_UNUSED = tg401*xformL02;
    
    CCTK_REAL csetemp46 CCTK_ATTRIBUTE_UNUSED = tg411*xformL12;
    
    CCTK_REAL csetemp47 CCTK_ATTRIBUTE_UNUSED = tg412*xformL22;
    
    CCTK_REAL csetemp48 CCTK_ATTRIBUTE_UNUSED = tg413*xformL32;
    
    CCTK_REAL csetemp49 CCTK_ATTRIBUTE_UNUSED = tg402*xformL02;
    
    CCTK_REAL csetemp50 CCTK_ATTRIBUTE_UNUSED = tg412*xformL12;
    
    CCTK_REAL csetemp51 CCTK_ATTRIBUTE_UNUSED = tg422*xformL22;
    
    CCTK_REAL csetemp52 CCTK_ATTRIBUTE_UNUSED = tg423*xformL32;
    
    CCTK_REAL csetemp53 CCTK_ATTRIBUTE_UNUSED = tg403*xformL02;
    
    CCTK_REAL csetemp54 CCTK_ATTRIBUTE_UNUSED = tg413*xformL12;
    
    CCTK_REAL csetemp55 CCTK_ATTRIBUTE_UNUSED = tg423*xformL22;
    
    CCTK_REAL csetemp56 CCTK_ATTRIBUTE_UNUSED = tg433*xformL32;
    
    CCTK_REAL g402 CCTK_ATTRIBUTE_UNUSED = (csetemp41 + csetemp42 + 
      csetemp43 + csetemp44)*xformL00 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL10 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL20 + (csetemp53 + csetemp54 + csetemp55 + 
      csetemp56)*xformL30;
    
    CCTK_REAL csetemp57 CCTK_ATTRIBUTE_UNUSED = tg400*xformL03;
    
    CCTK_REAL csetemp58 CCTK_ATTRIBUTE_UNUSED = tg401*xformL13;
    
    CCTK_REAL csetemp59 CCTK_ATTRIBUTE_UNUSED = tg402*xformL23;
    
    CCTK_REAL csetemp60 CCTK_ATTRIBUTE_UNUSED = tg403*xformL33;
    
    CCTK_REAL csetemp61 CCTK_ATTRIBUTE_UNUSED = tg401*xformL03;
    
    CCTK_REAL csetemp62 CCTK_ATTRIBUTE_UNUSED = tg411*xformL13;
    
    CCTK_REAL csetemp63 CCTK_ATTRIBUTE_UNUSED = tg412*xformL23;
    
    CCTK_REAL csetemp64 CCTK_ATTRIBUTE_UNUSED = tg413*xformL33;
    
    CCTK_REAL csetemp65 CCTK_ATTRIBUTE_UNUSED = tg402*xformL03;
    
    CCTK_REAL csetemp66 CCTK_ATTRIBUTE_UNUSED = tg412*xformL13;
    
    CCTK_REAL csetemp67 CCTK_ATTRIBUTE_UNUSED = tg422*xformL23;
    
    CCTK_REAL csetemp68 CCTK_ATTRIBUTE_UNUSED = tg423*xformL33;
    
    CCTK_REAL csetemp69 CCTK_ATTRIBUTE_UNUSED = tg403*xformL03;
    
    CCTK_REAL csetemp70 CCTK_ATTRIBUTE_UNUSED = tg413*xformL13;
    
    CCTK_REAL csetemp71 CCTK_ATTRIBUTE_UNUSED = tg423*xformL23;
    
    CCTK_REAL csetemp72 CCTK_ATTRIBUTE_UNUSED = tg433*xformL33;
    
    CCTK_REAL g403 CCTK_ATTRIBUTE_UNUSED = (csetemp57 + csetemp58 + 
      csetemp59 + csetemp60)*xformL00 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL10 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL20 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL30;
    
    CCTK_REAL g411 CCTK_ATTRIBUTE_UNUSED = (csetemp25 + csetemp26 + 
      csetemp27 + csetemp28)*xformL01 + (csetemp29 + csetemp30 + csetemp31 + 
      csetemp32)*xformL11 + (csetemp33 + csetemp34 + csetemp35 + 
      csetemp36)*xformL21 + (csetemp37 + csetemp38 + csetemp39 + 
      csetemp40)*xformL31;
    
    CCTK_REAL g412 CCTK_ATTRIBUTE_UNUSED = (csetemp41 + csetemp42 + 
      csetemp43 + csetemp44)*xformL01 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL11 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL21 + (csetemp53 + csetemp54 + csetemp55 + 
      csetemp56)*xformL31;
    
    CCTK_REAL g413 CCTK_ATTRIBUTE_UNUSED = (csetemp57 + csetemp58 + 
      csetemp59 + csetemp60)*xformL01 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL11 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL21 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL31;
    
    CCTK_REAL g422 CCTK_ATTRIBUTE_UNUSED = (csetemp41 + csetemp42 + 
      csetemp43 + csetemp44)*xformL02 + (csetemp45 + csetemp46 + csetemp47 + 
      csetemp48)*xformL12 + (csetemp49 + csetemp50 + csetemp51 + 
      csetemp52)*xformL22 + (csetemp53 + csetemp54 + csetemp55 + 
      csetemp56)*xformL32;
    
    CCTK_REAL g423 CCTK_ATTRIBUTE_UNUSED = (csetemp57 + csetemp58 + 
      csetemp59 + csetemp60)*xformL02 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL12 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL22 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL32;
    
    CCTK_REAL g433 CCTK_ATTRIBUTE_UNUSED = (csetemp57 + csetemp58 + 
      csetemp59 + csetemp60)*xformL03 + (csetemp61 + csetemp62 + csetemp63 + 
      csetemp64)*xformL13 + (csetemp65 + csetemp66 + csetemp67 + 
      csetemp68)*xformL23 + (csetemp69 + csetemp70 + csetemp71 + 
      csetemp72)*xformL33;
    
    CCTK_REAL csetemp73 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL00;
    
    CCTK_REAL csetemp74 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL10;
    
    CCTK_REAL csetemp75 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL20;
    
    CCTK_REAL csetemp76 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL30;
    
    CCTK_REAL csetemp77 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL00;
    
    CCTK_REAL csetemp78 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL10;
    
    CCTK_REAL csetemp79 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL20;
    
    CCTK_REAL csetemp80 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL30;
    
    CCTK_REAL csetemp81 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL00;
    
    CCTK_REAL csetemp82 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL10;
    
    CCTK_REAL csetemp83 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL20;
    
    CCTK_REAL csetemp84 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL30;
    
    CCTK_REAL csetemp85 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL00;
    
    CCTK_REAL csetemp86 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL10;
    
    CCTK_REAL csetemp87 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL20;
    
    CCTK_REAL csetemp88 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL30;
    
    CCTK_REAL csetemp89 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL00;
    
    CCTK_REAL csetemp90 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL10;
    
    CCTK_REAL csetemp91 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL20;
    
    CCTK_REAL csetemp92 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL30;
    
    CCTK_REAL csetemp93 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL00;
    
    CCTK_REAL csetemp94 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL10;
    
    CCTK_REAL csetemp95 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL20;
    
    CCTK_REAL csetemp96 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL30;
    
    CCTK_REAL csetemp97 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL00;
    
    CCTK_REAL csetemp98 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL10;
    
    CCTK_REAL csetemp99 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL20;
    
    CCTK_REAL csetemp100 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL30;
    
    CCTK_REAL csetemp101 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL00;
    
    CCTK_REAL csetemp102 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL10;
    
    CCTK_REAL csetemp103 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL20;
    
    CCTK_REAL csetemp104 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL30;
    
    CCTK_REAL csetemp105 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL00;
    
    CCTK_REAL csetemp106 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL10;
    
    CCTK_REAL csetemp107 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL20;
    
    CCTK_REAL csetemp108 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL30;
    
    CCTK_REAL csetemp109 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL00;
    
    CCTK_REAL csetemp110 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL10;
    
    CCTK_REAL csetemp111 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL20;
    
    CCTK_REAL csetemp112 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL30;
    
    CCTK_REAL dg4000 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp105 + csetemp106 + 
      csetemp107 + csetemp108)*xformL20*xformL30 + xformL00*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL10 + (csetemp81 + csetemp82 + 
      csetemp83 + csetemp84)*xformL20 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL30) + xformL10*((csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL20 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL30)) + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*pow(xformL00,2) + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*pow(xformL10,2) + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*pow(xformL20,2) + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*pow(xformL30,2);
    
    CCTK_REAL dg4010 CCTK_ATTRIBUTE_UNUSED = xformL20*((csetemp81 + 
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
    
    CCTK_REAL csetemp113 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL01;
    
    CCTK_REAL csetemp114 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL11;
    
    CCTK_REAL csetemp115 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL21;
    
    CCTK_REAL csetemp116 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL31;
    
    CCTK_REAL csetemp117 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL01;
    
    CCTK_REAL csetemp118 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL11;
    
    CCTK_REAL csetemp119 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL21;
    
    CCTK_REAL csetemp120 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL31;
    
    CCTK_REAL csetemp121 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL01;
    
    CCTK_REAL csetemp122 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL11;
    
    CCTK_REAL csetemp123 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL21;
    
    CCTK_REAL csetemp124 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL31;
    
    CCTK_REAL csetemp125 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL01;
    
    CCTK_REAL csetemp126 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL11;
    
    CCTK_REAL csetemp127 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL21;
    
    CCTK_REAL csetemp128 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL31;
    
    CCTK_REAL csetemp129 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL01;
    
    CCTK_REAL csetemp130 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL11;
    
    CCTK_REAL csetemp131 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL21;
    
    CCTK_REAL csetemp132 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL31;
    
    CCTK_REAL csetemp133 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL01;
    
    CCTK_REAL csetemp134 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL11;
    
    CCTK_REAL csetemp135 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL21;
    
    CCTK_REAL csetemp136 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL31;
    
    CCTK_REAL csetemp137 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL01;
    
    CCTK_REAL csetemp138 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL11;
    
    CCTK_REAL csetemp139 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL21;
    
    CCTK_REAL csetemp140 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL31;
    
    CCTK_REAL csetemp141 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL01;
    
    CCTK_REAL csetemp142 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL11;
    
    CCTK_REAL csetemp143 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL21;
    
    CCTK_REAL csetemp144 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL31;
    
    CCTK_REAL csetemp145 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL01;
    
    CCTK_REAL csetemp146 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL11;
    
    CCTK_REAL csetemp147 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL21;
    
    CCTK_REAL csetemp148 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL31;
    
    CCTK_REAL csetemp149 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL01;
    
    CCTK_REAL csetemp150 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL11;
    
    CCTK_REAL csetemp151 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL21;
    
    CCTK_REAL csetemp152 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL31;
    
    CCTK_REAL dg4011 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp113 + 
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
    
    CCTK_REAL csetemp153 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL02;
    
    CCTK_REAL csetemp154 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL12;
    
    CCTK_REAL csetemp155 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL22;
    
    CCTK_REAL csetemp156 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL32;
    
    CCTK_REAL csetemp157 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL02;
    
    CCTK_REAL csetemp158 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL12;
    
    CCTK_REAL csetemp159 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL22;
    
    CCTK_REAL csetemp160 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL32;
    
    CCTK_REAL csetemp161 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL02;
    
    CCTK_REAL csetemp162 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL12;
    
    CCTK_REAL csetemp163 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL22;
    
    CCTK_REAL csetemp164 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL32;
    
    CCTK_REAL csetemp165 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL02;
    
    CCTK_REAL csetemp166 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL12;
    
    CCTK_REAL csetemp167 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL22;
    
    CCTK_REAL csetemp168 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL32;
    
    CCTK_REAL csetemp169 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL02;
    
    CCTK_REAL csetemp170 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL12;
    
    CCTK_REAL csetemp171 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL22;
    
    CCTK_REAL csetemp172 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL32;
    
    CCTK_REAL csetemp173 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL02;
    
    CCTK_REAL csetemp174 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL12;
    
    CCTK_REAL csetemp175 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL22;
    
    CCTK_REAL csetemp176 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL32;
    
    CCTK_REAL csetemp177 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL02;
    
    CCTK_REAL csetemp178 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL12;
    
    CCTK_REAL csetemp179 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL22;
    
    CCTK_REAL csetemp180 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL32;
    
    CCTK_REAL csetemp181 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL02;
    
    CCTK_REAL csetemp182 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL12;
    
    CCTK_REAL csetemp183 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL22;
    
    CCTK_REAL csetemp184 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL32;
    
    CCTK_REAL csetemp185 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL02;
    
    CCTK_REAL csetemp186 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL12;
    
    CCTK_REAL csetemp187 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL22;
    
    CCTK_REAL csetemp188 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL32;
    
    CCTK_REAL csetemp189 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL02;
    
    CCTK_REAL csetemp190 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL12;
    
    CCTK_REAL csetemp191 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL22;
    
    CCTK_REAL csetemp192 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL32;
    
    CCTK_REAL dg4012 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp153 + 
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
    
    CCTK_REAL csetemp193 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL03;
    
    CCTK_REAL csetemp194 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL13;
    
    CCTK_REAL csetemp195 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL23;
    
    CCTK_REAL csetemp196 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL33;
    
    CCTK_REAL csetemp197 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL03;
    
    CCTK_REAL csetemp198 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL13;
    
    CCTK_REAL csetemp199 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL23;
    
    CCTK_REAL csetemp200 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL33;
    
    CCTK_REAL csetemp201 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL03;
    
    CCTK_REAL csetemp202 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL13;
    
    CCTK_REAL csetemp203 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL23;
    
    CCTK_REAL csetemp204 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL33;
    
    CCTK_REAL csetemp205 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL03;
    
    CCTK_REAL csetemp206 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL13;
    
    CCTK_REAL csetemp207 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL23;
    
    CCTK_REAL csetemp208 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL33;
    
    CCTK_REAL csetemp209 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL03;
    
    CCTK_REAL csetemp210 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL13;
    
    CCTK_REAL csetemp211 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL23;
    
    CCTK_REAL csetemp212 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL33;
    
    CCTK_REAL csetemp213 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL03;
    
    CCTK_REAL csetemp214 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL13;
    
    CCTK_REAL csetemp215 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL23;
    
    CCTK_REAL csetemp216 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL33;
    
    CCTK_REAL csetemp217 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL03;
    
    CCTK_REAL csetemp218 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL13;
    
    CCTK_REAL csetemp219 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL23;
    
    CCTK_REAL csetemp220 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL33;
    
    CCTK_REAL csetemp221 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL03;
    
    CCTK_REAL csetemp222 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL13;
    
    CCTK_REAL csetemp223 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL23;
    
    CCTK_REAL csetemp224 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL33;
    
    CCTK_REAL csetemp225 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL03;
    
    CCTK_REAL csetemp226 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL13;
    
    CCTK_REAL csetemp227 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL23;
    
    CCTK_REAL csetemp228 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL33;
    
    CCTK_REAL csetemp229 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL03;
    
    CCTK_REAL csetemp230 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL13;
    
    CCTK_REAL csetemp231 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL23;
    
    CCTK_REAL csetemp232 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL33;
    
    CCTK_REAL dg4013 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp193 + 
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
    
    CCTK_REAL dg4020 CCTK_ATTRIBUTE_UNUSED = xformL20*((csetemp81 + 
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
    
    CCTK_REAL dg4021 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp113 + 
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
    
    CCTK_REAL dg4022 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp153 + 
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
    
    CCTK_REAL dg4023 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp193 + 
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
    
    CCTK_REAL dg4030 CCTK_ATTRIBUTE_UNUSED = xformL20*((csetemp81 + 
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
    
    CCTK_REAL dg4031 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp113 + 
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
    
    CCTK_REAL dg4032 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp153 + 
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
    
    CCTK_REAL dg4033 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp193 + 
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
    
    CCTK_REAL dg4110 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp105 + csetemp106 + 
      csetemp107 + csetemp108)*xformL21*xformL31 + xformL01*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL11 + (csetemp81 + csetemp82 + 
      csetemp83 + csetemp84)*xformL21 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL31) + xformL11*((csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL21 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL31)) + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*pow(xformL01,2) + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*pow(xformL11,2) + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*pow(xformL21,2) + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*pow(xformL31,2);
    
    CCTK_REAL dg4111 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp145 + csetemp146 + 
      csetemp147 + csetemp148)*xformL21*xformL31 + xformL01*((csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL11 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL21 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL31) + xformL11*((csetemp133 
      + csetemp134 + csetemp135 + csetemp136)*xformL21 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL31)) + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*pow(xformL01,2) + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*pow(xformL11,2) + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*pow(xformL21,2) + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*pow(xformL31,2);
    
    CCTK_REAL dg4112 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp185 + csetemp186 + 
      csetemp187 + csetemp188)*xformL21*xformL31 + xformL01*((csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL11 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL21 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL31) + xformL11*((csetemp173 
      + csetemp174 + csetemp175 + csetemp176)*xformL21 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL31)) + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*pow(xformL01,2) + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*pow(xformL11,2) + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*pow(xformL21,2) + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*pow(xformL31,2);
    
    CCTK_REAL dg4113 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp225 + csetemp226 + 
      csetemp227 + csetemp228)*xformL21*xformL31 + xformL01*((csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL11 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL21 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL31) + xformL11*((csetemp213 
      + csetemp214 + csetemp215 + csetemp216)*xformL21 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL31)) + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*pow(xformL01,2) + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*pow(xformL11,2) + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*pow(xformL21,2) + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*pow(xformL31,2);
    
    CCTK_REAL dg4120 CCTK_ATTRIBUTE_UNUSED = xformL21*((csetemp81 + 
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
    
    CCTK_REAL dg4121 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp113 + 
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
    
    CCTK_REAL dg4122 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp153 + 
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
    
    CCTK_REAL dg4123 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp193 + 
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
    
    CCTK_REAL dg4130 CCTK_ATTRIBUTE_UNUSED = xformL21*((csetemp81 + 
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
    
    CCTK_REAL dg4131 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp113 + 
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
    
    CCTK_REAL dg4132 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp153 + 
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
    
    CCTK_REAL dg4133 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp193 + 
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
    
    CCTK_REAL dg4220 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp105 + csetemp106 + 
      csetemp107 + csetemp108)*xformL22*xformL32 + xformL02*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL12 + (csetemp81 + csetemp82 + 
      csetemp83 + csetemp84)*xformL22 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL32) + xformL12*((csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL22 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL32)) + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*pow(xformL02,2) + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*pow(xformL12,2) + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*pow(xformL22,2) + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*pow(xformL32,2);
    
    CCTK_REAL dg4221 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp145 + csetemp146 + 
      csetemp147 + csetemp148)*xformL22*xformL32 + xformL02*((csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL12 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL22 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL32) + xformL12*((csetemp133 
      + csetemp134 + csetemp135 + csetemp136)*xformL22 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL32)) + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*pow(xformL02,2) + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*pow(xformL12,2) + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*pow(xformL22,2) + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*pow(xformL32,2);
    
    CCTK_REAL dg4222 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp185 + csetemp186 + 
      csetemp187 + csetemp188)*xformL22*xformL32 + xformL02*((csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL12 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL22 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL32) + xformL12*((csetemp173 
      + csetemp174 + csetemp175 + csetemp176)*xformL22 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL32)) + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*pow(xformL02,2) + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*pow(xformL12,2) + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*pow(xformL22,2) + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*pow(xformL32,2);
    
    CCTK_REAL dg4223 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp225 + csetemp226 + 
      csetemp227 + csetemp228)*xformL22*xformL32 + xformL02*((csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL12 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL22 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL32) + xformL12*((csetemp213 
      + csetemp214 + csetemp215 + csetemp216)*xformL22 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL32)) + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*pow(xformL02,2) + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*pow(xformL12,2) + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*pow(xformL22,2) + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*pow(xformL32,2);
    
    CCTK_REAL dg4230 CCTK_ATTRIBUTE_UNUSED = xformL22*((csetemp81 + 
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
    
    CCTK_REAL dg4231 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp113 + 
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
    
    CCTK_REAL dg4232 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp153 + 
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
    
    CCTK_REAL dg4233 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp193 + 
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
    
    CCTK_REAL dg4330 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp105 + csetemp106 + 
      csetemp107 + csetemp108)*xformL23*xformL33 + xformL03*((csetemp77 + 
      csetemp78 + csetemp79 + csetemp80)*xformL13 + (csetemp81 + csetemp82 + 
      csetemp83 + csetemp84)*xformL23 + (csetemp85 + csetemp86 + csetemp87 + 
      csetemp88)*xformL33) + xformL13*((csetemp93 + csetemp94 + csetemp95 + 
      csetemp96)*xformL23 + (csetemp100 + csetemp97 + csetemp98 + 
      csetemp99)*xformL33)) + (csetemp73 + csetemp74 + csetemp75 + 
      csetemp76)*pow(xformL03,2) + (csetemp89 + csetemp90 + csetemp91 + 
      csetemp92)*pow(xformL13,2) + (csetemp101 + csetemp102 + csetemp103 + 
      csetemp104)*pow(xformL23,2) + (csetemp109 + csetemp110 + csetemp111 + 
      csetemp112)*pow(xformL33,2);
    
    CCTK_REAL dg4331 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp145 + csetemp146 + 
      csetemp147 + csetemp148)*xformL23*xformL33 + xformL03*((csetemp117 + 
      csetemp118 + csetemp119 + csetemp120)*xformL13 + (csetemp121 + 
      csetemp122 + csetemp123 + csetemp124)*xformL23 + (csetemp125 + 
      csetemp126 + csetemp127 + csetemp128)*xformL33) + xformL13*((csetemp133 
      + csetemp134 + csetemp135 + csetemp136)*xformL23 + (csetemp137 + 
      csetemp138 + csetemp139 + csetemp140)*xformL33)) + (csetemp113 + 
      csetemp114 + csetemp115 + csetemp116)*pow(xformL03,2) + (csetemp129 + 
      csetemp130 + csetemp131 + csetemp132)*pow(xformL13,2) + (csetemp141 + 
      csetemp142 + csetemp143 + csetemp144)*pow(xformL23,2) + (csetemp149 + 
      csetemp150 + csetemp151 + csetemp152)*pow(xformL33,2);
    
    CCTK_REAL dg4332 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp185 + csetemp186 + 
      csetemp187 + csetemp188)*xformL23*xformL33 + xformL03*((csetemp157 + 
      csetemp158 + csetemp159 + csetemp160)*xformL13 + (csetemp161 + 
      csetemp162 + csetemp163 + csetemp164)*xformL23 + (csetemp165 + 
      csetemp166 + csetemp167 + csetemp168)*xformL33) + xformL13*((csetemp173 
      + csetemp174 + csetemp175 + csetemp176)*xformL23 + (csetemp177 + 
      csetemp178 + csetemp179 + csetemp180)*xformL33)) + (csetemp153 + 
      csetemp154 + csetemp155 + csetemp156)*pow(xformL03,2) + (csetemp169 + 
      csetemp170 + csetemp171 + csetemp172)*pow(xformL13,2) + (csetemp181 + 
      csetemp182 + csetemp183 + csetemp184)*pow(xformL23,2) + (csetemp189 + 
      csetemp190 + csetemp191 + csetemp192)*pow(xformL33,2);
    
    CCTK_REAL dg4333 CCTK_ATTRIBUTE_UNUSED = 2*((csetemp225 + csetemp226 + 
      csetemp227 + csetemp228)*xformL23*xformL33 + xformL03*((csetemp197 + 
      csetemp198 + csetemp199 + csetemp200)*xformL13 + (csetemp201 + 
      csetemp202 + csetemp203 + csetemp204)*xformL23 + (csetemp205 + 
      csetemp206 + csetemp207 + csetemp208)*xformL33) + xformL13*((csetemp213 
      + csetemp214 + csetemp215 + csetemp216)*xformL23 + (csetemp217 + 
      csetemp218 + csetemp219 + csetemp220)*xformL33)) + (csetemp193 + 
      csetemp194 + csetemp195 + csetemp196)*pow(xformL03,2) + (csetemp209 + 
      csetemp210 + csetemp211 + csetemp212)*pow(xformL13,2) + (csetemp221 + 
      csetemp222 + csetemp223 + csetemp224)*pow(xformL23,2) + (csetemp229 + 
      csetemp230 + csetemp231 + csetemp232)*pow(xformL33,2);
    
    CCTK_REAL betal1 CCTK_ATTRIBUTE_UNUSED = g401;
    
    CCTK_REAL betal2 CCTK_ATTRIBUTE_UNUSED = g402;
    
    CCTK_REAL betal3 CCTK_ATTRIBUTE_UNUSED = g403;
    
    gExact11L = g411;
    
    gExact12L = g412;
    
    gExact13L = g413;
    
    gExact22L = g422;
    
    gExact23L = g423;
    
    gExact33L = g433;
    
    CCTK_REAL csetemp233 CCTK_ATTRIBUTE_UNUSED = 2*gxyL*gxzL*gyzL;
    
    CCTK_REAL csetemp234 CCTK_ATTRIBUTE_UNUSED = gxxL*gyyL*gzzL;
    
    CCTK_REAL csetemp235 CCTK_ATTRIBUTE_UNUSED = pow(gxzL,2);
    
    CCTK_REAL csetemp236 CCTK_ATTRIBUTE_UNUSED = pow(gyzL,2);
    
    CCTK_REAL csetemp237 CCTK_ATTRIBUTE_UNUSED = pow(gxyL,2);
    
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = csetemp233 + csetemp234 - 
      gyyL*csetemp235 - gxxL*csetemp236 - gzzL*csetemp237;
    
    CCTK_REAL csetemp238 CCTK_ATTRIBUTE_UNUSED = 
      2*gExact12L*gExact13L*gExact23L;
    
    CCTK_REAL csetemp239 CCTK_ATTRIBUTE_UNUSED = 
      gExact11L*gExact22L*gExact33L;
    
    CCTK_REAL csetemp240 CCTK_ATTRIBUTE_UNUSED = pow(detg,-1);
    
    CCTK_REAL csetemp241 CCTK_ATTRIBUTE_UNUSED = pow(gExact23L,2);
    
    CCTK_REAL csetemp242 CCTK_ATTRIBUTE_UNUSED = pow(gExact13L,2);
    
    CCTK_REAL csetemp243 CCTK_ATTRIBUTE_UNUSED = pow(gExact12L,2);
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = (csetemp233 + csetemp234 - 
      gyyL*csetemp235 - gxxL*csetemp236 - 
      gzzL*csetemp237)*csetemp240*(-(gExact22L*gExact33L) + 
      csetemp241)*pow(-csetemp238 - csetemp239 + gExact11L*csetemp241 + 
      gExact22L*csetemp242 + gExact33L*csetemp243,-1);
    
    CCTK_REAL gu12 CCTK_ATTRIBUTE_UNUSED = (gExact13L*gExact23L - 
      gExact12L*gExact33L)*(csetemp233 + csetemp234 - gyyL*csetemp235 - 
      gxxL*csetemp236 - gzzL*csetemp237)*csetemp240*pow(csetemp238 + 
      csetemp239 - gExact11L*csetemp241 - gExact22L*csetemp242 - 
      gExact33L*csetemp243,-1);
    
    CCTK_REAL gu13 CCTK_ATTRIBUTE_UNUSED = (gExact13L*gExact22L - 
      gExact12L*gExact23L)*(csetemp233 + csetemp234 - gyyL*csetemp235 - 
      gxxL*csetemp236 - gzzL*csetemp237)*csetemp240*pow(-csetemp238 - 
      csetemp239 + gExact11L*csetemp241 + gExact22L*csetemp242 + 
      gExact33L*csetemp243,-1);
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = (csetemp233 + csetemp234 - 
      gyyL*csetemp235 - gxxL*csetemp236 - 
      gzzL*csetemp237)*csetemp240*(-(gExact11L*gExact33L) + 
      csetemp242)*pow(-csetemp238 - csetemp239 + gExact11L*csetemp241 + 
      gExact22L*csetemp242 + gExact33L*csetemp243,-1);
    
    CCTK_REAL gu23 CCTK_ATTRIBUTE_UNUSED = (gExact12L*gExact13L - 
      gExact11L*gExact23L)*(csetemp233 + csetemp234 - gyyL*csetemp235 - 
      gxxL*csetemp236 - gzzL*csetemp237)*csetemp240*pow(csetemp238 + 
      csetemp239 - gExact11L*csetemp241 - gExact22L*csetemp242 - 
      gExact33L*csetemp243,-1);
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = (csetemp233 + csetemp234 - 
      gyyL*csetemp235 - gxxL*csetemp236 - 
      gzzL*csetemp237)*csetemp240*(-(gExact11L*gExact22L) + 
      csetemp243)*pow(-csetemp238 - csetemp239 + gExact11L*csetemp241 + 
      gExact22L*csetemp242 + gExact33L*csetemp243,-1);
    
    betaExact1L = betal1*gu11 + betal2*gu12 + betal3*gu13;
    
    betaExact2L = betal1*gu12 + betal2*gu22 + betal3*gu23;
    
    betaExact3L = betal1*gu13 + betal2*gu23 + betal3*gu33;
    
    CCTK_REAL betasq CCTK_ATTRIBUTE_UNUSED = betaExact1L*betal1 + 
      betaExact2L*betal2 + betaExact3L*betal3;
    
    alpExactL = pow(betasq - g400,0.5);
    
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
    
    CCTK_REAL csetemp244 CCTK_ATTRIBUTE_UNUSED = dtg11*gu11;
    
    CCTK_REAL csetemp245 CCTK_ATTRIBUTE_UNUSED = dtg12*gu12;
    
    CCTK_REAL csetemp246 CCTK_ATTRIBUTE_UNUSED = dtg13*gu13;
    
    CCTK_REAL csetemp247 CCTK_ATTRIBUTE_UNUSED = dtg12*gu11;
    
    CCTK_REAL csetemp248 CCTK_ATTRIBUTE_UNUSED = dtg22*gu12;
    
    CCTK_REAL csetemp249 CCTK_ATTRIBUTE_UNUSED = dtg23*gu13;
    
    CCTK_REAL csetemp250 CCTK_ATTRIBUTE_UNUSED = dtg13*gu11;
    
    CCTK_REAL csetemp251 CCTK_ATTRIBUTE_UNUSED = dtg23*gu12;
    
    CCTK_REAL csetemp252 CCTK_ATTRIBUTE_UNUSED = dtg33*gu13;
    
    CCTK_REAL dtgu11 CCTK_ATTRIBUTE_UNUSED = -((csetemp244 + csetemp245 + 
      csetemp246)*gu11) - (csetemp247 + csetemp248 + csetemp249)*gu12 - 
      (csetemp250 + csetemp251 + csetemp252)*gu13;
    
    CCTK_REAL dtgu12 CCTK_ATTRIBUTE_UNUSED = -((csetemp244 + csetemp245 + 
      csetemp246)*gu12) - (csetemp247 + csetemp248 + csetemp249)*gu22 - 
      (csetemp250 + csetemp251 + csetemp252)*gu23;
    
    CCTK_REAL dtgu13 CCTK_ATTRIBUTE_UNUSED = -((csetemp244 + csetemp245 + 
      csetemp246)*gu13) - (csetemp247 + csetemp248 + csetemp249)*gu23 - 
      (csetemp250 + csetemp251 + csetemp252)*gu33;
    
    CCTK_REAL csetemp253 CCTK_ATTRIBUTE_UNUSED = dtg11*gu12;
    
    CCTK_REAL csetemp254 CCTK_ATTRIBUTE_UNUSED = dtg12*gu22;
    
    CCTK_REAL csetemp255 CCTK_ATTRIBUTE_UNUSED = dtg13*gu23;
    
    CCTK_REAL csetemp256 CCTK_ATTRIBUTE_UNUSED = dtg22*gu22;
    
    CCTK_REAL csetemp257 CCTK_ATTRIBUTE_UNUSED = dtg23*gu23;
    
    CCTK_REAL csetemp258 CCTK_ATTRIBUTE_UNUSED = dtg13*gu12;
    
    CCTK_REAL csetemp259 CCTK_ATTRIBUTE_UNUSED = dtg23*gu22;
    
    CCTK_REAL csetemp260 CCTK_ATTRIBUTE_UNUSED = dtg33*gu23;
    
    CCTK_REAL dtgu22 CCTK_ATTRIBUTE_UNUSED = -((csetemp253 + csetemp254 + 
      csetemp255)*gu12) - (csetemp245 + csetemp256 + csetemp257)*gu22 - 
      (csetemp258 + csetemp259 + csetemp260)*gu23;
    
    CCTK_REAL dtgu23 CCTK_ATTRIBUTE_UNUSED = -((csetemp253 + csetemp254 + 
      csetemp255)*gu13) - (csetemp245 + csetemp256 + csetemp257)*gu23 - 
      (csetemp258 + csetemp259 + csetemp260)*gu33;
    
    CCTK_REAL dtgu33 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dtg11*gu13 + 
      dtg12*gu23 + dtg13*gu33)) - gu23*(dtg12*gu13 + dtg22*gu23 + dtg23*gu33) 
      - gu33*(csetemp246 + csetemp257 + dtg33*gu33);
    
    CCTK_REAL csetemp261 CCTK_ATTRIBUTE_UNUSED = dg111*gu11;
    
    CCTK_REAL csetemp262 CCTK_ATTRIBUTE_UNUSED = dg121*gu12;
    
    CCTK_REAL csetemp263 CCTK_ATTRIBUTE_UNUSED = dg131*gu13;
    
    CCTK_REAL csetemp264 CCTK_ATTRIBUTE_UNUSED = dg121*gu11;
    
    CCTK_REAL csetemp265 CCTK_ATTRIBUTE_UNUSED = dg221*gu12;
    
    CCTK_REAL csetemp266 CCTK_ATTRIBUTE_UNUSED = dg231*gu13;
    
    CCTK_REAL csetemp267 CCTK_ATTRIBUTE_UNUSED = dg131*gu11;
    
    CCTK_REAL csetemp268 CCTK_ATTRIBUTE_UNUSED = dg231*gu12;
    
    CCTK_REAL csetemp269 CCTK_ATTRIBUTE_UNUSED = dg331*gu13;
    
    CCTK_REAL dgu111 CCTK_ATTRIBUTE_UNUSED = -((csetemp261 + csetemp262 + 
      csetemp263)*gu11) - (csetemp264 + csetemp265 + csetemp266)*gu12 - 
      (csetemp267 + csetemp268 + csetemp269)*gu13;
    
    CCTK_REAL dgu121 CCTK_ATTRIBUTE_UNUSED = -((csetemp261 + csetemp262 + 
      csetemp263)*gu12) - (csetemp264 + csetemp265 + csetemp266)*gu22 - 
      (csetemp267 + csetemp268 + csetemp269)*gu23;
    
    CCTK_REAL dgu131 CCTK_ATTRIBUTE_UNUSED = -((csetemp261 + csetemp262 + 
      csetemp263)*gu13) - (csetemp264 + csetemp265 + csetemp266)*gu23 - 
      (csetemp267 + csetemp268 + csetemp269)*gu33;
    
    CCTK_REAL csetemp270 CCTK_ATTRIBUTE_UNUSED = dg111*gu12;
    
    CCTK_REAL csetemp271 CCTK_ATTRIBUTE_UNUSED = dg121*gu22;
    
    CCTK_REAL csetemp272 CCTK_ATTRIBUTE_UNUSED = dg131*gu23;
    
    CCTK_REAL csetemp273 CCTK_ATTRIBUTE_UNUSED = dg221*gu22;
    
    CCTK_REAL csetemp274 CCTK_ATTRIBUTE_UNUSED = dg231*gu23;
    
    CCTK_REAL csetemp275 CCTK_ATTRIBUTE_UNUSED = dg131*gu12;
    
    CCTK_REAL csetemp276 CCTK_ATTRIBUTE_UNUSED = dg231*gu22;
    
    CCTK_REAL csetemp277 CCTK_ATTRIBUTE_UNUSED = dg331*gu23;
    
    CCTK_REAL dgu221 CCTK_ATTRIBUTE_UNUSED = -((csetemp270 + csetemp271 + 
      csetemp272)*gu12) - (csetemp262 + csetemp273 + csetemp274)*gu22 - 
      (csetemp275 + csetemp276 + csetemp277)*gu23;
    
    CCTK_REAL dgu231 CCTK_ATTRIBUTE_UNUSED = -((csetemp270 + csetemp271 + 
      csetemp272)*gu13) - (csetemp262 + csetemp273 + csetemp274)*gu23 - 
      (csetemp275 + csetemp276 + csetemp277)*gu33;
    
    CCTK_REAL dgu331 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dg111*gu13 + 
      dg121*gu23 + dg131*gu33)) - gu23*(dg121*gu13 + dg221*gu23 + dg231*gu33) 
      - gu33*(csetemp263 + csetemp274 + dg331*gu33);
    
    CCTK_REAL csetemp278 CCTK_ATTRIBUTE_UNUSED = dg112*gu11;
    
    CCTK_REAL csetemp279 CCTK_ATTRIBUTE_UNUSED = dg122*gu12;
    
    CCTK_REAL csetemp280 CCTK_ATTRIBUTE_UNUSED = dg132*gu13;
    
    CCTK_REAL csetemp281 CCTK_ATTRIBUTE_UNUSED = dg122*gu11;
    
    CCTK_REAL csetemp282 CCTK_ATTRIBUTE_UNUSED = dg222*gu12;
    
    CCTK_REAL csetemp283 CCTK_ATTRIBUTE_UNUSED = dg232*gu13;
    
    CCTK_REAL csetemp284 CCTK_ATTRIBUTE_UNUSED = dg132*gu11;
    
    CCTK_REAL csetemp285 CCTK_ATTRIBUTE_UNUSED = dg232*gu12;
    
    CCTK_REAL csetemp286 CCTK_ATTRIBUTE_UNUSED = dg332*gu13;
    
    CCTK_REAL dgu112 CCTK_ATTRIBUTE_UNUSED = -((csetemp278 + csetemp279 + 
      csetemp280)*gu11) - (csetemp281 + csetemp282 + csetemp283)*gu12 - 
      (csetemp284 + csetemp285 + csetemp286)*gu13;
    
    CCTK_REAL dgu122 CCTK_ATTRIBUTE_UNUSED = -((csetemp278 + csetemp279 + 
      csetemp280)*gu12) - (csetemp281 + csetemp282 + csetemp283)*gu22 - 
      (csetemp284 + csetemp285 + csetemp286)*gu23;
    
    CCTK_REAL dgu132 CCTK_ATTRIBUTE_UNUSED = -((csetemp278 + csetemp279 + 
      csetemp280)*gu13) - (csetemp281 + csetemp282 + csetemp283)*gu23 - 
      (csetemp284 + csetemp285 + csetemp286)*gu33;
    
    CCTK_REAL csetemp287 CCTK_ATTRIBUTE_UNUSED = dg112*gu12;
    
    CCTK_REAL csetemp288 CCTK_ATTRIBUTE_UNUSED = dg122*gu22;
    
    CCTK_REAL csetemp289 CCTK_ATTRIBUTE_UNUSED = dg132*gu23;
    
    CCTK_REAL csetemp290 CCTK_ATTRIBUTE_UNUSED = dg222*gu22;
    
    CCTK_REAL csetemp291 CCTK_ATTRIBUTE_UNUSED = dg232*gu23;
    
    CCTK_REAL csetemp292 CCTK_ATTRIBUTE_UNUSED = dg132*gu12;
    
    CCTK_REAL csetemp293 CCTK_ATTRIBUTE_UNUSED = dg232*gu22;
    
    CCTK_REAL csetemp294 CCTK_ATTRIBUTE_UNUSED = dg332*gu23;
    
    CCTK_REAL dgu222 CCTK_ATTRIBUTE_UNUSED = -((csetemp287 + csetemp288 + 
      csetemp289)*gu12) - (csetemp279 + csetemp290 + csetemp291)*gu22 - 
      (csetemp292 + csetemp293 + csetemp294)*gu23;
    
    CCTK_REAL dgu232 CCTK_ATTRIBUTE_UNUSED = -((csetemp287 + csetemp288 + 
      csetemp289)*gu13) - (csetemp279 + csetemp290 + csetemp291)*gu23 - 
      (csetemp292 + csetemp293 + csetemp294)*gu33;
    
    CCTK_REAL dgu332 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dg112*gu13 + 
      dg122*gu23 + dg132*gu33)) - gu23*(dg122*gu13 + dg222*gu23 + dg232*gu33) 
      - gu33*(csetemp280 + csetemp291 + dg332*gu33);
    
    CCTK_REAL csetemp295 CCTK_ATTRIBUTE_UNUSED = dg113*gu11;
    
    CCTK_REAL csetemp296 CCTK_ATTRIBUTE_UNUSED = dg123*gu12;
    
    CCTK_REAL csetemp297 CCTK_ATTRIBUTE_UNUSED = dg133*gu13;
    
    CCTK_REAL csetemp298 CCTK_ATTRIBUTE_UNUSED = dg123*gu11;
    
    CCTK_REAL csetemp299 CCTK_ATTRIBUTE_UNUSED = dg223*gu12;
    
    CCTK_REAL csetemp300 CCTK_ATTRIBUTE_UNUSED = dg233*gu13;
    
    CCTK_REAL csetemp301 CCTK_ATTRIBUTE_UNUSED = dg133*gu11;
    
    CCTK_REAL csetemp302 CCTK_ATTRIBUTE_UNUSED = dg233*gu12;
    
    CCTK_REAL csetemp303 CCTK_ATTRIBUTE_UNUSED = dg333*gu13;
    
    CCTK_REAL dgu113 CCTK_ATTRIBUTE_UNUSED = -((csetemp295 + csetemp296 + 
      csetemp297)*gu11) - (csetemp298 + csetemp299 + csetemp300)*gu12 - 
      (csetemp301 + csetemp302 + csetemp303)*gu13;
    
    CCTK_REAL dgu123 CCTK_ATTRIBUTE_UNUSED = -((csetemp295 + csetemp296 + 
      csetemp297)*gu12) - (csetemp298 + csetemp299 + csetemp300)*gu22 - 
      (csetemp301 + csetemp302 + csetemp303)*gu23;
    
    CCTK_REAL dgu133 CCTK_ATTRIBUTE_UNUSED = -((csetemp295 + csetemp296 + 
      csetemp297)*gu13) - (csetemp298 + csetemp299 + csetemp300)*gu23 - 
      (csetemp301 + csetemp302 + csetemp303)*gu33;
    
    CCTK_REAL csetemp304 CCTK_ATTRIBUTE_UNUSED = dg113*gu12;
    
    CCTK_REAL csetemp305 CCTK_ATTRIBUTE_UNUSED = dg123*gu22;
    
    CCTK_REAL csetemp306 CCTK_ATTRIBUTE_UNUSED = dg133*gu23;
    
    CCTK_REAL csetemp307 CCTK_ATTRIBUTE_UNUSED = dg223*gu22;
    
    CCTK_REAL csetemp308 CCTK_ATTRIBUTE_UNUSED = dg233*gu23;
    
    CCTK_REAL csetemp309 CCTK_ATTRIBUTE_UNUSED = dg133*gu12;
    
    CCTK_REAL csetemp310 CCTK_ATTRIBUTE_UNUSED = dg233*gu22;
    
    CCTK_REAL csetemp311 CCTK_ATTRIBUTE_UNUSED = dg333*gu23;
    
    CCTK_REAL dgu223 CCTK_ATTRIBUTE_UNUSED = -((csetemp304 + csetemp305 + 
      csetemp306)*gu12) - (csetemp296 + csetemp307 + csetemp308)*gu22 - 
      (csetemp309 + csetemp310 + csetemp311)*gu23;
    
    CCTK_REAL dgu233 CCTK_ATTRIBUTE_UNUSED = -((csetemp304 + csetemp305 + 
      csetemp306)*gu13) - (csetemp296 + csetemp307 + csetemp308)*gu23 - 
      (csetemp309 + csetemp310 + csetemp311)*gu33;
    
    CCTK_REAL dgu333 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dg113*gu13 + 
      dg123*gu23 + dg133*gu33)) - gu23*(dg123*gu13 + dg223*gu23 + dg233*gu33) 
      - gu33*(csetemp297 + csetemp308 + dg333*gu33);
    
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
    
    dtbetaExact1L = betal1*dtgu11 + betal2*dtgu12 + betal3*dtgu13 + 
      dtbetal1*gu11 + dtbetal2*gu12 + dtbetal3*gu13;
    
    dtbetaExact2L = betal1*dtgu12 + betal2*dtgu22 + betal3*dtgu23 + 
      dtbetal1*gu12 + dtbetal2*gu22 + dtbetal3*gu23;
    
    dtbetaExact3L = betal1*dtgu13 + betal2*dtgu23 + betal3*dtgu33 + 
      dtbetal1*gu13 + dtbetal2*gu23 + dtbetal3*gu33;
    
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
    
    CCTK_REAL dtbetasq CCTK_ATTRIBUTE_UNUSED = dtbetaExact1L*betal1 + 
      dtbetaExact2L*betal2 + dtbetaExact3L*betal3 + betaExact1L*dtbetal1 + 
      betaExact2L*dtbetal2 + betaExact3L*dtbetal3;
    
    CCTK_REAL csetemp312 CCTK_ATTRIBUTE_UNUSED = pow(alpExactL,-1);
    
    CCTK_REAL dtalpExactL CCTK_ATTRIBUTE_UNUSED = 0.5*csetemp312*(-dg4000 
      + dtbetasq);
    
    CCTK_REAL kExact11L CCTK_ATTRIBUTE_UNUSED = 
      0.5*csetemp312*(2*(gExact11L*dbeta11 + gExact12L*dbeta21 + 
      gExact13L*dbeta31) + betaExact1L*dg111 + betaExact2L*dg112 + 
      betaExact3L*dg113 - dtg11);
    
    CCTK_REAL kExact12L CCTK_ATTRIBUTE_UNUSED = 
      0.5*csetemp312*(gExact11L*dbeta12 + gExact22L*dbeta21 + 
      gExact12L*(dbeta11 + dbeta22) + gExact23L*dbeta31 + gExact13L*dbeta32 + 
      betaExact1L*dg121 + betaExact2L*dg122 + betaExact3L*dg123 - dtg12);
    
    CCTK_REAL kExact13L CCTK_ATTRIBUTE_UNUSED = 
      0.5*csetemp312*(gExact11L*dbeta13 + gExact23L*dbeta21 + 
      gExact12L*dbeta23 + gExact33L*dbeta31 + gExact13L*(dbeta11 + dbeta33) + 
      betaExact1L*dg131 + betaExact2L*dg132 + betaExact3L*dg133 - dtg13);
    
    CCTK_REAL kExact22L CCTK_ATTRIBUTE_UNUSED = 
      0.5*csetemp312*(2*(gExact12L*dbeta12 + gExact22L*dbeta22 + 
      gExact23L*dbeta32) + betaExact1L*dg221 + betaExact2L*dg222 + 
      betaExact3L*dg223 - dtg22);
    
    CCTK_REAL kExact23L CCTK_ATTRIBUTE_UNUSED = 
      0.5*csetemp312*(gExact13L*dbeta12 + gExact12L*dbeta13 + 
      gExact22L*dbeta23 + gExact33L*dbeta32 + gExact23L*(dbeta22 + dbeta33) + 
      betaExact1L*dg231 + betaExact2L*dg232 + betaExact3L*dg233 - dtg23);
    
    CCTK_REAL kExact33L CCTK_ATTRIBUTE_UNUSED = 
      0.5*csetemp312*(2*(gExact13L*dbeta13 + gExact23L*dbeta23 + 
      gExact33L*dbeta33) + betaExact1L*dg331 + betaExact2L*dg332 + 
      betaExact3L*dg333 - dtg33);
    /* Copy local copies back to grid functions */
    alpExact[index] = alpExactL;
    betaExact1[index] = betaExact1L;
    betaExact2[index] = betaExact2L;
    betaExact3[index] = betaExact3L;
    dtalpExact[index] = dtalpExactL;
    dtbetaExact1[index] = dtbetaExact1L;
    dtbetaExact2[index] = dtbetaExact2L;
    dtbetaExact3[index] = dtbetaExact3L;
    gExact11[index] = gExact11L;
    gExact12[index] = gExact12L;
    gExact13[index] = gExact13L;
    gExact22[index] = gExact22L;
    gExact23[index] = gExact23L;
    gExact33[index] = gExact33L;
    kExact11[index] = kExact11L;
    kExact12[index] = kExact12L;
    kExact13[index] = kExact13L;
    kExact22[index] = kExact22L;
    kExact23[index] = kExact23L;
    kExact33[index] = kExact33L;
  }
  CCTK_ENDLOOP3(Vaidya2_exact);
}
extern "C" void Vaidya2_exact(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_Vaidya2_exact
  DECLARE_CCTK_ARGUMENTS_CHECKED(Vaidya2_exact);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering Vaidya2_exact_Body");
  }
  if (cctk_iteration % Vaidya2_exact_calc_every != Vaidya2_exact_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "admbase::metric",
    "Vaidya2::curv_exact",
    "Vaidya2::dtlapse_exact",
    "Vaidya2::dtshift_exact",
    "grid::coordinates",
    "Vaidya2::lapse_exact",
    "Vaidya2::metric_exact",
    "Vaidya2::shift_exact"};
  AssertGroupStorage(cctkGH, "Vaidya2_exact", 8, groups);
  
  
  LoopOverEverything(cctkGH, Vaidya2_exact_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving Vaidya2_exact_Body");
  }
}

} // namespace Vaidya2
