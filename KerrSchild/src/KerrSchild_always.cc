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
#include "cctk_Arguments_KerrSchild.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "Differencing.h"
#include "loopcontrol.h"

namespace KerrSchild {


static void KerrSchild_always_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS_KerrSchild_always;
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
  CCTK_LOOP3(KerrSchild_always,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
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
    
    CCTK_REAL xform2L11 CCTK_ATTRIBUTE_UNUSED = (csetemp6 + (-1 + csetemp7 
      + csetemp8)*(1 + pow(1 - csetemp6 - csetemp7 - csetemp8,0.5)))*pow(-1 + 
      csetemp6 + csetemp7 + csetemp8,-1)*pow(1 + pow(1 - csetemp6 - csetemp7 
      - csetemp8,0.5),-1);
    
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
    
    CCTK_REAL xform2L22 CCTK_ATTRIBUTE_UNUSED = (csetemp7 + (-1 + csetemp6 
      + csetemp8)*(1 + pow(1 - csetemp6 - csetemp7 - csetemp8,0.5)))*pow(-1 + 
      csetemp6 + csetemp7 + csetemp8,-1)*pow(1 + pow(1 - csetemp6 - csetemp7 
      - csetemp8,0.5),-1);
    
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
    
    CCTK_REAL txx1 CCTK_ATTRIBUTE_UNUSED = xformL10*xx0 + xformL11*xx1 + 
      xformL12*xx2 + xformL13*xx3;
    
    CCTK_REAL txx2 CCTK_ATTRIBUTE_UNUSED = xformL20*xx0 + xformL21*xx1 + 
      xformL22*xx2 + xformL23*xx3;
    
    CCTK_REAL txx3 CCTK_ATTRIBUTE_UNUSED = xformL30*xx0 + xformL31*xx1 + 
      xformL32*xx2 + xformL33*xx3;
    
    CCTK_REAL X CCTK_ATTRIBUTE_UNUSED = txx1;
    
    CCTK_REAL Y CCTK_ATTRIBUTE_UNUSED = txx2;
    
    CCTK_REAL Z CCTK_ATTRIBUTE_UNUSED = txx3;
    
    CCTK_REAL csetemp10 CCTK_ATTRIBUTE_UNUSED = pow(a,2);
    
    CCTK_REAL csetemp11 CCTK_ATTRIBUTE_UNUSED = pow(X,2);
    
    CCTK_REAL csetemp12 CCTK_ATTRIBUTE_UNUSED = pow(Y,2);
    
    CCTK_REAL csetemp13 CCTK_ATTRIBUTE_UNUSED = pow(Z,2);
    
    CCTK_REAL rXYZ CCTK_ATTRIBUTE_UNUSED = 
      0.707106781186547524400844362105*pow(-csetemp10 + csetemp11 + csetemp12 
      + csetemp13 + pow(4*csetemp10*csetemp13 + pow(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13,2),0.5),0.5);
    
    CCTK_REAL csetemp14 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,3);
    
    CCTK_REAL csetemp15 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,4);
    
    CCTK_REAL tg400 CCTK_ATTRIBUTE_UNUSED = -1 + 
      2*csetemp14*M*pow(csetemp10*csetemp13 + csetemp15,-1);
    
    CCTK_REAL csetemp16 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,2);
    
    CCTK_REAL csetemp17 CCTK_ATTRIBUTE_UNUSED = rXYZ*X;
    
    CCTK_REAL csetemp18 CCTK_ATTRIBUTE_UNUSED = a*Y;
    
    CCTK_REAL tg401 CCTK_ATTRIBUTE_UNUSED = 2*csetemp14*(csetemp17 + 
      csetemp18)*M*pow(csetemp10*csetemp13 + csetemp15,-1)*pow(csetemp10 + 
      csetemp16,-1);
    
    CCTK_REAL csetemp19 CCTK_ATTRIBUTE_UNUSED = a*X;
    
    CCTK_REAL csetemp20 CCTK_ATTRIBUTE_UNUSED = -csetemp19;
    
    CCTK_REAL csetemp21 CCTK_ATTRIBUTE_UNUSED = rXYZ*Y;
    
    CCTK_REAL tg402 CCTK_ATTRIBUTE_UNUSED = 2*csetemp14*(csetemp20 + 
      csetemp21)*M*pow(csetemp10*csetemp13 + csetemp15,-1)*pow(csetemp10 + 
      csetemp16,-1);
    
    CCTK_REAL tg403 CCTK_ATTRIBUTE_UNUSED = 
      2*csetemp16*M*Z*pow(csetemp10*csetemp13 + csetemp15,-1);
    
    CCTK_REAL tg411 CCTK_ATTRIBUTE_UNUSED = 1 + 
      2*csetemp14*M*pow(csetemp10*csetemp13 + csetemp15,-1)*pow(csetemp10 + 
      csetemp16,-2)*pow(csetemp17 + csetemp18,2);
    
    CCTK_REAL tg412 CCTK_ATTRIBUTE_UNUSED = 2*csetemp14*(csetemp17 + 
      csetemp18)*(csetemp20 + csetemp21)*M*pow(csetemp10*csetemp13 + 
      csetemp15,-1)*pow(csetemp10 + csetemp16,-2);
    
    CCTK_REAL tg413 CCTK_ATTRIBUTE_UNUSED = 2*csetemp16*(csetemp17 + 
      csetemp18)*M*Z*pow(csetemp10*csetemp13 + csetemp15,-1)*pow(csetemp10 + 
      csetemp16,-1);
    
    CCTK_REAL tg422 CCTK_ATTRIBUTE_UNUSED = 1 + 
      2*csetemp14*M*pow(csetemp10*csetemp13 + csetemp15,-1)*pow(csetemp10 + 
      csetemp16,-2)*pow(csetemp20 + csetemp21,2);
    
    CCTK_REAL tg423 CCTK_ATTRIBUTE_UNUSED = 2*csetemp16*(csetemp20 + 
      csetemp21)*M*Z*pow(csetemp10*csetemp13 + csetemp15,-1)*pow(csetemp10 + 
      csetemp16,-1);
    
    CCTK_REAL tg433 CCTK_ATTRIBUTE_UNUSED = 1 + 
      2*csetemp13*M*rXYZ*pow(csetemp10*csetemp13 + csetemp15,-1);
    
    CCTK_REAL tdg4000 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL csetemp22 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,7);
    
    CCTK_REAL tdg4001 CCTK_ATTRIBUTE_UNUSED = 
      2*(3*csetemp10*csetemp13*csetemp14 - 
      csetemp22)*M*X*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(4*csetemp10*csetemp13 + pow(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4002 CCTK_ATTRIBUTE_UNUSED = 
      2*(3*csetemp10*csetemp13*csetemp14 - 
      csetemp22)*M*Y*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(4*csetemp10*csetemp13 + pow(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL csetemp23 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,6);
    
    CCTK_REAL csetemp24 CCTK_ATTRIBUTE_UNUSED = pow(a,4);
    
    CCTK_REAL tdg4003 CCTK_ATTRIBUTE_UNUSED = 2*(csetemp10*(-5*csetemp15 + 
      (-2*csetemp10 + 2*(csetemp11 + csetemp12) + 5*csetemp13)*csetemp16) - 
      csetemp23 + 3*csetemp13*csetemp24)*M*rXYZ*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(4*csetemp10*csetemp13 + pow(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4010 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL csetemp25 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,9);
    
    CCTK_REAL csetemp26 CCTK_ATTRIBUTE_UNUSED = pow(a,3);
    
    CCTK_REAL csetemp27 CCTK_ATTRIBUTE_UNUSED = pow(a,5);
    
    CCTK_REAL csetemp28 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,5);
    
    CCTK_REAL tdg4011 CCTK_ATTRIBUTE_UNUSED = 2*csetemp14*M*((3*csetemp10 
      - 3*csetemp11 - csetemp12 - csetemp13)*csetemp22 + 2*csetemp25 + 
      csetemp10*((3*csetemp10 + csetemp11 - csetemp12 - 
      csetemp13)*csetemp13*csetemp14 + (csetemp10 - csetemp11 - csetemp12 + 
      csetemp13)*csetemp28) + (-3*a*csetemp23 - csetemp15*csetemp26 + 
      3*csetemp13*csetemp27)*X*Y + csetemp13*((csetemp10 + 3*csetemp11 - 
      csetemp12 - csetemp13)*csetemp24*rXYZ + 
      csetemp16*csetemp26*X*Y))*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-2)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL csetemp29 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,8);
    
    CCTK_REAL tdg4012 CCTK_ATTRIBUTE_UNUSED = 2*csetemp14*M*(((csetemp10 - 
      csetemp11 - 2*csetemp12 + csetemp13)*csetemp15 + (3*csetemp10 - 
      csetemp11 - csetemp13)*csetemp13*csetemp16)*csetemp26 + a*((3*csetemp10 
      - csetemp11 - 4*csetemp12 - csetemp13)*csetemp23 + 2*csetemp29) + 
      (-2*csetemp22 + 4*csetemp13*csetemp24*rXYZ)*X*Y + csetemp13*((csetemp10 
      - csetemp11 + 2*csetemp12 - csetemp13)*csetemp27 + 
      2*csetemp10*csetemp14*X*Y))*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-2)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4013 CCTK_ATTRIBUTE_UNUSED = 2*M*rXYZ*((-2*csetemp22 + 
      csetemp10*(2*(-csetemp10 + csetemp11 + csetemp12 + 
      2*csetemp13)*csetemp14 - 4*csetemp28) + 4*csetemp13*csetemp24*rXYZ)*X + 
      (-3*a*csetemp23 + (-5*csetemp15 + (-2*csetemp10 + 2*(csetemp11 + 
      csetemp12) + 3*csetemp13)*csetemp16)*csetemp26 + 
      3*csetemp13*csetemp27)*Y)*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-1)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4020 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4021 CCTK_ATTRIBUTE_UNUSED = -2*csetemp14*M*(((csetemp10 
      - 2*csetemp11 - csetemp12 + csetemp13)*csetemp15 + (3*csetemp10 - 
      csetemp12 - csetemp13)*csetemp13*csetemp16)*csetemp26 + a*((3*csetemp10 
      - 4*csetemp11 - csetemp12 - csetemp13)*csetemp23 + 2*csetemp29) + 
      (2*csetemp22 - 4*csetemp13*csetemp24*rXYZ)*X*Y + csetemp13*((csetemp10 
      + 2*csetemp11 - csetemp12 - csetemp13)*csetemp27 - 
      2*csetemp10*csetemp14*X*Y))*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-2)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4022 CCTK_ATTRIBUTE_UNUSED = 2*csetemp14*M*((3*csetemp10 
      - csetemp11 - 3*csetemp12 - csetemp13)*csetemp22 + 2*csetemp25 + 
      csetemp10*((3*csetemp10 - csetemp11 + csetemp12 - 
      csetemp13)*csetemp13*csetemp14 + (csetemp10 - csetemp11 - csetemp12 + 
      csetemp13)*csetemp28) + (3*a*csetemp23 + csetemp15*csetemp26 - 
      3*csetemp13*csetemp27)*X*Y + csetemp13*((csetemp10 - csetemp11 + 
      3*csetemp12 - csetemp13)*csetemp24*rXYZ - 
      csetemp16*csetemp26*X*Y))*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-2)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4023 CCTK_ATTRIBUTE_UNUSED = 2*M*rXYZ*((3*a*csetemp23 + 
      (5*csetemp15 + (2*csetemp10 - 2*(csetemp11 + csetemp12) - 
      3*csetemp13)*csetemp16)*csetemp26)*X + (-2*csetemp22 + 
      csetemp10*(2*(-csetemp10 + csetemp11 + csetemp12 + 
      2*csetemp13)*csetemp14 - 4*csetemp28))*Y + csetemp13*(-3*csetemp27*X + 
      4*csetemp24*rXYZ*Y))*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-1)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4030 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4031 CCTK_ATTRIBUTE_UNUSED = 
      4*(csetemp10*csetemp13*csetemp16 - 
      csetemp23)*M*X*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(4*csetemp10*csetemp13 + pow(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4032 CCTK_ATTRIBUTE_UNUSED = 
      4*(csetemp10*csetemp13*csetemp16 - 
      csetemp23)*M*Y*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(4*csetemp10*csetemp13 + pow(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4033 CCTK_ATTRIBUTE_UNUSED = 2*(csetemp10*csetemp13 - 
      csetemp15)*(2*csetemp10*csetemp13 - 2*csetemp15 + (-csetemp10 + 
      csetemp11 + csetemp12 + 
      3*csetemp13)*csetemp16)*M*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(4*csetemp10*csetemp13 + pow(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4110 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4111 CCTK_ATTRIBUTE_UNUSED = 2*csetemp14*(csetemp17 + 
      csetemp18)*M*((6*csetemp10 - 5*csetemp11 - 2*(csetemp12 + 
      csetemp13))*csetemp22 + 4*csetemp25 + csetemp10*(csetemp13*(6*csetemp10 
      - csetemp11 - 2*(csetemp12 + csetemp13))*csetemp14 + (-csetemp11 - 
      2*csetemp12 + 2*(csetemp10 + csetemp13))*csetemp28) + (-5*a*csetemp23 - 
      csetemp15*csetemp26 + 3*csetemp13*csetemp27)*X*Y + 
      csetemp13*((2*csetemp10 + 3*csetemp11 - 2*(csetemp12 + 
      csetemp13))*csetemp24*rXYZ - 
      csetemp16*csetemp26*X*Y))*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-3)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4112 CCTK_ATTRIBUTE_UNUSED = 2*csetemp14*(csetemp17 + 
      csetemp18)*M*(((-2*csetemp11 - 3*csetemp12 + 2*(csetemp10 + 
      csetemp13))*csetemp15 + csetemp13*(6*csetemp10 - 3*csetemp12 - 
      2*(csetemp11 + csetemp13))*csetemp16)*csetemp26 + a*((6*csetemp10 - 
      7*csetemp12 - 2*(csetemp11 + csetemp13))*csetemp23 + 4*csetemp29) + 
      (-3*csetemp22 + csetemp10*csetemp28 + 5*csetemp13*csetemp24*rXYZ)*X*Y + 
      csetemp13*((2*csetemp10 + csetemp12 - 2*(csetemp11 + 
      csetemp13))*csetemp27 + 
      csetemp10*csetemp14*X*Y))*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-3)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4113 CCTK_ATTRIBUTE_UNUSED = 2*(csetemp17 + 
      csetemp18)*M*rXYZ*((-3*csetemp22 + csetemp10*((-2*csetemp10 + 
      2*(csetemp11 + csetemp12) + 3*csetemp13)*csetemp14 - 3*csetemp28) + 
      5*csetemp13*csetemp24*rXYZ)*X + (-5*a*csetemp23 + (-5*csetemp15 + 
      (-2*csetemp10 + 2*(csetemp11 + csetemp12) + 
      csetemp13)*csetemp16)*csetemp26 + 
      3*csetemp13*csetemp27)*Y)*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-2)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4120 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL csetemp30 CCTK_ATTRIBUTE_UNUSED = pow(rXYZ,10);
    
    CCTK_REAL csetemp31 CCTK_ATTRIBUTE_UNUSED = pow(a,6);
    
    CCTK_REAL tdg4121 CCTK_ATTRIBUTE_UNUSED = 
      2*csetemp14*M*((a*(2*(-3*csetemp10 + 3*csetemp11 - csetemp12 + 
      csetemp13)*csetemp22 - 4*csetemp25) + 
      csetemp26*(2*csetemp13*(-3*csetemp10 + csetemp11 + csetemp12 + 
      csetemp13)*csetemp14 - 2*(csetemp10 - csetemp11 - csetemp12 + 
      csetemp13)*csetemp28))*X + (csetemp10*(2*(-csetemp10 + 3*csetemp11 + 
      csetemp13)*csetemp23 - csetemp15*(csetemp13*(csetemp12 + csetemp13) - 
      csetemp10*(csetemp12 + 2*(csetemp11 + csetemp13)) + csetemp24)) - 
      (-csetemp10 + 4*csetemp11 + csetemp12 + csetemp13)*csetemp29 + 
      2*csetemp30 + csetemp13*(-csetemp10 - 2*csetemp11 + csetemp12 + 
      csetemp13)*csetemp31)*Y - 2*csetemp13*((csetemp10 + csetemp11 - 
      3*csetemp12 - csetemp13)*csetemp27*rXYZ*X + (csetemp10 - 
      3*csetemp11)*csetemp16*csetemp24*Y))*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-3)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4122 CCTK_ATTRIBUTE_UNUSED = 
      -2*csetemp14*M*((csetemp10*csetemp15*(csetemp13*(csetemp11 + csetemp13) 
      - csetemp10*(csetemp11 + 2*(csetemp12 + csetemp13)) + csetemp24) + 
      (-csetemp10 + csetemp11 + 4*csetemp12 + csetemp13)*csetemp29 - 
      2*csetemp30 + (csetemp10 - csetemp11 + 2*csetemp12 - 
      csetemp13)*csetemp13*csetemp31)*X + (-4*a*csetemp25 - 2*((csetemp10 - 
      csetemp11 - csetemp12 + csetemp13)*csetemp26*csetemp28 + (csetemp10 - 
      3*csetemp11 + csetemp12 - csetemp13)*csetemp13*csetemp27*rXYZ))*Y + 
      2*((csetemp10*(csetemp10 - 3*csetemp12 - csetemp13)*csetemp23 + 
      (csetemp10 - 3*csetemp12)*csetemp13*csetemp16*csetemp24)*X + 
      (a*(-3*csetemp10 - csetemp11 + 3*csetemp12 + csetemp13)*csetemp22 + 
      csetemp13*(-3*csetemp10 + csetemp11 + csetemp12 + 
      csetemp13)*csetemp14*csetemp26)*Y))*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-3)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4123 CCTK_ATTRIBUTE_UNUSED = 2*M*rXYZ*((csetemp11 - 
      csetemp12)*(4*a*csetemp22 + csetemp26*(-2*(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13)*csetemp14 + 4*csetemp28)) + 
      (csetemp10*((2*(csetemp11 + csetemp12) + 3*(csetemp10 + 
      csetemp13))*csetemp15 + 2*csetemp23) + 2*(csetemp10 - csetemp11 - 
      csetemp12 + 2*csetemp13)*csetemp16*csetemp24 - 3*csetemp29)*X*Y + 
      csetemp13*(4*(-csetemp11 + csetemp12)*csetemp27*rXYZ - 
      3*csetemp31*X*Y))*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-2)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4130 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4131 CCTK_ATTRIBUTE_UNUSED = 
      -2*csetemp16*M*((-3*csetemp10 + 4*csetemp11 + csetemp12 + 
      csetemp13)*csetemp22 + csetemp10*(csetemp13*(-3*csetemp10 + csetemp12 + 
      csetemp13)*csetemp14 - (csetemp10 - 2*csetemp11 - csetemp12 + 
      csetemp13)*csetemp28) + csetemp13*(-csetemp10 - 2*csetemp11 + csetemp12 
      + csetemp13)*csetemp24*rXYZ + (4*a*csetemp23 + 
      2*csetemp15*csetemp26)*X*Y - 2*(csetemp25 + 
      csetemp13*csetemp27*X*Y))*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-2)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4132 CCTK_ATTRIBUTE_UNUSED = 
      -2*csetemp16*M*((-((csetemp10 - csetemp11 - 3*csetemp12 + 
      csetemp13)*csetemp15) + csetemp13*(-3*csetemp10 + csetemp11 + csetemp12 
      + csetemp13)*csetemp16)*csetemp26 + a*((-3*csetemp10 + csetemp11 + 
      5*csetemp12 + csetemp13)*csetemp23 - 2*csetemp29) + (3*csetemp22 + 
      csetemp10*csetemp28 - 3*csetemp13*csetemp24*rXYZ)*X*Y + 
      csetemp13*((-csetemp10 + csetemp11 - csetemp12 + csetemp13)*csetemp27 - 
      csetemp10*csetemp14*X*Y))*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-2)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL csetemp32 CCTK_ATTRIBUTE_UNUSED = pow(Z,4);
    
    CCTK_REAL tdg4133 CCTK_ATTRIBUTE_UNUSED = 2*M*((-((-csetemp10 + 
      csetemp11 + csetemp12 + 4*csetemp13)*csetemp22) + 2*csetemp25 + 
      csetemp10*csetemp13*((-csetemp10 + csetemp11 + csetemp12 + 
      2*csetemp13)*csetemp14 - 3*csetemp28) + 3*csetemp24*csetemp32*rXYZ)*X + 
      (a*(csetemp10 - csetemp11 - csetemp12 - 5*csetemp13)*csetemp23 + 
      csetemp13*(-4*csetemp15 + (-csetemp10 + csetemp11 + csetemp12 + 
      csetemp13)*csetemp16)*csetemp26 + 2*(a*csetemp29 + 
      csetemp27*csetemp32))*Y)*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-1)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4220 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL csetemp33 CCTK_ATTRIBUTE_UNUSED = -csetemp21;
    
    CCTK_REAL tdg4221 CCTK_ATTRIBUTE_UNUSED = 2*csetemp14*(csetemp19 + 
      csetemp33)*M*(((-3*csetemp11 - 2*csetemp12 + 2*(csetemp10 + 
      csetemp13))*csetemp15 + csetemp13*(6*csetemp10 - 3*csetemp11 - 
      2*(csetemp12 + csetemp13))*csetemp16)*csetemp26 + a*((6*csetemp10 - 
      7*csetemp11 - 2*(csetemp12 + csetemp13))*csetemp23 + 4*csetemp29) + 
      (3*csetemp22 - csetemp10*csetemp28 - 5*csetemp13*csetemp24*rXYZ)*X*Y + 
      csetemp13*((2*csetemp10 + csetemp11 - 2*(csetemp12 + 
      csetemp13))*csetemp27 - 
      csetemp10*csetemp14*X*Y))*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-3)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4222 CCTK_ATTRIBUTE_UNUSED = 2*csetemp14*(csetemp20 + 
      csetemp21)*M*((6*csetemp10 - 5*csetemp12 - 2*(csetemp11 + 
      csetemp13))*csetemp22 + 4*csetemp25 + csetemp10*(csetemp13*(6*csetemp10 
      - csetemp12 - 2*(csetemp11 + csetemp13))*csetemp14 + (-2*csetemp11 - 
      csetemp12 + 2*(csetemp10 + csetemp13))*csetemp28) + (5*a*csetemp23 + 
      csetemp15*csetemp26 - 3*csetemp13*csetemp27)*X*Y + 
      csetemp13*((2*csetemp10 + 3*csetemp12 - 2*(csetemp11 + 
      csetemp13))*csetemp24*rXYZ + 
      csetemp16*csetemp26*X*Y))*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-3)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4223 CCTK_ATTRIBUTE_UNUSED = 2*(csetemp19 + 
      csetemp33)*M*rXYZ*((-5*csetemp15 + (-2*csetemp10 + 2*(csetemp11 + 
      csetemp12) + csetemp13)*csetemp16)*csetemp26*X + 
      csetemp10*((2*csetemp10 - 2*(csetemp11 + csetemp12) - 
      3*csetemp13)*csetemp14 + 3*csetemp28)*Y + 3*(csetemp13*csetemp27*X + 
      csetemp22*Y) - 5*(a*csetemp23*X + 
      csetemp13*csetemp24*rXYZ*Y))*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-2)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4230 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4231 CCTK_ATTRIBUTE_UNUSED = 2*csetemp16*M*((-((csetemp10 
      - 3*csetemp11 - csetemp12 + csetemp13)*csetemp15) + 
      csetemp13*(-3*csetemp10 + csetemp11 + csetemp12 + 
      csetemp13)*csetemp16)*csetemp26 + a*((-3*csetemp10 + 5*csetemp11 + 
      csetemp12 + csetemp13)*csetemp23 - 2*csetemp29) + (-3*csetemp22 - 
      csetemp10*csetemp28 + 3*csetemp13*csetemp24*rXYZ)*X*Y + 
      csetemp13*((-csetemp10 - csetemp11 + csetemp12 + csetemp13)*csetemp27 + 
      csetemp10*csetemp14*X*Y))*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-2)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4232 CCTK_ATTRIBUTE_UNUSED = 
      -2*csetemp16*M*((-3*csetemp10 + csetemp11 + 4*csetemp12 + 
      csetemp13)*csetemp22 + csetemp10*(csetemp13*(-3*csetemp10 + csetemp11 + 
      csetemp13)*csetemp14 - (csetemp10 - csetemp11 - 2*csetemp12 + 
      csetemp13)*csetemp28) - 4*a*csetemp23*X*Y - 2*(csetemp25 + 
      csetemp15*csetemp26*X*Y) + csetemp13*((-csetemp10 + csetemp11 - 
      2*csetemp12 + csetemp13)*csetemp24*rXYZ + 
      2*csetemp27*X*Y))*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-2)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4233 CCTK_ATTRIBUTE_UNUSED = -2*M*(a*((csetemp10 - 
      csetemp11 - csetemp12 - 5*csetemp13)*csetemp23 + 2*csetemp29)*X + 
      ((-csetemp10 + csetemp11 + csetemp12 + 4*csetemp13)*csetemp22 - 
      2*csetemp25 + 3*csetemp10*csetemp13*csetemp28)*Y + 
      csetemp13*((-4*csetemp15 + (-csetemp10 + csetemp11 + csetemp12 + 
      csetemp13)*csetemp16)*csetemp26*X + csetemp10*(csetemp10 - csetemp11 - 
      csetemp12 - 2*csetemp13)*csetemp14*Y) + csetemp32*(2*csetemp27*X - 
      3*csetemp24*rXYZ*Y))*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(csetemp10 + csetemp16,-1)*pow(4*csetemp10*csetemp13 + 
      pow(-csetemp10 + csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4330 CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL tdg4331 CCTK_ATTRIBUTE_UNUSED = 2*csetemp13*M*(-3*csetemp28 
      + csetemp10*csetemp13*rXYZ)*X*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(4*csetemp10*csetemp13 + pow(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4332 CCTK_ATTRIBUTE_UNUSED = 2*csetemp13*M*(-3*csetemp28 
      + csetemp10*csetemp13*rXYZ)*Y*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(4*csetemp10*csetemp13 + pow(-csetemp10 + csetemp11 + 
      csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL tdg4333 CCTK_ATTRIBUTE_UNUSED = 2*((2*csetemp10 - 
      2*(csetemp11 + csetemp12) - 5*csetemp13)*csetemp23 + 4*csetemp29 + 
      csetemp24*csetemp32 + csetemp10*(-3*csetemp13*csetemp15 + 
      csetemp16*csetemp32))*M*Z*pow(csetemp10*csetemp13 + 
      csetemp15,-2)*pow(rXYZ,-1)*pow(4*csetemp10*csetemp13 + pow(-csetemp10 + 
      csetemp11 + csetemp12 + csetemp13,2),-0.5);
    
    CCTK_REAL g400 CCTK_ATTRIBUTE_UNUSED = 2*(tg423*xformL20*xformL30 + 
      xformL00*(tg401*xformL10 + tg402*xformL20 + tg403*xformL30) + 
      xformL10*(tg412*xformL20 + tg413*xformL30)) + tg400*pow(xformL00,2) + 
      tg411*pow(xformL10,2) + tg422*pow(xformL20,2) + tg433*pow(xformL30,2);
    
    CCTK_REAL csetemp34 CCTK_ATTRIBUTE_UNUSED = tg400*xformL01;
    
    CCTK_REAL csetemp35 CCTK_ATTRIBUTE_UNUSED = tg401*xformL11;
    
    CCTK_REAL csetemp36 CCTK_ATTRIBUTE_UNUSED = tg402*xformL21;
    
    CCTK_REAL csetemp37 CCTK_ATTRIBUTE_UNUSED = tg403*xformL31;
    
    CCTK_REAL csetemp38 CCTK_ATTRIBUTE_UNUSED = tg401*xformL01;
    
    CCTK_REAL csetemp39 CCTK_ATTRIBUTE_UNUSED = tg411*xformL11;
    
    CCTK_REAL csetemp40 CCTK_ATTRIBUTE_UNUSED = tg412*xformL21;
    
    CCTK_REAL csetemp41 CCTK_ATTRIBUTE_UNUSED = tg413*xformL31;
    
    CCTK_REAL csetemp42 CCTK_ATTRIBUTE_UNUSED = tg402*xformL01;
    
    CCTK_REAL csetemp43 CCTK_ATTRIBUTE_UNUSED = tg412*xformL11;
    
    CCTK_REAL csetemp44 CCTK_ATTRIBUTE_UNUSED = tg422*xformL21;
    
    CCTK_REAL csetemp45 CCTK_ATTRIBUTE_UNUSED = tg423*xformL31;
    
    CCTK_REAL csetemp46 CCTK_ATTRIBUTE_UNUSED = tg403*xformL01;
    
    CCTK_REAL csetemp47 CCTK_ATTRIBUTE_UNUSED = tg413*xformL11;
    
    CCTK_REAL csetemp48 CCTK_ATTRIBUTE_UNUSED = tg423*xformL21;
    
    CCTK_REAL csetemp49 CCTK_ATTRIBUTE_UNUSED = tg433*xformL31;
    
    CCTK_REAL g401 CCTK_ATTRIBUTE_UNUSED = (csetemp34 + csetemp35 + 
      csetemp36 + csetemp37)*xformL00 + (csetemp38 + csetemp39 + csetemp40 + 
      csetemp41)*xformL10 + (csetemp42 + csetemp43 + csetemp44 + 
      csetemp45)*xformL20 + (csetemp46 + csetemp47 + csetemp48 + 
      csetemp49)*xformL30;
    
    CCTK_REAL csetemp50 CCTK_ATTRIBUTE_UNUSED = tg400*xformL02;
    
    CCTK_REAL csetemp51 CCTK_ATTRIBUTE_UNUSED = tg401*xformL12;
    
    CCTK_REAL csetemp52 CCTK_ATTRIBUTE_UNUSED = tg402*xformL22;
    
    CCTK_REAL csetemp53 CCTK_ATTRIBUTE_UNUSED = tg403*xformL32;
    
    CCTK_REAL csetemp54 CCTK_ATTRIBUTE_UNUSED = tg401*xformL02;
    
    CCTK_REAL csetemp55 CCTK_ATTRIBUTE_UNUSED = tg411*xformL12;
    
    CCTK_REAL csetemp56 CCTK_ATTRIBUTE_UNUSED = tg412*xformL22;
    
    CCTK_REAL csetemp57 CCTK_ATTRIBUTE_UNUSED = tg413*xformL32;
    
    CCTK_REAL csetemp58 CCTK_ATTRIBUTE_UNUSED = tg402*xformL02;
    
    CCTK_REAL csetemp59 CCTK_ATTRIBUTE_UNUSED = tg412*xformL12;
    
    CCTK_REAL csetemp60 CCTK_ATTRIBUTE_UNUSED = tg422*xformL22;
    
    CCTK_REAL csetemp61 CCTK_ATTRIBUTE_UNUSED = tg423*xformL32;
    
    CCTK_REAL csetemp62 CCTK_ATTRIBUTE_UNUSED = tg403*xformL02;
    
    CCTK_REAL csetemp63 CCTK_ATTRIBUTE_UNUSED = tg413*xformL12;
    
    CCTK_REAL csetemp64 CCTK_ATTRIBUTE_UNUSED = tg423*xformL22;
    
    CCTK_REAL csetemp65 CCTK_ATTRIBUTE_UNUSED = tg433*xformL32;
    
    CCTK_REAL g402 CCTK_ATTRIBUTE_UNUSED = (csetemp50 + csetemp51 + 
      csetemp52 + csetemp53)*xformL00 + (csetemp54 + csetemp55 + csetemp56 + 
      csetemp57)*xformL10 + (csetemp58 + csetemp59 + csetemp60 + 
      csetemp61)*xformL20 + (csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL30;
    
    CCTK_REAL csetemp66 CCTK_ATTRIBUTE_UNUSED = tg400*xformL03;
    
    CCTK_REAL csetemp67 CCTK_ATTRIBUTE_UNUSED = tg401*xformL13;
    
    CCTK_REAL csetemp68 CCTK_ATTRIBUTE_UNUSED = tg402*xformL23;
    
    CCTK_REAL csetemp69 CCTK_ATTRIBUTE_UNUSED = tg403*xformL33;
    
    CCTK_REAL csetemp70 CCTK_ATTRIBUTE_UNUSED = tg401*xformL03;
    
    CCTK_REAL csetemp71 CCTK_ATTRIBUTE_UNUSED = tg411*xformL13;
    
    CCTK_REAL csetemp72 CCTK_ATTRIBUTE_UNUSED = tg412*xformL23;
    
    CCTK_REAL csetemp73 CCTK_ATTRIBUTE_UNUSED = tg413*xformL33;
    
    CCTK_REAL csetemp74 CCTK_ATTRIBUTE_UNUSED = tg402*xformL03;
    
    CCTK_REAL csetemp75 CCTK_ATTRIBUTE_UNUSED = tg412*xformL13;
    
    CCTK_REAL csetemp76 CCTK_ATTRIBUTE_UNUSED = tg422*xformL23;
    
    CCTK_REAL csetemp77 CCTK_ATTRIBUTE_UNUSED = tg423*xformL33;
    
    CCTK_REAL csetemp78 CCTK_ATTRIBUTE_UNUSED = tg403*xformL03;
    
    CCTK_REAL csetemp79 CCTK_ATTRIBUTE_UNUSED = tg413*xformL13;
    
    CCTK_REAL csetemp80 CCTK_ATTRIBUTE_UNUSED = tg423*xformL23;
    
    CCTK_REAL csetemp81 CCTK_ATTRIBUTE_UNUSED = tg433*xformL33;
    
    CCTK_REAL g403 CCTK_ATTRIBUTE_UNUSED = (csetemp66 + csetemp67 + 
      csetemp68 + csetemp69)*xformL00 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL10 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL20 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL30;
    
    CCTK_REAL g411 CCTK_ATTRIBUTE_UNUSED = (csetemp34 + csetemp35 + 
      csetemp36 + csetemp37)*xformL01 + (csetemp38 + csetemp39 + csetemp40 + 
      csetemp41)*xformL11 + (csetemp42 + csetemp43 + csetemp44 + 
      csetemp45)*xformL21 + (csetemp46 + csetemp47 + csetemp48 + 
      csetemp49)*xformL31;
    
    CCTK_REAL g412 CCTK_ATTRIBUTE_UNUSED = (csetemp50 + csetemp51 + 
      csetemp52 + csetemp53)*xformL01 + (csetemp54 + csetemp55 + csetemp56 + 
      csetemp57)*xformL11 + (csetemp58 + csetemp59 + csetemp60 + 
      csetemp61)*xformL21 + (csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL31;
    
    CCTK_REAL g413 CCTK_ATTRIBUTE_UNUSED = (csetemp66 + csetemp67 + 
      csetemp68 + csetemp69)*xformL01 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL11 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL21 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL31;
    
    CCTK_REAL g422 CCTK_ATTRIBUTE_UNUSED = (csetemp50 + csetemp51 + 
      csetemp52 + csetemp53)*xformL02 + (csetemp54 + csetemp55 + csetemp56 + 
      csetemp57)*xformL12 + (csetemp58 + csetemp59 + csetemp60 + 
      csetemp61)*xformL22 + (csetemp62 + csetemp63 + csetemp64 + 
      csetemp65)*xformL32;
    
    CCTK_REAL g423 CCTK_ATTRIBUTE_UNUSED = (csetemp66 + csetemp67 + 
      csetemp68 + csetemp69)*xformL02 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL12 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL22 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL32;
    
    CCTK_REAL g433 CCTK_ATTRIBUTE_UNUSED = (csetemp66 + csetemp67 + 
      csetemp68 + csetemp69)*xformL03 + (csetemp70 + csetemp71 + csetemp72 + 
      csetemp73)*xformL13 + (csetemp74 + csetemp75 + csetemp76 + 
      csetemp77)*xformL23 + (csetemp78 + csetemp79 + csetemp80 + 
      csetemp81)*xformL33;
    
    CCTK_REAL csetemp82 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL00;
    
    CCTK_REAL csetemp83 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL10;
    
    CCTK_REAL csetemp84 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL20;
    
    CCTK_REAL csetemp85 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL30;
    
    CCTK_REAL csetemp86 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL00;
    
    CCTK_REAL csetemp87 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL10;
    
    CCTK_REAL csetemp88 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL20;
    
    CCTK_REAL csetemp89 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL30;
    
    CCTK_REAL csetemp90 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL00;
    
    CCTK_REAL csetemp91 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL10;
    
    CCTK_REAL csetemp92 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL20;
    
    CCTK_REAL csetemp93 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL30;
    
    CCTK_REAL csetemp94 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL00;
    
    CCTK_REAL csetemp95 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL10;
    
    CCTK_REAL csetemp96 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL20;
    
    CCTK_REAL csetemp97 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL30;
    
    CCTK_REAL csetemp98 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL00;
    
    CCTK_REAL csetemp99 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL10;
    
    CCTK_REAL csetemp100 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL20;
    
    CCTK_REAL csetemp101 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL30;
    
    CCTK_REAL csetemp102 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL00;
    
    CCTK_REAL csetemp103 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL10;
    
    CCTK_REAL csetemp104 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL20;
    
    CCTK_REAL csetemp105 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL30;
    
    CCTK_REAL csetemp106 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL00;
    
    CCTK_REAL csetemp107 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL10;
    
    CCTK_REAL csetemp108 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL20;
    
    CCTK_REAL csetemp109 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL30;
    
    CCTK_REAL csetemp110 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL00;
    
    CCTK_REAL csetemp111 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL10;
    
    CCTK_REAL csetemp112 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL20;
    
    CCTK_REAL csetemp113 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL30;
    
    CCTK_REAL csetemp114 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL00;
    
    CCTK_REAL csetemp115 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL10;
    
    CCTK_REAL csetemp116 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL20;
    
    CCTK_REAL csetemp117 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL30;
    
    CCTK_REAL csetemp118 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL00;
    
    CCTK_REAL csetemp119 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL10;
    
    CCTK_REAL csetemp120 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL20;
    
    CCTK_REAL csetemp121 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL30;
    
    CCTK_REAL dg4000 CCTK_ATTRIBUTE_UNUSED = xformL10*((csetemp86 + 
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
    
    CCTK_REAL dg4010 CCTK_ATTRIBUTE_UNUSED = xformL10*((csetemp86 + 
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
    
    CCTK_REAL csetemp122 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL01;
    
    CCTK_REAL csetemp123 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL11;
    
    CCTK_REAL csetemp124 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL21;
    
    CCTK_REAL csetemp125 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL31;
    
    CCTK_REAL csetemp126 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL01;
    
    CCTK_REAL csetemp127 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL11;
    
    CCTK_REAL csetemp128 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL21;
    
    CCTK_REAL csetemp129 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL31;
    
    CCTK_REAL csetemp130 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL01;
    
    CCTK_REAL csetemp131 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL11;
    
    CCTK_REAL csetemp132 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL21;
    
    CCTK_REAL csetemp133 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL31;
    
    CCTK_REAL csetemp134 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL01;
    
    CCTK_REAL csetemp135 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL11;
    
    CCTK_REAL csetemp136 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL21;
    
    CCTK_REAL csetemp137 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL31;
    
    CCTK_REAL csetemp138 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL01;
    
    CCTK_REAL csetemp139 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL11;
    
    CCTK_REAL csetemp140 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL21;
    
    CCTK_REAL csetemp141 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL31;
    
    CCTK_REAL csetemp142 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL01;
    
    CCTK_REAL csetemp143 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL11;
    
    CCTK_REAL csetemp144 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL21;
    
    CCTK_REAL csetemp145 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL31;
    
    CCTK_REAL csetemp146 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL01;
    
    CCTK_REAL csetemp147 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL11;
    
    CCTK_REAL csetemp148 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL21;
    
    CCTK_REAL csetemp149 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL31;
    
    CCTK_REAL csetemp150 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL01;
    
    CCTK_REAL csetemp151 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL11;
    
    CCTK_REAL csetemp152 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL21;
    
    CCTK_REAL csetemp153 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL31;
    
    CCTK_REAL csetemp154 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL01;
    
    CCTK_REAL csetemp155 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL11;
    
    CCTK_REAL csetemp156 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL21;
    
    CCTK_REAL csetemp157 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL31;
    
    CCTK_REAL csetemp158 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL01;
    
    CCTK_REAL csetemp159 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL11;
    
    CCTK_REAL csetemp160 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL21;
    
    CCTK_REAL csetemp161 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL31;
    
    CCTK_REAL dg4011 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp122 + 
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
    
    CCTK_REAL csetemp162 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL02;
    
    CCTK_REAL csetemp163 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL12;
    
    CCTK_REAL csetemp164 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL22;
    
    CCTK_REAL csetemp165 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL32;
    
    CCTK_REAL csetemp166 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL02;
    
    CCTK_REAL csetemp167 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL12;
    
    CCTK_REAL csetemp168 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL22;
    
    CCTK_REAL csetemp169 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL32;
    
    CCTK_REAL csetemp170 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL02;
    
    CCTK_REAL csetemp171 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL12;
    
    CCTK_REAL csetemp172 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL22;
    
    CCTK_REAL csetemp173 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL32;
    
    CCTK_REAL csetemp174 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL02;
    
    CCTK_REAL csetemp175 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL12;
    
    CCTK_REAL csetemp176 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL22;
    
    CCTK_REAL csetemp177 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL32;
    
    CCTK_REAL csetemp178 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL02;
    
    CCTK_REAL csetemp179 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL12;
    
    CCTK_REAL csetemp180 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL22;
    
    CCTK_REAL csetemp181 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL32;
    
    CCTK_REAL csetemp182 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL02;
    
    CCTK_REAL csetemp183 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL12;
    
    CCTK_REAL csetemp184 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL22;
    
    CCTK_REAL csetemp185 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL32;
    
    CCTK_REAL csetemp186 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL02;
    
    CCTK_REAL csetemp187 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL12;
    
    CCTK_REAL csetemp188 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL22;
    
    CCTK_REAL csetemp189 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL32;
    
    CCTK_REAL csetemp190 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL02;
    
    CCTK_REAL csetemp191 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL12;
    
    CCTK_REAL csetemp192 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL22;
    
    CCTK_REAL csetemp193 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL32;
    
    CCTK_REAL csetemp194 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL02;
    
    CCTK_REAL csetemp195 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL12;
    
    CCTK_REAL csetemp196 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL22;
    
    CCTK_REAL csetemp197 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL32;
    
    CCTK_REAL csetemp198 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL02;
    
    CCTK_REAL csetemp199 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL12;
    
    CCTK_REAL csetemp200 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL22;
    
    CCTK_REAL csetemp201 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL32;
    
    CCTK_REAL dg4012 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp162 + 
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
    
    CCTK_REAL csetemp202 CCTK_ATTRIBUTE_UNUSED = tdg4000*xformL03;
    
    CCTK_REAL csetemp203 CCTK_ATTRIBUTE_UNUSED = tdg4001*xformL13;
    
    CCTK_REAL csetemp204 CCTK_ATTRIBUTE_UNUSED = tdg4002*xformL23;
    
    CCTK_REAL csetemp205 CCTK_ATTRIBUTE_UNUSED = tdg4003*xformL33;
    
    CCTK_REAL csetemp206 CCTK_ATTRIBUTE_UNUSED = tdg4010*xformL03;
    
    CCTK_REAL csetemp207 CCTK_ATTRIBUTE_UNUSED = tdg4011*xformL13;
    
    CCTK_REAL csetemp208 CCTK_ATTRIBUTE_UNUSED = tdg4012*xformL23;
    
    CCTK_REAL csetemp209 CCTK_ATTRIBUTE_UNUSED = tdg4013*xformL33;
    
    CCTK_REAL csetemp210 CCTK_ATTRIBUTE_UNUSED = tdg4020*xformL03;
    
    CCTK_REAL csetemp211 CCTK_ATTRIBUTE_UNUSED = tdg4021*xformL13;
    
    CCTK_REAL csetemp212 CCTK_ATTRIBUTE_UNUSED = tdg4022*xformL23;
    
    CCTK_REAL csetemp213 CCTK_ATTRIBUTE_UNUSED = tdg4023*xformL33;
    
    CCTK_REAL csetemp214 CCTK_ATTRIBUTE_UNUSED = tdg4030*xformL03;
    
    CCTK_REAL csetemp215 CCTK_ATTRIBUTE_UNUSED = tdg4031*xformL13;
    
    CCTK_REAL csetemp216 CCTK_ATTRIBUTE_UNUSED = tdg4032*xformL23;
    
    CCTK_REAL csetemp217 CCTK_ATTRIBUTE_UNUSED = tdg4033*xformL33;
    
    CCTK_REAL csetemp218 CCTK_ATTRIBUTE_UNUSED = tdg4110*xformL03;
    
    CCTK_REAL csetemp219 CCTK_ATTRIBUTE_UNUSED = tdg4111*xformL13;
    
    CCTK_REAL csetemp220 CCTK_ATTRIBUTE_UNUSED = tdg4112*xformL23;
    
    CCTK_REAL csetemp221 CCTK_ATTRIBUTE_UNUSED = tdg4113*xformL33;
    
    CCTK_REAL csetemp222 CCTK_ATTRIBUTE_UNUSED = tdg4120*xformL03;
    
    CCTK_REAL csetemp223 CCTK_ATTRIBUTE_UNUSED = tdg4121*xformL13;
    
    CCTK_REAL csetemp224 CCTK_ATTRIBUTE_UNUSED = tdg4122*xformL23;
    
    CCTK_REAL csetemp225 CCTK_ATTRIBUTE_UNUSED = tdg4123*xformL33;
    
    CCTK_REAL csetemp226 CCTK_ATTRIBUTE_UNUSED = tdg4130*xformL03;
    
    CCTK_REAL csetemp227 CCTK_ATTRIBUTE_UNUSED = tdg4131*xformL13;
    
    CCTK_REAL csetemp228 CCTK_ATTRIBUTE_UNUSED = tdg4132*xformL23;
    
    CCTK_REAL csetemp229 CCTK_ATTRIBUTE_UNUSED = tdg4133*xformL33;
    
    CCTK_REAL csetemp230 CCTK_ATTRIBUTE_UNUSED = tdg4220*xformL03;
    
    CCTK_REAL csetemp231 CCTK_ATTRIBUTE_UNUSED = tdg4221*xformL13;
    
    CCTK_REAL csetemp232 CCTK_ATTRIBUTE_UNUSED = tdg4222*xformL23;
    
    CCTK_REAL csetemp233 CCTK_ATTRIBUTE_UNUSED = tdg4223*xformL33;
    
    CCTK_REAL csetemp234 CCTK_ATTRIBUTE_UNUSED = tdg4230*xformL03;
    
    CCTK_REAL csetemp235 CCTK_ATTRIBUTE_UNUSED = tdg4231*xformL13;
    
    CCTK_REAL csetemp236 CCTK_ATTRIBUTE_UNUSED = tdg4232*xformL23;
    
    CCTK_REAL csetemp237 CCTK_ATTRIBUTE_UNUSED = tdg4233*xformL33;
    
    CCTK_REAL csetemp238 CCTK_ATTRIBUTE_UNUSED = tdg4330*xformL03;
    
    CCTK_REAL csetemp239 CCTK_ATTRIBUTE_UNUSED = tdg4331*xformL13;
    
    CCTK_REAL csetemp240 CCTK_ATTRIBUTE_UNUSED = tdg4332*xformL23;
    
    CCTK_REAL csetemp241 CCTK_ATTRIBUTE_UNUSED = tdg4333*xformL33;
    
    CCTK_REAL dg4013 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp202 + 
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
    
    CCTK_REAL dg4020 CCTK_ATTRIBUTE_UNUSED = xformL10*((csetemp86 + 
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
    
    CCTK_REAL dg4021 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp122 + 
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
    
    CCTK_REAL dg4022 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp162 + 
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
    
    CCTK_REAL dg4023 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp202 + 
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
    
    CCTK_REAL dg4030 CCTK_ATTRIBUTE_UNUSED = xformL10*((csetemp86 + 
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
    
    CCTK_REAL dg4031 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp122 + 
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
    
    CCTK_REAL dg4032 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp162 + 
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
    
    CCTK_REAL dg4033 CCTK_ATTRIBUTE_UNUSED = xformL00*((csetemp202 + 
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
    
    CCTK_REAL dg4110 CCTK_ATTRIBUTE_UNUSED = xformL11*((csetemp86 + 
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
    
    CCTK_REAL dg4111 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp122 + 
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
    
    CCTK_REAL dg4112 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp162 + 
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
    
    CCTK_REAL dg4113 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp202 + 
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
    
    CCTK_REAL dg4120 CCTK_ATTRIBUTE_UNUSED = xformL11*((csetemp86 + 
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
    
    CCTK_REAL dg4121 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp122 + 
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
    
    CCTK_REAL dg4122 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp162 + 
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
    
    CCTK_REAL dg4123 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp202 + 
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
    
    CCTK_REAL dg4130 CCTK_ATTRIBUTE_UNUSED = xformL11*((csetemp86 + 
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
    
    CCTK_REAL dg4131 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp122 + 
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
    
    CCTK_REAL dg4132 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp162 + 
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
    
    CCTK_REAL dg4133 CCTK_ATTRIBUTE_UNUSED = xformL01*((csetemp202 + 
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
    
    CCTK_REAL dg4220 CCTK_ATTRIBUTE_UNUSED = xformL12*((csetemp86 + 
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
    
    CCTK_REAL dg4221 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp122 + 
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
    
    CCTK_REAL dg4222 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp162 + 
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
    
    CCTK_REAL dg4223 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp202 + 
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
    
    CCTK_REAL dg4230 CCTK_ATTRIBUTE_UNUSED = xformL12*((csetemp86 + 
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
    
    CCTK_REAL dg4231 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp122 + 
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
    
    CCTK_REAL dg4232 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp162 + 
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
    
    CCTK_REAL dg4233 CCTK_ATTRIBUTE_UNUSED = xformL02*((csetemp202 + 
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
    
    CCTK_REAL dg4330 CCTK_ATTRIBUTE_UNUSED = xformL13*((csetemp86 + 
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
    
    CCTK_REAL dg4331 CCTK_ATTRIBUTE_UNUSED = xformL03*((csetemp122 + 
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
    
    CCTK_REAL dg4332 CCTK_ATTRIBUTE_UNUSED = xformL03*((csetemp162 + 
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
    
    CCTK_REAL dg4333 CCTK_ATTRIBUTE_UNUSED = xformL03*((csetemp202 + 
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
    
    CCTK_REAL betal1 CCTK_ATTRIBUTE_UNUSED = g401;
    
    CCTK_REAL betal2 CCTK_ATTRIBUTE_UNUSED = g402;
    
    CCTK_REAL betal3 CCTK_ATTRIBUTE_UNUSED = g403;
    
    gxxL = g411;
    
    gxyL = g412;
    
    gxzL = g413;
    
    gyyL = g422;
    
    gyzL = g423;
    
    gzzL = g433;
    
    CCTK_REAL csetemp242 CCTK_ATTRIBUTE_UNUSED = pow(gxzL,2);
    
    CCTK_REAL csetemp243 CCTK_ATTRIBUTE_UNUSED = pow(gyzL,2);
    
    CCTK_REAL csetemp244 CCTK_ATTRIBUTE_UNUSED = pow(gxyL,2);
    
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = 2*gxyL*gxzL*gyzL + 
      gyyL*(gxxL*gzzL - csetemp242) - gxxL*csetemp243 - gzzL*csetemp244;
    
    CCTK_REAL csetemp245 CCTK_ATTRIBUTE_UNUSED = pow(detg,-1);
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = (gyyL*gzzL - 
      csetemp243)*csetemp245;
    
    CCTK_REAL gu12 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyzL - 
      gxyL*gzzL)*csetemp245;
    
    CCTK_REAL gu13 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyyL) + 
      gxyL*gyzL)*csetemp245;
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = (gxxL*gzzL - 
      csetemp242)*csetemp245;
    
    CCTK_REAL gu23 CCTK_ATTRIBUTE_UNUSED = (gxyL*gxzL - 
      gxxL*gyzL)*csetemp245;
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = (gxxL*gyyL - 
      csetemp244)*csetemp245;
    
    betaxL = betal1*gu11 + betal2*gu12 + betal3*gu13;
    
    betayL = betal1*gu12 + betal2*gu22 + betal3*gu23;
    
    betazL = betal1*gu13 + betal2*gu23 + betal3*gu33;
    
    CCTK_REAL betasq CCTK_ATTRIBUTE_UNUSED = betaxL*betal1 + betayL*betal2 
      + betazL*betal3;
    
    alpL = pow(betasq - g400,0.5);
    
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
    
    CCTK_REAL csetemp246 CCTK_ATTRIBUTE_UNUSED = dtg11*gu11;
    
    CCTK_REAL csetemp247 CCTK_ATTRIBUTE_UNUSED = dtg12*gu12;
    
    CCTK_REAL csetemp248 CCTK_ATTRIBUTE_UNUSED = dtg13*gu13;
    
    CCTK_REAL csetemp249 CCTK_ATTRIBUTE_UNUSED = dtg12*gu11;
    
    CCTK_REAL csetemp250 CCTK_ATTRIBUTE_UNUSED = dtg22*gu12;
    
    CCTK_REAL csetemp251 CCTK_ATTRIBUTE_UNUSED = dtg23*gu13;
    
    CCTK_REAL csetemp252 CCTK_ATTRIBUTE_UNUSED = dtg13*gu11;
    
    CCTK_REAL csetemp253 CCTK_ATTRIBUTE_UNUSED = dtg23*gu12;
    
    CCTK_REAL csetemp254 CCTK_ATTRIBUTE_UNUSED = dtg33*gu13;
    
    CCTK_REAL dtgu11 CCTK_ATTRIBUTE_UNUSED = -((csetemp246 + csetemp247 + 
      csetemp248)*gu11) - (csetemp249 + csetemp250 + csetemp251)*gu12 - 
      (csetemp252 + csetemp253 + csetemp254)*gu13;
    
    CCTK_REAL dtgu12 CCTK_ATTRIBUTE_UNUSED = -((csetemp246 + csetemp247 + 
      csetemp248)*gu12) - (csetemp249 + csetemp250 + csetemp251)*gu22 - 
      (csetemp252 + csetemp253 + csetemp254)*gu23;
    
    CCTK_REAL dtgu13 CCTK_ATTRIBUTE_UNUSED = -((csetemp246 + csetemp247 + 
      csetemp248)*gu13) - (csetemp249 + csetemp250 + csetemp251)*gu23 - 
      (csetemp252 + csetemp253 + csetemp254)*gu33;
    
    CCTK_REAL csetemp255 CCTK_ATTRIBUTE_UNUSED = dtg11*gu12;
    
    CCTK_REAL csetemp256 CCTK_ATTRIBUTE_UNUSED = dtg12*gu22;
    
    CCTK_REAL csetemp257 CCTK_ATTRIBUTE_UNUSED = dtg13*gu23;
    
    CCTK_REAL csetemp258 CCTK_ATTRIBUTE_UNUSED = dtg22*gu22;
    
    CCTK_REAL csetemp259 CCTK_ATTRIBUTE_UNUSED = dtg23*gu23;
    
    CCTK_REAL csetemp260 CCTK_ATTRIBUTE_UNUSED = dtg13*gu12;
    
    CCTK_REAL csetemp261 CCTK_ATTRIBUTE_UNUSED = dtg23*gu22;
    
    CCTK_REAL csetemp262 CCTK_ATTRIBUTE_UNUSED = dtg33*gu23;
    
    CCTK_REAL dtgu22 CCTK_ATTRIBUTE_UNUSED = -((csetemp255 + csetemp256 + 
      csetemp257)*gu12) - (csetemp247 + csetemp258 + csetemp259)*gu22 - 
      (csetemp260 + csetemp261 + csetemp262)*gu23;
    
    CCTK_REAL dtgu23 CCTK_ATTRIBUTE_UNUSED = -((csetemp255 + csetemp256 + 
      csetemp257)*gu13) - (csetemp247 + csetemp258 + csetemp259)*gu23 - 
      (csetemp260 + csetemp261 + csetemp262)*gu33;
    
    CCTK_REAL dtgu33 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dtg11*gu13 + 
      dtg12*gu23 + dtg13*gu33)) - gu23*(dtg12*gu13 + dtg22*gu23 + dtg23*gu33) 
      - gu33*(csetemp248 + csetemp259 + dtg33*gu33);
    
    CCTK_REAL csetemp263 CCTK_ATTRIBUTE_UNUSED = dg111*gu11;
    
    CCTK_REAL csetemp264 CCTK_ATTRIBUTE_UNUSED = dg121*gu12;
    
    CCTK_REAL csetemp265 CCTK_ATTRIBUTE_UNUSED = dg131*gu13;
    
    CCTK_REAL csetemp266 CCTK_ATTRIBUTE_UNUSED = dg121*gu11;
    
    CCTK_REAL csetemp267 CCTK_ATTRIBUTE_UNUSED = dg221*gu12;
    
    CCTK_REAL csetemp268 CCTK_ATTRIBUTE_UNUSED = dg231*gu13;
    
    CCTK_REAL csetemp269 CCTK_ATTRIBUTE_UNUSED = dg131*gu11;
    
    CCTK_REAL csetemp270 CCTK_ATTRIBUTE_UNUSED = dg231*gu12;
    
    CCTK_REAL csetemp271 CCTK_ATTRIBUTE_UNUSED = dg331*gu13;
    
    CCTK_REAL dgu111 CCTK_ATTRIBUTE_UNUSED = -((csetemp263 + csetemp264 + 
      csetemp265)*gu11) - (csetemp266 + csetemp267 + csetemp268)*gu12 - 
      (csetemp269 + csetemp270 + csetemp271)*gu13;
    
    CCTK_REAL dgu121 CCTK_ATTRIBUTE_UNUSED = -((csetemp263 + csetemp264 + 
      csetemp265)*gu12) - (csetemp266 + csetemp267 + csetemp268)*gu22 - 
      (csetemp269 + csetemp270 + csetemp271)*gu23;
    
    CCTK_REAL dgu131 CCTK_ATTRIBUTE_UNUSED = -((csetemp263 + csetemp264 + 
      csetemp265)*gu13) - (csetemp266 + csetemp267 + csetemp268)*gu23 - 
      (csetemp269 + csetemp270 + csetemp271)*gu33;
    
    CCTK_REAL csetemp272 CCTK_ATTRIBUTE_UNUSED = dg111*gu12;
    
    CCTK_REAL csetemp273 CCTK_ATTRIBUTE_UNUSED = dg121*gu22;
    
    CCTK_REAL csetemp274 CCTK_ATTRIBUTE_UNUSED = dg131*gu23;
    
    CCTK_REAL csetemp275 CCTK_ATTRIBUTE_UNUSED = dg221*gu22;
    
    CCTK_REAL csetemp276 CCTK_ATTRIBUTE_UNUSED = dg231*gu23;
    
    CCTK_REAL csetemp277 CCTK_ATTRIBUTE_UNUSED = dg131*gu12;
    
    CCTK_REAL csetemp278 CCTK_ATTRIBUTE_UNUSED = dg231*gu22;
    
    CCTK_REAL csetemp279 CCTK_ATTRIBUTE_UNUSED = dg331*gu23;
    
    CCTK_REAL dgu221 CCTK_ATTRIBUTE_UNUSED = -((csetemp272 + csetemp273 + 
      csetemp274)*gu12) - (csetemp264 + csetemp275 + csetemp276)*gu22 - 
      (csetemp277 + csetemp278 + csetemp279)*gu23;
    
    CCTK_REAL dgu231 CCTK_ATTRIBUTE_UNUSED = -((csetemp272 + csetemp273 + 
      csetemp274)*gu13) - (csetemp264 + csetemp275 + csetemp276)*gu23 - 
      (csetemp277 + csetemp278 + csetemp279)*gu33;
    
    CCTK_REAL dgu331 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dg111*gu13 + 
      dg121*gu23 + dg131*gu33)) - gu23*(dg121*gu13 + dg221*gu23 + dg231*gu33) 
      - gu33*(csetemp265 + csetemp276 + dg331*gu33);
    
    CCTK_REAL csetemp280 CCTK_ATTRIBUTE_UNUSED = dg112*gu11;
    
    CCTK_REAL csetemp281 CCTK_ATTRIBUTE_UNUSED = dg122*gu12;
    
    CCTK_REAL csetemp282 CCTK_ATTRIBUTE_UNUSED = dg132*gu13;
    
    CCTK_REAL csetemp283 CCTK_ATTRIBUTE_UNUSED = dg122*gu11;
    
    CCTK_REAL csetemp284 CCTK_ATTRIBUTE_UNUSED = dg222*gu12;
    
    CCTK_REAL csetemp285 CCTK_ATTRIBUTE_UNUSED = dg232*gu13;
    
    CCTK_REAL csetemp286 CCTK_ATTRIBUTE_UNUSED = dg132*gu11;
    
    CCTK_REAL csetemp287 CCTK_ATTRIBUTE_UNUSED = dg232*gu12;
    
    CCTK_REAL csetemp288 CCTK_ATTRIBUTE_UNUSED = dg332*gu13;
    
    CCTK_REAL dgu112 CCTK_ATTRIBUTE_UNUSED = -((csetemp280 + csetemp281 + 
      csetemp282)*gu11) - (csetemp283 + csetemp284 + csetemp285)*gu12 - 
      (csetemp286 + csetemp287 + csetemp288)*gu13;
    
    CCTK_REAL dgu122 CCTK_ATTRIBUTE_UNUSED = -((csetemp280 + csetemp281 + 
      csetemp282)*gu12) - (csetemp283 + csetemp284 + csetemp285)*gu22 - 
      (csetemp286 + csetemp287 + csetemp288)*gu23;
    
    CCTK_REAL dgu132 CCTK_ATTRIBUTE_UNUSED = -((csetemp280 + csetemp281 + 
      csetemp282)*gu13) - (csetemp283 + csetemp284 + csetemp285)*gu23 - 
      (csetemp286 + csetemp287 + csetemp288)*gu33;
    
    CCTK_REAL csetemp289 CCTK_ATTRIBUTE_UNUSED = dg112*gu12;
    
    CCTK_REAL csetemp290 CCTK_ATTRIBUTE_UNUSED = dg122*gu22;
    
    CCTK_REAL csetemp291 CCTK_ATTRIBUTE_UNUSED = dg132*gu23;
    
    CCTK_REAL csetemp292 CCTK_ATTRIBUTE_UNUSED = dg222*gu22;
    
    CCTK_REAL csetemp293 CCTK_ATTRIBUTE_UNUSED = dg232*gu23;
    
    CCTK_REAL csetemp294 CCTK_ATTRIBUTE_UNUSED = dg132*gu12;
    
    CCTK_REAL csetemp295 CCTK_ATTRIBUTE_UNUSED = dg232*gu22;
    
    CCTK_REAL csetemp296 CCTK_ATTRIBUTE_UNUSED = dg332*gu23;
    
    CCTK_REAL dgu222 CCTK_ATTRIBUTE_UNUSED = -((csetemp289 + csetemp290 + 
      csetemp291)*gu12) - (csetemp281 + csetemp292 + csetemp293)*gu22 - 
      (csetemp294 + csetemp295 + csetemp296)*gu23;
    
    CCTK_REAL dgu232 CCTK_ATTRIBUTE_UNUSED = -((csetemp289 + csetemp290 + 
      csetemp291)*gu13) - (csetemp281 + csetemp292 + csetemp293)*gu23 - 
      (csetemp294 + csetemp295 + csetemp296)*gu33;
    
    CCTK_REAL dgu332 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dg112*gu13 + 
      dg122*gu23 + dg132*gu33)) - gu23*(dg122*gu13 + dg222*gu23 + dg232*gu33) 
      - gu33*(csetemp282 + csetemp293 + dg332*gu33);
    
    CCTK_REAL csetemp297 CCTK_ATTRIBUTE_UNUSED = dg113*gu11;
    
    CCTK_REAL csetemp298 CCTK_ATTRIBUTE_UNUSED = dg123*gu12;
    
    CCTK_REAL csetemp299 CCTK_ATTRIBUTE_UNUSED = dg133*gu13;
    
    CCTK_REAL csetemp300 CCTK_ATTRIBUTE_UNUSED = dg123*gu11;
    
    CCTK_REAL csetemp301 CCTK_ATTRIBUTE_UNUSED = dg223*gu12;
    
    CCTK_REAL csetemp302 CCTK_ATTRIBUTE_UNUSED = dg233*gu13;
    
    CCTK_REAL csetemp303 CCTK_ATTRIBUTE_UNUSED = dg133*gu11;
    
    CCTK_REAL csetemp304 CCTK_ATTRIBUTE_UNUSED = dg233*gu12;
    
    CCTK_REAL csetemp305 CCTK_ATTRIBUTE_UNUSED = dg333*gu13;
    
    CCTK_REAL dgu113 CCTK_ATTRIBUTE_UNUSED = -((csetemp297 + csetemp298 + 
      csetemp299)*gu11) - (csetemp300 + csetemp301 + csetemp302)*gu12 - 
      (csetemp303 + csetemp304 + csetemp305)*gu13;
    
    CCTK_REAL dgu123 CCTK_ATTRIBUTE_UNUSED = -((csetemp297 + csetemp298 + 
      csetemp299)*gu12) - (csetemp300 + csetemp301 + csetemp302)*gu22 - 
      (csetemp303 + csetemp304 + csetemp305)*gu23;
    
    CCTK_REAL dgu133 CCTK_ATTRIBUTE_UNUSED = -((csetemp297 + csetemp298 + 
      csetemp299)*gu13) - (csetemp300 + csetemp301 + csetemp302)*gu23 - 
      (csetemp303 + csetemp304 + csetemp305)*gu33;
    
    CCTK_REAL csetemp306 CCTK_ATTRIBUTE_UNUSED = dg113*gu12;
    
    CCTK_REAL csetemp307 CCTK_ATTRIBUTE_UNUSED = dg123*gu22;
    
    CCTK_REAL csetemp308 CCTK_ATTRIBUTE_UNUSED = dg133*gu23;
    
    CCTK_REAL csetemp309 CCTK_ATTRIBUTE_UNUSED = dg223*gu22;
    
    CCTK_REAL csetemp310 CCTK_ATTRIBUTE_UNUSED = dg233*gu23;
    
    CCTK_REAL csetemp311 CCTK_ATTRIBUTE_UNUSED = dg133*gu12;
    
    CCTK_REAL csetemp312 CCTK_ATTRIBUTE_UNUSED = dg233*gu22;
    
    CCTK_REAL csetemp313 CCTK_ATTRIBUTE_UNUSED = dg333*gu23;
    
    CCTK_REAL dgu223 CCTK_ATTRIBUTE_UNUSED = -((csetemp306 + csetemp307 + 
      csetemp308)*gu12) - (csetemp298 + csetemp309 + csetemp310)*gu22 - 
      (csetemp311 + csetemp312 + csetemp313)*gu23;
    
    CCTK_REAL dgu233 CCTK_ATTRIBUTE_UNUSED = -((csetemp306 + csetemp307 + 
      csetemp308)*gu13) - (csetemp298 + csetemp309 + csetemp310)*gu23 - 
      (csetemp311 + csetemp312 + csetemp313)*gu33;
    
    CCTK_REAL dgu333 CCTK_ATTRIBUTE_UNUSED = -(gu13*(dg113*gu13 + 
      dg123*gu23 + dg133*gu33)) - gu23*(dg123*gu13 + dg223*gu23 + dg233*gu33) 
      - gu33*(csetemp299 + csetemp310 + dg333*gu33);
    
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
    
    dtbetaxL = betal1*dtgu11 + betal2*dtgu12 + betal3*dtgu13 + 
      dtbetal1*gu11 + dtbetal2*gu12 + dtbetal3*gu13;
    
    dtbetayL = betal1*dtgu12 + betal2*dtgu22 + betal3*dtgu23 + 
      dtbetal1*gu12 + dtbetal2*gu22 + dtbetal3*gu23;
    
    dtbetazL = betal1*dtgu13 + betal2*dtgu23 + betal3*dtgu33 + 
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
    
    CCTK_REAL dtbetasq CCTK_ATTRIBUTE_UNUSED = dtbetaxL*betal1 + 
      dtbetayL*betal2 + dtbetazL*betal3 + betaxL*dtbetal1 + betayL*dtbetal2 + 
      betazL*dtbetal3;
    
    CCTK_REAL csetemp314 CCTK_ATTRIBUTE_UNUSED = pow(alpL,-1);
    
    CCTK_REAL dtalpL CCTK_ATTRIBUTE_UNUSED = 0.5*csetemp314*(-dg4000 + 
      dtbetasq);
    
    CCTK_REAL kxxL CCTK_ATTRIBUTE_UNUSED = 0.5*csetemp314*(2*(gxxL*dbeta11 
      + gxyL*dbeta21 + gxzL*dbeta31) + betaxL*dg111 + betayL*dg112 + 
      betazL*dg113 - dtg11);
    
    CCTK_REAL kxyL CCTK_ATTRIBUTE_UNUSED = 0.5*csetemp314*(gxxL*dbeta12 + 
      gyyL*dbeta21 + gxyL*(dbeta11 + dbeta22) + gyzL*dbeta31 + gxzL*dbeta32 + 
      betaxL*dg121 + betayL*dg122 + betazL*dg123 - dtg12);
    
    CCTK_REAL kxzL CCTK_ATTRIBUTE_UNUSED = 0.5*csetemp314*(gxxL*dbeta13 + 
      gyzL*dbeta21 + gxyL*dbeta23 + gzzL*dbeta31 + gxzL*(dbeta11 + dbeta33) + 
      betaxL*dg131 + betayL*dg132 + betazL*dg133 - dtg13);
    
    CCTK_REAL kyyL CCTK_ATTRIBUTE_UNUSED = 0.5*csetemp314*(2*(gxyL*dbeta12 
      + gyyL*dbeta22 + gyzL*dbeta32) + betaxL*dg221 + betayL*dg222 + 
      betazL*dg223 - dtg22);
    
    CCTK_REAL kyzL CCTK_ATTRIBUTE_UNUSED = 0.5*csetemp314*(gxzL*dbeta12 + 
      gxyL*dbeta13 + gyyL*dbeta23 + gzzL*dbeta32 + gyzL*(dbeta22 + dbeta33) + 
      betaxL*dg231 + betayL*dg232 + betazL*dg233 - dtg23);
    
    CCTK_REAL kzzL CCTK_ATTRIBUTE_UNUSED = 0.5*csetemp314*(2*(gxzL*dbeta13 
      + gyzL*dbeta23 + gzzL*dbeta33) + betaxL*dg331 + betayL*dg332 + 
      betazL*dg333 - dtg33);
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
  CCTK_ENDLOOP3(KerrSchild_always);
}
extern "C" void KerrSchild_always(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_KerrSchild_always;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering KerrSchild_always_Body");
  }
  if (cctk_iteration % KerrSchild_always_calc_every != KerrSchild_always_calc_offset)
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
  AssertGroupStorage(cctkGH, "KerrSchild_always", 7, groups);
  
  
  LoopOverEverything(cctkGH, KerrSchild_always_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving KerrSchild_always_Body");
  }
}

} // namespace KerrSchild
