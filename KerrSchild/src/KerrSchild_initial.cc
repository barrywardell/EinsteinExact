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
  CCTK_LOOP3(KerrSchild_initial,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL xL = x[index];
    CCTK_REAL yL = y[index];
    CCTK_REAL zL = z[index];
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL xx1 = xL;
    
    CCTK_REAL xx2 = yL;
    
    CCTK_REAL xx3 = zL;
    
    CCTK_REAL position1 = ToReal(positionx);
    
    CCTK_REAL position2 = ToReal(positiony);
    
    CCTK_REAL position3 = ToReal(positionz);
    
    CCTK_REAL shiftadd1 = ToReal(shiftaddx);
    
    CCTK_REAL shiftadd2 = ToReal(shiftaddy);
    
    CCTK_REAL shiftadd3 = ToReal(shiftaddz);
    
    CCTK_REAL csetemp0 = cos(ToReal(phi));
    
    CCTK_REAL csetemp1 = cos(ToReal(psi));
    
    CCTK_REAL csetemp2 = cos(ToReal(theta));
    
    CCTK_REAL csetemp3 = sin(ToReal(phi));
    
    CCTK_REAL csetemp4 = sin(ToReal(psi));
    
    CCTK_REAL Jac11 = csetemp0*csetemp1 - csetemp2*csetemp3*csetemp4;
    
    CCTK_REAL Jac12 = csetemp1*csetemp3 + csetemp0*csetemp2*csetemp4;
    
    CCTK_REAL csetemp5 = sin(ToReal(theta));
    
    CCTK_REAL Jac13 = csetemp4*csetemp5;
    
    CCTK_REAL Jac21 = -(csetemp1*csetemp2*csetemp3) - csetemp0*csetemp4;
    
    CCTK_REAL Jac22 = csetemp0*csetemp1*csetemp2 - csetemp3*csetemp4;
    
    CCTK_REAL Jac23 = csetemp1*csetemp5;
    
    CCTK_REAL Jac31 = csetemp3*csetemp5;
    
    CCTK_REAL Jac32 = -(csetemp0*csetemp5);
    
    CCTK_REAL Jac33 = csetemp2;
    
    CCTK_REAL InvJac11 = Jac11;
    
    CCTK_REAL InvJac12 = Jac21;
    
    CCTK_REAL InvJac13 = Jac31;
    
    CCTK_REAL InvJac21 = Jac12;
    
    CCTK_REAL InvJac22 = Jac22;
    
    CCTK_REAL InvJac23 = Jac32;
    
    CCTK_REAL InvJac31 = Jac13;
    
    CCTK_REAL InvJac32 = Jac23;
    
    CCTK_REAL InvJac33 = Jac33;
    
    CCTK_REAL T = ToReal(lapsefactor)*(t - ToReal(positiont));
    
    CCTK_REAL csetemp6 = -(shiftadd1*T);
    
    CCTK_REAL csetemp7 = -(shiftadd2*T);
    
    CCTK_REAL csetemp8 = -(shiftadd3*T);
    
    CCTK_REAL XX1 = Jac11*(csetemp6 - position1 + xx1) + Jac12*(csetemp7 - 
      position2 + xx2) + Jac13*(csetemp8 - position3 + xx3);
    
    CCTK_REAL XX2 = Jac21*(csetemp6 - position1 + xx1) + Jac22*(csetemp7 - 
      position2 + xx2) + Jac23*(csetemp8 - position3 + xx3);
    
    CCTK_REAL XX3 = Jac31*(csetemp6 - position1 + xx1) + Jac32*(csetemp7 - 
      position2 + xx2) + Jac33*(csetemp8 - position3 + xx3);
    
    CCTK_REAL X = XX1;
    
    CCTK_REAL Y = XX2;
    
    CCTK_REAL Z = XX3;
    
    CCTK_REAL csetemp9 = SQR(ToReal(a));
    
    CCTK_REAL csetemp10 = SQR(X);
    
    CCTK_REAL csetemp11 = SQR(Y);
    
    CCTK_REAL csetemp12 = SQR(Z);
    
    CCTK_REAL rXYZ = sqrt(INV(2)*(csetemp10 + csetemp11 + csetemp12 - 
      csetemp9 + sqrt(4*csetemp12*csetemp9 + SQR(csetemp10 + csetemp11 + 
      csetemp12 - csetemp9))));
    
    CCTK_REAL csetemp13 = CUB(rXYZ);
    
    CCTK_REAL csetemp14 = SQR(rXYZ);
    
    CCTK_REAL csetemp15 = rXYZ*X;
    
    CCTK_REAL csetemp16 = Y*ToReal(a);
    
    CCTK_REAL csetemp17 = QAD(rXYZ);
    
    CCTK_REAL G11 = 1 + 2*csetemp13*INV((csetemp17 + 
      csetemp12*csetemp9)*SQR(csetemp14 + csetemp9))*SQR(csetemp15 + 
      csetemp16)*ToReal(M);
    
    CCTK_REAL csetemp18 = X*ToReal(a);
    
    CCTK_REAL csetemp19 = -csetemp18;
    
    CCTK_REAL csetemp20 = rXYZ*Y;
    
    CCTK_REAL G21 = 2*csetemp13*(csetemp15 + csetemp16)*(csetemp19 + 
      csetemp20)*INV((csetemp17 + csetemp12*csetemp9)*SQR(csetemp14 + 
      csetemp9))*ToReal(M);
    
    CCTK_REAL G31 = 2*csetemp14*(csetemp15 + csetemp16)*Z*INV((csetemp14 + 
      csetemp9)*(csetemp17 + csetemp12*csetemp9))*ToReal(M);
    
    CCTK_REAL G22 = 1 + 2*csetemp13*INV((csetemp17 + 
      csetemp12*csetemp9)*SQR(csetemp14 + csetemp9))*SQR(csetemp19 + 
      csetemp20)*ToReal(M);
    
    CCTK_REAL G32 = 2*csetemp14*(csetemp19 + csetemp20)*Z*INV((csetemp14 + 
      csetemp9)*(csetemp17 + csetemp12*csetemp9))*ToReal(M);
    
    CCTK_REAL G33 = 1 + 2*csetemp12*rXYZ*INV(csetemp17 + 
      csetemp12*csetemp9)*ToReal(M);
    
    CCTK_REAL csetemp21 = 4*Y*ToReal(M);
    
    CCTK_REAL csetemp22 = csetemp21*X;
    
    CCTK_REAL csetemp23 = 2*X*Y*ToReal(M);
    
    CCTK_REAL csetemp24 = pow(rXYZ,8);
    
    CCTK_REAL csetemp25 = pow(rXYZ,5);
    
    CCTK_REAL csetemp26 = pow(rXYZ,7);
    
    CCTK_REAL csetemp27 = QAD(ToReal(a));
    
    CCTK_REAL csetemp28 = pow(rXYZ,6);
    
    CCTK_REAL csetemp29 = CUB(ToReal(a));
    
    CCTK_REAL csetemp30 = pow(ToReal(a),5);
    
    CCTK_REAL csetemp31 = pow(rXYZ,9);
    
    CCTK_REAL csetemp32 = pow(ToReal(a),6);
    
    CCTK_REAL K11 = -2*csetemp13*INV(QAD(csetemp14 + 
      csetemp9)*SQR(csetemp17 + csetemp12*csetemp9)*SQR(csetemp17 + 
      csetemp12*csetemp9 + 2*csetemp13*ToReal(M))*sqrt(4*csetemp12*csetemp9 + 
      SQR(csetemp10 + csetemp11 + csetemp12 - csetemp9)))*sqrt(1 + 
      2*csetemp13*INV(csetemp17 + 
      csetemp12*csetemp9)*ToReal(M))*ToReal(M)*((-(csetemp12*(csetemp15 + 
      csetemp16)*SQR(csetemp14 + csetemp9)*(-3*csetemp26*X - 
      3*csetemp25*csetemp9*X + csetemp13*(2*csetemp10 + 2*csetemp11 + 
      3*csetemp12 - 2*csetemp9)*csetemp9*X + 5*csetemp12*csetemp27*rXYZ*X - 
      5*csetemp17*csetemp29*Y + 3*csetemp12*csetemp30*Y + 
      csetemp14*csetemp29*(2*csetemp10 + 2*csetemp11 + csetemp12 - 
      2*csetemp9)*Y - 5*csetemp28*Y*ToReal(a))) - csetemp13*SQR(csetemp15 + 
      csetemp16)*(4*csetemp31 + csetemp26*(-5*csetemp10 - 2*(csetemp11 + 
      csetemp12) + 6*csetemp9) + csetemp9*(-(csetemp12*csetemp13*(csetemp10 + 
      2*(csetemp11 + csetemp12 - 3*csetemp9))) + csetemp25*(-csetemp10 + 
      2*(-csetemp11 + csetemp12 + csetemp9))) + (-(csetemp17*csetemp29) + 
      3*csetemp12*csetemp30)*X*Y + csetemp12*(csetemp27*(3*csetemp10 - 
      2*(csetemp11 + csetemp12) + 2*csetemp9)*rXYZ - csetemp14*csetemp29*X*Y) 
      - 5*csetemp28*X*Y*ToReal(a)))*ToReal(M)*(csetemp17 + csetemp12*csetemp9 
      + 2*csetemp13*ToReal(M)) + (csetemp15 + 
      csetemp16)*ToReal(M)*(4*csetemp12*(csetemp17 + 
      csetemp12*csetemp9)*rXYZ*X*CUB(csetemp14 + csetemp9)*(csetemp17 - 
      csetemp12*csetemp9 + csetemp13*ToReal(M)) - csetemp13*(csetemp19 + 
      csetemp20)*(csetemp17*csetemp29*(-2*csetemp10 - 3*csetemp11 + 
      2*(csetemp12 + csetemp9)) + 
      csetemp12*(-(csetemp14*csetemp29*(2*csetemp10 + 3*csetemp11 + 
      2*csetemp12 - 6*csetemp9)) + csetemp30*(-2*csetemp10 + csetemp11 - 
      2*csetemp12 + 2*csetemp9)) - 3*csetemp26*X*Y + 
      csetemp12*csetemp13*csetemp9*X*Y + csetemp25*csetemp9*X*Y + 
      5*csetemp12*csetemp27*rXYZ*X*Y + (4*csetemp24 - csetemp28*(2*csetemp10 
      + 7*csetemp11 + 2*csetemp12 - 6*csetemp9))*ToReal(a))*(csetemp17 + 
      csetemp12*csetemp9 + 2*csetemp13*ToReal(M))) + (csetemp17 + 
      csetemp12*csetemp9)*(csetemp24 + 2*csetemp12*csetemp14*csetemp27 + 
      csetemp12*csetemp32 + 2*csetemp28*csetemp9 + 
      csetemp17*ToReal(a)*(csetemp22 + csetemp29 + csetemp12*ToReal(a)) + 
      2*csetemp10*csetemp25*ToReal(M) + 
      2*csetemp11*csetemp13*csetemp9*ToReal(M))*(-2*csetemp31 + 
      csetemp26*(3*csetemp10 + csetemp11 + csetemp12 - 3*csetemp9) + 
      csetemp12*csetemp27*(-3*csetemp10 + csetemp11 + csetemp12 - 
      csetemp9)*rXYZ - 3*csetemp12*csetemp30*X*Y + 
      csetemp12*(csetemp13*(-csetemp10 + csetemp11 + csetemp12 - 
      3*csetemp9)*csetemp9 - csetemp14*csetemp29*X*Y) + 
      csetemp25*ToReal(a)*(csetemp22 - csetemp29 + (csetemp10 + csetemp11 - 
      csetemp12)*ToReal(a)) - 4*csetemp24*ToReal(M) + 
      csetemp28*(3*X*Y*ToReal(a) + 2*(2*csetemp10 + csetemp11 + 
      csetemp12)*ToReal(M) - 6*csetemp9*ToReal(M)) + 
      csetemp17*csetemp9*(X*Y*ToReal(a) + 2*(csetemp11 + csetemp12)*ToReal(M) 
      - 2*csetemp9*ToReal(M))) + 2*csetemp13*(csetemp15 + 
      csetemp16)*(csetemp19 + csetemp20)*(csetemp17 + 
      csetemp12*csetemp9)*ToReal(M)*(csetemp12*(-(csetemp14*csetemp29*(csetemp11 
      + csetemp12 - 3*csetemp9)) - csetemp30*(-2*csetemp10 + csetemp11 + 
      csetemp12 - csetemp9)) - 4*csetemp12*csetemp27*rXYZ*X*Y + 
      2*csetemp24*ToReal(a) + csetemp28*(csetemp23 + 3*csetemp29 - 
      (4*csetemp10 + csetemp11 + csetemp12)*ToReal(a)) - 
      2*csetemp25*(3*csetemp10 + csetemp11 + csetemp12 - 
      3*csetemp9)*ToReal(a)*ToReal(M) + csetemp17*csetemp9*(csetemp29 + 
      (-2*csetemp10 - csetemp11 + csetemp12)*ToReal(a) - 2*X*Y*ToReal(M)) + 
      2*csetemp26*(X*Y + 2*ToReal(a)*ToReal(M)) + 
      2*csetemp13*csetemp9*(-(csetemp12*X*Y) + (csetemp29 - (csetemp10 + 
      csetemp11 + csetemp12)*ToReal(a))*ToReal(M))));
    
    CCTK_REAL csetemp33 = pow(rXYZ,19);
    
    CCTK_REAL csetemp34 = pow(rXYZ,18);
    
    CCTK_REAL csetemp35 = pow(rXYZ,17);
    
    CCTK_REAL csetemp36 = SQR(ToReal(M));
    
    CCTK_REAL csetemp37 = pow(ToReal(a),11);
    
    CCTK_REAL csetemp38 = pow(Z,6);
    
    CCTK_REAL csetemp39 = pow(ToReal(a),10);
    
    CCTK_REAL csetemp40 = QAD(Z);
    
    CCTK_REAL csetemp41 = pow(rXYZ,16);
    
    CCTK_REAL csetemp42 = pow(rXYZ,15);
    
    CCTK_REAL csetemp43 = pow(rXYZ,14);
    
    CCTK_REAL csetemp44 = pow(rXYZ,10);
    
    CCTK_REAL csetemp45 = pow(rXYZ,13);
    
    CCTK_REAL csetemp46 = pow(rXYZ,12);
    
    CCTK_REAL csetemp47 = pow(ToReal(a),9);
    
    CCTK_REAL csetemp48 = pow(ToReal(a),8);
    
    CCTK_REAL csetemp49 = pow(ToReal(a),7);
    
    CCTK_REAL csetemp50 = pow(rXYZ,11);
    
    CCTK_REAL K21 = -(csetemp13*INV(QAD(csetemp14 + 
      csetemp9)*SQR(csetemp17 + csetemp12*csetemp9)*SQR(csetemp17 + 
      csetemp12*csetemp9 + 2*csetemp13*ToReal(M))*sqrt(4*csetemp12*csetemp9 + 
      SQR(csetemp10 + csetemp11 + csetemp12 - csetemp9)))*sqrt(1 + 
      2*csetemp13*INV(csetemp17 + 
      csetemp12*csetemp9)*ToReal(M))*ToReal(M)*(-3*csetemp34*(csetemp22 + 
      (csetemp10 - csetemp11)*ToReal(a)) + 
      3*csetemp38*pow(ToReal(a),12)*(csetemp23 + (csetemp10 - 
      csetemp11)*ToReal(a)) + csetemp41*((csetemp10 - 
      csetemp11)*(-7*csetemp29 + 32*csetemp36*ToReal(a)) + (22*(csetemp10 + 
      csetemp11 + csetemp12) - 4*csetemp9)*X*Y*ToReal(M)) + 
      csetemp43*ToReal(a)*((csetemp10 - csetemp11)*(-5*csetemp27 - 
      36*(csetemp10 + csetemp11 + csetemp12)*csetemp36 + (-5*csetemp12 + 
      48*csetemp36)*csetemp9) + X*Y*(12*csetemp29 + 
      8*csetemp12*ToReal(a))*ToReal(M)) + 
      csetemp40*(csetemp17*csetemp48*((csetemp10 - csetemp11)*(5*csetemp29 + 
      (5*csetemp12 + 12*csetemp36)*ToReal(a)) + (-2*(5*(csetemp10 + 
      csetemp11) + 24*csetemp12) - 4*csetemp9)*X*Y*ToReal(M)) + 
      csetemp14*csetemp39*(7*(csetemp10 - csetemp11)*csetemp12*ToReal(a) + 
      4*(csetemp10 + csetemp11 - csetemp12 - csetemp9)*X*Y*ToReal(M))) + 
      4*(csetemp33*X*Y + csetemp37*csetemp38*rXYZ*(-2*X*Y*ToReal(a) + 
      (3*csetemp10 - 3*csetemp11)*ToReal(M)) + csetemp35*(2*(-4*csetemp36 + 
      csetemp9)*X*Y + 3*(csetemp10 - csetemp11)*ToReal(a)*ToReal(M)) + 
      csetemp42*((csetemp27 + 9*(csetemp10 + csetemp11 + csetemp12)*csetemp36 
      + (csetemp12 - 4*csetemp36)*csetemp9)*X*Y + (csetemp10 - 
      csetemp11)*(4*csetemp29 - 5*(csetemp10 + csetemp11 + 
      csetemp12)*ToReal(a))*ToReal(M))) - 
      2*(csetemp12*csetemp25*csetemp49*(X*Y*(4*csetemp29*(2*csetemp12 + 
      csetemp36) - 4*((csetemp10 + csetemp11 + 2*csetemp12)*csetemp36 - 
      2*csetemp40)*ToReal(a)) + (csetemp10 - 
      csetemp11)*csetemp12*(-9*(csetemp10 + csetemp11) - 10*csetemp12 + 
      2*csetemp9)*ToReal(M)) + csetemp13*csetemp40*csetemp47*(-((csetemp10 - 
      csetemp11)*(5*(csetemp10 + csetemp11) + 16*csetemp12)*ToReal(M)) + 
      2*((5*csetemp12 - 3*csetemp36)*X*Y*ToReal(a) + (csetemp10 - 
      csetemp11)*csetemp9*ToReal(M)))) + 
      csetemp12*(-(csetemp28*csetemp32*((csetemp10 - 
      csetemp11)*(csetemp29*(-9*csetemp12 + 8*csetemp36) - 
      (8*csetemp10*csetemp36 + 8*csetemp11*csetemp36 + 28*csetemp12*csetemp36 
      + csetemp40)*ToReal(a)) + 4*csetemp12*(9*csetemp10 + 9*csetemp11 + 
      11*csetemp12)*X*Y*ToReal(M) + 4*csetemp27*X*Y*ToReal(M) - 4*(csetemp10 
      + csetemp11 - csetemp12)*csetemp9*X*Y*ToReal(M))) + 
      csetemp24*csetemp27*((csetemp10 - csetemp11)*csetemp30 + (csetemp10 - 
      csetemp11)*csetemp29*(3*csetemp12 - 20*csetemp36) + 4*(csetemp10 - 
      csetemp11)*(3*csetemp10 + 3*csetemp11 + csetemp12)*csetemp36*ToReal(a) 
      - 6*csetemp12*(csetemp10 + csetemp11 + csetemp12)*X*Y*ToReal(M) - 
      6*csetemp27*X*Y*ToReal(M) - 4*(4*csetemp10 + 4*csetemp11 + 
      11*csetemp12)*csetemp9*X*Y*ToReal(M)) - 
      csetemp29*csetemp44*(3*(csetemp10 - csetemp11)*csetemp27 + 
      12*(csetemp10 - csetemp11)*(csetemp10 + csetemp11 + 
      csetemp12)*csetemp36 + (csetemp10 - csetemp11)*(csetemp12 + 
      20*csetemp36)*csetemp9 - 12*csetemp29*X*Y*ToReal(M) + 12*(3*csetemp10 + 
      3*csetemp11 + 2*csetemp12)*X*Y*ToReal(a)*ToReal(M)) - 
      4*(csetemp26*csetemp32*((-((csetemp10 + csetemp11 - 
      6*csetemp12)*csetemp36) + csetemp40)*X*Y + 3*(3*csetemp12 + 
      csetemp36)*csetemp9*X*Y + (csetemp10 - csetemp11)*csetemp29*ToReal(M) - 
      2*(csetemp10 - csetemp11)*(csetemp10 + csetemp11 + 
      csetemp12)*ToReal(a)*ToReal(M)) + csetemp27*csetemp31*(2*csetemp27*X*Y 
      + 2*(3*csetemp10 + 3*csetemp11 + 2*csetemp12)*csetemp36*X*Y + 
      6*csetemp12*csetemp9*X*Y + (csetemp10 - csetemp11)*(3*csetemp29 - 
      (csetemp10 + csetemp11 - 3*csetemp12)*ToReal(a))*ToReal(M)))) + 
      csetemp9*(csetemp46*((-csetemp10 + csetemp11)*csetemp30 + (csetemp10 - 
      csetemp11)*(csetemp29*(-9*csetemp12 + 16*csetemp36) - 4*(5*(csetemp10 + 
      csetemp11) + 11*csetemp12)*csetemp36*ToReal(a)) + 
      (16*csetemp12*(csetemp10 + csetemp11 + csetemp12) + 4*csetemp27 - 
      6*(csetemp10 + csetemp11)*csetemp9)*X*Y*ToReal(M)) + 
      2*(csetemp45*ToReal(M)*((csetemp10 - csetemp11)*(2*csetemp29 - 
      (7*(csetemp10 + csetemp11) + 10*csetemp12)*ToReal(a)) + 
      (4*(-2*(csetemp10 + csetemp11) + csetemp12) + 
      16*csetemp9)*X*Y*ToReal(M)) + csetemp50*((6*csetemp12*(csetemp10 + 
      csetemp11 + csetemp12)*csetemp36 + 2*csetemp27*(-3*csetemp12 + 
      4*csetemp36) - 2*((5*(csetemp10 + csetemp11) + 6*csetemp12)*csetemp36 + 
      csetemp40)*csetemp9)*X*Y + (csetemp10 - csetemp11)*(-((csetemp10 + 
      csetemp11 + 4*csetemp12)*csetemp29) - 10*csetemp12*(csetemp10 + 
      csetemp11 + csetemp12)*ToReal(a))*ToReal(M))))));
    
    CCTK_REAL csetemp51 = 4*X*ToReal(M);
    
    CCTK_REAL csetemp52 = pow(ToReal(a),13);
    
    CCTK_REAL csetemp53 = QAD(X);
    
    CCTK_REAL csetemp54 = QAD(Y);
    
    CCTK_REAL K31 = rXYZ*Z*INV(CUB(csetemp14 + csetemp9)*SQR(csetemp17 + 
      csetemp12*csetemp9)*SQR(csetemp17 + csetemp12*csetemp9 + 
      2*csetemp13*ToReal(M))*sqrt(4*csetemp12*csetemp9 + SQR(csetemp10 + 
      csetemp11 + csetemp12 - csetemp9)))*sqrt(1 + 2*csetemp13*INV(csetemp17 
      + csetemp12*csetemp9)*ToReal(M))*ToReal(M)*(-4*csetemp33*X + 
      3*csetemp34*(csetemp51 - Y*ToReal(a)) + csetemp38*(3*csetemp52*Y + 
      2*csetemp37*rXYZ*(csetemp21 + 3*X*ToReal(a))) - 
      csetemp41*(11*csetemp29*Y - 32*csetemp36*Y*ToReal(a) + 22*(csetemp10 + 
      csetemp11 + csetemp12)*X*ToReal(M) - 24*csetemp9*X*ToReal(M)) + 
      csetemp12*csetemp24*csetemp27*(-15*csetemp30*Y + csetemp29*(8*csetemp10 
      + 8*csetemp11 + 7*csetemp12 - 4*csetemp36)*Y - 4*(csetemp10 + csetemp11 
      + 3*csetemp12)*csetemp36*Y*ToReal(a) + 2*(4*csetemp10*csetemp11 + 
      7*csetemp10*csetemp12 + 7*csetemp11*csetemp12 + 5*csetemp40 + 
      2*csetemp53 + 2*csetemp54)*X*ToReal(M) - 4*(csetemp10 + csetemp11 + 
      3*csetemp12)*csetemp9*X*ToReal(M)) + 
      csetemp12*csetemp28*csetemp32*((4*csetemp10 + 4*csetemp11 + 
      3*csetemp12)*csetemp29*Y - 4*csetemp30*Y + csetemp12*(2*csetemp10 + 
      2*csetemp11 + 3*csetemp12 + 4*csetemp36)*Y*ToReal(a) + 
      2*(4*csetemp10*csetemp11 + 17*csetemp10*csetemp12 + 
      17*csetemp11*csetemp12 + 19*csetemp40 + 2*csetemp53 + 
      2*csetemp54)*X*ToReal(M) - 4*(csetemp10 + csetemp11 + 
      3*csetemp12)*csetemp9*X*ToReal(M)) + 
      2*csetemp13*csetemp40*csetemp47*(X*(-csetemp29 + (csetemp10 + csetemp11 
      + 9*csetemp12)*ToReal(a)) + (5*csetemp10 + 5*csetemp11 + 12*csetemp12 - 
      2*csetemp9)*Y*ToReal(M)) + 
      2*csetemp12*csetemp25*csetemp49*(csetemp12*csetemp29*X + 
      csetemp12*(2*csetemp10 + 2*csetemp11 + 9*csetemp12 + 
      2*csetemp36)*X*ToReal(a) + (11*csetemp11*csetemp12 + 
      csetemp10*(4*csetemp11 + 11*csetemp12) + 2*(5*csetemp40 + csetemp53 + 
      csetemp54))*Y*ToReal(M) - 2*(csetemp10 + csetemp11 + 
      4*csetemp12)*csetemp9*Y*ToReal(M)) + 
      2*csetemp12*csetemp26*csetemp30*((2*csetemp10 + 2*csetemp11 + 
      9*csetemp12)*csetemp29*X - 2*csetemp30*X + csetemp12*(csetemp10 + 
      csetemp11 + 3*csetemp12 + 2*csetemp36)*X*ToReal(a) - 4*(csetemp10 + 
      csetemp11 + 3*csetemp12)*csetemp9*Y*ToReal(M) + 2*Y*SQR(csetemp10 + 
      csetemp11 + csetemp12)*ToReal(M)) - csetemp43*ToReal(a)*(36*(csetemp10 
      + csetemp11 + csetemp12)*csetemp36*Y + (15*csetemp27 - (2*csetemp10 + 
      2*csetemp11 - 3*csetemp12 + 64*csetemp36)*csetemp9)*Y - 
      16*csetemp29*X*ToReal(M) + 2*(17*csetemp10 + 17*csetemp11 + 
      21*csetemp12)*X*ToReal(a)*ToReal(M)) - 
      csetemp29*csetemp44*(12*csetemp12*(csetemp10 + csetemp11 + 
      csetemp12)*csetemp36*Y + (2*csetemp32 - csetemp27*(2*csetemp10 + 
      2*csetemp11 - 21*csetemp12 + 8*csetemp36))*Y + 
      (-(csetemp12*(4*csetemp10 + 4*csetemp11 + 3*csetemp12)) + 
      4*(2*csetemp10 + 2*csetemp11 + 7*csetemp12)*csetemp36)*csetemp9*Y + 
      2*(2*csetemp10 + 2*csetemp11 + 5*csetemp12)*csetemp29*X*ToReal(M) - 
      4*(2*csetemp10*csetemp11 - 4*csetemp40 + csetemp53 + 
      csetemp54)*X*ToReal(a)*ToReal(M)) + 2*csetemp35*(16*csetemp36*X - 
      7*csetemp9*X + 6*Y*ToReal(a)*ToReal(M)) - 2*csetemp42*(18*(csetemp10 + 
      csetemp11 + csetemp12)*csetemp36*X + (9*csetemp27 - (csetemp10 + 
      csetemp11 - csetemp12 + 32*csetemp36)*csetemp9)*X - 
      12*csetemp29*Y*ToReal(M) + 10*(csetemp10 + csetemp11 + 
      csetemp12)*Y*ToReal(a)*ToReal(M)) - 
      2*csetemp27*csetemp31*(2*csetemp12*(csetemp10 + csetemp11 + 
      3*csetemp12)*csetemp36*X + csetemp12*(7*csetemp27 - (4*csetemp10 + 
      4*csetemp11 + 11*csetemp12 - 2*csetemp36)*csetemp9)*X + 2*(csetemp10 + 
      csetemp11 + 2*csetemp12)*csetemp29*Y*ToReal(M) - 
      2*(2*csetemp10*(csetemp11 - csetemp12) - 2*csetemp11*csetemp12 - 
      7*csetemp40 + csetemp53 + csetemp54)*Y*ToReal(a)*ToReal(M)) - 
      2*csetemp50*csetemp9*(6*csetemp12*(csetemp10 + csetemp11 + 
      csetemp12)*csetemp36*X + (csetemp32 - csetemp27*(csetemp10 + csetemp11 
      - 9*csetemp12 + 4*csetemp36))*X + 2*(-(csetemp12*(csetemp10 + csetemp11 
      + 2*csetemp12)) + (2*csetemp10 + 2*csetemp11 + 
      7*csetemp12)*csetemp36)*csetemp9*X + 3*(3*csetemp10 + 3*csetemp11 + 
      4*csetemp12)*csetemp29*Y*ToReal(M) - 2*csetemp30*Y*ToReal(M) - 
      2*(2*csetemp10*csetemp11 - 3*csetemp10*csetemp12 - 
      3*csetemp11*csetemp12 - 4*csetemp40 + csetemp53 + 
      csetemp54)*Y*ToReal(a)*ToReal(M)) + 
      csetemp40*(csetemp14*csetemp39*(-2*csetemp29*Y + (2*csetemp10 + 
      2*csetemp11 + 9*csetemp12)*Y*ToReal(a) + 14*csetemp12*X*ToReal(M)) + 
      csetemp17*csetemp48*(-3*csetemp29*Y + (4*csetemp10 + 4*csetemp11 + 
      9*csetemp12 + 4*csetemp36)*Y*ToReal(a) + 2*(8*csetemp10 + 8*csetemp11 + 
      21*csetemp12)*X*ToReal(M) - 4*csetemp9*X*ToReal(M))) + 
      csetemp9*(-(csetemp46*((9*csetemp30 - csetemp29*(4*csetemp10 + 
      4*csetemp11 - 13*csetemp12 + 40*csetemp36))*Y + 12*(3*csetemp10 + 
      3*csetemp11 + 5*csetemp12)*csetemp36*Y*ToReal(a) - 
      4*csetemp27*X*ToReal(M) - 4*(2*csetemp10*(csetemp11 - csetemp12) - 
      2*csetemp11*csetemp12 - 3*csetemp40 + csetemp53 + 
      csetemp54)*X*ToReal(M) + 10*(2*csetemp10 + 2*csetemp11 + 
      3*csetemp12)*csetemp9*X*ToReal(M))) - 2*csetemp45*(6*(3*csetemp10 + 
      3*csetemp11 + 5*csetemp12)*csetemp36*X + (5*csetemp27 - (2*csetemp10 + 
      2*csetemp11 - 5*csetemp12 + 20*csetemp36)*csetemp9)*X - 
      8*csetemp29*Y*ToReal(M) + 3*(5*csetemp10 + 5*csetemp11 + 
      6*csetemp12)*Y*ToReal(a)*ToReal(M))));
    
    CCTK_REAL csetemp55 = -csetemp20;
    
    CCTK_REAL K22 = -2*csetemp13*INV(QAD(csetemp14 + 
      csetemp9)*SQR(csetemp17 + csetemp12*csetemp9)*SQR(csetemp17 + 
      csetemp12*csetemp9 + 2*csetemp13*ToReal(M))*sqrt(4*csetemp12*csetemp9 + 
      SQR(csetemp10 + csetemp11 + csetemp12 - csetemp9)))*sqrt(1 + 
      2*csetemp13*INV(csetemp17 + 
      csetemp12*csetemp9)*ToReal(M))*ToReal(M)*((-(csetemp12*(csetemp18 + 
      csetemp55)*SQR(csetemp14 + csetemp9)*(-5*csetemp17*csetemp29*X + 
      3*csetemp12*csetemp30*X + csetemp14*csetemp29*(2*csetemp10 + 
      2*csetemp11 + csetemp12 - 2*csetemp9)*X + 3*csetemp26*Y + 
      3*csetemp25*csetemp9*Y + csetemp13*csetemp9*(-2*csetemp10 - 2*csetemp11 
      - 3*csetemp12 + 2*csetemp9)*Y - 5*csetemp12*csetemp27*rXYZ*Y - 
      5*csetemp28*X*ToReal(a))) - csetemp13*SQR(csetemp18 + 
      csetemp55)*(4*csetemp31 - csetemp26*(2*csetemp10 + 5*csetemp11 + 
      2*csetemp12 - 6*csetemp9) + 
      csetemp9*(-(csetemp12*csetemp13*(2*csetemp10 + csetemp11 + 2*csetemp12 
      - 6*csetemp9)) + csetemp25*(-2*csetemp10 - csetemp11 + 2*(csetemp12 + 
      csetemp9))) + csetemp12*csetemp27*(-2*csetemp10 + 3*csetemp11 - 
      2*csetemp12 + 2*csetemp9)*rXYZ + csetemp12*csetemp14*csetemp29*X*Y + 
      csetemp17*csetemp29*X*Y - 3*csetemp12*csetemp30*X*Y + 
      5*csetemp28*X*Y*ToReal(a)))*ToReal(M)*(csetemp17 + csetemp12*csetemp9 + 
      2*csetemp13*ToReal(M)) + (csetemp19 + 
      csetemp20)*(4*csetemp12*(csetemp17 + 
      csetemp12*csetemp9)*rXYZ*Y*CUB(csetemp14 + 
      csetemp9)*ToReal(M)*(csetemp17 - csetemp12*csetemp9 + 
      csetemp13*ToReal(M)) + csetemp13*(csetemp15 + 
      csetemp16)*(csetemp12*csetemp14*csetemp29*(-3*csetemp10 - 2*(csetemp11 
      + csetemp12) + 6*csetemp9) + csetemp17*csetemp29*(-3*csetemp10 + 
      2*(-csetemp11 + csetemp12 + csetemp9)) + (3*csetemp26 - 
      csetemp25*csetemp9)*X*Y - 5*csetemp12*csetemp27*rXYZ*X*Y + 
      csetemp12*(csetemp30*(csetemp10 - 2*(csetemp11 + csetemp12) + 
      2*csetemp9) - csetemp13*csetemp9*X*Y) + 4*csetemp24*ToReal(a) + 
      csetemp28*(-7*csetemp10 - 2*(csetemp11 + csetemp12) + 
      6*csetemp9)*ToReal(a))*ToReal(M)*(csetemp17 + csetemp12*csetemp9 + 
      2*csetemp13*ToReal(M))) + (csetemp17 + 
      csetemp12*csetemp9)*(-2*csetemp31 + csetemp26*(csetemp10 + 3*csetemp11 
      + csetemp12 - 3*csetemp9) + csetemp12*csetemp13*(csetemp10 - csetemp11 
      + csetemp12 - 3*csetemp9)*csetemp9 + csetemp12*csetemp27*(csetemp10 - 
      3*csetemp11 + csetemp12 - csetemp9)*rXYZ + 
      csetemp12*csetemp14*csetemp29*X*Y + 3*csetemp12*csetemp30*X*Y - 
      csetemp25*ToReal(a)*(csetemp22 + csetemp29 - (csetemp10 + csetemp11 - 
      csetemp12)*ToReal(a)) - 4*csetemp24*ToReal(M) + 
      csetemp28*(-3*X*Y*ToReal(a) + 2*(csetemp10 + 2*csetemp11 + 
      csetemp12)*ToReal(M) - 6*csetemp9*ToReal(M)) + 
      csetemp17*csetemp9*(-(X*Y*ToReal(a)) + 2*(csetemp10 + 
      csetemp12)*ToReal(M) - 2*csetemp9*ToReal(M)))*(csetemp24 + 
      2*csetemp12*csetemp14*csetemp27 + csetemp12*csetemp32 + 
      2*csetemp28*csetemp9 + 2*csetemp11*csetemp25*ToReal(M) + 
      2*csetemp10*csetemp13*csetemp9*ToReal(M) + 
      csetemp17*ToReal(a)*(csetemp29 + csetemp12*ToReal(a) - 
      4*X*Y*ToReal(M))) + 2*csetemp13*(csetemp15 + csetemp16)*(csetemp19 + 
      csetemp20)*(csetemp17 + 
      csetemp12*csetemp9)*ToReal(M)*(csetemp12*csetemp14*csetemp29*(csetemp10 
      + csetemp12 - 3*csetemp9) + csetemp12*csetemp30*(csetemp10 - 
      2*csetemp11 + csetemp12 - csetemp9) - 4*csetemp12*csetemp27*rXYZ*X*Y - 
      2*csetemp24*ToReal(a) + csetemp28*(csetemp23 - 3*csetemp29 + (csetemp10 
      + 4*csetemp11 + csetemp12)*ToReal(a)) + 2*csetemp25*(csetemp10 + 
      3*csetemp11 + csetemp12 - 3*csetemp9)*ToReal(a)*ToReal(M) + 
      csetemp26*(2*X*Y - 4*ToReal(a)*ToReal(M)) + 
      csetemp9*(-(csetemp17*(csetemp23 + csetemp29 + (-csetemp10 - 
      2*csetemp11 + csetemp12)*ToReal(a))) + 2*csetemp13*(-(csetemp12*X*Y) + 
      (-csetemp29 + (csetemp10 + csetemp11 + 
      csetemp12)*ToReal(a))*ToReal(M)))));
    
    CCTK_REAL K32 = -(rXYZ*Z*INV(CUB(csetemp14 + csetemp9)*SQR(csetemp17 + 
      csetemp12*csetemp9)*SQR(csetemp17 + csetemp12*csetemp9 + 
      2*csetemp13*ToReal(M))*sqrt(4*csetemp12*csetemp9 + SQR(csetemp10 + 
      csetemp11 + csetemp12 - csetemp9)))*sqrt(1 + 2*csetemp13*INV(csetemp17 
      + csetemp12*csetemp9)*ToReal(M))*ToReal(M)*(-3*(csetemp18 + 
      csetemp21)*csetemp34 + 4*csetemp33*Y + csetemp38*(3*csetemp52*X + 
      2*csetemp37*rXYZ*(csetemp51 - 3*Y*ToReal(a))) + 
      csetemp41*(X*(-11*csetemp29 + 32*csetemp36*ToReal(a)) + (22*(csetemp10 
      + csetemp11 + csetemp12) - 24*csetemp9)*Y*ToReal(M)) + 
      csetemp43*ToReal(a)*((-15*csetemp27 - 36*(csetemp10 + csetemp11 + 
      csetemp12)*csetemp36 + (2*(csetemp10 + csetemp11) - 3*csetemp12 + 
      64*csetemp36)*csetemp9)*X + Y*(-16*csetemp29 + 2*(17*(csetemp10 + 
      csetemp11) + 21*csetemp12)*ToReal(a))*ToReal(M)) + 
      csetemp29*csetemp44*((-2*csetemp32 - 12*csetemp12*(csetemp10 + 
      csetemp11 + csetemp12)*csetemp36 + csetemp27*(2*(csetemp10 + csetemp11) 
      - 21*csetemp12 + 8*csetemp36) + (csetemp12*(4*(csetemp10 + csetemp11) + 
      3*csetemp12) - 4*(2*(csetemp10 + csetemp11) + 
      7*csetemp12)*csetemp36)*csetemp9)*X + Y*(2*(2*(csetemp10 + csetemp11) + 
      5*csetemp12)*csetemp29 - 4*(2*csetemp10*csetemp11 - 4*csetemp40 + 
      csetemp53 + csetemp54)*ToReal(a))*ToReal(M)) + 
      csetemp12*(-2*csetemp25*csetemp49*(csetemp12*Y*(csetemp29 + 
      (9*csetemp12 + 2*(csetemp10 + csetemp11 + csetemp36))*ToReal(a)) + 
      (-11*csetemp11*csetemp12 - csetemp10*(4*csetemp11 + 11*csetemp12) - 
      2*(5*csetemp40 + csetemp53 + csetemp54) + 2*(csetemp10 + csetemp11 + 
      4*csetemp12)*csetemp9)*X*ToReal(M)) - 
      csetemp24*csetemp27*(X*(15*csetemp30 - csetemp29*(8*(csetemp10 + 
      csetemp11) + 7*csetemp12 - 4*csetemp36) + 4*(csetemp10 + csetemp11 + 
      3*csetemp12)*csetemp36*ToReal(a)) + (2*(7*csetemp11*csetemp12 + 
      csetemp10*(4*csetemp11 + 7*csetemp12) + 5*csetemp40 + 2*(csetemp53 + 
      csetemp54)) - 4*(csetemp10 + csetemp11 + 
      3*csetemp12)*csetemp9)*Y*ToReal(M)) + 
      csetemp28*csetemp32*(X*((4*(csetemp10 + csetemp11) + 
      3*csetemp12)*csetemp29 - 4*csetemp30 + csetemp12*(2*(csetemp10 + 
      csetemp11) + 3*csetemp12 + 4*csetemp36)*ToReal(a)) + 
      (-2*(17*csetemp11*csetemp12 + csetemp10*(4*csetemp11 + 17*csetemp12) + 
      19*csetemp40 + 2*(csetemp53 + csetemp54)) + 4*(csetemp10 + csetemp11 + 
      3*csetemp12)*csetemp9)*Y*ToReal(M)) + 
      2*csetemp26*csetemp30*(Y*(-((2*csetemp10 + 2*csetemp11 + 
      9*csetemp12)*csetemp29) + 2*csetemp30 - csetemp12*(csetemp10 + 
      csetemp11 + 3*csetemp12 + 2*csetemp36)*ToReal(a)) + X*(-4*(csetemp10 + 
      csetemp11 + 3*csetemp12)*csetemp9 + 2*SQR(csetemp10 + csetemp11 + 
      csetemp12))*ToReal(M))) + csetemp9*(csetemp46*(X*(-9*csetemp30 + 
      csetemp29*(4*(csetemp10 + csetemp11) - 13*csetemp12 + 40*csetemp36) - 
      12*(3*(csetemp10 + csetemp11) + 5*csetemp12)*csetemp36*ToReal(a)) + 
      (-4*(2*csetemp10*(csetemp11 - csetemp12) - 2*csetemp11*csetemp12 + 
      csetemp27 - 3*csetemp40 + csetemp53 + csetemp54) + 10*(2*(csetemp10 + 
      csetemp11) + 3*csetemp12)*csetemp9)*Y*ToReal(M)) + 
      2*csetemp45*((5*csetemp27 + 6*(3*(csetemp10 + csetemp11) + 
      5*csetemp12)*csetemp36 - (2*(csetemp10 + csetemp11) - 5*csetemp12 + 
      20*csetemp36)*csetemp9)*Y + X*(8*csetemp29 - 3*(5*(csetemp10 + 
      csetemp11) + 6*csetemp12)*ToReal(a))*ToReal(M))) + 
      2*(csetemp35*((-16*csetemp36 + 7*csetemp9)*Y + 6*X*ToReal(a)*ToReal(M)) 
      + csetemp42*((9*csetemp27 + 18*(csetemp10 + csetemp11 + 
      csetemp12)*csetemp36 - (csetemp10 + csetemp11 - csetemp12 + 
      32*csetemp36)*csetemp9)*Y + X*(12*csetemp29 - 10*(csetemp10 + csetemp11 
      + csetemp12)*ToReal(a))*ToReal(M)) + 
      csetemp27*csetemp31*(csetemp12*(7*csetemp27 + 2*(csetemp10 + csetemp11 
      + 3*csetemp12)*csetemp36 - (4*(csetemp10 + csetemp11) + 11*csetemp12 - 
      2*csetemp36)*csetemp9)*Y + X*(-2*(csetemp10 + csetemp11 + 
      2*csetemp12)*csetemp29 + 2*(2*csetemp10*(csetemp11 - csetemp12) - 
      2*csetemp11*csetemp12 - 7*csetemp40 + csetemp53 + 
      csetemp54)*ToReal(a))*ToReal(M)) + csetemp50*csetemp9*((csetemp32 + 
      6*csetemp12*(csetemp10 + csetemp11 + csetemp12)*csetemp36 - 
      csetemp27*(csetemp10 + csetemp11 - 9*csetemp12 + 4*csetemp36) + 
      2*(-(csetemp12*(csetemp10 + csetemp11 + 2*csetemp12)) + (2*(csetemp10 + 
      csetemp11) + 7*csetemp12)*csetemp36)*csetemp9)*Y + X*(-3*(3*(csetemp10 
      + csetemp11) + 4*csetemp12)*csetemp29 + 2*(csetemp30 + 
      (csetemp10*(2*csetemp11 - 3*csetemp12) - 3*csetemp11*csetemp12 - 
      4*csetemp40 + csetemp53 + csetemp54)*ToReal(a)))*ToReal(M))) + 
      csetemp40*(2*csetemp13*csetemp47*(Y*(csetemp29 - (csetemp10 + csetemp11 
      + 9*csetemp12)*ToReal(a)) + (5*(csetemp10 + csetemp11) + 12*csetemp12 - 
      2*csetemp9)*X*ToReal(M)) + csetemp17*csetemp48*(X*(-3*csetemp29 + 
      (9*csetemp12 + 4*(csetemp10 + csetemp11 + csetemp36))*ToReal(a)) - 
      2*(8*(csetemp10 + csetemp11) + 21*csetemp12 - 2*csetemp9)*Y*ToReal(M)) 
      + csetemp14*csetemp39*((2*(csetemp10 + csetemp11) + 
      9*csetemp12)*X*ToReal(a) - 2*(csetemp29*X + 
      7*csetemp12*Y*ToReal(M))))));
    
    CCTK_REAL csetemp56 = pow(Z,8);
    
    CCTK_REAL K33 = -2*INV((csetemp14 + csetemp9)*SQR(csetemp17 + 
      csetemp12*csetemp9)*SQR(csetemp17 + csetemp12*csetemp9 + 
      2*csetemp13*ToReal(M))*sqrt(4*csetemp12*csetemp9 + SQR(csetemp10 + 
      csetemp11 + csetemp12 - csetemp9)))*sqrt(1 + 2*csetemp13*INV(csetemp17 
      + csetemp12*csetemp9)*ToReal(M))*ToReal(M)*(-2*(csetemp34 + 
      csetemp39*csetemp56) + csetemp38*(-(csetemp14*csetemp48*(csetemp10 + 
      csetemp11 + 5*csetemp12 - csetemp9)) - csetemp17*csetemp32*(csetemp10 + 
      csetemp11 + 3*csetemp12 - 2*csetemp36 - csetemp9)) + 
      csetemp41*(csetemp10 + csetemp11 + 3*csetemp12 - 3*csetemp9) + 
      csetemp43*(-csetemp27 - 16*csetemp12*csetemp36 + (csetemp10 + csetemp11 
      + 3*csetemp12)*csetemp9) + csetemp12*csetemp46*(-csetemp27 + 
      18*(csetemp10 + csetemp11 + csetemp12)*csetemp36 + (csetemp10 + 
      csetemp11 + 3*(csetemp12 - 8*csetemp36))*csetemp9) + 
      csetemp12*csetemp44*csetemp9*(-csetemp27 + 4*(2*csetemp10 + 2*csetemp11 
      + 5*csetemp12)*csetemp36 + (csetemp10 + csetemp11 + 7*csetemp12 - 
      8*csetemp36)*csetemp9) + csetemp27*csetemp28*csetemp40*(csetemp27 + 
      4*(csetemp10 + csetemp11 + 2*csetemp12)*csetemp36 - (csetemp10 + 
      csetemp11 + 3*csetemp12 + 4*csetemp36)*csetemp9) - 
      4*csetemp35*ToReal(M) + 2*csetemp42*(csetemp10 + csetemp11 - 
      2*csetemp12 - 3*csetemp9)*ToReal(M) + 
      2*csetemp26*csetemp27*csetemp40*(4*csetemp10 + 4*csetemp11 + 
      7*csetemp12 - 2*csetemp9)*ToReal(M) + 
      2*csetemp12*csetemp50*(9*csetemp10 + 9*csetemp11 + 11*csetemp12 - 
      7*csetemp9)*csetemp9*ToReal(M) + 
      4*csetemp13*csetemp32*csetemp38*(-2*csetemp10 - 2*csetemp11 - 
      3*csetemp12 + csetemp9)*ToReal(M) + csetemp45*(11*csetemp12*(csetemp10 
      + csetemp11 + csetemp12) - 2*csetemp27 + 2*(csetemp10 + csetemp11 - 
      7*csetemp12)*csetemp9)*ToReal(M) + 
      csetemp12*csetemp31*csetemp9*(-4*csetemp27 - 4*(2*csetemp10*csetemp11 - 
      csetemp40 + csetemp53 + csetemp54) + (8*csetemp10 + 8*csetemp11 + 
      9*csetemp12)*csetemp9)*ToReal(M) - 5*csetemp48*csetemp56*rXYZ*ToReal(M) 
      + csetemp40*(csetemp24*csetemp9*(5*csetemp27 + 6*(csetemp10 + csetemp11 
      + csetemp12)*csetemp36 - (csetemp10 + csetemp11 + 3*csetemp12 + 
      2*csetemp36)*csetemp9) - csetemp25*csetemp27*(11*csetemp11*csetemp12 + 
      csetemp10*(8*csetemp11 + 11*csetemp12) + 2*csetemp27 + 7*csetemp40 + 
      4*csetemp53 + 4*csetemp54 - 2*(3*csetemp10 + 3*csetemp11 + 
      7*csetemp12)*csetemp9)*ToReal(M)));
    
    CCTK_REAL alpp = INV(sqrt(1 + 2*csetemp13*INV(csetemp17 + 
      csetemp12*csetemp9)*ToReal(M)));
    
    CCTK_REAL dtalpp = 0;
    
    CCTK_REAL betap1 = 2*csetemp13*(csetemp15 + csetemp16)*INV((csetemp14 
      + csetemp9)*(csetemp17 + csetemp12*csetemp9 + 
      2*csetemp13*ToReal(M)))*ToReal(M);
    
    CCTK_REAL betap2 = 2*csetemp13*(csetemp19 + csetemp20)*INV((csetemp14 
      + csetemp9)*(csetemp17 + csetemp12*csetemp9 + 
      2*csetemp13*ToReal(M)))*ToReal(M);
    
    CCTK_REAL betap3 = 2*csetemp14*Z*INV(csetemp17 + csetemp12*csetemp9 + 
      2*csetemp13*ToReal(M))*ToReal(M);
    
    CCTK_REAL dtbetap1 = 0;
    
    CCTK_REAL dtbetap2 = 0;
    
    CCTK_REAL dtbetap3 = 0;
    
    CCTK_REAL csetemp57 = SQR(Jac11);
    
    CCTK_REAL csetemp58 = SQR(Jac21);
    
    CCTK_REAL csetemp59 = SQR(Jac31);
    
    CCTK_REAL gxxL = csetemp57*G11 + csetemp58*G22 + csetemp59*G33 + 
      2*(G32*Jac21*Jac31 + Jac11*(G21*Jac21 + G31*Jac31));
    
    CCTK_REAL gxyL = Jac12*(G11*Jac11 + G21*Jac21 + G31*Jac31) + 
      Jac22*(G21*Jac11 + G22*Jac21 + G32*Jac31) + (G31*Jac11 + G32*Jac21 + 
      G33*Jac31)*Jac32;
    
    CCTK_REAL gxzL = Jac13*(G11*Jac11 + G21*Jac21 + G31*Jac31) + 
      Jac23*(G21*Jac11 + G22*Jac21 + G32*Jac31) + (G31*Jac11 + G32*Jac21 + 
      G33*Jac31)*Jac33;
    
    CCTK_REAL csetemp60 = SQR(Jac12);
    
    CCTK_REAL csetemp61 = SQR(Jac22);
    
    CCTK_REAL csetemp62 = SQR(Jac32);
    
    CCTK_REAL gyyL = csetemp60*G11 + csetemp61*G22 + csetemp62*G33 + 
      2*(G32*Jac22*Jac32 + Jac12*(G21*Jac22 + G31*Jac32));
    
    CCTK_REAL gyzL = Jac13*(G11*Jac12 + G21*Jac22 + G31*Jac32) + 
      Jac23*(G21*Jac12 + G22*Jac22 + G32*Jac32) + (G31*Jac12 + G32*Jac22 + 
      G33*Jac32)*Jac33;
    
    CCTK_REAL csetemp63 = SQR(Jac13);
    
    CCTK_REAL csetemp64 = SQR(Jac23);
    
    CCTK_REAL csetemp65 = SQR(Jac33);
    
    CCTK_REAL gzzL = csetemp63*G11 + csetemp64*G22 + csetemp65*G33 + 
      2*(G32*Jac23*Jac33 + Jac13*(G21*Jac23 + G31*Jac33));
    
    CCTK_REAL kxxL = csetemp57*K11 + csetemp58*K22 + 2*(Jac11*(Jac21*K21 + 
      Jac31*K31) + Jac21*Jac31*K32) + csetemp59*K33;
    
    CCTK_REAL kxyL = Jac11*(Jac12*K11 + Jac22*K21 + Jac32*K31) + 
      Jac21*(Jac12*K21 + Jac22*K22 + Jac32*K32) + Jac31*(Jac12*K31 + 
      Jac22*K32 + Jac32*K33);
    
    CCTK_REAL kxzL = Jac11*(Jac13*K11 + Jac23*K21 + Jac33*K31) + 
      Jac21*(Jac13*K21 + Jac23*K22 + Jac33*K32) + Jac31*(Jac13*K31 + 
      Jac23*K32 + Jac33*K33);
    
    CCTK_REAL kyyL = csetemp60*K11 + csetemp61*K22 + 2*(Jac12*(Jac22*K21 + 
      Jac32*K31) + Jac22*Jac32*K32) + csetemp62*K33;
    
    CCTK_REAL kyzL = Jac12*(Jac13*K11 + Jac23*K21 + Jac33*K31) + 
      Jac22*(Jac13*K21 + Jac23*K22 + Jac33*K32) + Jac32*(Jac13*K31 + 
      Jac23*K32 + Jac33*K33);
    
    CCTK_REAL kzzL = csetemp63*K11 + csetemp64*K22 + 2*(Jac13*(Jac23*K21 + 
      Jac33*K31) + Jac23*Jac33*K32) + csetemp65*K33;
    
    CCTK_REAL alpL = alpp*ToReal(lapsefactor);
    
    CCTK_REAL dtalpL = dtalpp*ToReal(lapsefactor);
    
    CCTK_REAL betaxL = betap1*InvJac11 + betap2*InvJac12 + betap3*InvJac13 
      + shiftadd1;
    
    CCTK_REAL betayL = betap1*InvJac21 + betap2*InvJac22 + betap3*InvJac23 
      + shiftadd2;
    
    CCTK_REAL betazL = betap1*InvJac31 + betap2*InvJac32 + betap3*InvJac33 
      + shiftadd3;
    
    CCTK_REAL dtbetaxL = dtbetap1*InvJac11 + dtbetap2*InvJac12 + 
      dtbetap3*InvJac13;
    
    CCTK_REAL dtbetayL = dtbetap1*InvJac21 + dtbetap2*InvJac22 + 
      dtbetap3*InvJac23;
    
    CCTK_REAL dtbetazL = dtbetap1*InvJac31 + dtbetap2*InvJac32 + 
      dtbetap3*InvJac33;
    
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
