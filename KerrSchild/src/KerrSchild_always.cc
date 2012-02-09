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
#define QAD(x) (SQR(SQR(x)))
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))

static void KerrSchild_always_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_LOOP3(KerrSchild_always,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2])
  {
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL = alp[index];
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
    
    CCTK_REAL Jac11 = cos(ToReal(phi))*cos(ToReal(psi)) - 
      cos(ToReal(theta))*sin(ToReal(phi))*sin(ToReal(psi));
    
    CCTK_REAL Jac12 = cos(ToReal(psi))*sin(ToReal(phi)) + 
      cos(ToReal(phi))*cos(ToReal(theta))*sin(ToReal(psi));
    
    CCTK_REAL Jac13 = sin(ToReal(psi))*sin(ToReal(theta));
    
    CCTK_REAL Jac21 = 
      -(cos(ToReal(psi))*cos(ToReal(theta))*sin(ToReal(phi))) - 
      cos(ToReal(phi))*sin(ToReal(psi));
    
    CCTK_REAL Jac22 = cos(ToReal(phi))*cos(ToReal(psi))*cos(ToReal(theta)) 
      - sin(ToReal(phi))*sin(ToReal(psi));
    
    CCTK_REAL Jac23 = cos(ToReal(psi))*sin(ToReal(theta));
    
    CCTK_REAL Jac31 = sin(ToReal(phi))*sin(ToReal(theta));
    
    CCTK_REAL Jac32 = -(cos(ToReal(phi))*sin(ToReal(theta)));
    
    CCTK_REAL Jac33 = cos(ToReal(theta));
    
    CCTK_REAL InvJac11 = Jac11;
    
    CCTK_REAL InvJac12 = Jac21;
    
    CCTK_REAL InvJac13 = Jac31;
    
    CCTK_REAL InvJac21 = Jac12;
    
    CCTK_REAL InvJac22 = Jac22;
    
    CCTK_REAL InvJac23 = Jac32;
    
    CCTK_REAL InvJac31 = Jac13;
    
    CCTK_REAL InvJac32 = Jac23;
    
    CCTK_REAL InvJac33 = Jac33;
    
    CCTK_REAL T = t - ToReal(positiont);
    
    CCTK_REAL XX1 = -(Jac11*(position1 + shiftadd1*T - xx1)) - 
      Jac12*(position2 + shiftadd2*T - xx2) - Jac13*(position3 + shiftadd3*T 
      - xx3);
    
    CCTK_REAL XX2 = -(Jac21*(position1 + shiftadd1*T - xx1)) - 
      Jac22*(position2 + shiftadd2*T - xx2) - Jac23*(position3 + shiftadd3*T 
      - xx3);
    
    CCTK_REAL XX3 = -(Jac31*(position1 + shiftadd1*T - xx1)) - 
      Jac32*(position2 + shiftadd2*T - xx2) - Jac33*(position3 + shiftadd3*T 
      - xx3);
    
    CCTK_REAL X = XX1;
    
    CCTK_REAL Y = XX2;
    
    CCTK_REAL Z = XX3;
    
    CCTK_REAL rXYZ = INV(sqrt(2))*sqrt(SQR(X) + SQR(Y) + SQR(Z) - 
      SQR(ToReal(a)) + sqrt(SQR(SQR(X) + SQR(Y) + SQR(Z) - SQR(ToReal(a))) + 
      4*SQR(Z*ToReal(a))));
    
    alpL = INV(sqrt(1 + 2*CUB(rXYZ)*INV(QAD(rXYZ) + 
      SQR(Z*ToReal(a)))*ToReal(M)));
    
    CCTK_REAL dtalpL = 0;
    
    CCTK_REAL G11 = 1 + 2*CUB(rXYZ)*INV(SQR(SQR(rXYZ) + 
      SQR(ToReal(a)))*(QAD(rXYZ) + SQR(Z*ToReal(a))))*SQR(rXYZ*X + 
      Y*ToReal(a))*ToReal(M);
    
    CCTK_REAL G21 = 2*CUB(rXYZ)*INV(SQR(SQR(rXYZ) + 
      SQR(ToReal(a)))*(QAD(rXYZ) + SQR(Z*ToReal(a))))*(rXYZ*Y - 
      X*ToReal(a))*(rXYZ*X + Y*ToReal(a))*ToReal(M);
    
    CCTK_REAL G31 = 2*Z*INV((SQR(rXYZ) + SQR(ToReal(a)))*(QAD(rXYZ) + 
      SQR(Z*ToReal(a))))*SQR(rXYZ)*(rXYZ*X + Y*ToReal(a))*ToReal(M);
    
    CCTK_REAL G22 = 1 + 2*CUB(rXYZ)*INV(SQR(SQR(rXYZ) + 
      SQR(ToReal(a)))*(QAD(rXYZ) + SQR(Z*ToReal(a))))*SQR(-(rXYZ*Y) + 
      X*ToReal(a))*ToReal(M);
    
    CCTK_REAL G32 = 2*Z*INV((SQR(rXYZ) + SQR(ToReal(a)))*(QAD(rXYZ) + 
      SQR(Z*ToReal(a))))*SQR(rXYZ)*(rXYZ*Y - X*ToReal(a))*ToReal(M);
    
    CCTK_REAL G33 = 1 + 2*rXYZ*INV(QAD(rXYZ) + 
      SQR(Z*ToReal(a)))*SQR(Z)*ToReal(M);
    
    CCTK_REAL K11 = -2*CUB(rXYZ)*INV(alpL*QAD(SQR(rXYZ) + 
      SQR(ToReal(a)))*SQR(QAD(rXYZ) + SQR(Z*ToReal(a)))*SQR(QAD(rXYZ) + 
      SQR(Z*ToReal(a)) + 2*CUB(rXYZ)*ToReal(M))*sqrt(SQR(SQR(X) + SQR(Y) + 
      SQR(Z) - SQR(ToReal(a))) + 
      4*SQR(Z*ToReal(a))))*ToReal(M)*(CUB(rXYZ)*(-((rXYZ*Y - 
      X*ToReal(a))*(rXYZ*X + Y*ToReal(a))*(-3*X*Y*pow(rXYZ,7) + 
      2*pow(ToReal(a),7)*SQR(Z) + 5*rXYZ*X*Y*QAD(ToReal(a))*SQR(Z) + 
      CUB(ToReal(a))*SQR(rXYZ)*(6*QAD(rXYZ) + SQR(rXYZ)*(-2*SQR(X) - 3*SQR(Y) 
      + 2*SQR(Z)) - SQR(Z)*(2*SQR(X) + 3*SQR(Y) + 2*SQR(Z))) + 
      pow(ToReal(a),5)*(2*QAD(rXYZ) + (-2*SQR(X) + SQR(Y) - 2*SQR(Z))*SQR(Z) 
      + 6*SQR(rXYZ*Z)) + X*Y*CUB(rXYZ)*(SQR(rXYZ) + SQR(Z))*SQR(ToReal(a)) + 
      pow(rXYZ,6)*(4*SQR(rXYZ) - 2*SQR(X) - 7*SQR(Y) - 2*SQR(Z))*ToReal(a))) 
      - SQR(rXYZ*X + Y*ToReal(a))*(3*X*Y*pow(ToReal(a),5)*SQR(Z) + 
      2*rXYZ*pow(ToReal(a),6)*SQR(Z) + pow(rXYZ,7)*(4*SQR(rXYZ) - 5*SQR(X) - 
      2*(SQR(Y) + SQR(Z))) + rXYZ*QAD(ToReal(a))*(2*QAD(rXYZ) + 
      SQR(Z)*(3*SQR(X) - 2*(SQR(Y) + SQR(Z))) + 6*SQR(rXYZ*Z)) + 
      CUB(rXYZ)*(6*QAD(rXYZ) - SQR(rXYZ)*(SQR(X) + 2*SQR(Y) - 2*SQR(Z)) - 
      SQR(Z)*(SQR(X) + 2*(SQR(Y) + SQR(Z))))*SQR(ToReal(a)) + 
      X*Y*(-(CUB(ToReal(a))*SQR(rXYZ)*(SQR(rXYZ) + SQR(Z))) - 
      5*pow(rXYZ,6)*ToReal(a))))*ToReal(M)*(QAD(rXYZ) + SQR(Z*ToReal(a)) + 
      2*CUB(rXYZ)*ToReal(M)) + (rXYZ*X + 
      Y*ToReal(a))*ToReal(M)*(4*rXYZ*X*CUB(SQR(rXYZ) + 
      SQR(ToReal(a)))*SQR(Z)*(QAD(rXYZ) + SQR(Z*ToReal(a)))*(QAD(rXYZ) - 
      SQR(Z*ToReal(a)) + CUB(rXYZ)*ToReal(M)) - SQR(Z*(SQR(rXYZ) + 
      SQR(ToReal(a))))*(-3*X*pow(rXYZ,7) + 
      Y*CUB(ToReal(a))*SQR(rXYZ)*(-5*SQR(rXYZ) + 2*SQR(X) + 2*SQR(Y) + 
      SQR(Z)) + QAD(ToReal(a))*(-2*X*CUB(rXYZ) + 5*rXYZ*X*SQR(Z)) + 
      pow(ToReal(a),5)*(-2*Y*SQR(rXYZ) + 3*Y*SQR(Z)) + 
      X*CUB(rXYZ)*(-3*SQR(rXYZ) + 2*SQR(X) + 2*SQR(Y) + 
      3*SQR(Z))*SQR(ToReal(a)) - 5*Y*pow(rXYZ,6)*ToReal(a))*(QAD(rXYZ) + 
      SQR(Z*ToReal(a)) + 2*CUB(rXYZ)*ToReal(M))) + (QAD(rXYZ) + 
      SQR(Z*ToReal(a)))*(pow(rXYZ,8) + pow(ToReal(a),6)*SQR(Z) + 
      QAD(ToReal(a))*(QAD(rXYZ) + 2*SQR(rXYZ*Z)) + 
      2*pow(rXYZ,5)*SQR(X)*ToReal(M) + 4*X*Y*QAD(rXYZ)*ToReal(a)*ToReal(M) + 
      CUB(rXYZ)*SQR(ToReal(a))*(2*CUB(rXYZ) + rXYZ*SQR(Z) + 
      2*SQR(Y)*ToReal(M)))*(-2*pow(rXYZ,9) + X*Y*(-3*pow(ToReal(a),5)*SQR(Z) 
      - CUB(ToReal(a))*SQR(rXYZ*Z)) + pow(rXYZ,7)*(3*SQR(X) + SQR(Y) + SQR(Z) 
      - 3*SQR(ToReal(a))) + rXYZ*QAD(ToReal(a))*SQR(Z)*(-3*SQR(X) + SQR(Y) + 
      SQR(Z) - SQR(ToReal(a))) + CUB(rXYZ)*(-SQR(X) + SQR(Y) + SQR(Z) - 
      3*SQR(ToReal(a)))*SQR(Z*ToReal(a)) - 4*pow(rXYZ,8)*ToReal(M) + 
      pow(rXYZ,5)*ToReal(a)*(-CUB(ToReal(a)) + (SQR(X) + SQR(Y) - 
      SQR(Z))*ToReal(a) + 4*X*Y*ToReal(M)) + pow(rXYZ,6)*(3*X*Y*ToReal(a) + 
      2*(2*SQR(X) + SQR(Y) + SQR(Z))*ToReal(M) - 6*SQR(ToReal(a))*ToReal(M)) 
      + QAD(rXYZ)*SQR(ToReal(a))*(X*Y*ToReal(a) + 2*(SQR(Y) + 
      SQR(Z))*ToReal(M) - 2*SQR(ToReal(a))*ToReal(M))) + 
      2*CUB(rXYZ)*(QAD(rXYZ) + SQR(Z*ToReal(a)))*(rXYZ*Y - 
      X*ToReal(a))*(rXYZ*X + 
      Y*ToReal(a))*ToReal(M)*(-4*rXYZ*X*Y*QAD(ToReal(a))*SQR(Z) + 
      pow(ToReal(a),5)*SQR(Z)*(2*SQR(X) - SQR(Y) - SQR(Z) + SQR(ToReal(a))) + 
      CUB(ToReal(a))*SQR(rXYZ*Z)*(-SQR(Y) - SQR(Z) + 3*SQR(ToReal(a))) + 
      2*pow(rXYZ,8)*ToReal(a) + 2*pow(rXYZ,5)*(-3*SQR(X) - SQR(Y) - SQR(Z) + 
      3*SQR(ToReal(a)))*ToReal(a)*ToReal(M) + 
      QAD(rXYZ)*SQR(ToReal(a))*(CUB(ToReal(a)) + (-2*SQR(X) - SQR(Y) + 
      SQR(Z))*ToReal(a) - 2*X*Y*ToReal(M)) + pow(rXYZ,6)*(3*CUB(ToReal(a)) - 
      (4*SQR(X) + SQR(Y) + SQR(Z))*ToReal(a) + 2*X*Y*ToReal(M)) + 
      2*pow(rXYZ,7)*(X*Y + 2*ToReal(a)*ToReal(M)) + 
      2*CUB(rXYZ)*SQR(ToReal(a))*(-(X*Y*SQR(Z)) + (CUB(ToReal(a)) - (SQR(X) + 
      SQR(Y) + SQR(Z))*ToReal(a))*ToReal(M))));
    
    CCTK_REAL K21 = CUB(rXYZ)*INV(alpL*QAD(SQR(rXYZ) + 
      SQR(ToReal(a)))*SQR(QAD(rXYZ) + SQR(Z*ToReal(a)))*SQR(QAD(rXYZ) + 
      SQR(Z*ToReal(a)) + 2*CUB(rXYZ)*ToReal(M))*sqrt(SQR(SQR(X) + SQR(Y) + 
      SQR(Z) - SQR(ToReal(a))) + 
      4*SQR(Z*ToReal(a))))*ToReal(M)*(3*pow(Z,6)*pow(ToReal(a),13)*(-SQR(X) + 
      SQR(Y)) + QAD(Z)*(2*X*Y*pow(ToReal(a),12)*(SQR(Z)*(4*rXYZ - 
      3*ToReal(M)) + 2*SQR(rXYZ)*ToReal(M)) + rXYZ*pow(ToReal(a),11)*(SQR(X) 
      - SQR(Y))*(-5*CUB(rXYZ) - 7*rXYZ*SQR(Z) + 4*(SQR(rXYZ) - 
      3*SQR(Z))*ToReal(M))) + (SQR(X) - 
      SQR(Y))*(pow(rXYZ,14)*ToReal(a)*(3*QAD(rXYZ) - 4*(8*SQR(rXYZ) - 
      9*(SQR(X) + SQR(Y) + SQR(Z)))*SQR(ToReal(M)) - 4*rXYZ*(3*SQR(rXYZ) - 
      5*(SQR(X) + SQR(Y) + SQR(Z)))*ToReal(M)) - 
      pow(rXYZ,8)*pow(ToReal(a),5)*(-(SQR(rXYZ)*(5*QAD(rXYZ) + QAD(Z) + 
      9*SQR(rXYZ*Z))) + 4*(4*QAD(rXYZ) + SQR(Z)*(3*SQR(X) + 3*SQR(Y) + 
      SQR(Z)) - 5*SQR(rXYZ*Z))*SQR(ToReal(M)) + 2*rXYZ*(2*QAD(rXYZ) + 
      2*(SQR(X) + SQR(Y) - 3*SQR(Z))*SQR(Z) - SQR(rXYZ)*(SQR(X) + SQR(Y) + 
      4*SQR(Z)))*ToReal(M)) + CUB(ToReal(a))*pow(rXYZ,10)*(7*pow(rXYZ,6) + 
      5*QAD(rXYZ)*SQR(Z) + 4*(-12*QAD(rXYZ) + 3*SQR(Z)*(SQR(X) + SQR(Y) + 
      SQR(Z)) + SQR(rXYZ)*(5*SQR(X) + 5*SQR(Y) + 11*SQR(Z)))*SQR(ToReal(M)) + 
      2*rXYZ*(-8*QAD(rXYZ) + 10*SQR(Z)*(SQR(X) + SQR(Y) + SQR(Z)) + 
      SQR(rXYZ)*(7*SQR(X) + 7*SQR(Y) + 10*SQR(Z)))*ToReal(M)) - 
      CUB(rXYZ)*pow(ToReal(a),9)*SQR(Z)*(pow(rXYZ,5) + 5*rXYZ*QAD(Z) + 
      9*CUB(rXYZ)*SQR(Z) + (-8*CUB(rXYZ) + 12*rXYZ*SQR(Z))*SQR(ToReal(M)) + 
      (-4*QAD(rXYZ) + 2*SQR(Z)*(5*SQR(X) + 5*SQR(Y) + 16*SQR(Z)) - 
      4*SQR(rXYZ*Z))*ToReal(M)) - pow(rXYZ,5)*pow(ToReal(a),7)*(-pow(rXYZ,7) 
      - 3*pow(rXYZ,5)*SQR(Z) + rXYZ*SQR(Z)*(QAD(Z) + 4*(2*SQR(X) + 2*SQR(Y) + 
      7*SQR(Z))*SQR(ToReal(M))) + CUB(rXYZ)*(3*QAD(Z) - 20*SQR(Z*ToReal(M))) 
      - 12*QAD(rXYZ)*SQR(Z)*ToReal(M) + 2*QAD(Z)*(9*SQR(X) + 9*SQR(Y) + 
      10*SQR(Z))*ToReal(M) + 8*(SQR(X) + SQR(Y) + 
      SQR(Z))*SQR(rXYZ*Z)*ToReal(M))) + X*Y*(-2*pow(rXYZ,15)*(2*QAD(rXYZ) - 
      2*(8*SQR(rXYZ) - 9*(SQR(X) + SQR(Y) + SQR(Z)))*SQR(ToReal(M)) + 
      rXYZ*(-6*SQR(rXYZ) + 11*(SQR(X) + SQR(Y) + SQR(Z)))*ToReal(M)) + 
      4*pow(ToReal(a),10)*SQR(rXYZ*Z)*(5*rXYZ*QAD(Z) + 4*CUB(rXYZ)*SQR(Z) + 
      (2*CUB(rXYZ) - 3*rXYZ*SQR(Z))*SQR(ToReal(M)) + (QAD(rXYZ) + 
      SQR(Z)*(-SQR(X) - SQR(Y) + SQR(Z)) + SQR(rXYZ*Z))*ToReal(M)) + 
      2*(pow(ToReal(a),8)*QAD(rXYZ)*SQR(Z)*(4*pow(rXYZ,5) + 
      18*CUB(rXYZ)*SQR(Z) + rXYZ*(8*QAD(Z) + 2*(3*SQR(rXYZ) - 2*(SQR(X) + 
      SQR(Y) + 2*SQR(Z)))*SQR(ToReal(M))) + (3*QAD(rXYZ) - 
      2*SQR(rXYZ)*(SQR(X) + SQR(Y) - SQR(Z)) + SQR(Z)*(5*(SQR(X) + SQR(Y)) + 
      24*SQR(Z)))*ToReal(M)) + pow(rXYZ,8)*QAD(ToReal(a))*(-2*pow(rXYZ,7) + 
      2*(CUB(rXYZ)*QAD(Z) + rXYZ*(-8*QAD(rXYZ) + 2*SQR(Z)*(3*(SQR(X) + 
      SQR(Y)) + 2*SQR(Z)) + SQR(rXYZ)*(5*(SQR(X) + SQR(Y)) + 
      6*SQR(Z)))*SQR(ToReal(M))) + 3*(-2*pow(rXYZ,6) + QAD(rXYZ)*(SQR(X) + 
      SQR(Y)) + QAD(Z)*(SQR(X) + SQR(Y) + SQR(Z)) + 2*(3*(SQR(X) + SQR(Y)) + 
      2*SQR(Z))*SQR(rXYZ*Z))*ToReal(M))) - 
      4*(pow(rXYZ,11)*SQR(ToReal(a))*(QAD(rXYZ)*(2*SQR(rXYZ) + SQR(Z)) + 
      (-4*QAD(rXYZ) + 3*SQR(Z)*(SQR(X) + SQR(Y) + SQR(Z)) + 
      SQR(rXYZ)*(-4*(SQR(X) + SQR(Y)) + 2*SQR(Z)))*SQR(ToReal(M)) + 
      rXYZ*(-QAD(rXYZ) + 4*SQR(Z)*(SQR(X) + SQR(Y) + SQR(Z)) + 
      2*SQR(rXYZ*Z))*ToReal(M)) + 
      pow(rXYZ,6)*pow(ToReal(a),6)*(rXYZ*(-(SQR(Z)*(3*QAD(rXYZ) + QAD(Z) + 
      6*SQR(rXYZ*Z))) + (4*QAD(rXYZ) + (SQR(X) + SQR(Y) - 
      6*SQR(Z))*SQR(Z))*SQR(ToReal(M))) + (pow(rXYZ,6) + 3*QAD(rXYZ)*SQR(Z) - 
      QAD(Z)*(9*SQR(X) + 9*SQR(Y) + 11*SQR(Z)) - (4*SQR(X) + 4*SQR(Y) + 
      11*SQR(Z))*SQR(rXYZ*Z))*ToReal(M)))));
    
    CCTK_REAL K31 = rXYZ*Z*INV(alpL*CUB(SQR(rXYZ) + 
      SQR(ToReal(a)))*SQR(QAD(rXYZ) + SQR(Z*ToReal(a)))*SQR(QAD(rXYZ) + 
      SQR(Z*ToReal(a)) + 2*CUB(rXYZ)*ToReal(M))*sqrt(SQR(SQR(X) + SQR(Y) + 
      SQR(Z) - SQR(ToReal(a))) + 
      4*SQR(Z*ToReal(a))))*ToReal(M)*(Y*(pow(ToReal(a),13)*(3*pow(Z,6) - 
      2*QAD(Z)*SQR(rXYZ)) - pow(rXYZ,14)*ToReal(a)*(3*QAD(rXYZ) - 
      4*(8*SQR(rXYZ) - 9*(SQR(X) + SQR(Y) + SQR(Z)))*SQR(ToReal(M)) - 
      4*rXYZ*(3*SQR(rXYZ) - 5*(SQR(X) + SQR(Y) + SQR(Z)))*ToReal(M)) + 
      rXYZ*pow(ToReal(a),11)*SQR(Z)*(-4*pow(rXYZ,5) - 3*CUB(rXYZ)*SQR(Z) + 
      rXYZ*SQR(Z)*(2*SQR(X) + 2*SQR(Y) + 9*SQR(Z)) + 8*QAD(Z)*ToReal(M) - 
      4*SQR(rXYZ*Z)*ToReal(M)) + CUB(rXYZ)*pow(ToReal(a),9)*(-2*pow(rXYZ,7) - 
      15*pow(rXYZ,5)*SQR(Z) + CUB(rXYZ)*SQR(Z)*(4*SQR(X) + 4*SQR(Y) + 
      3*SQR(Z)) + rXYZ*QAD(Z)*(4*SQR(X) + 4*SQR(Y) + 9*SQR(Z) + 
      4*SQR(ToReal(M))) + 2*QAD(Z)*(5*SQR(X) + 5*SQR(Y) + 
      12*SQR(Z))*ToReal(M) - 4*(SQR(X) + SQR(Y) + 
      4*SQR(Z))*SQR(rXYZ*Z)*ToReal(M)) - 
      CUB(ToReal(a))*pow(rXYZ,10)*(QAD(rXYZ)*(11*SQR(rXYZ) - 2*SQR(X) - 
      2*SQR(Y) + 3*SQR(Z)) + 4*(-16*QAD(rXYZ) + 3*SQR(Z)*(SQR(X) + SQR(Y) + 
      SQR(Z)) + 3*SQR(rXYZ)*(3*SQR(X) + 3*SQR(Y) + 5*SQR(Z)))*SQR(ToReal(M)) 
      + 2*rXYZ*(-12*QAD(rXYZ) + 3*SQR(rXYZ)*(5*SQR(X) + 5*SQR(Y) + 6*SQR(Z)) 
      - 2*(QAD(X) + QAD(Y) - 4*QAD(Z) + SQR(X)*(2*SQR(Y) - 3*SQR(Z)) - 
      3*SQR(Y*Z)))*ToReal(M)) + 
      pow(rXYZ,5)*pow(ToReal(a),7)*(rXYZ*(-9*pow(rXYZ,6) + 
      QAD(rXYZ)*(2*SQR(X) + 2*SQR(Y) - 21*SQR(Z)) + QAD(Z)*(2*SQR(X) + 
      2*SQR(Y) + 3*SQR(Z)) + (8*SQR(X) + 8*SQR(Y) + 7*SQR(Z))*SQR(rXYZ*Z)) + 
      4*rXYZ*(2*QAD(rXYZ) + QAD(Z) - SQR(rXYZ*Z))*SQR(ToReal(M)) + 
      (4*pow(rXYZ,6) - 4*QAD(rXYZ)*(SQR(X) + SQR(Y) + 2*SQR(Z)) - 8*(SQR(X) + 
      SQR(Y) + 3*SQR(Z))*SQR(rXYZ*Z) + 2*SQR(Z)*(2*QAD(X) + 2*QAD(Y) + 
      10*QAD(Z) + SQR(X)*(4*SQR(Y) + 11*SQR(Z)) + 11*SQR(Y*Z)))*ToReal(M)) + 
      pow(rXYZ,7)*pow(ToReal(a),5)*(CUB(rXYZ)*(-15*QAD(rXYZ) + 
      SQR(rXYZ)*(4*SQR(X) + 4*SQR(Y) - 13*SQR(Z)) + SQR(Z)*(4*SQR(X) + 
      4*SQR(Y) + 3*SQR(Z))) + 4*rXYZ*(10*QAD(rXYZ) - SQR(Z)*(SQR(X) + SQR(Y) 
      + 3*SQR(Z)) - SQR(rXYZ)*(2*SQR(X) + 2*SQR(Y) + 
      7*SQR(Z)))*SQR(ToReal(M)) + 2*(8*pow(rXYZ,6) - 3*QAD(rXYZ)*(3*SQR(X) + 
      3*SQR(Y) + 4*SQR(Z)) + 2*SQR(rXYZ)*(QAD(X) + QAD(Y) - 7*QAD(Z) + 
      2*SQR(X)*(SQR(Y) - SQR(Z)) - 2*SQR(Y*Z)) + 2*SQR(Z*(SQR(X) + SQR(Y) + 
      SQR(Z))))*ToReal(M))) + X*(-2*(pow(ToReal(a),12)*QAD(Z)*(CUB(rXYZ) - 
      3*rXYZ*SQR(Z)) + pow(rXYZ,15)*(2*QAD(rXYZ) - 2*(8*SQR(rXYZ) - 9*(SQR(X) 
      + SQR(Y) + SQR(Z)))*SQR(ToReal(M)) + rXYZ*(-6*SQR(rXYZ) + 11*(SQR(X) + 
      SQR(Y) + SQR(Z)))*ToReal(M)) + 
      pow(rXYZ,11)*SQR(ToReal(a))*(QAD(rXYZ)*(7*SQR(rXYZ) - SQR(X) - SQR(Y) + 
      SQR(Z)) + (-32*QAD(rXYZ) + 6*(SQR(Z)*(SQR(X) + SQR(Y) + SQR(Z)) + 
      SQR(rXYZ)*(3*(SQR(X) + SQR(Y)) + 5*SQR(Z))))*SQR(ToReal(M)) + 
      rXYZ*(-12*QAD(rXYZ) + SQR(rXYZ)*(17*(SQR(X) + SQR(Y)) + 21*SQR(Z)) - 
      2*(QAD(X) + QAD(Y) - 3*QAD(Z) + 2*SQR(X)*(SQR(Y) - SQR(Z)) - 
      2*SQR(Y*Z)))*ToReal(M)) + pow(ToReal(a),8)*QAD(rXYZ)*(pow(rXYZ,7) + 
      SQR(Z)*(7*pow(rXYZ,5) - CUB(rXYZ)*(2*SQR(X) + 2*SQR(Y) + 9*SQR(Z))) + 
      2*(SQR(X) + SQR(Y) + 3*SQR(Z))*SQR(rXYZ*Z)*ToReal(M) + 
      QAD(Z)*(-(rXYZ*(2*SQR(X) + 2*SQR(Y) + 9*SQR(Z) + 2*SQR(ToReal(M)))) - 
      (8*SQR(X) + 8*SQR(Y) + 21*SQR(Z))*ToReal(M)))) + 
      2*(pow(rXYZ,8)*QAD(ToReal(a))*(-(CUB(rXYZ)*(SQR(rXYZ) + 
      SQR(Z))*(9*SQR(rXYZ) - 2*(SQR(X) + SQR(Y) + 2*SQR(Z)))) + 
      2*rXYZ*(10*QAD(rXYZ) - SQR(Z)*(SQR(X) + SQR(Y) + 3*SQR(Z)) - 
      SQR(rXYZ)*(2*SQR(X) + 2*SQR(Y) + 7*SQR(Z)))*SQR(ToReal(M)) + 
      (8*pow(rXYZ,6) + 2*SQR(rXYZ)*(QAD(X) + QAD(Y) - 4*QAD(Z) + 2*SQR(X*Y)) 
      - 5*QAD(rXYZ)*(2*(SQR(X) + SQR(Y)) + 3*SQR(Z)) + SQR(Z)*(2*(QAD(X) + 
      QAD(Y)) + 5*QAD(Z) + SQR(X)*(4*SQR(Y) + 7*SQR(Z)) + 
      7*SQR(Y*Z)))*ToReal(M)) + 
      pow(rXYZ,6)*pow(ToReal(a),6)*(rXYZ*(-5*pow(rXYZ,6) + QAD(rXYZ)*(SQR(X) 
      + SQR(Y) - 9*SQR(Z)) + QAD(Z)*(SQR(X) + SQR(Y) + 3*SQR(Z)) + (4*(SQR(X) 
      + SQR(Y)) + 11*SQR(Z))*SQR(rXYZ*Z) + 2*(2*QAD(rXYZ) + QAD(Z) - 
      SQR(rXYZ*Z))*SQR(ToReal(M))) + (2*pow(rXYZ,6) - QAD(rXYZ)*(2*(SQR(X) + 
      SQR(Y)) + 5*SQR(Z)) - 2*(SQR(X) + SQR(Y) + 3*SQR(Z))*SQR(rXYZ*Z) + 
      SQR(Z)*(2*(QAD(X) + QAD(Y)) + 19*QAD(Z) + SQR(X)*(4*SQR(Y) + 17*SQR(Z)) 
      + 17*SQR(Y*Z)))*ToReal(M)) + 
      pow(ToReal(a),10)*SQR(rXYZ*Z)*(SQR(Z)*(CUB(rXYZ) + rXYZ*(SQR(X) + 
      SQR(Y) + 9*SQR(Z))) + 7*QAD(Z)*ToReal(M) - 2*(pow(rXYZ,5) + 
      SQR(rXYZ*Z)*ToReal(M))))));
    
    CCTK_REAL K22 = -2*CUB(rXYZ)*INV(alpL*QAD(SQR(rXYZ) + 
      SQR(ToReal(a)))*SQR(QAD(rXYZ) + SQR(Z*ToReal(a)))*SQR(QAD(rXYZ) + 
      SQR(Z*ToReal(a)) + 2*CUB(rXYZ)*ToReal(M))*sqrt(SQR(SQR(X) + SQR(Y) + 
      SQR(Z) - SQR(ToReal(a))) + 4*SQR(Z*ToReal(a))))*ToReal(M)*((QAD(rXYZ) + 
      SQR(Z*ToReal(a)))*(-2*pow(rXYZ,9) + X*Y*(3*pow(ToReal(a),5)*SQR(Z) + 
      CUB(ToReal(a))*SQR(rXYZ*Z)) + pow(rXYZ,7)*(SQR(X) + 3*SQR(Y) + SQR(Z) - 
      3*SQR(ToReal(a))) + rXYZ*QAD(ToReal(a))*SQR(Z)*(SQR(X) - 3*SQR(Y) + 
      SQR(Z) - SQR(ToReal(a))) + CUB(rXYZ)*(SQR(X) - SQR(Y) + SQR(Z) - 
      3*SQR(ToReal(a)))*SQR(Z*ToReal(a)) - 4*pow(rXYZ,8)*ToReal(M) - 
      pow(rXYZ,5)*ToReal(a)*(CUB(ToReal(a)) - (SQR(X) + SQR(Y) - 
      SQR(Z))*ToReal(a) + 4*X*Y*ToReal(M)) + pow(rXYZ,6)*(-3*X*Y*ToReal(a) + 
      (2*(SQR(X) + 2*SQR(Y) + SQR(Z)) - 6*SQR(ToReal(a)))*ToReal(M)) + 
      QAD(rXYZ)*SQR(ToReal(a))*(-(X*Y*ToReal(a)) + (2*(SQR(X) + SQR(Z)) - 
      2*SQR(ToReal(a)))*ToReal(M)))*(pow(rXYZ,8) + pow(ToReal(a),6)*SQR(Z) + 
      QAD(ToReal(a))*(QAD(rXYZ) + 2*SQR(rXYZ*Z)) + (2*pow(rXYZ,5)*SQR(Y) - 
      4*X*Y*QAD(rXYZ)*ToReal(a))*ToReal(M) + 
      CUB(rXYZ)*SQR(ToReal(a))*(rXYZ*SQR(Z) + 2*(CUB(rXYZ) + 
      SQR(X)*ToReal(M)))) + ToReal(M)*((-(SQR(Z*(SQR(rXYZ) + 
      SQR(ToReal(a))))*(-(rXYZ*Y) + X*ToReal(a))*(Y*(3*pow(rXYZ,7) + 
      rXYZ*QAD(ToReal(a))*(2*SQR(rXYZ) - 5*SQR(Z)) + CUB(rXYZ)*(3*SQR(rXYZ) - 
      2*(SQR(X) + SQR(Y)) - 3*SQR(Z))*SQR(ToReal(a))) + 
      X*(CUB(ToReal(a))*SQR(rXYZ)*(-5*SQR(rXYZ) + 2*(SQR(X) + SQR(Y)) + 
      SQR(Z)) + pow(ToReal(a),5)*(-2*SQR(rXYZ) + 3*SQR(Z)) - 
      5*pow(rXYZ,6)*ToReal(a)))) + CUB(rXYZ)*((rXYZ*Y - X*ToReal(a))*(rXYZ*X 
      + Y*ToReal(a))*(pow(ToReal(a),5)*SQR(Z)*(SQR(X) - 2*(SQR(Y) + SQR(Z)) + 
      2*SQR(ToReal(a))) + CUB(ToReal(a))*(SQR(rXYZ*Z)*(-3*SQR(X) - 2*(SQR(Y) 
      + SQR(Z)) + 6*SQR(ToReal(a))) + QAD(rXYZ)*(-3*SQR(X) - 2*SQR(Y) + 
      2*(SQR(Z) + SQR(ToReal(a))))) + X*Y*(3*pow(rXYZ,7) - 
      5*rXYZ*QAD(ToReal(a))*SQR(Z) - pow(rXYZ,5)*SQR(ToReal(a)) - 
      CUB(rXYZ)*SQR(Z*ToReal(a))) + (4*pow(rXYZ,8) + pow(rXYZ,6)*(-7*SQR(X) - 
      2*(SQR(Y) + SQR(Z)) + 6*SQR(ToReal(a))))*ToReal(a)) - SQR(-(rXYZ*Y) + 
      X*ToReal(a))*(pow(rXYZ,7)*(4*SQR(rXYZ) - 5*SQR(Y) - 2*(SQR(X) + 
      SQR(Z))) + rXYZ*(2*pow(ToReal(a),6)*SQR(Z) + 
      QAD(ToReal(a))*(2*QAD(rXYZ) + SQR(Z)*(3*SQR(Y) - 2*(SQR(X) + SQR(Z))) + 
      6*SQR(rXYZ*Z))) + CUB(rXYZ)*(6*QAD(rXYZ) - SQR(rXYZ)*(2*SQR(X) + SQR(Y) 
      - 2*SQR(Z)) - SQR(Z)*(2*SQR(X) + SQR(Y) + 2*SQR(Z)))*SQR(ToReal(a)) + 
      X*Y*(-3*pow(ToReal(a),5)*SQR(Z) + CUB(ToReal(a))*SQR(rXYZ)*(SQR(rXYZ) + 
      SQR(Z)) + 5*pow(rXYZ,6)*ToReal(a)))))*(QAD(rXYZ) + SQR(Z*ToReal(a)) + 
      2*CUB(rXYZ)*ToReal(M)) + (QAD(rXYZ) + SQR(Z*ToReal(a)))*(rXYZ*Y - 
      X*ToReal(a))*(4*rXYZ*Y*CUB(SQR(rXYZ) + 
      SQR(ToReal(a)))*SQR(Z)*(QAD(rXYZ) - SQR(Z*ToReal(a)) + 
      CUB(rXYZ)*ToReal(M)) + 2*CUB(rXYZ)*(rXYZ*X + 
      Y*ToReal(a))*((-pow(ToReal(a),7) - 4*rXYZ*X*Y*QAD(ToReal(a)))*SQR(Z) + 
      pow(ToReal(a),5)*(-QAD(rXYZ) + SQR(Z)*(SQR(X) - 2*SQR(Y) + SQR(Z)) - 
      3*SQR(rXYZ*Z) - 2*CUB(rXYZ)*ToReal(M)) + 
      CUB(ToReal(a))*SQR(rXYZ)*(-3*QAD(rXYZ) + SQR(rXYZ)*(SQR(X) + 2*SQR(Y) - 
      SQR(Z)) + SQR(Z)*(SQR(X) + SQR(Z)) + 2*rXYZ*(-3*SQR(rXYZ) + SQR(X) + 
      SQR(Y) + SQR(Z))*ToReal(M)) + pow(rXYZ,5)*ToReal(a)*(rXYZ*(-2*SQR(rXYZ) 
      + SQR(X) + 4*SQR(Y) + SQR(Z)) + 2*(-2*SQR(rXYZ) + SQR(X) + 3*SQR(Y) + 
      SQR(Z))*ToReal(M)) + X*Y*(2*pow(rXYZ,6)*(rXYZ + ToReal(M)) - 
      2*CUB(rXYZ)*SQR(ToReal(a))*(SQR(Z) + rXYZ*ToReal(M)))))));
    
    CCTK_REAL K32 = rXYZ*Z*INV(alpL*CUB(SQR(rXYZ) + 
      SQR(ToReal(a)))*SQR(QAD(rXYZ) + SQR(Z*ToReal(a)))*SQR(QAD(rXYZ) + 
      SQR(Z*ToReal(a)) + 2*CUB(rXYZ)*ToReal(M))*sqrt(SQR(SQR(X) + SQR(Y) + 
      SQR(Z) - SQR(ToReal(a))) + 
      4*SQR(Z*ToReal(a))))*ToReal(M)*(QAD(Z)*(X*pow(ToReal(a),13)*(2*SQR(rXYZ) 
      - 3*SQR(Z)) - 2*Y*pow(ToReal(a),12)*(CUB(rXYZ) - 3*rXYZ*SQR(Z))) + 
      Y*(-2*pow(rXYZ,15)*(2*QAD(rXYZ) - 2*(8*SQR(rXYZ) - 9*(SQR(X) + SQR(Y) + 
      SQR(Z)))*SQR(ToReal(M)) + rXYZ*(-6*SQR(rXYZ) + 11*(SQR(X) + SQR(Y) + 
      SQR(Z)))*ToReal(M)) + 2*pow(ToReal(a),10)*SQR(rXYZ*Z)*(-2*pow(rXYZ,5) + 
      CUB(rXYZ)*SQR(Z) + rXYZ*SQR(Z)*(SQR(X) + SQR(Y) + 9*SQR(Z)) + 
      7*QAD(Z)*ToReal(M) - 2*SQR(rXYZ*Z)*ToReal(M))) + 
      X*(pow(rXYZ,14)*ToReal(a)*(3*QAD(rXYZ) - 4*(8*SQR(rXYZ) - 9*(SQR(X) + 
      SQR(Y) + SQR(Z)))*SQR(ToReal(M)) - 4*rXYZ*(3*SQR(rXYZ) - 5*(SQR(X) + 
      SQR(Y) + SQR(Z)))*ToReal(M)) + 
      rXYZ*pow(ToReal(a),11)*SQR(Z)*(4*pow(rXYZ,5) + SQR(Z)*(3*CUB(rXYZ) - 
      rXYZ*(2*SQR(X) + 2*SQR(Y) + 9*SQR(Z))) - 8*QAD(Z)*ToReal(M) + 
      4*SQR(rXYZ*Z)*ToReal(M))) + 
      2*(Y*pow(rXYZ,8)*QAD(ToReal(a))*(-(CUB(rXYZ)*(SQR(rXYZ) + 
      SQR(Z))*(9*SQR(rXYZ) - 2*(SQR(X) + SQR(Y) + 2*SQR(Z)))) + 
      2*rXYZ*(10*QAD(rXYZ) - SQR(Z)*(SQR(X) + SQR(Y) + 3*SQR(Z)) - 
      SQR(rXYZ)*(2*SQR(X) + 2*SQR(Y) + 7*SQR(Z)))*SQR(ToReal(M)) + 
      (8*pow(rXYZ,6) + 2*SQR(rXYZ)*(QAD(X) + QAD(Y) - 4*QAD(Z) + 2*SQR(X*Y)) 
      - 5*QAD(rXYZ)*(2*SQR(X) + 2*SQR(Y) + 3*SQR(Z)) + SQR(Z)*(2*QAD(X) + 
      2*QAD(Y) + 5*QAD(Z) + SQR(X)*(4*SQR(Y) + 7*SQR(Z)) + 
      7*SQR(Y*Z)))*ToReal(M)) + 
      Y*pow(rXYZ,6)*pow(ToReal(a),6)*(rXYZ*(-5*pow(rXYZ,6) + 
      QAD(rXYZ)*(SQR(X) + SQR(Y) - 9*SQR(Z)) + QAD(Z)*(SQR(X) + SQR(Y) + 
      3*SQR(Z)) + (4*SQR(X) + 4*SQR(Y) + 11*SQR(Z))*SQR(rXYZ*Z)) + 
      2*rXYZ*(2*QAD(rXYZ) + QAD(Z) - SQR(rXYZ*Z))*SQR(ToReal(M)) + 
      (2*pow(rXYZ,6) - QAD(rXYZ)*(2*SQR(X) + 2*SQR(Y) + 5*SQR(Z)) - 2*(SQR(X) 
      + SQR(Y) + 3*SQR(Z))*SQR(rXYZ*Z) + SQR(Z)*(2*QAD(X) + 2*QAD(Y) + 
      19*QAD(Z) + SQR(X)*(4*SQR(Y) + 17*SQR(Z)) + 17*SQR(Y*Z)))*ToReal(M))) + 
      X*(-(pow(rXYZ,5)*pow(ToReal(a),7)*(rXYZ*(-9*pow(rXYZ,6) + 
      QAD(rXYZ)*(2*SQR(X) + 2*SQR(Y) - 21*SQR(Z)) + QAD(Z)*(2*SQR(X) + 
      2*SQR(Y) + 3*SQR(Z)) + (8*SQR(X) + 8*SQR(Y) + 7*SQR(Z))*SQR(rXYZ*Z)) + 
      4*rXYZ*(2*QAD(rXYZ) + QAD(Z) - SQR(rXYZ*Z))*SQR(ToReal(M)) + 
      (4*pow(rXYZ,6) - 4*QAD(rXYZ)*(SQR(X) + SQR(Y) + 2*SQR(Z)) - 8*(SQR(X) + 
      SQR(Y) + 3*SQR(Z))*SQR(rXYZ*Z) + 2*SQR(Z)*(2*QAD(X) + 2*QAD(Y) + 
      10*QAD(Z) + SQR(X)*(4*SQR(Y) + 11*SQR(Z)) + 11*SQR(Y*Z)))*ToReal(M))) - 
      pow(rXYZ,7)*pow(ToReal(a),5)*(CUB(rXYZ)*(-15*QAD(rXYZ) + 
      SQR(rXYZ)*(4*SQR(X) + 4*SQR(Y) - 13*SQR(Z)) + SQR(Z)*(4*SQR(X) + 
      4*SQR(Y) + 3*SQR(Z))) + 4*rXYZ*(10*QAD(rXYZ) - SQR(Z)*(SQR(X) + SQR(Y) 
      + 3*SQR(Z)) - SQR(rXYZ)*(2*SQR(X) + 2*SQR(Y) + 
      7*SQR(Z)))*SQR(ToReal(M)) + 2*(8*pow(rXYZ,6) - 3*QAD(rXYZ)*(3*SQR(X) + 
      3*SQR(Y) + 4*SQR(Z)) + 2*SQR(rXYZ)*(QAD(X) + QAD(Y) - 7*QAD(Z) + 
      2*SQR(X)*(SQR(Y) - SQR(Z)) - 2*SQR(Y*Z)) + 2*SQR(Z*(SQR(X) + SQR(Y) + 
      SQR(Z))))*ToReal(M))) + 
      X*(CUB(ToReal(a))*pow(rXYZ,10)*(QAD(rXYZ)*(11*SQR(rXYZ) - 2*SQR(X) - 
      2*SQR(Y) + 3*SQR(Z)) + 4*(-16*QAD(rXYZ) + 3*SQR(Z)*(SQR(X) + SQR(Y) + 
      SQR(Z)) + 3*SQR(rXYZ)*(3*SQR(X) + 3*SQR(Y) + 5*SQR(Z)))*SQR(ToReal(M)) 
      + 2*rXYZ*(-12*QAD(rXYZ) + 3*SQR(rXYZ)*(5*SQR(X) + 5*SQR(Y) + 6*SQR(Z)) 
      - 2*(QAD(X) + QAD(Y) - 4*QAD(Z) + SQR(X)*(2*SQR(Y) - 3*SQR(Z)) - 
      3*SQR(Y*Z)))*ToReal(M)) + CUB(rXYZ)*pow(ToReal(a),9)*(2*pow(rXYZ,7) + 
      SQR(Z)*(15*pow(rXYZ,5) - CUB(rXYZ)*(4*SQR(X) + 4*SQR(Y) + 3*SQR(Z))) + 
      4*(SQR(X) + SQR(Y) + 4*SQR(Z))*SQR(rXYZ*Z)*ToReal(M) + 
      QAD(Z)*(-(rXYZ*(4*SQR(X) + 4*SQR(Y) + 9*SQR(Z) + 4*SQR(ToReal(M)))) - 
      2*(5*SQR(X) + 5*SQR(Y) + 12*SQR(Z))*ToReal(M)))) - 
      2*(Y*pow(rXYZ,11)*SQR(ToReal(a))*(QAD(rXYZ)*(7*SQR(rXYZ) - SQR(X) - 
      SQR(Y) + SQR(Z)) + (-32*QAD(rXYZ) + 6*SQR(Z)*(SQR(X) + SQR(Y) + SQR(Z)) 
      + 6*SQR(rXYZ)*(3*SQR(X) + 3*SQR(Y) + 5*SQR(Z)))*SQR(ToReal(M)) + 
      rXYZ*(-12*QAD(rXYZ) + SQR(rXYZ)*(17*SQR(X) + 17*SQR(Y) + 21*SQR(Z)) - 
      2*(QAD(X) + QAD(Y) - 3*QAD(Z) + 2*SQR(X)*(SQR(Y) - SQR(Z)) - 
      2*SQR(Y*Z)))*ToReal(M)) + Y*pow(ToReal(a),8)*QAD(rXYZ)*(pow(rXYZ,7) + 
      SQR(Z)*(7*pow(rXYZ,5) - CUB(rXYZ)*(2*SQR(X) + 2*SQR(Y) + 9*SQR(Z))) + 
      2*(SQR(X) + SQR(Y) + 3*SQR(Z))*SQR(rXYZ*Z)*ToReal(M) + 
      QAD(Z)*(-(rXYZ*(2*SQR(X) + 2*SQR(Y) + 9*SQR(Z) + 2*SQR(ToReal(M)))) - 
      (8*SQR(X) + 8*SQR(Y) + 21*SQR(Z))*ToReal(M)))));
    
    CCTK_REAL K33 = -2*INV(alpL*SQR(QAD(rXYZ) + 
      SQR(Z*ToReal(a)))*(SQR(rXYZ) + SQR(ToReal(a)))*SQR(QAD(rXYZ) + 
      SQR(Z*ToReal(a)) + 2*CUB(rXYZ)*ToReal(M))*sqrt(SQR(SQR(X) + SQR(Y) + 
      SQR(Z) - SQR(ToReal(a))) + 
      4*SQR(Z*ToReal(a))))*ToReal(M)*(-2*(pow(rXYZ,18) + 
      pow(Z,8)*pow(ToReal(a),10)) + pow(rXYZ,16)*(SQR(X) + SQR(Y) + 3*SQR(Z) 
      - 3*SQR(ToReal(a))) + pow(rXYZ,12)*SQR(Z)*(-QAD(ToReal(a)) + 
      SQR(ToReal(a))*(SQR(X) + SQR(Y) + 3*SQR(Z) - 24*SQR(ToReal(M))) + 
      18*(SQR(X) + SQR(Y) + SQR(Z))*SQR(ToReal(M))) + 
      pow(rXYZ,8)*QAD(Z)*SQR(ToReal(a))*(5*QAD(ToReal(a)) + 6*(SQR(X) + 
      SQR(Y) + SQR(Z))*SQR(ToReal(M)) - SQR(ToReal(a))*(SQR(X) + SQR(Y) + 
      3*SQR(Z) + 2*SQR(ToReal(M)))) + pow(rXYZ,14)*(-QAD(ToReal(a)) + (SQR(X) 
      + SQR(Y) + 3*SQR(Z))*SQR(ToReal(a)) - 16*SQR(Z*ToReal(M))) + 
      (-4*pow(rXYZ,17) - 5*rXYZ*pow(Z,8)*pow(ToReal(a),8) + 
      2*pow(rXYZ,15)*(SQR(X) + SQR(Y) - 2*SQR(Z) - 3*SQR(ToReal(a))) + 
      pow(rXYZ,13)*(-2*QAD(ToReal(a)) + 11*SQR(Z)*(SQR(X) + SQR(Y) + SQR(Z)) 
      + 2*(SQR(X) + SQR(Y) - 7*SQR(Z))*SQR(ToReal(a))) - 
      pow(rXYZ,5)*QAD(Z*ToReal(a))*(4*(QAD(X) + QAD(Y)) + 7*QAD(Z) + 
      2*QAD(ToReal(a)) + SQR(X)*(8*SQR(Y) + 11*SQR(Z)) + 11*SQR(Y*Z) - 
      2*(3*(SQR(X) + SQR(Y)) + 7*SQR(Z))*SQR(ToReal(a))) + 
      pow(rXYZ,9)*(-4*(QAD(X) + QAD(Y) - QAD(Z) + QAD(ToReal(a)) + 
      2*SQR(X*Y)) + (8*(SQR(X) + SQR(Y)) + 
      9*SQR(Z))*SQR(ToReal(a)))*SQR(Z*ToReal(a)))*ToReal(M) + 
      SQR(Z*ToReal(a))*(pow(rXYZ,10)*(-QAD(ToReal(a)) + 
      SQR(ToReal(a))*(SQR(X) + SQR(Y) + 7*SQR(Z) - 8*SQR(ToReal(M))) + 
      4*(2*(SQR(X) + SQR(Y)) + 5*SQR(Z))*SQR(ToReal(M))) + 
      2*pow(rXYZ,11)*(9*(SQR(X) + SQR(Y)) + 11*SQR(Z) - 
      7*SQR(ToReal(a)))*ToReal(M)) + 
      QAD(Z*ToReal(a))*(pow(rXYZ,6)*(QAD(ToReal(a)) + 4*(SQR(X) + SQR(Y) + 
      2*SQR(Z))*SQR(ToReal(M)) - SQR(ToReal(a))*(SQR(X) + SQR(Y) + 3*SQR(Z) + 
      4*SQR(ToReal(M)))) + 2*pow(rXYZ,7)*(4*(SQR(X) + SQR(Y)) + 7*SQR(Z) - 
      2*SQR(ToReal(a)))*ToReal(M)) + 
      pow(Z,6)*(pow(ToReal(a),8)*SQR(rXYZ)*(-SQR(X) - SQR(Y) - 5*SQR(Z) + 
      SQR(ToReal(a))) + pow(ToReal(a),6)*(QAD(rXYZ)*(-SQR(X) - SQR(Y) - 
      3*SQR(Z) + SQR(ToReal(a)) + 2*SQR(ToReal(M))) + 4*CUB(rXYZ)*(-2*(SQR(X) 
      + SQR(Y)) - 3*SQR(Z) + SQR(ToReal(a)))*ToReal(M))));
    
    CCTK_REAL betap1 = 2*CUB(rXYZ)*INV((SQR(rXYZ) + 
      SQR(ToReal(a)))*(QAD(rXYZ) + SQR(Z*ToReal(a)) + 
      2*CUB(rXYZ)*ToReal(M)))*(rXYZ*X + Y*ToReal(a))*ToReal(M);
    
    CCTK_REAL betap2 = 2*CUB(rXYZ)*INV((SQR(rXYZ) + 
      SQR(ToReal(a)))*(QAD(rXYZ) + SQR(Z*ToReal(a)) + 
      2*CUB(rXYZ)*ToReal(M)))*(rXYZ*Y - X*ToReal(a))*ToReal(M);
    
    CCTK_REAL betap3 = 2*Z*INV(QAD(rXYZ) + SQR(Z*ToReal(a)) + 
      2*CUB(rXYZ)*ToReal(M))*SQR(rXYZ)*ToReal(M);
    
    CCTK_REAL dtbetap1 = 0;
    
    CCTK_REAL dtbetap2 = 0;
    
    CCTK_REAL dtbetap3 = 0;
    
    CCTK_REAL gxxL = 2*(G32*Jac21*Jac31 + Jac11*(G21*Jac21 + G31*Jac31)) + 
      G11*SQR(Jac11) + G22*SQR(Jac21) + G33*SQR(Jac31);
    
    CCTK_REAL gxyL = Jac12*(G11*Jac11 + G21*Jac21 + G31*Jac31) + 
      Jac22*(G21*Jac11 + G22*Jac21 + G32*Jac31) + (G31*Jac11 + G32*Jac21 + 
      G33*Jac31)*Jac32;
    
    CCTK_REAL gxzL = Jac13*(G11*Jac11 + G21*Jac21 + G31*Jac31) + 
      Jac23*(G21*Jac11 + G22*Jac21 + G32*Jac31) + (G31*Jac11 + G32*Jac21 + 
      G33*Jac31)*Jac33;
    
    CCTK_REAL gyyL = 2*(G32*Jac22*Jac32 + Jac12*(G21*Jac22 + G31*Jac32)) + 
      G11*SQR(Jac12) + G22*SQR(Jac22) + G33*SQR(Jac32);
    
    CCTK_REAL gyzL = Jac13*(G11*Jac12 + G21*Jac22 + G31*Jac32) + 
      Jac23*(G21*Jac12 + G22*Jac22 + G32*Jac32) + (G31*Jac12 + G32*Jac22 + 
      G33*Jac32)*Jac33;
    
    CCTK_REAL gzzL = 2*(G32*Jac23*Jac33 + Jac13*(G21*Jac23 + G31*Jac33)) + 
      G11*SQR(Jac13) + G22*SQR(Jac23) + G33*SQR(Jac33);
    
    CCTK_REAL kxxL = 2*(Jac11*(Jac21*K21 + Jac31*K31) + Jac21*Jac31*K32) + 
      K11*SQR(Jac11) + K22*SQR(Jac21) + K33*SQR(Jac31);
    
    CCTK_REAL kxyL = Jac11*(Jac12*K11 + Jac22*K21 + Jac32*K31) + 
      Jac21*(Jac12*K21 + Jac22*K22 + Jac32*K32) + Jac31*(Jac12*K31 + 
      Jac22*K32 + Jac32*K33);
    
    CCTK_REAL kxzL = Jac11*(Jac13*K11 + Jac23*K21 + Jac33*K31) + 
      Jac21*(Jac13*K21 + Jac23*K22 + Jac33*K32) + Jac31*(Jac13*K31 + 
      Jac23*K32 + Jac33*K33);
    
    CCTK_REAL kyyL = 2*(Jac12*(Jac22*K21 + Jac32*K31) + Jac22*Jac32*K32) + 
      K11*SQR(Jac12) + K22*SQR(Jac22) + K33*SQR(Jac32);
    
    CCTK_REAL kyzL = Jac12*(Jac13*K11 + Jac23*K21 + Jac33*K31) + 
      Jac22*(Jac13*K21 + Jac23*K22 + Jac33*K32) + Jac32*(Jac13*K31 + 
      Jac23*K32 + Jac33*K33);
    
    CCTK_REAL kzzL = 2*(Jac13*(Jac23*K21 + Jac33*K31) + Jac23*Jac33*K32) + 
      K11*SQR(Jac13) + K22*SQR(Jac23) + K33*SQR(Jac33);
    
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
  CCTK_ENDLOOP3(KerrSchild_always);
}

extern "C" void KerrSchild_always(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering KerrSchild_always_Body");
  }
  
  if (cctk_iteration % KerrSchild_always_calc_every != KerrSchild_always_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"admbase::curv","admbase::dtlapse","admbase::dtshift","admbase::lapse","admbase::metric","admbase::shift","grid::coordinates"};
  GenericFD_AssertGroupStorage(cctkGH, "KerrSchild_always", 7, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, KerrSchild_always_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving KerrSchild_always_Body");
  }
}
