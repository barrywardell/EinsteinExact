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

/* Define macros used in calculations */
#define INITVALUE (42)
#define QAD(x) (SQR(SQR(x)))
#define INV(x) ((1.0) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))

static void KerrSchild_initial_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const min[3], int const max[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  /* Declare the variables used for looping over grid points */
  CCTK_INT i, j, k;
  // CCTK_INT index = INITVALUE;
  
  /* Declare finite differencing variables */
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering KerrSchild_initial_Body");
  }
  
  if (cctk_iteration % KerrSchild_initial_calc_every != KerrSchild_initial_calc_offset)
  {
    return;
  }
  
  const char *groups[] = {"admbase::curv","admbase::dtlapse","admbase::dtshift","admbase::lapse","admbase::metric","admbase::shift","grid::coordinates"};
  GenericFD_AssertGroupStorage(cctkGH, "KerrSchild_initial", 7, groups);
  
  
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
  
  /* Loop over the grid points */
  for (k = min[2]; k < max[2]; k++)
  {
    for (j = min[1]; j < max[1]; j++)
    {
      for (i = min[0]; i < max[0]; i++)
      {
         int  const  index  =  CCTK_GFINDEX3D(cctkGH,i,j,k) ;
        
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
        
        CCTK_REAL Jac11 = Cos(ToReal(phi))*Cos(ToReal(psi)) - 
          Cos(ToReal(theta))*Sin(ToReal(phi))*Sin(ToReal(psi));
        
        CCTK_REAL Jac12 = Cos(ToReal(psi))*Sin(ToReal(phi)) + 
          Cos(ToReal(phi))*Cos(ToReal(theta))*Sin(ToReal(psi));
        
        CCTK_REAL Jac13 = Sin(ToReal(psi))*Sin(ToReal(theta));
        
        CCTK_REAL Jac21 = 
          -(Cos(ToReal(psi))*Cos(ToReal(theta))*Sin(ToReal(phi))) - 
          Cos(ToReal(phi))*Sin(ToReal(psi));
        
        CCTK_REAL Jac22 = Cos(ToReal(phi))*Cos(ToReal(psi))*Cos(ToReal(theta)) 
          - Sin(ToReal(phi))*Sin(ToReal(psi));
        
        CCTK_REAL Jac23 = Cos(ToReal(psi))*Sin(ToReal(theta));
        
        CCTK_REAL Jac31 = Sin(ToReal(phi))*Sin(ToReal(theta));
        
        CCTK_REAL Jac32 = -(Cos(ToReal(phi))*Sin(ToReal(theta)));
        
        CCTK_REAL Jac33 = Cos(ToReal(theta));
        
        CCTK_REAL XX1 = Jac11*xx1 + Jac12*xx2 + Jac13*xx3;
        
        CCTK_REAL XX2 = Jac21*xx1 + Jac22*xx2 + Jac23*xx3;
        
        CCTK_REAL XX3 = Jac31*xx1 + Jac32*xx2 + Jac33*xx3;
        
        CCTK_REAL X = XX1;
        
        CCTK_REAL Y = XX2;
        
        CCTK_REAL Z = XX3;
        
        CCTK_REAL rXYZ = INV(sqrt(2))*sqrt(pow(4*SQR(Z)*SQR(ToReal(a)) + 
          SQR(Power(X,2) + Power(Y,2) + Power(Z,2) - Power(ToReal(a),2)),0.5) + 
          SQR(X) + SQR(Y) + SQR(Z) - SQR(ToReal(a)));
        
        CCTK_REAL d100rXYZ = X*INV(sqrt(2))*INV(sqrt(SQR(pow(X,2) + pow(Y,2) + 
          pow(Z,2) - pow(ToReal(a),2)) + 
          4*SQR(Z)*SQR(ToReal(a))))*sqrt(pow(4*SQR(Z)*SQR(ToReal(a)) + 
          SQR(Power(X,2) + Power(Y,2) + Power(Z,2) - Power(ToReal(a),2)),0.5) + 
          SQR(X) + SQR(Y) + SQR(Z) - SQR(ToReal(a)));
        
        CCTK_REAL d010rXYZ = Y*INV(sqrt(2))*INV(sqrt(SQR(pow(X,2) + pow(Y,2) + 
          pow(Z,2) - pow(ToReal(a),2)) + 
          4*SQR(Z)*SQR(ToReal(a))))*sqrt(pow(4*SQR(Z)*SQR(ToReal(a)) + 
          SQR(Power(X,2) + Power(Y,2) + Power(Z,2) - Power(ToReal(a),2)),0.5) + 
          SQR(X) + SQR(Y) + SQR(Z) - SQR(ToReal(a)));
        
        CCTK_REAL d001rXYZ = Z*INV(sqrt(2))*INV(sqrt(SQR(X) + SQR(Y) + SQR(Z) 
          - SQR(ToReal(a)) + sqrt(SQR(pow(X,2) + pow(Y,2) + pow(Z,2) - 
          pow(ToReal(a),2)) + 4*SQR(Z)*SQR(ToReal(a)))))*(1 + 
          INV(sqrt(SQR(pow(X,2) + pow(Y,2) + pow(Z,2) - pow(ToReal(a),2)) + 
          4*SQR(Z)*SQR(ToReal(a))))*(SQR(X) + SQR(Y) + SQR(Z) + SQR(ToReal(a))));
        
        CCTK_REAL G11 = 1 + 2*CUB(rXYZ)*INV(SQR(SQR(rXYZ) + 
          SQR(ToReal(a))))*INV(QAD(rXYZ) + SQR(Z)*SQR(ToReal(a)))*SQR(rXYZ*X + 
          Y*ToReal(a))*ToReal(M);
        
        CCTK_REAL G21 = 2*CUB(rXYZ)*INV(SQR(SQR(rXYZ) + 
          SQR(ToReal(a))))*INV(QAD(rXYZ) + SQR(Z)*SQR(ToReal(a)))*(rXYZ*Y - 
          X*ToReal(a))*(rXYZ*X + Y*ToReal(a))*ToReal(M);
        
        CCTK_REAL G31 = 2*Z*INV(SQR(rXYZ) + SQR(ToReal(a)))*INV(QAD(rXYZ) + 
          SQR(Z)*SQR(ToReal(a)))*SQR(rXYZ)*(rXYZ*X + Y*ToReal(a))*ToReal(M);
        
        CCTK_REAL G22 = 1 + 2*CUB(rXYZ)*INV(SQR(SQR(rXYZ) + 
          SQR(ToReal(a))))*INV(QAD(rXYZ) + SQR(Z)*SQR(ToReal(a)))*SQR(-(rXYZ*Y) + 
          X*ToReal(a))*ToReal(M);
        
        CCTK_REAL G32 = 2*Z*INV(SQR(rXYZ) + SQR(ToReal(a)))*INV(QAD(rXYZ) + 
          SQR(Z)*SQR(ToReal(a)))*SQR(rXYZ)*(rXYZ*Y - X*ToReal(a))*ToReal(M);
        
        CCTK_REAL G33 = 1 + 2*rXYZ*INV(QAD(rXYZ) + 
          SQR(Z)*SQR(ToReal(a)))*SQR(Z)*ToReal(M);
        
        CCTK_REAL K11 = 2*INV(SQR(QAD(rXYZ) + 
          SQR(Z)*SQR(ToReal(a))))*INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*pow(SQR(rXYZ) + 
          SQR(ToReal(a)),-3)*SQR(rXYZ)*sqrt(INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(-SQR(rXYZ) + SQR(X) + SQR(Y) + 
          SQR(Z))*ToReal(M)) + rXYZ*SQR(ToReal(a))*(CUB(rXYZ) - 
          2*SQR(rXYZ)*ToReal(M) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*(QAD(ToReal(a))*SQR(Z) + CUB(rXYZ)*(CUB(rXYZ) + 
          2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M)))))*ToReal(M)*((rXYZ*(rXYZ + 
          4*d100rXYZ*X)*pow(ToReal(a),10) + 
          3*d100rXYZ*Y*pow(ToReal(a),11))*QAD(Z) + 
          pow(ToReal(a),6)*QAD(rXYZ)*(pow(rXYZ,6) + Z*SQR(rXYZ)*(3*CUB(Z) + 
          (2*d100rXYZ*X*Z - d001rXYZ*SQR(Y))*ToReal(M)) + 
          rXYZ*(8*d100rXYZ*X*QAD(Z) - 6*SQR(X)*SQR(Z)*ToReal(M)) + 
          SQR(Z)*(8*d100rXYZ*X*CUB(rXYZ) + 6*QAD(rXYZ) + (d010rXYZ*(3*CUB(Y) - 
          8*Y*SQR(X)) + d001rXYZ*Z*(5*SQR(X) + 2*SQR(Y)) + d100rXYZ*X*(11*SQR(Y) 
          + 6*SQR(Z)))*ToReal(M))) + pow(rXYZ,12)*(QAD(rXYZ) + 
          X*(-2*d100rXYZ*CUB(rXYZ) - (3*X*(d010rXYZ*Y + d001rXYZ*Z) + 
          d100rXYZ*(SQR(X) - 2*(SQR(Y) + SQR(Z))))*ToReal(M))) + 
          pow(rXYZ,8)*SQR(ToReal(a))*(3*pow(rXYZ,6) + (4*CUB(rXYZ)*(-SQR(X) + 
          SQR(Y)) + SQR(rXYZ)*(d010rXYZ*(-5*CUB(Y) + 9*Y*SQR(X)) + 
          d001rXYZ*Z*(-2*SQR(X) - 5*SQR(Y)) + d100rXYZ*X*(3*SQR(X) - 11*SQR(Y) + 
          6*SQR(Z))))*ToReal(M) + SQR(Z)*(2*QAD(rXYZ) - 2*rXYZ*SQR(X)*ToReal(M)) 
          + X*(-4*d100rXYZ*pow(rXYZ,5) + SQR(Z)*(X*(d010rXYZ*Y + d001rXYZ*Z) + 
          d100rXYZ*(3*SQR(X) + 2*(SQR(Y) + SQR(Z))))*ToReal(M))) + 
          pow(rXYZ,6)*QAD(ToReal(a))*(-2*d100rXYZ*X*pow(rXYZ,5) + 3*pow(rXYZ,6) + 
          4*CUB(rXYZ)*(d100rXYZ*X*SQR(Z) + (-SQR(X) + SQR(Y))*ToReal(M)) + 
          SQR(rXYZ)*(QAD(Z) + (-(d010rXYZ*CUB(Y)) + d001rXYZ*Z*SQR(X) - 
          6*d001rXYZ*Z*SQR(Y) - d100rXYZ*X*(SQR(Y) - 6*SQR(Z)))*ToReal(M)) + 
          SQR(Z)*(6*QAD(rXYZ) + (d010rXYZ*(-CUB(Y) + 5*Y*SQR(X)) + 
          d001rXYZ*Z*(6*SQR(X) - SQR(Y)) + d100rXYZ*X*(7*SQR(X) + SQR(Y) + 
          6*SQR(Z)))*ToReal(M) + 2*rXYZ*(d100rXYZ*X*SQR(Z) + (-4*SQR(X) + 
          SQR(Y))*ToReal(M)))) + 
          SQR(Z)*(pow(ToReal(a),8)*SQR(rXYZ)*(4*d100rXYZ*X*CUB(rXYZ) + 
          2*QAD(rXYZ) + 3*SQR(rXYZ)*SQR(Z) + Z*(2*d100rXYZ*X*Z + 
          3*d001rXYZ*SQR(Y))*ToReal(M) - 2*rXYZ*(-5*d100rXYZ*X*SQR(Z) + 
          SQR(Y)*ToReal(M))) + d100rXYZ*rXYZ*Y*pow(ToReal(a),9)*(7*rXYZ*SQR(Z) + 
          2*(CUB(rXYZ) + SQR(Z)*ToReal(M)))) + 
          pow(rXYZ,5)*pow(ToReal(a),5)*(X*(d010rXYZ*SQR(rXYZ)*SQR(Y) + 
          (8*d001rXYZ*Y*Z + d010rXYZ*(-5*SQR(X) + 9*SQR(Y)))*SQR(Z) + 
          Y*(-4*CUB(rXYZ) - 8*rXYZ*SQR(Z)))*ToReal(M) + 
          d100rXYZ*Y*(-5*pow(rXYZ,5) + rXYZ*QAD(Z) - SQR(rXYZ)*(SQR(Y) - 
          6*SQR(Z))*ToReal(M) + SQR(Z)*(-2*CUB(rXYZ) + (15*SQR(X) + SQR(Y) + 
          6*SQR(Z))*ToReal(M)))) + pow(rXYZ,11)*ToReal(a)*(d010rXYZ*X*(3*SQR(X) - 
          8*SQR(Y))*ToReal(M) + Y*(4*X*(rXYZ - 2*d001rXYZ*Z)*ToReal(M) + 
          d100rXYZ*(-3*CUB(rXYZ) + (-9*SQR(X) + 2*(SQR(Y) + SQR(Z)))*ToReal(M)))) 
          + Y*CUB(rXYZ)*pow(ToReal(a),7)*(X*(-8*rXYZ - 3*d010rXYZ*Y + 
          8*d001rXYZ*Z)*SQR(Z)*ToReal(M) + d100rXYZ*(-pow(rXYZ,5) + 5*rXYZ*QAD(Z) 
          + SQR(Z)*(3*(SQR(Y) + 2*SQR(Z))*ToReal(M) + 2*(CUB(rXYZ) + 
          SQR(rXYZ)*ToReal(M))))) - 
          CUB(ToReal(a))*pow(rXYZ,7)*(d010rXYZ*X*(SQR(rXYZ)*(SQR(X) - 5*SQR(Y)) + 
          SQR(X)*SQR(Z))*ToReal(M) + Y*(8*d001rXYZ*X*Z*SQR(rXYZ)*ToReal(M) + 
          d100rXYZ*(7*pow(rXYZ,5) - 3*SQR(rXYZ)*(SQR(X) - SQR(Y) + 
          2*SQR(Z))*ToReal(M) + SQR(Z)*(2*CUB(rXYZ) - (3*SQR(X) + 2*(SQR(Y) + 
          SQR(Z)))*ToReal(M))))));
        
        CCTK_REAL K21 = -(INV(SQR(QAD(rXYZ) + 
          SQR(Z)*SQR(ToReal(a))))*INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*pow(SQR(rXYZ) + 
          SQR(ToReal(a)),-3)*SQR(rXYZ)*sqrt(INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(-SQR(rXYZ) + SQR(X) + SQR(Y) + 
          SQR(Z))*ToReal(M)) + rXYZ*SQR(ToReal(a))*(CUB(rXYZ) - 
          2*SQR(rXYZ)*ToReal(M) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*(QAD(ToReal(a))*SQR(Z) + CUB(rXYZ)*(CUB(rXYZ) + 
          2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M)))))*ToReal(M)*((-4*rXYZ*(d010rXYZ*X + 
          d100rXYZ*Y)*pow(ToReal(a),10) + 3*(d100rXYZ*X - 
          d010rXYZ*Y)*pow(ToReal(a),11))*QAD(Z) + 
          pow(rXYZ,11)*ToReal(a)*(3*d010rXYZ*Y*CUB(rXYZ) - 2*((-2*rXYZ + 
          4*d001rXYZ*Z)*(SQR(X) - SQR(Y)) + d010rXYZ*Y*(8*SQR(X) - 3*SQR(Y) + 
          SQR(Z)))*ToReal(M) + d100rXYZ*X*(-3*CUB(rXYZ) + 2*(-3*SQR(X) + 8*SQR(Y) 
          + SQR(Z))*ToReal(M))) + 
          pow(rXYZ,5)*pow(ToReal(a),5)*(d010rXYZ*rXYZ*Y*(5*QAD(rXYZ) - QAD(Z) + 
          2*SQR(rXYZ)*SQR(Z)) - 2*(-(d010rXYZ*Y*SQR(rXYZ)*(SQR(X) - 3*SQR(Z))) + 
          (SQR(X) - SQR(Y))*(2*CUB(rXYZ) + 4*rXYZ*SQR(Z)) + 
          SQR(Z)*(4*d001rXYZ*Z*(-SQR(X) + SQR(Y)) + d010rXYZ*(5*CUB(Y) + 
          Y*(-9*SQR(X) + 3*SQR(Z)))))*ToReal(M) + d100rXYZ*X*(-5*pow(rXYZ,5) + 
          rXYZ*QAD(Z) + 2*SQR(Z)*(5*SQR(X) - 9*SQR(Y) + 3*SQR(Z))*ToReal(M) - 
          2*(CUB(rXYZ)*SQR(Z) + SQR(rXYZ)*(SQR(Y) - 3*SQR(Z))*ToReal(M)))) + 
          CUB(ToReal(a))*pow(rXYZ,7)*(8*d001rXYZ*Z*SQR(rXYZ)*(-SQR(X) + 
          SQR(Y))*ToReal(M) + d010rXYZ*Y*(7*pow(rXYZ,5) - 2*SQR(Z)*(SQR(Y) + 
          SQR(Z))*ToReal(M) + 2*(CUB(rXYZ)*SQR(Z) + SQR(rXYZ)*(5*SQR(X) - SQR(Y) 
          - 3*SQR(Z))*ToReal(M))) + d100rXYZ*X*(-7*pow(rXYZ,5) + 
          2*SQR(rXYZ)*(SQR(X) - 5*SQR(Y) + 3*SQR(Z))*ToReal(M) + 
          SQR(Z)*(-2*CUB(rXYZ) + 2*(SQR(X) + SQR(Z))*ToReal(M)))) + 
          CUB(rXYZ)*pow(ToReal(a),7)*(8*(-rXYZ + d001rXYZ*Z)*(SQR(X) - 
          SQR(Y))*SQR(Z)*ToReal(M) + d010rXYZ*Y*(pow(rXYZ,5) - 5*rXYZ*QAD(Z) + 
          SQR(Z)*(-6*(SQR(X) + SQR(Z))*ToReal(M) - 2*(CUB(rXYZ) + 
          SQR(rXYZ)*ToReal(M)))) + d100rXYZ*X*(-pow(rXYZ,5) + 5*rXYZ*QAD(Z) + 
          SQR(Z)*(6*(SQR(Y) + SQR(Z))*ToReal(M) + 2*(CUB(rXYZ) + 
          SQR(rXYZ)*ToReal(M))))) - 
          2*Z*pow(ToReal(a),6)*QAD(rXYZ)*(d001rXYZ*X*Y*(SQR(rXYZ) + 
          3*SQR(Z))*ToReal(M) + Z*(d010rXYZ*X*(4*rXYZ*(SQR(rXYZ) + SQR(Z)) + 
          (SQR(rXYZ) + 4*SQR(X) - 7*SQR(Y) + 3*SQR(Z))*ToReal(M)) + 
          Y*(-6*rXYZ*X*ToReal(M) + d100rXYZ*(4*rXYZ*(SQR(rXYZ) + SQR(Z)) + 
          (SQR(rXYZ) - 7*SQR(X) + 4*SQR(Y) + 3*SQR(Z))*ToReal(M))))) + 
          2*(pow(rXYZ,12)*(d010rXYZ*X*(CUB(rXYZ) - (SQR(X) - 2*SQR(Y) + 
          SQR(Z))*ToReal(M)) + Y*(3*d001rXYZ*X*Z*ToReal(M) + d100rXYZ*(CUB(rXYZ) 
          + (2*SQR(X) - SQR(Y) - SQR(Z))*ToReal(M)))) + 
          pow(rXYZ,6)*QAD(ToReal(a))*(-(d010rXYZ*X*(-pow(rXYZ,5) + rXYZ*QAD(Z) + 
          2*CUB(rXYZ)*SQR(Z) + SQR(rXYZ)*(SQR(Y) + 3*SQR(Z))*ToReal(M) + 
          SQR(Z)*(SQR(X) + 7*SQR(Y) + 3*SQR(Z))*ToReal(M))) - Y*(X*(-8*CUB(rXYZ) 
          + 7*d001rXYZ*CUB(Z) + 7*d001rXYZ*Z*SQR(rXYZ) - 
          10*rXYZ*SQR(Z))*ToReal(M) + d100rXYZ*(-pow(rXYZ,5) + rXYZ*QAD(Z) + 
          2*CUB(rXYZ)*SQR(Z) + SQR(rXYZ)*(SQR(X) + 3*SQR(Z))*ToReal(M) + 
          SQR(Z)*(7*SQR(X) + SQR(Y) + 3*SQR(Z))*ToReal(M)))) + 
          pow(rXYZ,8)*SQR(ToReal(a))*(d010rXYZ*X*(2*pow(rXYZ,5) + 
          (SQR(rXYZ)*(3*SQR(X) - 11*SQR(Y) - 3*SQR(Z)) - SQR(Z)*(SQR(X) + 
          2*SQR(Y) + SQR(Z)))*ToReal(M)) + Y*(X*(8*CUB(rXYZ) + d001rXYZ*(-CUB(Z) 
          - 3*Z*SQR(rXYZ)) + 2*rXYZ*SQR(Z))*ToReal(M) + d100rXYZ*(2*pow(rXYZ,5) + 
          (SQR(rXYZ)*(-11*SQR(X) + 3*SQR(Y) - 3*SQR(Z)) - SQR(Z)*(2*SQR(X) + 
          SQR(Y) + SQR(Z)))*ToReal(M))))) + SQR(Z)*(rXYZ*(d100rXYZ*X - 
          d010rXYZ*Y)*pow(ToReal(a),9)*(7*rXYZ*SQR(Z) + 2*(CUB(rXYZ) + 
          SQR(Z)*ToReal(M))) - 
          2*pow(ToReal(a),8)*SQR(rXYZ)*(d010rXYZ*X*(2*CUB(rXYZ) + SQR(Z)*(5*rXYZ 
          + ToReal(M))) + Y*(X*(2*rXYZ - 3*d001rXYZ*Z)*ToReal(M) + 
          d100rXYZ*(2*CUB(rXYZ) + SQR(Z)*(5*rXYZ + ToReal(M))))))));
        
        CCTK_REAL K31 = -(rXYZ*INV(SQR(SQR(rXYZ) + 
          SQR(ToReal(a))))*INV(SQR(QAD(rXYZ) + 
          SQR(Z)*SQR(ToReal(a))))*INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*sqrt(INV(QAD(ToReal(a))*SQR(Z) + CUB(rXYZ)*(CUB(rXYZ) + 
          2*(-SQR(rXYZ) + SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) - 2*SQR(rXYZ)*ToReal(M) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*(QAD(ToReal(a))*SQR(Z) + CUB(rXYZ)*(CUB(rXYZ) + 
          2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M)))))*ToReal(M)*(-2*d100rXYZ*pow(Z,5)*pow(ToReal(a),10) + 
          pow(rXYZ,10)*ToReal(a)*(2*Z*((-2*rXYZ + 7*d100rXYZ*X)*Y + 
          d010rXYZ*(-3*SQR(X) + 4*SQR(Y)))*ToReal(M) + d001rXYZ*Y*(3*CUB(rXYZ) - 
          2*(SQR(X) + SQR(Y) - 3*SQR(Z))*ToReal(M))) + rXYZ*CUB(Z)*(Y*(2*rXYZ - 
          3*d001rXYZ*Z)*pow(ToReal(a),9) + 
          2*pow(ToReal(a),8)*(-(rXYZ*Z*(2*d001rXYZ*X + 3*d100rXYZ*Z)) + 
          X*SQR(rXYZ) - d100rXYZ*SQR(Z)*ToReal(M))) + 
          pow(ToReal(a),5)*QAD(rXYZ)*(d001rXYZ*Y*(pow(rXYZ,5) + QAD(Z)*(-rXYZ - 
          8*ToReal(M)) + 2*SQR(rXYZ)*SQR(Z)*ToReal(M)) + 2*Z*((-5*d100rXYZ*X*Y + 
          d010rXYZ*(3*SQR(X) - 2*SQR(Y)))*SQR(Z)*ToReal(M) + Y*(SQR(rXYZ)*SQR(Z) 
          + 2*(QAD(rXYZ) + rXYZ*SQR(Z)*ToReal(M))))) + 
          2*(pow(rXYZ,11)*(d001rXYZ*X*(CUB(rXYZ) - (SQR(X) + SQR(Y) - 
          2*SQR(Z))*ToReal(M)) + Z*(3*d010rXYZ*X*Y*ToReal(M) + 
          d100rXYZ*(CUB(rXYZ) + (2*SQR(X) - SQR(Y) - SQR(Z))*ToReal(M)))) + 
          CUB(ToReal(a))*pow(rXYZ,6)*(-(d001rXYZ*Y*(-2*pow(rXYZ,5) + 
          SQR(Z)*(-CUB(rXYZ) + (-4*SQR(rXYZ) + SQR(X) + SQR(Y) + 
          SQR(Z))*ToReal(M)))) + Z*(d010rXYZ*(-(SQR(rXYZ)*(SQR(X) - 2*SQR(Y))) + 
          SQR(X)*SQR(Z))*ToReal(M) + Y*(QAD(rXYZ) - (2*CUB(rXYZ) + 
          d100rXYZ*X*(-3*SQR(rXYZ) + SQR(Z)))*ToReal(M)))) + 
          pow(rXYZ,7)*SQR(ToReal(a))*(d001rXYZ*X*(pow(rXYZ,5) + (2*SQR(rXYZ) - 
          SQR(X) - SQR(Y) - 2*SQR(Z))*SQR(Z)*ToReal(M)) - 
          Z*(d100rXYZ*(-3*pow(rXYZ,5) + (3*SQR(rXYZ)*(-SQR(Y) + SQR(Z)) + 
          SQR(Z)*(2*SQR(X) + SQR(Y) + SQR(Z)))*ToReal(M)) - X*(QAD(rXYZ) + 
          (d010rXYZ*Y*(-3*SQR(rXYZ) - SQR(Z)) + 2*(CUB(rXYZ) + 
          rXYZ*SQR(Z)))*ToReal(M)))) + 
          Z*(pow(rXYZ,5)*QAD(ToReal(a))*(X*(rXYZ*(rXYZ - d001rXYZ*Z)*(2*SQR(rXYZ) 
          + SQR(Z)) + (2*CUB(rXYZ) - 2*d010rXYZ*Y*SQR(rXYZ) + (6*rXYZ - 
          3*(d010rXYZ*Y + 2*d001rXYZ*Z))*SQR(Z))*ToReal(M)) + 
          d100rXYZ*(3*pow(rXYZ,5) - rXYZ*QAD(Z) + (SQR(rXYZ)*(2*SQR(Y) - 
          3*SQR(Z)) - SQR(Z)*(4*SQR(X) + SQR(Y) + 3*SQR(Z)))*ToReal(M))) + 
          Y*pow(ToReal(a),7)*SQR(rXYZ)*(QAD(rXYZ) + d001rXYZ*(-(Z*CUB(rXYZ)) - 
          3*CUB(Z)*ToReal(M)) + 2*SQR(Z)*(SQR(rXYZ) + rXYZ*(-(d001rXYZ*Z) + 
          ToReal(M)))) + CUB(rXYZ)*pow(ToReal(a),6)*(d100rXYZ*(pow(rXYZ,5) + 
          (-SQR(rXYZ) - 2*SQR(Y))*SQR(Z)*ToReal(M) - 3*QAD(Z)*(rXYZ + ToReal(M))) 
          + X*(-2*d001rXYZ*Z*CUB(rXYZ) + QAD(rXYZ) + SQR(Z)*(rXYZ*(-3*d001rXYZ*Z 
          + 4*ToReal(M)) + 2*(SQR(rXYZ) + (d010rXYZ*Y - 
          2*d001rXYZ*Z)*ToReal(M)))))))));
        
        CCTK_REAL K22 = 2*INV(SQR(QAD(rXYZ) + 
          SQR(Z)*SQR(ToReal(a))))*INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*pow(SQR(rXYZ) + 
          SQR(ToReal(a)),-3)*SQR(rXYZ)*sqrt(INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(-SQR(rXYZ) + SQR(X) + SQR(Y) + 
          SQR(Z))*ToReal(M)) + rXYZ*SQR(ToReal(a))*(CUB(rXYZ) - 
          2*SQR(rXYZ)*ToReal(M) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*(QAD(ToReal(a))*SQR(Z) + CUB(rXYZ)*(CUB(rXYZ) + 
          2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M)))))*ToReal(M)*((rXYZ*(rXYZ + 
          4*d010rXYZ*Y)*pow(ToReal(a),10) - 
          3*d010rXYZ*X*pow(ToReal(a),11))*QAD(Z) + 
          pow(ToReal(a),6)*QAD(rXYZ)*(pow(rXYZ,6) + Z*SQR(rXYZ)*(3*CUB(Z) + 
          (2*d010rXYZ*Y*Z - d001rXYZ*SQR(X))*ToReal(M)) + 
          rXYZ*(8*d010rXYZ*Y*QAD(Z) - 6*SQR(Y)*SQR(Z)*ToReal(M)) + 
          SQR(Z)*(8*d010rXYZ*Y*CUB(rXYZ) + 6*QAD(rXYZ) + ((11*d010rXYZ*Y + 
          2*d001rXYZ*Z)*SQR(X) + 5*d001rXYZ*Z*SQR(Y) + d100rXYZ*(3*CUB(X) - 
          8*X*SQR(Y)) + 6*d010rXYZ*Y*SQR(Z))*ToReal(M))) + 
          pow(rXYZ,11)*ToReal(a)*(Y*(X*(-4*rXYZ + 8*d001rXYZ*Z) + 
          d100rXYZ*(8*SQR(X) - 3*SQR(Y)))*ToReal(M) + d010rXYZ*X*(3*CUB(rXYZ) + 
          (9*SQR(Y) - 2*(SQR(X) + SQR(Z)))*ToReal(M))) + pow(rXYZ,12)*(QAD(rXYZ) 
          + Y*(-2*d010rXYZ*CUB(rXYZ) - (3*Y*(d100rXYZ*X + d001rXYZ*Z) + 
          d010rXYZ*(SQR(Y) - 2*(SQR(X) + SQR(Z))))*ToReal(M))) + 
          pow(rXYZ,8)*SQR(ToReal(a))*(3*pow(rXYZ,6) + (4*CUB(rXYZ)*(SQR(X) - 
          SQR(Y)) + SQR(rXYZ)*(-5*(d100rXYZ*CUB(X) + d001rXYZ*Z*SQR(X)) + 
          (9*d100rXYZ*X - 2*d001rXYZ*Z)*SQR(Y) + d010rXYZ*(3*CUB(Y) + 
          Y*(-11*SQR(X) + 6*SQR(Z)))))*ToReal(M) + SQR(Z)*(2*QAD(rXYZ) - 
          2*rXYZ*SQR(Y)*ToReal(M)) + Y*(-4*d010rXYZ*pow(rXYZ,5) + 
          SQR(Z)*(Y*(d100rXYZ*X + d001rXYZ*Z) + d010rXYZ*(3*SQR(Y) + 2*(SQR(X) + 
          SQR(Z))))*ToReal(M))) + 
          pow(rXYZ,6)*QAD(ToReal(a))*(-2*d010rXYZ*Y*pow(rXYZ,5) + 3*pow(rXYZ,6) + 
          4*CUB(rXYZ)*(d010rXYZ*Y*SQR(Z) + (SQR(X) - SQR(Y))*ToReal(M)) + 
          SQR(rXYZ)*(QAD(Z) + (-(d100rXYZ*CUB(X)) - 6*d001rXYZ*Z*SQR(X) + 
          d001rXYZ*Z*SQR(Y) - d010rXYZ*Y*(SQR(X) - 6*SQR(Z)))*ToReal(M)) + 
          SQR(Z)*(6*QAD(rXYZ) + (-(d001rXYZ*Z*(SQR(X) - 6*SQR(Y))) - 
          d100rXYZ*(CUB(X) - 5*X*SQR(Y)) + d010rXYZ*Y*(SQR(X) + 7*SQR(Y) + 
          6*SQR(Z)))*ToReal(M) + 2*rXYZ*(d010rXYZ*Y*SQR(Z) + (SQR(X) - 
          4*SQR(Y))*ToReal(M)))) + 
          SQR(Z)*(pow(ToReal(a),8)*SQR(rXYZ)*(4*d010rXYZ*Y*CUB(rXYZ) + 
          2*QAD(rXYZ) + 3*SQR(rXYZ)*SQR(Z) + Z*(2*d010rXYZ*Y*Z + 
          3*d001rXYZ*SQR(X))*ToReal(M) - 2*rXYZ*(-5*d010rXYZ*Y*SQR(Z) + 
          SQR(X)*ToReal(M))) - d010rXYZ*rXYZ*X*pow(ToReal(a),9)*(7*rXYZ*SQR(Z) + 
          2*(CUB(rXYZ) + SQR(Z)*ToReal(M)))) + 
          pow(rXYZ,5)*pow(ToReal(a),5)*(Y*(-(d100rXYZ*SQR(rXYZ)*SQR(X)) + 
          (-8*d001rXYZ*X*Z + d100rXYZ*(-9*SQR(X) + 5*SQR(Y)))*SQR(Z) + 
          X*(4*CUB(rXYZ) + 8*rXYZ*SQR(Z)))*ToReal(M) + d010rXYZ*X*(5*pow(rXYZ,5) 
          - rXYZ*QAD(Z) + SQR(rXYZ)*(SQR(X) - 6*SQR(Z))*ToReal(M) + 
          SQR(Z)*(2*CUB(rXYZ) - (SQR(X) + 15*SQR(Y) + 6*SQR(Z))*ToReal(M)))) + 
          CUB(ToReal(a))*pow(rXYZ,7)*(Y*(8*d001rXYZ*X*Z*SQR(rXYZ) + 
          d100rXYZ*(SQR(rXYZ)*(-5*SQR(X) + SQR(Y)) + SQR(Y)*SQR(Z)))*ToReal(M) + 
          d010rXYZ*X*(7*pow(rXYZ,5) + 3*SQR(rXYZ)*(SQR(X) - SQR(Y) - 
          2*SQR(Z))*ToReal(M) + SQR(Z)*(2*CUB(rXYZ) - (3*SQR(Y) + 2*(SQR(X) + 
          SQR(Z)))*ToReal(M)))) + X*CUB(rXYZ)*pow(ToReal(a),7)*(Y*(8*rXYZ + 
          3*d100rXYZ*X - 8*d001rXYZ*Z)*SQR(Z)*ToReal(M) + d010rXYZ*(pow(rXYZ,5) - 
          5*rXYZ*QAD(Z) + SQR(Z)*(-3*(SQR(X) + 2*SQR(Z))*ToReal(M) - 2*(CUB(rXYZ) 
          + SQR(rXYZ)*ToReal(M))))));
        
        CCTK_REAL K32 = -(rXYZ*INV(SQR(SQR(rXYZ) + 
          SQR(ToReal(a))))*INV(SQR(QAD(rXYZ) + 
          SQR(Z)*SQR(ToReal(a))))*INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*sqrt(INV(QAD(ToReal(a))*SQR(Z) + CUB(rXYZ)*(CUB(rXYZ) + 
          2*(-SQR(rXYZ) + SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) - 2*SQR(rXYZ)*ToReal(M) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*(QAD(ToReal(a))*SQR(Z) + CUB(rXYZ)*(CUB(rXYZ) + 
          2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M)))))*ToReal(M)*(pow(rXYZ,10)*ToReal(a)*(2*Z*(X*(2*rXYZ - 
          7*d010rXYZ*Y) + d100rXYZ*(-4*SQR(X) + 3*SQR(Y)))*ToReal(M) + 
          d001rXYZ*X*(-3*CUB(rXYZ) + 2*(SQR(X) + SQR(Y) - 3*SQR(Z))*ToReal(M))) + 
          rXYZ*CUB(Z)*(X*(-2*rXYZ + 3*d001rXYZ*Z)*pow(ToReal(a),9) + 
          2*pow(ToReal(a),8)*(-(rXYZ*Z*(2*d001rXYZ*Y + 3*d010rXYZ*Z)) + 
          Y*SQR(rXYZ) - d010rXYZ*SQR(Z)*ToReal(M))) - 
          pow(ToReal(a),5)*QAD(rXYZ)*(d001rXYZ*X*(pow(rXYZ,5) + QAD(Z)*(-rXYZ - 
          8*ToReal(M)) + 2*SQR(rXYZ)*SQR(Z)*ToReal(M)) + 2*Z*((-5*d010rXYZ*X*Y + 
          d100rXYZ*(-2*SQR(X) + 3*SQR(Y)))*SQR(Z)*ToReal(M) + X*(SQR(rXYZ)*SQR(Z) 
          + 2*(QAD(rXYZ) + rXYZ*SQR(Z)*ToReal(M))))) - 
          2*(d010rXYZ*pow(Z,5)*pow(ToReal(a),10) + 
          X*Z*pow(ToReal(a),7)*SQR(rXYZ)*(QAD(rXYZ) + d001rXYZ*(-(Z*CUB(rXYZ)) - 
          3*CUB(Z)*ToReal(M)) + 2*SQR(Z)*(SQR(rXYZ) + rXYZ*(-(d001rXYZ*Z) + 
          ToReal(M)))) + CUB(ToReal(a))*pow(rXYZ,6)*(Z*(X*QAD(rXYZ) + 
          (-2*X*CUB(rXYZ) + SQR(rXYZ)*(3*d010rXYZ*X*Y + d100rXYZ*(2*SQR(X) - 
          SQR(Y))) + Y*(-(d010rXYZ*X) + d100rXYZ*Y)*SQR(Z))*ToReal(M)) - 
          d001rXYZ*X*(-2*pow(rXYZ,5) + SQR(Z)*(-CUB(rXYZ) + (-4*SQR(rXYZ) + 
          SQR(X) + SQR(Y) + SQR(Z))*ToReal(M))))) + 
          2*(pow(rXYZ,11)*(d001rXYZ*Y*(CUB(rXYZ) - (SQR(X) + SQR(Y) - 
          2*SQR(Z))*ToReal(M)) + Z*(3*d100rXYZ*X*Y*ToReal(M) + 
          d010rXYZ*(CUB(rXYZ) - (SQR(X) - 2*SQR(Y) + SQR(Z))*ToReal(M)))) + 
          pow(rXYZ,7)*SQR(ToReal(a))*(d001rXYZ*Y*(pow(rXYZ,5) + (2*SQR(rXYZ) - 
          SQR(X) - SQR(Y) - 2*SQR(Z))*SQR(Z)*ToReal(M)) - 
          Z*(d010rXYZ*(-3*pow(rXYZ,5) + (3*SQR(rXYZ)*(-SQR(X) + SQR(Z)) + 
          SQR(Z)*(SQR(X) + 2*SQR(Y) + SQR(Z)))*ToReal(M)) - Y*(QAD(rXYZ) + 
          (d100rXYZ*X*(-3*SQR(rXYZ) - SQR(Z)) + 2*(CUB(rXYZ) + 
          rXYZ*SQR(Z)))*ToReal(M)))) + 
          Z*(pow(rXYZ,5)*QAD(ToReal(a))*(Y*(rXYZ*(rXYZ - d001rXYZ*Z)*(2*SQR(rXYZ) 
          + SQR(Z)) + (2*CUB(rXYZ) - 2*d100rXYZ*X*SQR(rXYZ) + (6*rXYZ - 
          3*(d100rXYZ*X + 2*d001rXYZ*Z))*SQR(Z))*ToReal(M)) + 
          d010rXYZ*(3*pow(rXYZ,5) - rXYZ*QAD(Z) + (SQR(rXYZ)*(2*SQR(X) - 
          3*SQR(Z)) - SQR(Z)*(SQR(X) + 4*SQR(Y) + 3*SQR(Z)))*ToReal(M))) + 
          CUB(rXYZ)*pow(ToReal(a),6)*(d010rXYZ*(pow(rXYZ,5) + (-SQR(rXYZ) - 
          2*SQR(X))*SQR(Z)*ToReal(M) - 3*QAD(Z)*(rXYZ + ToReal(M))) + 
          Y*(-2*d001rXYZ*Z*CUB(rXYZ) + QAD(rXYZ) + SQR(Z)*(rXYZ*(-3*d001rXYZ*Z + 
          4*ToReal(M)) + 2*(SQR(rXYZ) + (d100rXYZ*X - 
          2*d001rXYZ*Z)*ToReal(M)))))))));
        
        CCTK_REAL K33 = 2*rXYZ*INV(SQR(QAD(rXYZ) + 
          SQR(Z)*SQR(ToReal(a))))*INV(SQR(rXYZ) + 
          SQR(ToReal(a)))*INV(QAD(ToReal(a))*SQR(Z) + CUB(rXYZ)*(CUB(rXYZ) + 
          2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*sqrt(INV(QAD(ToReal(a))*SQR(Z) + CUB(rXYZ)*(CUB(rXYZ) + 
          2*(-SQR(rXYZ) + SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) - 2*SQR(rXYZ)*ToReal(M) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*(QAD(ToReal(a))*SQR(Z) + CUB(rXYZ)*(CUB(rXYZ) + 
          2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M)))))*ToReal(M)*((d010rXYZ*X - 
          d100rXYZ*Y)*SQR(Z)*(CUB(ToReal(a))*QAD(rXYZ)*(3*SQR(rXYZ) - SQR(Z)) + 
          3*pow(rXYZ,8)*ToReal(a))*ToReal(M) + 
          pow(rXYZ,5)*SQR(ToReal(a))*(-4*d001rXYZ*Z*pow(rXYZ,5) + 2*pow(rXYZ,6) + 
          (-2*rXYZ*QAD(Z) - (3*(d100rXYZ*X + d010rXYZ*Y) + 
          2*d001rXYZ*Z)*SQR(rXYZ)*SQR(Z) + CUB(Z)*((d100rXYZ*X + d010rXYZ*Y)*Z + 
          d001rXYZ*(2*(SQR(X) + SQR(Y)) + 3*SQR(Z))))*ToReal(M)) + QAD(Z)*((-rXYZ 
          + 2*d001rXYZ*Z)*pow(ToReal(a),8) + (-(d010rXYZ*X) + 
          d100rXYZ*Y)*pow(ToReal(a),5)*SQR(rXYZ)*ToReal(M) + 
          rXYZ*pow(ToReal(a),6)*(-2*SQR(rXYZ) + rXYZ*(4*d001rXYZ*Z - 2*ToReal(M)) 
          + 3*d001rXYZ*Z*ToReal(M))) + 
          CUB(rXYZ)*QAD(ToReal(a))*(-2*d001rXYZ*Z*pow(rXYZ,5) + pow(rXYZ,6) - 
          CUB(Z)*SQR(rXYZ)*(Z + d001rXYZ*ToReal(M)) + QAD(Z)*(2*rXYZ*(d001rXYZ*Z 
          - 2*ToReal(M)) + (d100rXYZ*X + d010rXYZ*Y + 6*d001rXYZ*Z)*ToReal(M))) + 
          pow(rXYZ,9)*(QAD(rXYZ) + Z*(-2*d001rXYZ*CUB(rXYZ) - (3*(d100rXYZ*X + 
          d010rXYZ*Y)*Z + d001rXYZ*(-2*(SQR(X) + SQR(Y)) + SQR(Z)))*ToReal(M))));
        
        CCTK_REAL betap1 = 2*CUB(rXYZ)*INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 2*ToReal(M))))*(rXYZ*X 
          + Y*ToReal(a))*ToReal(M);
        
        CCTK_REAL betap2 = 2*CUB(rXYZ)*INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 2*ToReal(M))))*(rXYZ*Y 
          - X*ToReal(a))*ToReal(M);
        
        CCTK_REAL betap3 = 2*Z*INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*SQR(rXYZ)*(SQR(rXYZ) + SQR(ToReal(a)))*ToReal(M);
        
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
        
        CCTK_REAL alpL = INV(sqrt(INV(QAD(ToReal(a))*SQR(Z) + 
          CUB(rXYZ)*(CUB(rXYZ) + 2*(-SQR(rXYZ) + SQR(X) + SQR(Y) + 
          SQR(Z))*ToReal(M)) + rXYZ*SQR(ToReal(a))*(CUB(rXYZ) - 
          2*SQR(rXYZ)*ToReal(M) + SQR(Z)*(rXYZ + 
          2*ToReal(M))))*(QAD(ToReal(a))*SQR(Z) + CUB(rXYZ)*(CUB(rXYZ) + 
          2*(SQR(X) + SQR(Y) + SQR(Z))*ToReal(M)) + 
          rXYZ*SQR(ToReal(a))*(CUB(rXYZ) + SQR(Z)*(rXYZ + 2*ToReal(M))))));
        
        CCTK_REAL betaxL = betap1*Jac11 + betap2*Jac21 + betap3*Jac31;
        
        CCTK_REAL betayL = betap1*Jac12 + betap2*Jac22 + betap3*Jac32;
        
        CCTK_REAL betazL = betap1*Jac13 + betap2*Jac23 + betap3*Jac33;
        
        CCTK_REAL dtalpL = 0;
        
        CCTK_REAL dtbetaxL = dtbetap1*Jac11 + dtbetap2*Jac21 + dtbetap3*Jac31;
        
        CCTK_REAL dtbetayL = dtbetap1*Jac12 + dtbetap2*Jac22 + dtbetap3*Jac32;
        
        CCTK_REAL dtbetazL = dtbetap1*Jac13 + dtbetap2*Jac23 + dtbetap3*Jac33;
        
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
    }
  }
}

extern "C" void KerrSchild_initial(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  GenericFD_LoopOverEverything(cctkGH, &KerrSchild_initial_Body);
}
