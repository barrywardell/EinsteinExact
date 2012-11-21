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
#include "vectors.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define ScalarINV(x) ((CCTK_REAL)1.0 / (x))
#define ScalarSQR(x) ((x) * (x))
#define ScalarCUB(x) ((x) * ScalarSQR(x))
#define ScalarQAD(x) (ScalarSQR(ScalarSQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))
#define QAD(x) (SQR(SQR(x)))

static void Vaidya2_initial_Body(cGH const * restrict const cctkGH, int const dir, int const face, CCTK_REAL const normal[3], CCTK_REAL const tangentA[3], CCTK_REAL const tangentB[3], int const imin[3], int const imax[3], int const n_subblock_gfs, CCTK_REAL * restrict const subblock_gfs[])
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
  CCTK_REAL_VEC const dx CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL_VEC const dy CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL_VEC const dz CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL_VEC const dt CCTK_ATTRIBUTE_UNUSED  = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL_VEC const t CCTK_ATTRIBUTE_UNUSED  = ToReal(cctk_time);
  CCTK_REAL_VEC const dxi CCTK_ATTRIBUTE_UNUSED  = INV(dx);
  CCTK_REAL_VEC const dyi CCTK_ATTRIBUTE_UNUSED  = INV(dy);
  CCTK_REAL_VEC const dzi CCTK_ATTRIBUTE_UNUSED  = INV(dz);
  CCTK_REAL_VEC const khalf CCTK_ATTRIBUTE_UNUSED  = ToReal(0.5);
  CCTK_REAL_VEC const kthird CCTK_ATTRIBUTE_UNUSED  = ToReal(1.0/3.0);
  CCTK_REAL_VEC const ktwothird CCTK_ATTRIBUTE_UNUSED  = ToReal(2.0/3.0);
  CCTK_REAL_VEC const kfourthird CCTK_ATTRIBUTE_UNUSED  = ToReal(4.0/3.0);
  CCTK_REAL_VEC const keightthird CCTK_ATTRIBUTE_UNUSED  = ToReal(8.0/3.0);
  CCTK_REAL_VEC const hdxi CCTK_ATTRIBUTE_UNUSED  = kmul(ToReal(0.5), dxi);
  CCTK_REAL_VEC const hdyi CCTK_ATTRIBUTE_UNUSED  = kmul(ToReal(0.5), dyi);
  CCTK_REAL_VEC const hdzi CCTK_ATTRIBUTE_UNUSED  = kmul(ToReal(0.5), dzi);
  
  /* Initialize predefined quantities */
  
  /* Assign local copies of arrays functions */
  
  
  
  /* Calculate temporaries and arrays functions */
  
  /* Copy local copies back to grid functions */
  
  /* Loop over the grid points */
  #pragma omp parallel
  LC_LOOP3VEC(Vaidya2_initial,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    CCTK_REAL_VEC_SIZE)
  {
    ptrdiff_t const index CCTK_ATTRIBUTE_UNUSED  = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC xL CCTK_ATTRIBUTE_UNUSED = vec_load(x[index]);
    CCTK_REAL_VEC yL CCTK_ATTRIBUTE_UNUSED = vec_load(y[index]);
    CCTK_REAL_VEC zL CCTK_ATTRIBUTE_UNUSED = vec_load(z[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED xx1 = xL;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED xx2 = yL;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED xx3 = zL;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED position1 = ToReal(positionx);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED position2 = ToReal(positiony);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED position3 = ToReal(positionz);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED shiftadd1 = ToReal(shiftaddx);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED shiftadd2 = ToReal(shiftaddy);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED shiftadd3 = ToReal(shiftaddz);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp0 = kcos(ToReal(phi));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp1 = kcos(ToReal(psi));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp2 = kcos(ToReal(theta));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp3 = ksin(ToReal(phi));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp4 = ksin(ToReal(psi));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED Jac11 = 
      kmsub(csetemp0,csetemp1,kmul(csetemp2,kmul(csetemp3,csetemp4)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED Jac12 = 
      kmadd(csetemp1,csetemp3,kmul(csetemp0,kmul(csetemp2,csetemp4)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp5 = ksin(ToReal(theta));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED Jac13 = kmul(csetemp4,csetemp5);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED Jac21 = 
      knmadd(csetemp1,kmul(csetemp2,csetemp3),kmul(csetemp0,csetemp4));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED Jac22 = 
      kmsub(csetemp0,kmul(csetemp1,csetemp2),kmul(csetemp3,csetemp4));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED Jac23 = kmul(csetemp1,csetemp5);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED Jac31 = kmul(csetemp3,csetemp5);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED Jac32 = 
      kneg(kmul(csetemp0,csetemp5));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED Jac33 = csetemp2;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED InvJac11 = Jac11;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED InvJac12 = Jac21;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED InvJac13 = Jac31;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED InvJac21 = Jac12;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED InvJac22 = Jac22;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED InvJac23 = Jac32;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED InvJac31 = Jac13;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED InvJac32 = Jac23;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED InvJac33 = Jac33;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED T = 
      kmul(ToReal(lapsefactor),ksub(t,ToReal(positiont)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp6 = 
      kneg(kmul(shiftadd1,T));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp7 = 
      kneg(kmul(shiftadd2,T));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp8 = 
      kneg(kmul(shiftadd3,T));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED XX1 = 
      kmadd(Jac11,kadd(csetemp6,ksub(xx1,position1)),kmadd(Jac12,kadd(csetemp7,ksub(xx2,position2)),kmul(Jac13,kadd(csetemp8,ksub(xx3,position3)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED XX2 = 
      kmadd(Jac21,kadd(csetemp6,ksub(xx1,position1)),kmadd(Jac22,kadd(csetemp7,ksub(xx2,position2)),kmul(Jac23,kadd(csetemp8,ksub(xx3,position3)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED XX3 = 
      kmadd(Jac31,kadd(csetemp6,ksub(xx1,position1)),kmadd(Jac32,kadd(csetemp7,ksub(xx2,position2)),kmul(Jac33,kadd(csetemp8,ksub(xx3,position3)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED X = XX1;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED Y = XX2;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED Z = XX3;
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp9 = kmul(X,X);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp10 = kmul(Y,Y);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp11 = kmul(Z,Z);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED rXYZ = 
      ksqrt(kadd(csetemp10,kadd(csetemp11,csetemp9)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp12 = ToReal(ScalarINV(M));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED mTXYZ = 
      kmul(kmadd(ktanh(kmul(csetemp12,kmul(kadd(T,ksqrt(kadd(csetemp10,kadd(csetemp11,csetemp9)))),ToReal(dM)))),ktanh(kmul(csetemp12,kmul(kadd(T,ksqrt(kadd(csetemp10,kadd(csetemp11,csetemp9)))),ToReal(dM)))),ToReal(1)),ToReal(M));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp13 = 
      kdiv(ToReal(1),kmul(rXYZ,kmul(rXYZ,rXYZ)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED G11 = 
      kmadd(csetemp13,kmul(csetemp9,kmul(mTXYZ,ToReal(2))),ToReal(1));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED G21 = 
      kmul(csetemp13,kmul(mTXYZ,kmul(X,kmul(Y,ToReal(2)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED G31 = 
      kmul(csetemp13,kmul(mTXYZ,kmul(X,kmul(Z,ToReal(2)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED G22 = 
      kmadd(csetemp10,kmul(csetemp13,kmul(mTXYZ,ToReal(2))),ToReal(1));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED G32 = 
      kmul(csetemp13,kmul(mTXYZ,kmul(Y,kmul(Z,ToReal(2)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED G33 = 
      kmadd(csetemp11,kmul(csetemp13,kmul(mTXYZ,ToReal(2))),ToReal(1));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp14 = 
      kdiv(ToReal(1),kmul(kmul(rXYZ,rXYZ),kmul(rXYZ,rXYZ)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp15 = 
      kmul(rXYZ,kmul(rXYZ,rXYZ));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp16 = 
      kmul(rXYZ,kmul(kmul(rXYZ,rXYZ),kmul(rXYZ,rXYZ)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp17 = kmul(mTXYZ,mTXYZ);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp18 = 
      kmul(kmul(rXYZ,rXYZ),kmul(rXYZ,rXYZ));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp19 = kadd(rXYZ,T);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED K11 = 
      kdiv(kmul(csetemp13,kmul(csetemp14,kmul(ToReal(-2),kmsub(csetemp9,kmsub(csetemp17,kadd(csetemp10,kadd(csetemp11,csetemp9)),kdiv(kmul(csetemp18,kmul(ktanh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM))))))),kmul(mTXYZ,kmadd(csetemp15,kmul(csetemp9,ToReal(-2)),csetemp16)))))),ksqrt(kmul(csetemp13,kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED K21 = 
      kdiv(kmul(csetemp13,kmul(csetemp14,kmul(X,kmul(Y,kmul(ToReal(-2),kmadd(csetemp17,kadd(csetemp10,kadd(csetemp11,csetemp9)),kmsub(csetemp15,kmul(mTXYZ,ToReal(2)),kdiv(kmul(csetemp18,kmul(ktanh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM))))))))))))),ksqrt(kmul(csetemp13,kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED K31 = 
      kdiv(kmul(csetemp13,kmul(csetemp14,kmul(X,kmul(Z,kmul(ToReal(-2),kmadd(csetemp17,kadd(csetemp10,kadd(csetemp11,csetemp9)),kmsub(csetemp15,kmul(mTXYZ,ToReal(2)),kdiv(kmul(csetemp18,kmul(ktanh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM))))))))))))),ksqrt(kmul(csetemp13,kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED K22 = 
      kdiv(kmul(csetemp13,kmul(csetemp14,kmul(ToReal(-2),kmadd(mTXYZ,kmsub(csetemp10,kmul(csetemp15,ToReal(2)),csetemp16),kmul(csetemp10,kmsub(csetemp17,kadd(csetemp10,kadd(csetemp11,csetemp9)),kdiv(kmul(csetemp18,kmul(ktanh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))))))))))),ksqrt(kmul(csetemp13,kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED K32 = 
      kdiv(kmul(csetemp13,kmul(csetemp14,kmul(Y,kmul(Z,kmul(ToReal(-2),kmadd(csetemp17,kadd(csetemp10,kadd(csetemp11,csetemp9)),kmsub(csetemp15,kmul(mTXYZ,ToReal(2)),kdiv(kmul(csetemp18,kmul(ktanh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM))))))))))))),ksqrt(kmul(csetemp13,kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED K33 = 
      kdiv(kmul(csetemp13,kmul(csetemp14,kmul(ToReal(-2),kmadd(mTXYZ,kmsub(csetemp11,kmul(csetemp15,ToReal(2)),csetemp16),kmul(csetemp11,kmsub(csetemp17,kadd(csetemp10,kadd(csetemp11,csetemp9)),kdiv(kmul(csetemp18,kmul(ktanh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))))))))))),ksqrt(kmul(csetemp13,kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED alpp = 
      kdiv(ToReal(1),ksqrt(kmul(csetemp13,kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED dtalpp = 
      kdiv(kmul(csetemp13,kmul(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(kpow(kmul(csetemp13,kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15)),-1.5),kmul(ktanh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),ToReal(-2*dM))))),kmul(kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM))))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED betap1 = 
      kdiv(kmul(mTXYZ,kmul(rXYZ,kmul(X,ToReal(2)))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED betap2 = 
      kdiv(kmul(mTXYZ,kmul(rXYZ,kmul(Y,ToReal(2)))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED betap3 = 
      kdiv(kmul(mTXYZ,kmul(rXYZ,kmul(Z,ToReal(2)))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED dtbetap1 = 
      kdiv(kmul(csetemp18,kmul(X,kmul(ktanh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),ToReal(4*dM)))),kmul(kmul(kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15)),kmul(kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED dtbetap2 = 
      kdiv(kmul(csetemp18,kmul(Y,kmul(ktanh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),ToReal(4*dM)))),kmul(kmul(kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15)),kmul(kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED dtbetap3 = 
      kdiv(kmul(csetemp18,kmul(Z,kmul(ktanh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),ToReal(4*dM)))),kmul(kmul(kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2)),csetemp15)),kmul(kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp19,ToReal(dM)))))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp20 = kmul(Jac11,Jac11);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp21 = kmul(Jac21,Jac21);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp22 = kmul(Jac31,Jac31);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED gxxL = 
      kmadd(csetemp20,G11,kmadd(csetemp21,G22,kmadd(csetemp22,G33,kmul(kmadd(G32,kmul(Jac21,Jac31),kmul(Jac11,kmadd(G21,Jac21,kmul(G31,Jac31)))),ToReal(2)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED gxyL = 
      kmadd(Jac12,kmadd(G11,Jac11,kmadd(G21,Jac21,kmul(G31,Jac31))),kmadd(Jac22,kmadd(G21,Jac11,kmadd(G22,Jac21,kmul(G32,Jac31))),kmul(kmadd(G31,Jac11,kmadd(G32,Jac21,kmul(G33,Jac31))),Jac32)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED gxzL = 
      kmadd(Jac13,kmadd(G11,Jac11,kmadd(G21,Jac21,kmul(G31,Jac31))),kmadd(Jac23,kmadd(G21,Jac11,kmadd(G22,Jac21,kmul(G32,Jac31))),kmul(kmadd(G31,Jac11,kmadd(G32,Jac21,kmul(G33,Jac31))),Jac33)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp23 = kmul(Jac12,Jac12);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp24 = kmul(Jac22,Jac22);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp25 = kmul(Jac32,Jac32);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED gyyL = 
      kmadd(csetemp23,G11,kmadd(csetemp24,G22,kmadd(csetemp25,G33,kmul(kmadd(G32,kmul(Jac22,Jac32),kmul(Jac12,kmadd(G21,Jac22,kmul(G31,Jac32)))),ToReal(2)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED gyzL = 
      kmadd(Jac13,kmadd(G11,Jac12,kmadd(G21,Jac22,kmul(G31,Jac32))),kmadd(Jac23,kmadd(G21,Jac12,kmadd(G22,Jac22,kmul(G32,Jac32))),kmul(kmadd(G31,Jac12,kmadd(G32,Jac22,kmul(G33,Jac32))),Jac33)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp26 = kmul(Jac13,Jac13);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp27 = kmul(Jac23,Jac23);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED csetemp28 = kmul(Jac33,Jac33);
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED gzzL = 
      kmadd(csetemp26,G11,kmadd(csetemp27,G22,kmadd(csetemp28,G33,kmul(kmadd(G32,kmul(Jac23,Jac33),kmul(Jac13,kmadd(G21,Jac23,kmul(G31,Jac33)))),ToReal(2)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED kxxL = 
      kmadd(csetemp20,K11,kmadd(csetemp21,K22,kmadd(csetemp22,K33,kmul(kmadd(Jac11,kmadd(Jac21,K21,kmul(Jac31,K31)),kmul(Jac21,kmul(Jac31,K32))),ToReal(2)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED kxyL = 
      kmadd(Jac11,kmadd(Jac12,K11,kmadd(Jac22,K21,kmul(Jac32,K31))),kmadd(Jac21,kmadd(Jac12,K21,kmadd(Jac22,K22,kmul(Jac32,K32))),kmul(Jac31,kmadd(Jac12,K31,kmadd(Jac22,K32,kmul(Jac32,K33))))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED kxzL = 
      kmadd(Jac11,kmadd(Jac13,K11,kmadd(Jac23,K21,kmul(Jac33,K31))),kmadd(Jac21,kmadd(Jac13,K21,kmadd(Jac23,K22,kmul(Jac33,K32))),kmul(Jac31,kmadd(Jac13,K31,kmadd(Jac23,K32,kmul(Jac33,K33))))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED kyyL = 
      kmadd(csetemp23,K11,kmadd(csetemp24,K22,kmadd(csetemp25,K33,kmul(kmadd(Jac12,kmadd(Jac22,K21,kmul(Jac32,K31)),kmul(Jac22,kmul(Jac32,K32))),ToReal(2)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED kyzL = 
      kmadd(Jac12,kmadd(Jac13,K11,kmadd(Jac23,K21,kmul(Jac33,K31))),kmadd(Jac22,kmadd(Jac13,K21,kmadd(Jac23,K22,kmul(Jac33,K32))),kmul(Jac32,kmadd(Jac13,K31,kmadd(Jac23,K32,kmul(Jac33,K33))))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED kzzL = 
      kmadd(csetemp26,K11,kmadd(csetemp27,K22,kmadd(csetemp28,K33,kmul(kmadd(Jac13,kmadd(Jac23,K21,kmul(Jac33,K31)),kmul(Jac23,kmul(Jac33,K32))),ToReal(2)))));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED alpL = 
      kmul(alpp,ToReal(lapsefactor));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED dtalpL = 
      kmul(dtalpp,ToReal(lapsefactor));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED betaxL = 
      kmadd(betap1,InvJac11,kmadd(betap2,InvJac12,kmadd(betap3,InvJac13,shiftadd1)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED betayL = 
      kmadd(betap1,InvJac21,kmadd(betap2,InvJac22,kmadd(betap3,InvJac23,shiftadd2)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED betazL = 
      kmadd(betap1,InvJac31,kmadd(betap2,InvJac32,kmadd(betap3,InvJac33,shiftadd3)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED dtbetaxL = 
      kmadd(dtbetap1,InvJac11,kmadd(dtbetap2,InvJac12,kmul(dtbetap3,InvJac13)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED dtbetayL = 
      kmadd(dtbetap1,InvJac21,kmadd(dtbetap2,InvJac22,kmul(dtbetap3,InvJac23)));
    
    CCTK_REAL_VEC CCTK_ATTRIBUTE_UNUSED dtbetazL = 
      kmadd(dtbetap1,InvJac31,kmadd(dtbetap2,InvJac32,kmul(dtbetap3,InvJac33)));
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,lc_imin,lc_imax);
    vec_store_nta_partial(alp[index],alpL);
    vec_store_nta_partial(betax[index],betaxL);
    vec_store_nta_partial(betay[index],betayL);
    vec_store_nta_partial(betaz[index],betazL);
    vec_store_nta_partial(dtalp[index],dtalpL);
    vec_store_nta_partial(dtbetax[index],dtbetaxL);
    vec_store_nta_partial(dtbetay[index],dtbetayL);
    vec_store_nta_partial(dtbetaz[index],dtbetazL);
    vec_store_nta_partial(gxx[index],gxxL);
    vec_store_nta_partial(gxy[index],gxyL);
    vec_store_nta_partial(gxz[index],gxzL);
    vec_store_nta_partial(gyy[index],gyyL);
    vec_store_nta_partial(gyz[index],gyzL);
    vec_store_nta_partial(gzz[index],gzzL);
    vec_store_nta_partial(kxx[index],kxxL);
    vec_store_nta_partial(kxy[index],kxyL);
    vec_store_nta_partial(kxz[index],kxzL);
    vec_store_nta_partial(kyy[index],kyyL);
    vec_store_nta_partial(kyz[index],kyzL);
    vec_store_nta_partial(kzz[index],kzzL);
  }
  LC_ENDLOOP3VEC(Vaidya2_initial);
}

extern "C" void Vaidya2_initial(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering Vaidya2_initial_Body");
  }
  
  if (cctk_iteration % Vaidya2_initial_calc_every != Vaidya2_initial_calc_offset)
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
  GenericFD_AssertGroupStorage(cctkGH, "Vaidya2_initial", 7, groups);
  
  
  GenericFD_LoopOverEverything(cctkGH, Vaidya2_initial_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving Vaidya2_initial_Body");
  }
}
