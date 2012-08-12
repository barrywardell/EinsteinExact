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
  ptrdiff_t const di = 1;
  ptrdiff_t const dj = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const dk = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  ptrdiff_t const cdi = sizeof(CCTK_REAL) * di;
  ptrdiff_t const cdj = sizeof(CCTK_REAL) * dj;
  ptrdiff_t const cdk = sizeof(CCTK_REAL) * dk;
  CCTK_REAL_VEC const dx = ToReal(CCTK_DELTA_SPACE(0));
  CCTK_REAL_VEC const dy = ToReal(CCTK_DELTA_SPACE(1));
  CCTK_REAL_VEC const dz = ToReal(CCTK_DELTA_SPACE(2));
  CCTK_REAL_VEC const dt = ToReal(CCTK_DELTA_TIME);
  CCTK_REAL_VEC const t = ToReal(cctk_time);
  CCTK_REAL_VEC const dxi = INV(dx);
  CCTK_REAL_VEC const dyi = INV(dy);
  CCTK_REAL_VEC const dzi = INV(dz);
  CCTK_REAL_VEC const khalf = ToReal(0.5);
  CCTK_REAL_VEC const kthird = ToReal(1.0/3.0);
  CCTK_REAL_VEC const ktwothird = ToReal(2.0/3.0);
  CCTK_REAL_VEC const kfourthird = ToReal(4.0/3.0);
  CCTK_REAL_VEC const keightthird = ToReal(8.0/3.0);
  CCTK_REAL_VEC const hdxi = kmul(ToReal(0.5), dxi);
  CCTK_REAL_VEC const hdyi = kmul(ToReal(0.5), dyi);
  CCTK_REAL_VEC const hdzi = kmul(ToReal(0.5), dzi);
  
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
    ptrdiff_t const index = di*i + dj*j + dk*k;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alpL = vec_load(alp[index]);
    CCTK_REAL_VEC xL = vec_load(x[index]);
    CCTK_REAL_VEC yL = vec_load(y[index]);
    CCTK_REAL_VEC zL = vec_load(z[index]);
    
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC xx1 = xL;
    
    CCTK_REAL_VEC xx2 = yL;
    
    CCTK_REAL_VEC xx3 = zL;
    
    CCTK_REAL_VEC position1 = ToReal(positionx);
    
    CCTK_REAL_VEC position2 = ToReal(positiony);
    
    CCTK_REAL_VEC position3 = ToReal(positionz);
    
    CCTK_REAL_VEC shiftadd1 = ToReal(shiftaddx);
    
    CCTK_REAL_VEC shiftadd2 = ToReal(shiftaddy);
    
    CCTK_REAL_VEC shiftadd3 = ToReal(shiftaddz);
    
    CCTK_REAL_VEC csetemp0 = kcos(ToReal(phi));
    
    CCTK_REAL_VEC csetemp1 = kcos(ToReal(psi));
    
    CCTK_REAL_VEC csetemp2 = kcos(ToReal(theta));
    
    CCTK_REAL_VEC csetemp3 = ksin(ToReal(phi));
    
    CCTK_REAL_VEC csetemp4 = ksin(ToReal(psi));
    
    CCTK_REAL_VEC Jac11 = 
      kmsub(csetemp0,csetemp1,kmul(csetemp2,kmul(csetemp3,csetemp4)));
    
    CCTK_REAL_VEC Jac12 = 
      kmadd(csetemp1,csetemp3,kmul(csetemp0,kmul(csetemp2,csetemp4)));
    
    CCTK_REAL_VEC csetemp5 = ksin(ToReal(theta));
    
    CCTK_REAL_VEC Jac13 = kmul(csetemp4,csetemp5);
    
    CCTK_REAL_VEC Jac21 = 
      knmadd(csetemp1,kmul(csetemp2,csetemp3),kmul(csetemp0,csetemp4));
    
    CCTK_REAL_VEC Jac22 = 
      kmsub(csetemp0,kmul(csetemp1,csetemp2),kmul(csetemp3,csetemp4));
    
    CCTK_REAL_VEC Jac23 = kmul(csetemp1,csetemp5);
    
    CCTK_REAL_VEC Jac31 = kmul(csetemp3,csetemp5);
    
    CCTK_REAL_VEC Jac32 = kneg(kmul(csetemp0,csetemp5));
    
    CCTK_REAL_VEC Jac33 = csetemp2;
    
    CCTK_REAL_VEC InvJac11 = Jac11;
    
    CCTK_REAL_VEC InvJac12 = Jac21;
    
    CCTK_REAL_VEC InvJac13 = Jac31;
    
    CCTK_REAL_VEC InvJac21 = Jac12;
    
    CCTK_REAL_VEC InvJac22 = Jac22;
    
    CCTK_REAL_VEC InvJac23 = Jac32;
    
    CCTK_REAL_VEC InvJac31 = Jac13;
    
    CCTK_REAL_VEC InvJac32 = Jac23;
    
    CCTK_REAL_VEC InvJac33 = Jac33;
    
    CCTK_REAL_VEC T = ksub(t,ToReal(positiont));
    
    CCTK_REAL_VEC csetemp6 = kneg(kmul(shiftadd1,T));
    
    CCTK_REAL_VEC csetemp7 = kneg(kmul(shiftadd2,T));
    
    CCTK_REAL_VEC csetemp8 = kneg(kmul(shiftadd3,T));
    
    CCTK_REAL_VEC XX1 = 
      kmadd(Jac11,kadd(csetemp6,ksub(xx1,position1)),kmadd(Jac12,kadd(csetemp7,ksub(xx2,position2)),kmul(Jac13,kadd(csetemp8,ksub(xx3,position3)))));
    
    CCTK_REAL_VEC XX2 = 
      kmadd(Jac21,kadd(csetemp6,ksub(xx1,position1)),kmadd(Jac22,kadd(csetemp7,ksub(xx2,position2)),kmul(Jac23,kadd(csetemp8,ksub(xx3,position3)))));
    
    CCTK_REAL_VEC XX3 = 
      kmadd(Jac31,kadd(csetemp6,ksub(xx1,position1)),kmadd(Jac32,kadd(csetemp7,ksub(xx2,position2)),kmul(Jac33,kadd(csetemp8,ksub(xx3,position3)))));
    
    CCTK_REAL_VEC X = XX1;
    
    CCTK_REAL_VEC Y = XX2;
    
    CCTK_REAL_VEC Z = XX3;
    
    CCTK_REAL_VEC csetemp9 = kmul(X,X);
    
    CCTK_REAL_VEC csetemp10 = kmul(Y,Y);
    
    CCTK_REAL_VEC csetemp11 = kmul(Z,Z);
    
    CCTK_REAL_VEC rXYZ = ksqrt(kadd(csetemp10,kadd(csetemp11,csetemp9)));
    
    CCTK_REAL_VEC csetemp12 = ToReal(ScalarINV(M));
    
    CCTK_REAL_VEC mTXYZ = 
      kmul(kmadd(ktanh(kmul(csetemp12,kmul(kadd(T,ksqrt(kadd(csetemp10,kadd(csetemp11,csetemp9)))),ToReal(dM)))),ktanh(kmul(csetemp12,kmul(kadd(T,ksqrt(kadd(csetemp10,kadd(csetemp11,csetemp9)))),ToReal(dM)))),ToReal(1.)),ToReal(M));
    
    CCTK_REAL_VEC csetemp13 = kdiv(ToReal(1.),kmul(rXYZ,kmul(rXYZ,rXYZ)));
    
    CCTK_REAL_VEC csetemp14 = kmul(rXYZ,kmul(rXYZ,rXYZ));
    
    alpL = 
      kdiv(ToReal(1.),ksqrt(kmul(csetemp13,kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14))));
    
    CCTK_REAL_VEC csetemp15 = kadd(rXYZ,T);
    
    CCTK_REAL_VEC dtalpL = 
      kdiv(kmul(csetemp13,kmul(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(kpow(kmul(csetemp13,kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14)),-1.5),kmul(ktanh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),ToReal(-2.*dM))))),kmul(kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM))))));
    
    CCTK_REAL_VEC G11 = 
      kmadd(csetemp13,kmul(csetemp9,kmul(mTXYZ,ToReal(2.))),ToReal(1.));
    
    CCTK_REAL_VEC G21 = 
      kmul(csetemp13,kmul(mTXYZ,kmul(X,kmul(Y,ToReal(2.)))));
    
    CCTK_REAL_VEC G31 = 
      kmul(csetemp13,kmul(mTXYZ,kmul(X,kmul(Z,ToReal(2.)))));
    
    CCTK_REAL_VEC G22 = 
      kmadd(csetemp10,kmul(csetemp13,kmul(mTXYZ,ToReal(2.))),ToReal(1.));
    
    CCTK_REAL_VEC G32 = 
      kmul(csetemp13,kmul(mTXYZ,kmul(Y,kmul(Z,ToReal(2.)))));
    
    CCTK_REAL_VEC G33 = 
      kmadd(csetemp11,kmul(csetemp13,kmul(mTXYZ,ToReal(2.))),ToReal(1.));
    
    CCTK_REAL_VEC csetemp16 = kdiv(ToReal(1.),alpL);
    
    CCTK_REAL_VEC csetemp17 = 
      kdiv(ToReal(1.),kmul(kmul(rXYZ,rXYZ),kmul(rXYZ,rXYZ)));
    
    CCTK_REAL_VEC csetemp18 = 
      kmul(rXYZ,kmul(kmul(rXYZ,rXYZ),kmul(rXYZ,rXYZ)));
    
    CCTK_REAL_VEC csetemp19 = kmul(mTXYZ,mTXYZ);
    
    CCTK_REAL_VEC csetemp20 = kmul(kmul(rXYZ,rXYZ),kmul(rXYZ,rXYZ));
    
    CCTK_REAL_VEC K11 = 
      kdiv(kmul(csetemp16,kmul(csetemp17,kmul(ToReal(-2.),kmsub(csetemp9,kmsub(csetemp19,kadd(csetemp10,kadd(csetemp11,csetemp9)),kdiv(kmul(csetemp20,kmul(ktanh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM))))))),kmul(mTXYZ,kmadd(csetemp14,kmul(csetemp9,ToReal(-2.)),csetemp18)))))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14));
    
    CCTK_REAL_VEC K21 = 
      kdiv(kmul(csetemp16,kmul(csetemp17,kmul(X,kmul(Y,kmul(ToReal(-2.),kmadd(csetemp19,kadd(csetemp10,kadd(csetemp11,csetemp9)),kmsub(csetemp14,kmul(mTXYZ,ToReal(2.)),kdiv(kmul(csetemp20,kmul(ktanh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM))))))))))))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14));
    
    CCTK_REAL_VEC K31 = 
      kdiv(kmul(csetemp16,kmul(csetemp17,kmul(X,kmul(Z,kmul(ToReal(-2.),kmadd(csetemp19,kadd(csetemp10,kadd(csetemp11,csetemp9)),kmsub(csetemp14,kmul(mTXYZ,ToReal(2.)),kdiv(kmul(csetemp20,kmul(ktanh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM))))))))))))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14));
    
    CCTK_REAL_VEC K22 = 
      kdiv(kmul(csetemp16,kmul(csetemp17,kmul(ToReal(-2.),kmadd(mTXYZ,kmsub(csetemp10,kmul(csetemp14,ToReal(2.)),csetemp18),kmul(csetemp10,kmsub(csetemp19,kadd(csetemp10,kadd(csetemp11,csetemp9)),kdiv(kmul(csetemp20,kmul(ktanh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))))))))))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14));
    
    CCTK_REAL_VEC K32 = 
      kdiv(kmul(csetemp16,kmul(csetemp17,kmul(Y,kmul(Z,kmul(ToReal(-2.),kmadd(csetemp19,kadd(csetemp10,kadd(csetemp11,csetemp9)),kmsub(csetemp14,kmul(mTXYZ,ToReal(2.)),kdiv(kmul(csetemp20,kmul(ktanh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM))))))))))))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14));
    
    CCTK_REAL_VEC K33 = 
      kdiv(kmul(csetemp16,kmul(csetemp17,kmul(ToReal(-2.),kmadd(mTXYZ,kmsub(csetemp11,kmul(csetemp14,ToReal(2.)),csetemp18),kmul(csetemp11,kmsub(csetemp19,kadd(csetemp10,kadd(csetemp11,csetemp9)),kdiv(kmul(csetemp20,kmul(ktanh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),ToReal(dM))),kmul(kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))))))))))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14));
    
    CCTK_REAL_VEC betap1 = 
      kdiv(kmul(mTXYZ,kmul(rXYZ,kmul(X,ToReal(2.)))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14));
    
    CCTK_REAL_VEC betap2 = 
      kdiv(kmul(mTXYZ,kmul(rXYZ,kmul(Y,ToReal(2.)))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14));
    
    CCTK_REAL_VEC betap3 = 
      kdiv(kmul(mTXYZ,kmul(rXYZ,kmul(Z,ToReal(2.)))),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14));
    
    CCTK_REAL_VEC dtbetap1 = 
      kdiv(kmul(csetemp20,kmul(X,kmul(ktanh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),ToReal(4.*dM)))),kmul(kmul(kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14)),kmul(kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))))));
    
    CCTK_REAL_VEC dtbetap2 = 
      kdiv(kmul(csetemp20,kmul(Y,kmul(ktanh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),ToReal(4.*dM)))),kmul(kmul(kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14)),kmul(kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))))));
    
    CCTK_REAL_VEC dtbetap3 = 
      kdiv(kmul(csetemp20,kmul(Z,kmul(ktanh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),ToReal(4.*dM)))),kmul(kmul(kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp9)),kmul(mTXYZ,ToReal(2.)),csetemp14)),kmul(kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))),kcosh(kmul(csetemp12,kmul(csetemp15,ToReal(dM)))))));
    
    CCTK_REAL_VEC csetemp21 = kmul(Jac11,Jac11);
    
    CCTK_REAL_VEC csetemp22 = kmul(Jac21,Jac21);
    
    CCTK_REAL_VEC csetemp23 = kmul(Jac31,Jac31);
    
    CCTK_REAL_VEC gxxL = 
      kmadd(csetemp21,G11,kmadd(csetemp22,G22,kmadd(csetemp23,G33,kmul(kmadd(G32,kmul(Jac21,Jac31),kmul(Jac11,kmadd(G21,Jac21,kmul(G31,Jac31)))),ToReal(2.)))));
    
    CCTK_REAL_VEC gxyL = 
      kmadd(Jac12,kmadd(G11,Jac11,kmadd(G21,Jac21,kmul(G31,Jac31))),kmadd(Jac22,kmadd(G21,Jac11,kmadd(G22,Jac21,kmul(G32,Jac31))),kmul(kmadd(G31,Jac11,kmadd(G32,Jac21,kmul(G33,Jac31))),Jac32)));
    
    CCTK_REAL_VEC gxzL = 
      kmadd(Jac13,kmadd(G11,Jac11,kmadd(G21,Jac21,kmul(G31,Jac31))),kmadd(Jac23,kmadd(G21,Jac11,kmadd(G22,Jac21,kmul(G32,Jac31))),kmul(kmadd(G31,Jac11,kmadd(G32,Jac21,kmul(G33,Jac31))),Jac33)));
    
    CCTK_REAL_VEC csetemp24 = kmul(Jac12,Jac12);
    
    CCTK_REAL_VEC csetemp25 = kmul(Jac22,Jac22);
    
    CCTK_REAL_VEC csetemp26 = kmul(Jac32,Jac32);
    
    CCTK_REAL_VEC gyyL = 
      kmadd(csetemp24,G11,kmadd(csetemp25,G22,kmadd(csetemp26,G33,kmul(kmadd(G32,kmul(Jac22,Jac32),kmul(Jac12,kmadd(G21,Jac22,kmul(G31,Jac32)))),ToReal(2.)))));
    
    CCTK_REAL_VEC gyzL = 
      kmadd(Jac13,kmadd(G11,Jac12,kmadd(G21,Jac22,kmul(G31,Jac32))),kmadd(Jac23,kmadd(G21,Jac12,kmadd(G22,Jac22,kmul(G32,Jac32))),kmul(kmadd(G31,Jac12,kmadd(G32,Jac22,kmul(G33,Jac32))),Jac33)));
    
    CCTK_REAL_VEC csetemp27 = kmul(Jac13,Jac13);
    
    CCTK_REAL_VEC csetemp28 = kmul(Jac23,Jac23);
    
    CCTK_REAL_VEC csetemp29 = kmul(Jac33,Jac33);
    
    CCTK_REAL_VEC gzzL = 
      kmadd(csetemp27,G11,kmadd(csetemp28,G22,kmadd(csetemp29,G33,kmul(kmadd(G32,kmul(Jac23,Jac33),kmul(Jac13,kmadd(G21,Jac23,kmul(G31,Jac33)))),ToReal(2.)))));
    
    CCTK_REAL_VEC kxxL = 
      kmadd(csetemp21,K11,kmadd(csetemp22,K22,kmadd(csetemp23,K33,kmul(kmadd(Jac11,kmadd(Jac21,K21,kmul(Jac31,K31)),kmul(Jac21,kmul(Jac31,K32))),ToReal(2.)))));
    
    CCTK_REAL_VEC kxyL = 
      kmadd(Jac11,kmadd(Jac12,K11,kmadd(Jac22,K21,kmul(Jac32,K31))),kmadd(Jac21,kmadd(Jac12,K21,kmadd(Jac22,K22,kmul(Jac32,K32))),kmul(Jac31,kmadd(Jac12,K31,kmadd(Jac22,K32,kmul(Jac32,K33))))));
    
    CCTK_REAL_VEC kxzL = 
      kmadd(Jac11,kmadd(Jac13,K11,kmadd(Jac23,K21,kmul(Jac33,K31))),kmadd(Jac21,kmadd(Jac13,K21,kmadd(Jac23,K22,kmul(Jac33,K32))),kmul(Jac31,kmadd(Jac13,K31,kmadd(Jac23,K32,kmul(Jac33,K33))))));
    
    CCTK_REAL_VEC kyyL = 
      kmadd(csetemp24,K11,kmadd(csetemp25,K22,kmadd(csetemp26,K33,kmul(kmadd(Jac12,kmadd(Jac22,K21,kmul(Jac32,K31)),kmul(Jac22,kmul(Jac32,K32))),ToReal(2.)))));
    
    CCTK_REAL_VEC kyzL = 
      kmadd(Jac12,kmadd(Jac13,K11,kmadd(Jac23,K21,kmul(Jac33,K31))),kmadd(Jac22,kmadd(Jac13,K21,kmadd(Jac23,K22,kmul(Jac33,K32))),kmul(Jac32,kmadd(Jac13,K31,kmadd(Jac23,K32,kmul(Jac33,K33))))));
    
    CCTK_REAL_VEC kzzL = 
      kmadd(csetemp27,K11,kmadd(csetemp28,K22,kmadd(csetemp29,K33,kmul(kmadd(Jac13,kmadd(Jac23,K21,kmul(Jac33,K31)),kmul(Jac23,kmul(Jac33,K32))),ToReal(2.)))));
    
    CCTK_REAL_VEC betaxL = 
      kmadd(betap1,InvJac11,kmadd(betap2,InvJac12,kmadd(betap3,InvJac13,shiftadd1)));
    
    CCTK_REAL_VEC betayL = 
      kmadd(betap1,InvJac21,kmadd(betap2,InvJac22,kmadd(betap3,InvJac23,shiftadd2)));
    
    CCTK_REAL_VEC betazL = 
      kmadd(betap1,InvJac31,kmadd(betap2,InvJac32,kmadd(betap3,InvJac33,shiftadd3)));
    
    CCTK_REAL_VEC dtbetaxL = 
      kmadd(dtbetap1,InvJac11,kmadd(dtbetap2,InvJac12,kmul(dtbetap3,InvJac13)));
    
    CCTK_REAL_VEC dtbetayL = 
      kmadd(dtbetap1,InvJac21,kmadd(dtbetap2,InvJac22,kmul(dtbetap3,InvJac23)));
    
    CCTK_REAL_VEC dtbetazL = 
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
