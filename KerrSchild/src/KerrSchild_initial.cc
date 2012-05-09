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
#define QAD(x) (SQR(SQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))

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
  LC_LOOP3VEC(KerrSchild_initial,
    i,j,k, imin[0],imin[1],imin[2], imax[0],imax[1],imax[2],
    cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],
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
    
    CCTK_REAL_VEC csetemp9 = SQR(ToReal(a));
    
    CCTK_REAL_VEC csetemp10 = SQR(X);
    
    CCTK_REAL_VEC csetemp11 = SQR(Y);
    
    CCTK_REAL_VEC csetemp12 = SQR(Z);
    
    CCTK_REAL_VEC rXYZ = 
      kmul(INV(ksqrt(ToReal(2))),ksqrt(kadd(csetemp10,kadd(csetemp11,kadd(csetemp12,ksub(ksqrt(kmadd(csetemp12,kmul(csetemp9,ToReal(4)),SQR(kadd(csetemp10,kadd(csetemp11,ksub(csetemp12,csetemp9)))))),csetemp9))))));
    
    CCTK_REAL_VEC csetemp13 = CUB(rXYZ);
    
    CCTK_REAL_VEC csetemp14 = QAD(rXYZ);
    
    alpL = 
      INV(ksqrt(kmadd(csetemp13,kmul(INV(kmadd(csetemp12,csetemp9,csetemp14)),kmul(ToReal(2),ToReal(M))),ToReal(1))));
    
    CCTK_REAL_VEC dtalpL = ToReal(0);
    
    CCTK_REAL_VEC csetemp15 = SQR(rXYZ);
    
    CCTK_REAL_VEC csetemp16 = kmul(rXYZ,X);
    
    CCTK_REAL_VEC csetemp17 = kmul(Y,ToReal(a));
    
    CCTK_REAL_VEC G11 = 
      kmadd(csetemp13,kmul(INV(kmul(kmadd(csetemp12,csetemp9,csetemp14),SQR(kadd(csetemp15,csetemp9)))),kmul(SQR(kadd(csetemp16,csetemp17)),kmul(ToReal(2),ToReal(M)))),ToReal(1));
    
    CCTK_REAL_VEC csetemp18 = kmul(X,ToReal(a));
    
    CCTK_REAL_VEC csetemp19 = kneg(csetemp18);
    
    CCTK_REAL_VEC csetemp20 = kmul(rXYZ,Y);
    
    CCTK_REAL_VEC G21 = 
      kmul(csetemp13,kmul(kadd(csetemp16,csetemp17),kmul(kadd(csetemp19,csetemp20),kmul(INV(kmul(kmadd(csetemp12,csetemp9,csetemp14),SQR(kadd(csetemp15,csetemp9)))),kmul(ToReal(2),ToReal(M))))));
    
    CCTK_REAL_VEC G31 = 
      kmul(csetemp15,kmul(kadd(csetemp16,csetemp17),kmul(Z,kmul(INV(kmul(kadd(csetemp15,csetemp9),kmadd(csetemp12,csetemp9,csetemp14))),kmul(ToReal(2),ToReal(M))))));
    
    CCTK_REAL_VEC G22 = 
      kmadd(csetemp13,kmul(INV(kmul(kmadd(csetemp12,csetemp9,csetemp14),SQR(kadd(csetemp15,csetemp9)))),kmul(SQR(kadd(csetemp19,csetemp20)),kmul(ToReal(2),ToReal(M)))),ToReal(1));
    
    CCTK_REAL_VEC G32 = 
      kmul(csetemp15,kmul(kadd(csetemp19,csetemp20),kmul(Z,kmul(INV(kmul(kadd(csetemp15,csetemp9),kmadd(csetemp12,csetemp9,csetemp14))),kmul(ToReal(2),ToReal(M))))));
    
    CCTK_REAL_VEC G33 = 
      kmadd(csetemp12,kmul(rXYZ,kmul(INV(kmadd(csetemp12,csetemp9,csetemp14)),kmul(ToReal(2),ToReal(M)))),ToReal(1));
    
    CCTK_REAL_VEC csetemp21 = kmul(Y,kmul(ToReal(4),ToReal(M)));
    
    CCTK_REAL_VEC csetemp22 = kmul(csetemp21,X);
    
    CCTK_REAL_VEC csetemp23 = kmul(X,kmul(Y,kmul(ToReal(2),ToReal(M))));
    
    CCTK_REAL_VEC csetemp24 = INV(alpL);
    
    CCTK_REAL_VEC csetemp25 = kpow(rXYZ,8);
    
    CCTK_REAL_VEC csetemp26 = kpow(rXYZ,5);
    
    CCTK_REAL_VEC csetemp27 = kpow(rXYZ,7);
    
    CCTK_REAL_VEC csetemp28 = QAD(ToReal(a));
    
    CCTK_REAL_VEC csetemp29 = kpow(rXYZ,6);
    
    CCTK_REAL_VEC csetemp30 = CUB(ToReal(a));
    
    CCTK_REAL_VEC csetemp31 = kpow(ToReal(a),5);
    
    CCTK_REAL_VEC csetemp32 = kpow(rXYZ,9);
    
    CCTK_REAL_VEC csetemp33 = kpow(ToReal(a),6);
    
    CCTK_REAL_VEC K11 = 
      kmul(csetemp13,kmul(csetemp24,kmul(INV(kmul(QAD(kadd(csetemp15,csetemp9)),kmul(SQR(kmadd(csetemp12,csetemp9,csetemp14)),kmul(SQR(kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M)))))),ksqrt(kmadd(csetemp12,kmul(csetemp9,ToReal(4)),SQR(kadd(csetemp10,kadd(csetemp11,ksub(csetemp12,csetemp9)))))))))),kmul(ToReal(-2),kmul(ToReal(M),knmsub(kmadd(csetemp12,kmul(kadd(csetemp16,csetemp17),kmul(SQR(kadd(csetemp15,csetemp9)),kmadd(csetemp14,kmul(csetemp30,kmul(Y,ToReal(-5))),kmadd(csetemp27,kmul(X,ToReal(-3)),kmadd(csetemp26,kmul(csetemp9,kmul(X,ToReal(-3))),kmadd(csetemp15,kmul(csetemp30,kmul(Y,kadd(csetemp12,kmadd(csetemp9,ToReal(-2),kmadd(csetemp10,ToReal(2),kmul(csetemp11,ToReal(2))))))),kmadd(csetemp12,kmul(csetemp31,kmul(Y,ToReal(3))),kmadd(csetemp13,kmul(csetemp9,kmul(X,kmadd(csetemp9,ToReal(-2),kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp12,ToReal(3))))))),kmadd(csetemp12,kmul(csetemp28,kmul(rXYZ,kmul(X,ToReal(5)))),kmul(csetemp29,kmul(Y,kmul(ToReal(-5),ToReal(a))))))))))))),kmul(csetemp13,kmul(SQR(kadd(csetemp16,csetemp17)),kmadd(csetemp9,kmsub(csetemp26,kmsub(kadd(csetemp12,ksub(csetemp9,csetemp11)),ToReal(2),csetemp10),kmul(csetemp12,kmul(csetemp13,kmadd(kadd(csetemp11,kmadd(csetemp9,ToReal(-3),csetemp12)),ToReal(2),csetemp10)))),kmadd(X,kmul(Y,kmsub(csetemp12,kmul(csetemp31,ToReal(3)),kmul(csetemp14,csetemp30))),kmadd(csetemp12,kmsub(csetemp28,kmul(rXYZ,kmadd(kadd(csetemp11,csetemp12),ToReal(-2),kmadd(csetemp9,ToReal(2),kmul(csetemp10,ToReal(3))))),kmul(csetemp15,kmul(csetemp30,kmul(X,Y)))),kmadd(csetemp32,ToReal(4),kmadd(csetemp27,kmadd(csetemp10,ToReal(-5),kmadd(kadd(csetemp11,csetemp12),ToReal(-2),kmul(csetemp9,ToReal(6)))),kmul(csetemp29,kmul(X,kmul(Y,kmul(ToReal(-5),ToReal(a))))))))))))),kmul(ToReal(M),kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M)))))),kmadd(kadd(csetemp16,csetemp17),kmul(ToReal(M),kmsub(csetemp12,kmul(kmadd(csetemp12,csetemp9,csetemp14),kmul(rXYZ,kmul(X,kmul(CUB(kadd(csetemp15,csetemp9)),kmul(ToReal(4),kadd(csetemp14,kmsub(csetemp13,ToReal(M),kmul(csetemp12,csetemp9)))))))),kmul(csetemp13,kmul(kadd(csetemp19,csetemp20),kmul(kmadd(csetemp12,kmul(csetemp13,kmul(csetemp9,kmul(X,Y))),kmadd(csetemp26,kmul(csetemp9,kmul(X,Y)),kmadd(csetemp27,kmul(X,kmul(Y,ToReal(-3))),kmadd(csetemp14,kmul(csetemp30,kmadd(csetemp11,ToReal(-3),kmadd(csetemp10,ToReal(-2),kmul(kadd(csetemp12,csetemp9),ToReal(2))))),kmadd(csetemp12,kmsub(csetemp31,kadd(csetemp11,kmadd(csetemp10,ToReal(-2),kmadd(csetemp12,ToReal(-2),kmul(csetemp9,ToReal(2))))),kmul(csetemp15,kmul(csetemp30,kmadd(csetemp9,ToReal(-6),kmadd(csetemp10,ToReal(2),kmadd(csetemp12,ToReal(2),kmul(csetemp11,ToReal(3)))))))),kmadd(csetemp12,kmul(csetemp28,kmul(rXYZ,kmul(X,kmul(Y,ToReal(5))))),kmul(kmsub(csetemp25,ToReal(4),kmul(csetemp29,kmadd(csetemp9,ToReal(-6),kmadd(csetemp10,ToReal(2),kmadd(csetemp12,ToReal(2),kmul(csetemp11,ToReal(7))))))),ToReal(a)))))))),kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M)))))))))),kmadd(kmadd(csetemp12,csetemp9,csetemp14),kmul(kadd(csetemp25,kmadd(csetemp12,csetemp33,kmadd(csetemp12,kmul(csetemp15,kmul(csetemp28,ToReal(2))),kmadd(csetemp29,kmul(csetemp9,ToReal(2)),kmadd(csetemp14,kmul(ToReal(a),kadd(csetemp22,kmadd(csetemp12,ToReal(a),csetemp30))),kmadd(csetemp10,kmul(csetemp26,kmul(ToReal(2),ToReal(M))),kmul(csetemp11,kmul(csetemp13,kmul(csetemp9,kmul(ToReal(2),ToReal(M))))))))))),kmadd(csetemp12,kmul(csetemp31,kmul(X,kmul(Y,ToReal(-3)))),kmadd(csetemp32,ToReal(-2),kmadd(csetemp12,kmul(csetemp28,kmul(rXYZ,kadd(csetemp11,kadd(csetemp12,kmsub(csetemp10,ToReal(-3),csetemp9))))),kmadd(csetemp12,kmsub(csetemp13,kmul(csetemp9,kadd(csetemp11,kadd(csetemp12,kmsub(csetemp9,ToReal(-3),csetemp10)))),kmul(csetemp15,kmul(csetemp30,kmul(X,Y)))),kmadd(csetemp27,kadd(csetemp11,kadd(csetemp12,kmadd(csetemp9,ToReal(-3),kmul(csetemp10,ToReal(3))))),kmadd(csetemp26,kmul(ToReal(a),kadd(csetemp22,kmsub(kadd(csetemp10,ksub(csetemp11,csetemp12)),ToReal(a),csetemp30))),kmadd(csetemp25,kmul(ToReal(-4),ToReal(M)),kmadd(csetemp14,kmul(csetemp9,kmadd(X,kmul(Y,ToReal(a)),kmadd(csetemp9,kmul(ToReal(-2),ToReal(M)),kmul(kadd(csetemp11,csetemp12),kmul(ToReal(2),ToReal(M)))))),kmul(csetemp29,kmadd(X,kmul(Y,kmul(ToReal(3),ToReal(a))),kmadd(csetemp9,kmul(ToReal(-6),ToReal(M)),kmul(ToReal(2),kmul(kadd(csetemp11,kmadd(csetemp10,ToReal(2),csetemp12)),ToReal(M))))))))))))))),kmul(csetemp13,kmul(kadd(csetemp16,csetemp17),kmul(kadd(csetemp19,csetemp20),kmul(kmadd(csetemp12,csetemp9,csetemp14),kmul(ToReal(2),kmul(ToReal(M),kmadd(csetemp12,kmul(csetemp28,kmul(rXYZ,kmul(X,kmul(Y,ToReal(-4))))),knmsub(csetemp12,kmadd(csetemp15,kmul(csetemp30,kadd(csetemp11,kmadd(csetemp9,ToReal(-3),csetemp12))),kmul(csetemp31,kadd(csetemp11,kadd(csetemp12,kmsub(csetemp10,ToReal(-2),csetemp9))))),kmadd(csetemp25,kmul(ToReal(2),ToReal(a)),kmadd(csetemp29,kadd(csetemp23,kmsub(csetemp30,ToReal(3),kmul(kadd(csetemp11,kmadd(csetemp10,ToReal(4),csetemp12)),ToReal(a)))),kmadd(csetemp26,kmul(ToReal(-2),kmul(kadd(csetemp11,kadd(csetemp12,kmadd(csetemp9,ToReal(-3),kmul(csetemp10,ToReal(3))))),kmul(ToReal(a),ToReal(M)))),kmadd(csetemp14,kmul(csetemp9,kadd(csetemp30,kmadd(kadd(csetemp12,kmsub(csetemp10,ToReal(-2),csetemp11)),ToReal(a),kmul(X,kmul(Y,kmul(ToReal(-2),ToReal(M))))))),kmadd(csetemp27,kmul(ToReal(2),kmadd(X,Y,kmul(ToReal(2),kmul(ToReal(a),ToReal(M))))),kmul(csetemp13,kmul(csetemp9,kmul(ToReal(2),kmsub(knmsub(kadd(csetemp10,kadd(csetemp11,csetemp12)),ToReal(a),csetemp30),ToReal(M),kmul(csetemp12,kmul(X,Y)))))))))))))))))))))))))));
    
    CCTK_REAL_VEC csetemp34 = kpow(rXYZ,19);
    
    CCTK_REAL_VEC csetemp35 = kpow(rXYZ,18);
    
    CCTK_REAL_VEC csetemp36 = kpow(rXYZ,17);
    
    CCTK_REAL_VEC csetemp37 = SQR(ToReal(M));
    
    CCTK_REAL_VEC csetemp38 = kpow(ToReal(a),11);
    
    CCTK_REAL_VEC csetemp39 = kpow(Z,6);
    
    CCTK_REAL_VEC csetemp40 = kpow(ToReal(a),10);
    
    CCTK_REAL_VEC csetemp41 = QAD(Z);
    
    CCTK_REAL_VEC csetemp42 = kpow(rXYZ,16);
    
    CCTK_REAL_VEC csetemp43 = kpow(rXYZ,15);
    
    CCTK_REAL_VEC csetemp44 = kpow(rXYZ,14);
    
    CCTK_REAL_VEC csetemp45 = kpow(rXYZ,10);
    
    CCTK_REAL_VEC csetemp46 = kpow(rXYZ,13);
    
    CCTK_REAL_VEC csetemp47 = kpow(rXYZ,12);
    
    CCTK_REAL_VEC csetemp48 = kpow(ToReal(a),9);
    
    CCTK_REAL_VEC csetemp49 = kpow(ToReal(a),8);
    
    CCTK_REAL_VEC csetemp50 = kpow(ToReal(a),7);
    
    CCTK_REAL_VEC csetemp51 = kpow(rXYZ,11);
    
    CCTK_REAL_VEC K21 = 
      kneg(kmul(csetemp13,kmul(csetemp24,kmul(INV(kmul(QAD(kadd(csetemp15,csetemp9)),kmul(SQR(kmadd(csetemp12,csetemp9,csetemp14)),kmul(SQR(kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M)))))),ksqrt(kmadd(csetemp12,kmul(csetemp9,ToReal(4)),SQR(kadd(csetemp10,kadd(csetemp11,ksub(csetemp12,csetemp9)))))))))),kmul(ToReal(M),kmadd(csetemp35,kmul(ToReal(-3),kmadd(ksub(csetemp10,csetemp11),ToReal(a),csetemp22)),kmadd(csetemp39,kmul(kpow(ToReal(a),12),kmul(ToReal(3),kmadd(ksub(csetemp10,csetemp11),ToReal(a),csetemp23))),kmadd(csetemp42,kmadd(ksub(csetemp10,csetemp11),kmadd(csetemp30,ToReal(-7),kmul(csetemp37,kmul(ToReal(32),ToReal(a)))),kmul(X,kmul(Y,kmul(kmadd(csetemp9,ToReal(-4),kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),ToReal(22))),ToReal(M))))),kmadd(csetemp44,kmul(ToReal(a),kmadd(ksub(csetemp10,csetemp11),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,ToReal(-36)),kmadd(csetemp28,ToReal(-5),kmul(csetemp9,kmadd(csetemp12,ToReal(-5),kmul(csetemp37,ToReal(48)))))),kmul(X,kmul(Y,kmul(kmadd(csetemp30,ToReal(12),kmul(csetemp12,kmul(ToReal(8),ToReal(a)))),ToReal(M)))))),kmadd(csetemp41,kmadd(csetemp15,kmul(csetemp40,kmadd(csetemp12,kmul(ksub(csetemp10,csetemp11),kmul(ToReal(7),ToReal(a))),kmul(X,kmul(Y,kmul(kadd(csetemp10,ksub(csetemp11,kadd(csetemp9,csetemp12))),kmul(ToReal(4),ToReal(M))))))),kmul(csetemp14,kmul(csetemp49,kmadd(ksub(csetemp10,csetemp11),kmadd(csetemp30,ToReal(5),kmul(kmadd(csetemp12,ToReal(5),kmul(csetemp37,ToReal(12))),ToReal(a))),kmul(X,kmul(Y,kmul(kmadd(csetemp9,ToReal(-4),kmul(ToReal(-2),kmadd(kadd(csetemp10,csetemp11),ToReal(5),kmul(csetemp12,ToReal(24))))),ToReal(M)))))))),kmadd(ToReal(4),kmadd(csetemp34,kmul(X,Y),kmadd(csetemp38,kmul(csetemp39,kmul(rXYZ,kmadd(X,kmul(Y,kmul(ToReal(-2),ToReal(a))),kmul(kmadd(csetemp11,ToReal(-3),kmul(csetemp10,ToReal(3))),ToReal(M))))),kmadd(csetemp36,kmadd(X,kmul(Y,kmul(kmadd(csetemp37,ToReal(-4),csetemp9),ToReal(2))),kmul(ksub(csetemp10,csetemp11),kmul(ToReal(3),kmul(ToReal(a),ToReal(M))))),kmul(csetemp43,kmadd(X,kmul(Y,kadd(csetemp28,kmadd(csetemp9,kmadd(csetemp37,ToReal(-4),csetemp12),kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,ToReal(9)))))),kmul(ksub(csetemp10,csetemp11),kmul(kmadd(csetemp30,ToReal(4),kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(ToReal(-5),ToReal(a)))),ToReal(M)))))))),kmadd(ToReal(-2),kmadd(csetemp12,kmul(csetemp26,kmul(csetemp50,kmadd(X,kmul(Y,kmadd(csetemp30,kmul(kmadd(csetemp12,ToReal(2),csetemp37),ToReal(4)),kmul(ToReal(-4),kmul(kmadd(csetemp41,ToReal(-2),kmul(csetemp37,kadd(csetemp10,kmadd(csetemp12,ToReal(2),csetemp11)))),ToReal(a))))),kmul(csetemp12,kmul(ksub(csetemp10,csetemp11),kmul(kmadd(csetemp12,ToReal(-10),kmadd(kadd(csetemp10,csetemp11),ToReal(-9),kmul(csetemp9,ToReal(2)))),ToReal(M))))))),kmul(csetemp13,kmul(csetemp41,kmul(csetemp48,kmsub(ToReal(2),kmadd(X,kmul(Y,kmul(kmadd(csetemp37,ToReal(-3),kmul(csetemp12,ToReal(5))),ToReal(a))),kmul(csetemp9,kmul(ksub(csetemp10,csetemp11),ToReal(M)))),kmul(ksub(csetemp10,csetemp11),kmul(kmadd(kadd(csetemp10,csetemp11),ToReal(5),kmul(csetemp12,ToReal(16))),ToReal(M)))))))),kmadd(csetemp9,kmadd(csetemp47,kmadd(csetemp31,ksub(csetemp11,csetemp10),kmadd(ksub(csetemp10,csetemp11),kmadd(csetemp30,kmadd(csetemp12,ToReal(-9),kmul(csetemp37,ToReal(16))),kmul(csetemp37,kmul(ToReal(-4),kmul(kmadd(kadd(csetemp10,csetemp11),ToReal(5),kmul(csetemp12,ToReal(11))),ToReal(a))))),kmul(X,kmul(Y,kmul(kmadd(kadd(csetemp10,csetemp11),kmul(csetemp9,ToReal(-6)),kmadd(csetemp28,ToReal(4),kmul(csetemp12,kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),ToReal(16))))),ToReal(M)))))),kmul(ToReal(2),kmadd(csetemp46,kmul(ToReal(M),kmadd(ksub(csetemp10,csetemp11),kmsub(csetemp30,ToReal(2),kmul(kmadd(kadd(csetemp10,csetemp11),ToReal(7),kmul(csetemp12,ToReal(10))),ToReal(a))),kmul(X,kmul(Y,kmul(kmadd(kmadd(kadd(csetemp10,csetemp11),ToReal(-2),csetemp12),ToReal(4),kmul(csetemp9,ToReal(16))),ToReal(M)))))),kmul(csetemp51,kmadd(X,kmul(Y,kmadd(csetemp28,kmul(ToReal(2),kmadd(csetemp12,ToReal(-3),kmul(csetemp37,ToReal(4)))),kmadd(csetemp12,kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,ToReal(6))),kmul(csetemp9,kmul(ToReal(-2),kmadd(csetemp37,kmadd(kadd(csetemp10,csetemp11),ToReal(5),kmul(csetemp12,ToReal(6))),csetemp41)))))),kmul(ksub(csetemp10,csetemp11),kmul(kmsub(csetemp12,kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(ToReal(-10),ToReal(a))),kmul(csetemp30,kadd(csetemp10,kmadd(csetemp12,ToReal(4),csetemp11)))),ToReal(M)))))))),kmul(csetemp12,kmadd(csetemp25,kmul(csetemp28,kmadd(csetemp31,ksub(csetemp10,csetemp11),kmadd(csetemp30,kmul(ksub(csetemp10,csetemp11),kmadd(csetemp37,ToReal(-20),kmul(csetemp12,ToReal(3)))),kmadd(csetemp37,kmul(ksub(csetemp10,csetemp11),kmul(kadd(csetemp12,kmadd(csetemp10,ToReal(3),kmul(csetemp11,ToReal(3)))),kmul(ToReal(4),ToReal(a)))),kmadd(csetemp12,kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(X,kmul(Y,kmul(ToReal(-6),ToReal(M))))),kmadd(csetemp28,kmul(X,kmul(Y,kmul(ToReal(-6),ToReal(M)))),kmul(csetemp9,kmul(X,kmul(Y,kmul(ToReal(-4),kmul(kmadd(csetemp10,ToReal(4),kmadd(csetemp11,ToReal(4),kmul(csetemp12,ToReal(11)))),ToReal(M)))))))))))),knmsub(csetemp29,kmul(csetemp33,kmadd(ksub(csetemp10,csetemp11),kmsub(csetemp30,kmadd(csetemp12,ToReal(-9),kmul(csetemp37,ToReal(8))),kmul(kadd(csetemp41,kmadd(csetemp10,kmul(csetemp37,ToReal(8)),kmadd(csetemp11,kmul(csetemp37,ToReal(8)),kmul(csetemp12,kmul(csetemp37,ToReal(28)))))),ToReal(a))),kmadd(csetemp9,kmul(X,kmul(Y,kmul(ToReal(-4),kmul(kadd(csetemp10,ksub(csetemp11,csetemp12)),ToReal(M))))),kmadd(csetemp28,kmul(X,kmul(Y,kmul(ToReal(4),ToReal(M)))),kmul(csetemp12,kmul(X,kmul(Y,kmul(ToReal(4),kmul(kmadd(csetemp10,ToReal(9),kmadd(csetemp11,ToReal(9),kmul(csetemp12,ToReal(11)))),ToReal(M)))))))))),kmsub(ToReal(-4),kmadd(csetemp27,kmul(csetemp33,kmadd(X,kmul(Y,knmsub(csetemp37,kadd(csetemp10,kmadd(csetemp12,ToReal(-6),csetemp11)),csetemp41)),kmadd(csetemp9,kmul(X,kmul(Y,kmul(ToReal(3),kmadd(csetemp12,ToReal(3),csetemp37)))),kmadd(csetemp30,kmul(ksub(csetemp10,csetemp11),ToReal(M)),kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(ToReal(-2),kmul(ksub(csetemp10,csetemp11),kmul(ToReal(a),ToReal(M))))))))),kmul(csetemp28,kmul(csetemp32,kmadd(csetemp28,kmul(X,kmul(Y,ToReal(2))),kmadd(csetemp37,kmul(X,kmul(Y,kmul(ToReal(2),kmadd(csetemp12,ToReal(2),kmadd(csetemp10,ToReal(3),kmul(csetemp11,ToReal(3))))))),kmadd(csetemp12,kmul(csetemp9,kmul(X,kmul(Y,ToReal(6)))),kmul(ksub(csetemp10,csetemp11),kmul(kmsub(csetemp30,ToReal(3),kmul(kadd(csetemp10,kmadd(csetemp12,ToReal(-3),csetemp11)),ToReal(a))),ToReal(M))))))))),kmul(csetemp30,kmul(csetemp45,kmadd(csetemp28,kmul(ksub(csetemp10,csetemp11),ToReal(3)),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,kmul(ksub(csetemp10,csetemp11),ToReal(12))),kmadd(csetemp9,kmul(ksub(csetemp10,csetemp11),kmadd(csetemp37,ToReal(20),csetemp12)),kmadd(csetemp30,kmul(X,kmul(Y,kmul(ToReal(-12),ToReal(M)))),kmul(X,kmul(Y,kmul(kmadd(csetemp12,ToReal(2),kmadd(csetemp10,ToReal(3),kmul(csetemp11,ToReal(3)))),kmul(ToReal(12),kmul(ToReal(a),ToReal(M)))))))))))))))))))))))))))));
    
    CCTK_REAL_VEC csetemp52 = kmul(X,kmul(ToReal(4),ToReal(M)));
    
    CCTK_REAL_VEC csetemp53 = kpow(ToReal(a),13);
    
    CCTK_REAL_VEC csetemp54 = QAD(X);
    
    CCTK_REAL_VEC csetemp55 = QAD(Y);
    
    CCTK_REAL_VEC K31 = 
      kmul(csetemp24,kmul(rXYZ,kmul(Z,kmul(INV(kmul(CUB(kadd(csetemp15,csetemp9)),kmul(SQR(kmadd(csetemp12,csetemp9,csetemp14)),kmul(SQR(kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M)))))),ksqrt(kmadd(csetemp12,kmul(csetemp9,ToReal(4)),SQR(kadd(csetemp10,kadd(csetemp11,ksub(csetemp12,csetemp9)))))))))),kmul(ToReal(M),kmadd(csetemp34,kmul(X,ToReal(-4)),kmadd(csetemp35,kmul(ToReal(3),knmsub(Y,ToReal(a),csetemp52)),kmadd(csetemp39,kmadd(csetemp53,kmul(Y,ToReal(3)),kmul(csetemp38,kmul(rXYZ,kmul(ToReal(2),kmadd(X,kmul(ToReal(3),ToReal(a)),csetemp21))))),kmadd(csetemp12,kmul(csetemp27,kmul(csetemp31,kmul(ToReal(2),kmadd(csetemp31,kmul(X,ToReal(-2)),kmadd(csetemp30,kmul(X,kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp12,ToReal(9))))),kmadd(csetemp12,kmul(X,kmul(kadd(csetemp10,kadd(csetemp11,kmadd(csetemp37,ToReal(2),kmul(csetemp12,ToReal(3))))),ToReal(a))),kmadd(Y,kmul(SQR(kadd(csetemp10,kadd(csetemp11,csetemp12))),kmul(ToReal(2),ToReal(M))),kmul(csetemp9,kmul(Y,kmul(ToReal(-4),kmul(kadd(csetemp10,kmadd(csetemp12,ToReal(3),csetemp11)),ToReal(M)))))))))))),kmadd(csetemp12,kmul(csetemp25,kmul(csetemp28,kmadd(csetemp31,kmul(Y,ToReal(-15)),kmadd(csetemp30,kmul(Y,kmadd(csetemp37,ToReal(-4),kmadd(csetemp12,ToReal(7),kmadd(csetemp10,ToReal(8),kmul(csetemp11,ToReal(8)))))),kmadd(csetemp37,kmul(Y,kmul(ToReal(-4),kmul(kadd(csetemp10,kmadd(csetemp12,ToReal(3),csetemp11)),ToReal(a)))),kmadd(csetemp9,kmul(X,kmul(ToReal(-4),kmul(kadd(csetemp10,kmadd(csetemp12,ToReal(3),csetemp11)),ToReal(M)))),kmul(X,kmul(ToReal(2),kmul(kmadd(csetemp54,ToReal(2),kmadd(csetemp55,ToReal(2),kmadd(csetemp10,kmul(csetemp11,ToReal(4)),kmadd(csetemp41,ToReal(5),kmadd(csetemp10,kmul(csetemp12,ToReal(7)),kmul(csetemp11,kmul(csetemp12,ToReal(7)))))))),ToReal(M)))))))))),kmadd(csetemp12,kmul(csetemp26,kmul(csetemp50,kmul(ToReal(2),kmadd(csetemp12,kmul(csetemp30,X),kmadd(csetemp12,kmul(X,kmul(kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmadd(csetemp37,ToReal(2),kmul(csetemp12,ToReal(9))))),ToReal(a))),kmadd(csetemp9,kmul(Y,kmul(ToReal(-2),kmul(kadd(csetemp10,kmadd(csetemp12,ToReal(4),csetemp11)),ToReal(M)))),kmul(Y,kmul(kmadd(ToReal(2),kadd(csetemp54,kmadd(csetemp41,ToReal(5),csetemp55)),kmadd(csetemp11,kmul(csetemp12,ToReal(11)),kmul(csetemp10,kmadd(csetemp11,ToReal(4),kmul(csetemp12,ToReal(11)))))),ToReal(M))))))))),kmadd(csetemp13,kmul(csetemp41,kmul(csetemp48,kmul(ToReal(2),kmadd(X,kmsub(kadd(csetemp10,kmadd(csetemp12,ToReal(9),csetemp11)),ToReal(a),csetemp30),kmul(Y,kmul(kmadd(csetemp9,ToReal(-2),kmadd(csetemp10,ToReal(5),kmadd(csetemp11,ToReal(5),kmul(csetemp12,ToReal(12))))),ToReal(M))))))),kmadd(csetemp12,kmul(csetemp29,kmul(csetemp33,kmadd(csetemp31,kmul(Y,ToReal(-4)),kmadd(csetemp30,kmul(Y,kmadd(csetemp12,ToReal(3),kmadd(csetemp10,ToReal(4),kmul(csetemp11,ToReal(4))))),kmadd(csetemp12,kmul(Y,kmul(kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmadd(csetemp12,ToReal(3),kmul(csetemp37,ToReal(4))))),ToReal(a))),kmadd(csetemp9,kmul(X,kmul(ToReal(-4),kmul(kadd(csetemp10,kmadd(csetemp12,ToReal(3),csetemp11)),ToReal(M)))),kmul(X,kmul(ToReal(2),kmul(kmadd(csetemp54,ToReal(2),kmadd(csetemp55,ToReal(2),kmadd(csetemp10,kmul(csetemp11,ToReal(4)),kmadd(csetemp10,kmul(csetemp12,ToReal(17)),kmadd(csetemp11,kmul(csetemp12,ToReal(17)),kmul(csetemp41,ToReal(19))))))),ToReal(M)))))))))),knmsub(csetemp42,kmadd(csetemp30,kmul(Y,ToReal(11)),kmadd(csetemp37,kmul(Y,kmul(ToReal(-32),ToReal(a))),kmadd(csetemp9,kmul(X,kmul(ToReal(-24),ToReal(M))),kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(X,kmul(ToReal(22),ToReal(M))))))),knmsub(csetemp30,kmul(csetemp45,kmadd(csetemp9,kmul(Y,kmsub(csetemp37,kmul(ToReal(4),kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp12,ToReal(7))))),kmul(csetemp12,kmadd(csetemp12,ToReal(3),kmadd(csetemp10,ToReal(4),kmul(csetemp11,ToReal(4))))))),kmadd(Y,kmsub(csetemp33,ToReal(2),kmul(csetemp28,kmadd(csetemp12,ToReal(-21),kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp37,ToReal(8))))))),kmadd(csetemp12,kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,kmul(Y,ToReal(12)))),kmadd(csetemp30,kmul(X,kmul(ToReal(2),kmul(kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp12,ToReal(5)))),ToReal(M)))),kmul(X,kmul(ToReal(-4),kmul(kadd(csetemp54,kadd(csetemp55,kmadd(csetemp41,ToReal(-4),kmul(csetemp10,kmul(csetemp11,ToReal(2)))))),kmul(ToReal(a),ToReal(M)))))))))),kmadd(csetemp51,kmul(csetemp9,kmul(ToReal(-2),kmadd(X,knmsub(csetemp28,kadd(csetemp10,kadd(csetemp11,kmadd(csetemp12,ToReal(-9),kmul(csetemp37,ToReal(4))))),csetemp33),kmadd(csetemp12,kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,kmul(X,ToReal(6)))),kmadd(csetemp9,kmul(X,kmul(ToReal(2),kmsub(csetemp37,kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp12,ToReal(7)))),kmul(csetemp12,kadd(csetemp10,kmadd(csetemp12,ToReal(2),csetemp11)))))),kmadd(csetemp31,kmul(Y,kmul(ToReal(-2),ToReal(M))),kmadd(csetemp30,kmul(Y,kmul(ToReal(3),kmul(kmadd(csetemp10,ToReal(3),kmadd(csetemp11,ToReal(3),kmul(csetemp12,ToReal(4)))),ToReal(M)))),kmul(Y,kmul(ToReal(-2),kmul(kadd(csetemp54,kadd(csetemp55,kmadd(csetemp41,ToReal(-4),kmadd(csetemp10,kmul(csetemp12,ToReal(-3)),kmadd(csetemp11,kmul(csetemp12,ToReal(-3)),kmul(csetemp10,kmul(csetemp11,ToReal(2)))))))),kmul(ToReal(a),ToReal(M)))))))))))),kmadd(csetemp28,kmul(csetemp32,kmul(ToReal(-2),kmadd(csetemp12,kmul(csetemp37,kmul(X,kmul(ToReal(2),kadd(csetemp10,kmadd(csetemp12,ToReal(3),csetemp11))))),kmadd(csetemp12,kmul(X,kmsub(csetemp28,ToReal(7),kmul(csetemp9,kmadd(csetemp37,ToReal(-2),kmadd(csetemp10,ToReal(4),kmadd(csetemp11,ToReal(4),kmul(csetemp12,ToReal(11)))))))),kmadd(csetemp30,kmul(Y,kmul(ToReal(2),kmul(kadd(csetemp10,kmadd(csetemp12,ToReal(2),csetemp11)),ToReal(M)))),kmul(Y,kmul(ToReal(-2),kmul(kadd(csetemp54,kadd(csetemp55,kmadd(csetemp41,ToReal(-7),kmadd(csetemp11,kmul(csetemp12,ToReal(-2)),kmul(csetemp10,kmul(ksub(csetemp11,csetemp12),ToReal(2))))))),kmul(ToReal(a),ToReal(M)))))))))),kmadd(csetemp36,kmul(ToReal(2),kmadd(csetemp9,kmul(X,ToReal(-7)),kmadd(csetemp37,kmul(X,ToReal(16)),kmul(Y,kmul(ToReal(6),kmul(ToReal(a),ToReal(M))))))),kmadd(csetemp43,kmul(ToReal(-2),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,kmul(X,ToReal(18))),kmadd(X,kmsub(csetemp28,ToReal(9),kmul(csetemp9,kadd(csetemp10,kadd(csetemp11,kmsub(csetemp37,ToReal(32),csetemp12))))),kmadd(csetemp30,kmul(Y,kmul(ToReal(-12),ToReal(M))),kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(Y,kmul(ToReal(10),kmul(ToReal(a),ToReal(M))))))))),knmsub(csetemp44,kmul(ToReal(a),kmadd(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,kmul(Y,ToReal(36))),kmadd(Y,kmsub(csetemp28,ToReal(15),kmul(csetemp9,kmadd(csetemp12,ToReal(-3),kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp37,ToReal(64))))))),kmadd(csetemp30,kmul(X,kmul(ToReal(-16),ToReal(M))),kmul(X,kmul(ToReal(2),kmul(kmadd(csetemp10,ToReal(17),kmadd(csetemp11,ToReal(17),kmul(csetemp12,ToReal(21)))),kmul(ToReal(a),ToReal(M))))))))),kmadd(csetemp41,kmadd(csetemp15,kmul(csetemp40,kmadd(csetemp30,kmul(Y,ToReal(-2)),kmadd(Y,kmul(kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp12,ToReal(9)))),ToReal(a)),kmul(csetemp12,kmul(X,kmul(ToReal(14),ToReal(M))))))),kmul(csetemp14,kmul(csetemp49,kmadd(csetemp30,kmul(Y,ToReal(-3)),kmadd(Y,kmul(kmadd(csetemp10,ToReal(4),kmadd(csetemp11,ToReal(4),kmadd(csetemp37,ToReal(4),kmul(csetemp12,ToReal(9))))),ToReal(a)),kmadd(csetemp9,kmul(X,kmul(ToReal(-4),ToReal(M))),kmul(X,kmul(ToReal(2),kmul(kmadd(csetemp10,ToReal(8),kmadd(csetemp11,ToReal(8),kmul(csetemp12,ToReal(21)))),ToReal(M)))))))))),kmul(csetemp9,kmsub(csetemp46,kmul(ToReal(-2),kmadd(csetemp37,kmul(X,kmul(kmadd(csetemp10,ToReal(3),kmadd(csetemp11,ToReal(3),kmul(csetemp12,ToReal(5)))),ToReal(6))),kmadd(X,kmsub(csetemp28,ToReal(5),kmul(csetemp9,kmadd(csetemp12,ToReal(-5),kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp37,ToReal(20))))))),kmadd(csetemp30,kmul(Y,kmul(ToReal(-8),ToReal(M))),kmul(Y,kmul(ToReal(3),kmul(kmadd(csetemp10,ToReal(5),kmadd(csetemp11,ToReal(5),kmul(csetemp12,ToReal(6)))),kmul(ToReal(a),ToReal(M))))))))),kmul(csetemp47,kmadd(Y,kmsub(csetemp31,ToReal(9),kmul(csetemp30,kmadd(csetemp12,ToReal(-13),kmadd(csetemp10,ToReal(4),kmadd(csetemp11,ToReal(4),kmul(csetemp37,ToReal(40))))))),kmadd(csetemp37,kmul(Y,kmul(kmadd(csetemp10,ToReal(3),kmadd(csetemp11,ToReal(3),kmul(csetemp12,ToReal(5)))),kmul(ToReal(12),ToReal(a)))),kmadd(csetemp28,kmul(X,kmul(ToReal(-4),ToReal(M))),kmadd(X,kmul(ToReal(-4),kmul(kadd(csetemp54,kadd(csetemp55,kmadd(csetemp41,ToReal(-3),kmadd(csetemp11,kmul(csetemp12,ToReal(-2)),kmul(csetemp10,kmul(ksub(csetemp11,csetemp12),ToReal(2))))))),ToReal(M))),kmul(csetemp9,kmul(X,kmul(kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp12,ToReal(3)))),kmul(ToReal(10),ToReal(M)))))))))))))))))))))))))))))))));
    
    CCTK_REAL_VEC csetemp56 = kneg(csetemp20);
    
    CCTK_REAL_VEC K22 = 
      kmul(csetemp13,kmul(csetemp24,kmul(INV(kmul(QAD(kadd(csetemp15,csetemp9)),kmul(SQR(kmadd(csetemp12,csetemp9,csetemp14)),kmul(SQR(kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M)))))),ksqrt(kmadd(csetemp12,kmul(csetemp9,ToReal(4)),SQR(kadd(csetemp10,kadd(csetemp11,ksub(csetemp12,csetemp9)))))))))),kmul(ToReal(-2),kmul(ToReal(M),knmsub(kmadd(csetemp12,kmul(kadd(csetemp18,csetemp56),kmul(SQR(kadd(csetemp15,csetemp9)),kmadd(csetemp14,kmul(csetemp30,kmul(X,ToReal(-5))),kmadd(csetemp12,kmul(csetemp28,kmul(rXYZ,kmul(Y,ToReal(-5)))),kmadd(csetemp15,kmul(csetemp30,kmul(X,kadd(csetemp12,kmadd(csetemp9,ToReal(-2),kmadd(csetemp10,ToReal(2),kmul(csetemp11,ToReal(2))))))),kmadd(csetemp13,kmul(csetemp9,kmul(Y,kmadd(csetemp12,ToReal(-3),kmadd(csetemp10,ToReal(-2),kmadd(csetemp11,ToReal(-2),kmul(csetemp9,ToReal(2))))))),kmadd(csetemp12,kmul(csetemp31,kmul(X,ToReal(3))),kmadd(csetemp27,kmul(Y,ToReal(3)),kmadd(csetemp26,kmul(csetemp9,kmul(Y,ToReal(3))),kmul(csetemp29,kmul(X,kmul(ToReal(-5),ToReal(a))))))))))))),kmul(csetemp13,kmul(SQR(kadd(csetemp18,csetemp56)),kmadd(csetemp14,kmul(csetemp30,kmul(X,Y)),kmadd(csetemp12,kmul(csetemp15,kmul(csetemp30,kmul(X,Y))),kmadd(csetemp12,kmul(csetemp31,kmul(X,kmul(Y,ToReal(-3)))),kmadd(csetemp9,kmsub(csetemp26,kmadd(csetemp10,ToReal(-2),kmsub(kadd(csetemp12,csetemp9),ToReal(2),csetemp11)),kmul(csetemp12,kmul(csetemp13,kadd(csetemp11,kmadd(csetemp9,ToReal(-6),kmadd(csetemp10,ToReal(2),kmul(csetemp12,ToReal(2)))))))),kmadd(csetemp12,kmul(csetemp28,kmul(rXYZ,kmadd(csetemp10,ToReal(-2),kmadd(csetemp12,ToReal(-2),kmadd(csetemp9,ToReal(2),kmul(csetemp11,ToReal(3))))))),kmadd(csetemp32,ToReal(4),kmsub(csetemp29,kmul(X,kmul(Y,kmul(ToReal(5),ToReal(a)))),kmul(csetemp27,kmadd(csetemp9,ToReal(-6),kmadd(csetemp10,ToReal(2),kmadd(csetemp12,ToReal(2),kmul(csetemp11,ToReal(5)))))))))))))))),kmul(ToReal(M),kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M)))))),kmadd(kadd(csetemp19,csetemp20),kmadd(csetemp12,kmul(kmadd(csetemp12,csetemp9,csetemp14),kmul(rXYZ,kmul(Y,kmul(CUB(kadd(csetemp15,csetemp9)),kmul(ToReal(4),kmul(ToReal(M),kadd(csetemp14,kmsub(csetemp13,ToReal(M),kmul(csetemp12,csetemp9))))))))),kmul(csetemp13,kmul(kadd(csetemp16,csetemp17),kmul(kmadd(csetemp12,kmul(csetemp28,kmul(rXYZ,kmul(X,kmul(Y,ToReal(-5))))),kmadd(csetemp14,kmul(csetemp30,kmadd(csetemp10,ToReal(-3),kmul(kadd(csetemp12,ksub(csetemp9,csetemp11)),ToReal(2)))),kmadd(csetemp12,kmsub(csetemp31,kadd(csetemp10,kmadd(kadd(csetemp11,csetemp12),ToReal(-2),kmul(csetemp9,ToReal(2)))),kmul(csetemp13,kmul(csetemp9,kmul(X,Y)))),kmadd(X,kmul(Y,kmsub(csetemp27,ToReal(3),kmul(csetemp26,csetemp9))),kmadd(csetemp12,kmul(csetemp15,kmul(csetemp30,kmadd(csetemp10,ToReal(-3),kmadd(kadd(csetemp11,csetemp12),ToReal(-2),kmul(csetemp9,ToReal(6)))))),kmadd(csetemp25,kmul(ToReal(4),ToReal(a)),kmul(csetemp29,kmul(kmadd(csetemp10,ToReal(-7),kmadd(kadd(csetemp11,csetemp12),ToReal(-2),kmul(csetemp9,ToReal(6)))),ToReal(a))))))))),kmul(ToReal(M),kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M)))))))))),kmadd(kmadd(csetemp12,csetemp9,csetemp14),kmul(kadd(csetemp25,kmadd(csetemp12,csetemp33,kmadd(csetemp12,kmul(csetemp15,kmul(csetemp28,ToReal(2))),kmadd(csetemp29,kmul(csetemp9,ToReal(2)),kmadd(csetemp11,kmul(csetemp26,kmul(ToReal(2),ToReal(M))),kmadd(csetemp10,kmul(csetemp13,kmul(csetemp9,kmul(ToReal(2),ToReal(M)))),kmul(csetemp14,kmul(ToReal(a),kadd(csetemp30,kmadd(csetemp12,ToReal(a),kmul(X,kmul(Y,kmul(ToReal(-4),ToReal(M)))))))))))))),kmadd(csetemp12,kmul(csetemp15,kmul(csetemp30,kmul(X,Y))),kmadd(csetemp32,ToReal(-2),kmadd(csetemp12,kmul(csetemp13,kmul(csetemp9,kadd(csetemp10,kadd(csetemp12,kmsub(csetemp9,ToReal(-3),csetemp11))))),kmadd(csetemp12,kmul(csetemp28,kmul(rXYZ,kadd(csetemp10,kadd(csetemp12,kmsub(csetemp11,ToReal(-3),csetemp9))))),kmadd(csetemp12,kmul(csetemp31,kmul(X,kmul(Y,ToReal(3)))),kmadd(csetemp27,kadd(csetemp10,kadd(csetemp12,kmadd(csetemp9,ToReal(-3),kmul(csetemp11,ToReal(3))))),knmsub(csetemp26,kmul(ToReal(a),kadd(csetemp22,knmsub(kadd(csetemp10,ksub(csetemp11,csetemp12)),ToReal(a),csetemp30))),kmadd(csetemp25,kmul(ToReal(-4),ToReal(M)),kmadd(csetemp14,kmul(csetemp9,knmsub(X,kmul(Y,ToReal(a)),kmadd(csetemp9,kmul(ToReal(-2),ToReal(M)),kmul(kadd(csetemp10,csetemp12),kmul(ToReal(2),ToReal(M)))))),kmul(csetemp29,kmadd(X,kmul(Y,kmul(ToReal(-3),ToReal(a))),kmadd(csetemp9,kmul(ToReal(-6),ToReal(M)),kmul(ToReal(2),kmul(kadd(csetemp10,kmadd(csetemp11,ToReal(2),csetemp12)),ToReal(M)))))))))))))))),kmul(csetemp13,kmul(kadd(csetemp16,csetemp17),kmul(kadd(csetemp19,csetemp20),kmul(kmadd(csetemp12,csetemp9,csetemp14),kmul(ToReal(2),kmul(ToReal(M),kmadd(csetemp12,kmul(csetemp28,kmul(rXYZ,kmul(X,kmul(Y,ToReal(-4))))),kmadd(csetemp12,kmul(csetemp15,kmul(csetemp30,kadd(csetemp10,kmadd(csetemp9,ToReal(-3),csetemp12)))),kmadd(csetemp12,kmul(csetemp31,kadd(csetemp10,kadd(csetemp12,kmsub(csetemp11,ToReal(-2),csetemp9)))),kmadd(csetemp25,kmul(ToReal(-2),ToReal(a)),kmadd(csetemp29,kadd(csetemp23,kmadd(csetemp30,ToReal(-3),kmul(kadd(csetemp10,kmadd(csetemp11,ToReal(4),csetemp12)),ToReal(a)))),kmadd(csetemp26,kmul(ToReal(2),kmul(kadd(csetemp10,kadd(csetemp12,kmadd(csetemp9,ToReal(-3),kmul(csetemp11,ToReal(3))))),kmul(ToReal(a),ToReal(M)))),kmadd(csetemp27,kmadd(X,kmul(Y,ToReal(2)),kmul(ToReal(-4),kmul(ToReal(a),ToReal(M)))),kmul(csetemp9,kmsub(csetemp13,kmul(ToReal(2),kmsub(kmsub(kadd(csetemp10,kadd(csetemp11,csetemp12)),ToReal(a),csetemp30),ToReal(M),kmul(csetemp12,kmul(X,Y)))),kmul(csetemp14,kadd(csetemp23,kmadd(kadd(csetemp12,kmsub(csetemp11,ToReal(-2),csetemp10)),ToReal(a),csetemp30))))))))))))))))))))))))));
    
    CCTK_REAL_VEC K32 = 
      kneg(kmul(csetemp24,kmul(rXYZ,kmul(Z,kmul(INV(kmul(CUB(kadd(csetemp15,csetemp9)),kmul(SQR(kmadd(csetemp12,csetemp9,csetemp14)),kmul(SQR(kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M)))))),ksqrt(kmadd(csetemp12,kmul(csetemp9,ToReal(4)),SQR(kadd(csetemp10,kadd(csetemp11,ksub(csetemp12,csetemp9)))))))))),kmul(ToReal(M),kmadd(kadd(csetemp18,csetemp21),kmul(csetemp35,ToReal(-3)),kmadd(csetemp34,kmul(Y,ToReal(4)),kmadd(csetemp39,kmadd(csetemp53,kmul(X,ToReal(3)),kmul(csetemp38,kmul(rXYZ,kmul(ToReal(2),kmadd(Y,kmul(ToReal(-3),ToReal(a)),csetemp52))))),kmadd(csetemp42,kmadd(X,kmadd(csetemp30,ToReal(-11),kmul(csetemp37,kmul(ToReal(32),ToReal(a)))),kmul(Y,kmul(kmadd(csetemp9,ToReal(-24),kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),ToReal(22))),ToReal(M)))),kmadd(csetemp30,kmul(csetemp45,kmadd(X,kmadd(csetemp12,kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,ToReal(-12))),kmadd(csetemp33,ToReal(-2),kmadd(csetemp9,kmadd(csetemp12,kmadd(csetemp12,ToReal(3),kmul(kadd(csetemp10,csetemp11),ToReal(4))),kmul(csetemp37,kmul(ToReal(-4),kmadd(kadd(csetemp10,csetemp11),ToReal(2),kmul(csetemp12,ToReal(7)))))),kmul(csetemp28,kmadd(csetemp12,ToReal(-21),kmadd(kadd(csetemp10,csetemp11),ToReal(2),kmul(csetemp37,ToReal(8)))))))),kmul(Y,kmul(kmadd(csetemp30,kmul(ToReal(2),kmadd(kadd(csetemp10,csetemp11),ToReal(2),kmul(csetemp12,ToReal(5)))),kmul(ToReal(-4),kmul(kadd(csetemp54,kadd(csetemp55,kmadd(csetemp41,ToReal(-4),kmul(csetemp10,kmul(csetemp11,ToReal(2)))))),ToReal(a)))),ToReal(M))))),kmadd(csetemp44,kmul(ToReal(a),kmadd(X,kmadd(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,ToReal(-36)),kmadd(csetemp28,ToReal(-15),kmul(csetemp9,kmadd(csetemp12,ToReal(-3),kmadd(kadd(csetemp10,csetemp11),ToReal(2),kmul(csetemp37,ToReal(64))))))),kmul(Y,kmul(kmadd(csetemp30,ToReal(-16),kmul(ToReal(2),kmul(kmadd(kadd(csetemp10,csetemp11),ToReal(17),kmul(csetemp12,ToReal(21))),ToReal(a)))),ToReal(M))))),kmadd(csetemp12,kmadd(csetemp27,kmul(csetemp31,kmul(ToReal(2),kmadd(Y,kmsub(csetemp31,ToReal(2),kmadd(csetemp12,kmul(kadd(csetemp10,kadd(csetemp11,kmadd(csetemp37,ToReal(2),kmul(csetemp12,ToReal(3))))),ToReal(a)),kmul(csetemp30,kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp12,ToReal(9))))))),kmul(X,kmul(kmadd(SQR(kadd(csetemp10,kadd(csetemp11,csetemp12))),ToReal(2),kmul(csetemp9,kmul(ToReal(-4),kadd(csetemp10,kmadd(csetemp12,ToReal(3),csetemp11))))),ToReal(M)))))),knmsub(csetemp25,kmul(csetemp28,kmadd(X,knmsub(csetemp30,kmadd(csetemp37,ToReal(-4),kmadd(csetemp12,ToReal(7),kmul(kadd(csetemp10,csetemp11),ToReal(8)))),kmadd(csetemp31,ToReal(15),kmul(csetemp37,kmul(kadd(csetemp10,kmadd(csetemp12,ToReal(3),csetemp11)),kmul(ToReal(4),ToReal(a)))))),kmul(Y,kmul(kmadd(csetemp9,kmul(ToReal(-4),kadd(csetemp10,kmadd(csetemp12,ToReal(3),csetemp11))),kmul(ToReal(2),kmadd(kadd(csetemp54,csetemp55),ToReal(2),kmadd(csetemp41,ToReal(5),kmadd(csetemp11,kmul(csetemp12,ToReal(7)),kmul(csetemp10,kmadd(csetemp11,ToReal(4),kmul(csetemp12,ToReal(7))))))))),ToReal(M))))),kmadd(csetemp26,kmul(csetemp50,kmul(ToReal(-2),kmadd(csetemp12,kmul(Y,kmadd(kmadd(kadd(csetemp10,kadd(csetemp11,csetemp37)),ToReal(2),kmul(csetemp12,ToReal(9))),ToReal(a),csetemp30)),kmul(X,kmul(kmadd(csetemp11,kmul(csetemp12,ToReal(-11)),kmadd(csetemp9,kmul(ToReal(2),kadd(csetemp10,kmadd(csetemp12,ToReal(4),csetemp11))),kmsub(ToReal(-2),kadd(csetemp54,kmadd(csetemp41,ToReal(5),csetemp55)),kmul(csetemp10,kmadd(csetemp11,ToReal(4),kmul(csetemp12,ToReal(11))))))),ToReal(M)))))),kmul(csetemp29,kmul(csetemp33,kmadd(X,kmadd(csetemp31,ToReal(-4),kmadd(csetemp30,kmadd(csetemp12,ToReal(3),kmul(kadd(csetemp10,csetemp11),ToReal(4))),kmul(csetemp12,kmul(kmadd(kadd(csetemp10,csetemp11),ToReal(2),kmadd(csetemp12,ToReal(3),kmul(csetemp37,ToReal(4)))),ToReal(a))))),kmul(Y,kmul(kmadd(csetemp9,kmul(kadd(csetemp10,kmadd(csetemp12,ToReal(3),csetemp11)),ToReal(4)),kmul(ToReal(-2),kmadd(kadd(csetemp54,csetemp55),ToReal(2),kmadd(csetemp11,kmul(csetemp12,ToReal(17)),kmadd(csetemp10,kmadd(csetemp11,ToReal(4),kmul(csetemp12,ToReal(17))),kmul(csetemp41,ToReal(19))))))),ToReal(M))))))))),kmadd(csetemp9,kmadd(csetemp47,kmadd(X,kmadd(csetemp31,ToReal(-9),kmadd(csetemp30,kmadd(csetemp12,ToReal(-13),kmadd(kadd(csetemp10,csetemp11),ToReal(4),kmul(csetemp37,ToReal(40)))),kmul(csetemp37,kmul(ToReal(-12),kmul(kmadd(kadd(csetemp10,csetemp11),ToReal(3),kmul(csetemp12,ToReal(5))),ToReal(a)))))),kmul(Y,kmul(kmadd(ToReal(-4),kadd(csetemp28,kadd(csetemp54,kadd(csetemp55,kmadd(csetemp41,ToReal(-3),kmadd(csetemp11,kmul(csetemp12,ToReal(-2)),kmul(csetemp10,kmul(ksub(csetemp11,csetemp12),ToReal(2)))))))),kmul(csetemp9,kmul(kmadd(kadd(csetemp10,csetemp11),ToReal(2),kmul(csetemp12,ToReal(3))),ToReal(10)))),ToReal(M)))),kmul(csetemp46,kmul(ToReal(2),kmadd(Y,kmadd(csetemp28,ToReal(5),kmsub(csetemp37,kmul(kmadd(kadd(csetemp10,csetemp11),ToReal(3),kmul(csetemp12,ToReal(5))),ToReal(6)),kmul(csetemp9,kmadd(csetemp12,ToReal(-5),kmadd(kadd(csetemp10,csetemp11),ToReal(2),kmul(csetemp37,ToReal(20))))))),kmul(X,kmul(kmadd(csetemp30,ToReal(8),kmul(ToReal(-3),kmul(kmadd(kadd(csetemp10,csetemp11),ToReal(5),kmul(csetemp12,ToReal(6))),ToReal(a)))),ToReal(M))))))),kmadd(ToReal(2),kmadd(csetemp36,kmadd(Y,kmadd(csetemp37,ToReal(-16),kmul(csetemp9,ToReal(7))),kmul(X,kmul(ToReal(6),kmul(ToReal(a),ToReal(M))))),kmadd(csetemp43,kmadd(Y,kmadd(csetemp28,ToReal(9),kmsub(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,ToReal(18)),kmul(csetemp9,kadd(csetemp10,kadd(csetemp11,kmsub(csetemp37,ToReal(32),csetemp12)))))),kmul(X,kmul(kmadd(csetemp30,ToReal(12),kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(ToReal(-10),ToReal(a)))),ToReal(M)))),kmadd(csetemp28,kmul(csetemp32,kmadd(csetemp12,kmul(Y,kmadd(csetemp37,kmul(ToReal(2),kadd(csetemp10,kmadd(csetemp12,ToReal(3),csetemp11))),kmsub(csetemp28,ToReal(7),kmul(csetemp9,kmadd(csetemp37,ToReal(-2),kmadd(kadd(csetemp10,csetemp11),ToReal(4),kmul(csetemp12,ToReal(11)))))))),kmul(X,kmul(kmadd(csetemp30,kmul(ToReal(-2),kadd(csetemp10,kmadd(csetemp12,ToReal(2),csetemp11))),kmul(ToReal(2),kmul(kadd(csetemp54,kadd(csetemp55,kmadd(csetemp41,ToReal(-7),kmadd(csetemp11,kmul(csetemp12,ToReal(-2)),kmul(csetemp10,kmul(ksub(csetemp11,csetemp12),ToReal(2))))))),ToReal(a)))),ToReal(M))))),kmul(csetemp51,kmul(csetemp9,kmadd(Y,kadd(csetemp33,knmsub(csetemp28,kadd(csetemp10,kadd(csetemp11,kmadd(csetemp12,ToReal(-9),kmul(csetemp37,ToReal(4))))),kmadd(csetemp12,kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,ToReal(6))),kmul(csetemp9,kmul(ToReal(2),kmsub(csetemp37,kmadd(kadd(csetemp10,csetemp11),ToReal(2),kmul(csetemp12,ToReal(7))),kmul(csetemp12,kadd(csetemp10,kmadd(csetemp12,ToReal(2),csetemp11))))))))),kmul(X,kmul(kmadd(csetemp30,kmul(ToReal(-3),kmadd(kadd(csetemp10,csetemp11),ToReal(3),kmul(csetemp12,ToReal(4)))),kmul(ToReal(2),kmadd(kadd(csetemp54,kadd(csetemp55,kmadd(csetemp41,ToReal(-4),kmadd(csetemp11,kmul(csetemp12,ToReal(-3)),kmul(csetemp10,kmadd(csetemp12,ToReal(-3),kmul(csetemp11,ToReal(2)))))))),ToReal(a),csetemp31))),ToReal(M))))))))),kmul(csetemp41,kmadd(csetemp13,kmul(csetemp48,kmul(ToReal(2),kmadd(Y,knmsub(kadd(csetemp10,kmadd(csetemp12,ToReal(9),csetemp11)),ToReal(a),csetemp30),kmul(X,kmul(kmadd(csetemp9,ToReal(-2),kmadd(kadd(csetemp10,csetemp11),ToReal(5),kmul(csetemp12,ToReal(12)))),ToReal(M)))))),kmadd(csetemp14,kmul(csetemp49,kmadd(X,kmadd(csetemp30,ToReal(-3),kmul(kmadd(kadd(csetemp10,kadd(csetemp11,csetemp37)),ToReal(4),kmul(csetemp12,ToReal(9))),ToReal(a))),kmul(Y,kmul(ToReal(-2),kmul(kmadd(csetemp9,ToReal(-2),kmadd(kadd(csetemp10,csetemp11),ToReal(8),kmul(csetemp12,ToReal(21)))),ToReal(M)))))),kmul(csetemp15,kmul(csetemp40,kmadd(X,kmul(kmadd(kadd(csetemp10,csetemp11),ToReal(2),kmul(csetemp12,ToReal(9))),ToReal(a)),kmul(ToReal(-2),kmadd(csetemp30,X,kmul(csetemp12,kmul(Y,kmul(ToReal(7),ToReal(M)))))))))))))))))))))))))));
    
    CCTK_REAL_VEC csetemp57 = kpow(Z,8);
    
    CCTK_REAL_VEC K33 = 
      kmul(csetemp24,kmul(INV(kmul(kadd(csetemp15,csetemp9),kmul(SQR(kmadd(csetemp12,csetemp9,csetemp14)),kmul(SQR(kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M)))))),ksqrt(kmadd(csetemp12,kmul(csetemp9,ToReal(4)),SQR(kadd(csetemp10,kadd(csetemp11,ksub(csetemp12,csetemp9)))))))))),kmul(ToReal(-2),kmul(ToReal(M),kmadd(kmadd(csetemp40,csetemp57,csetemp35),ToReal(-2),kmadd(csetemp42,kadd(csetemp10,kadd(csetemp11,kmadd(csetemp9,ToReal(-3),kmul(csetemp12,ToReal(3))))),kmadd(csetemp44,kmadd(csetemp12,kmul(csetemp37,ToReal(-16)),kmsub(csetemp9,kadd(csetemp10,kmadd(csetemp12,ToReal(3),csetemp11)),csetemp28)),kmadd(csetemp28,kmul(csetemp29,kmul(csetemp41,kadd(csetemp28,kmsub(csetemp37,kmul(kadd(csetemp10,kmadd(csetemp12,ToReal(2),csetemp11)),ToReal(4)),kmul(csetemp9,kadd(csetemp10,kadd(csetemp11,kmadd(csetemp12,ToReal(3),kmul(csetemp37,ToReal(4)))))))))),knmsub(csetemp39,kmadd(csetemp15,kmul(csetemp49,kadd(csetemp10,kadd(csetemp11,kmsub(csetemp12,ToReal(5),csetemp9)))),kmul(csetemp14,kmul(csetemp33,kadd(csetemp10,kadd(csetemp11,kmadd(csetemp37,ToReal(-2),kmsub(csetemp12,ToReal(3),csetemp9))))))),kmadd(csetemp12,kmul(csetemp45,kmul(csetemp9,ksub(kmadd(csetemp37,kmul(ToReal(4),kmadd(csetemp10,ToReal(2),kmadd(csetemp11,ToReal(2),kmul(csetemp12,ToReal(5))))),kmul(csetemp9,kadd(csetemp10,kadd(csetemp11,kmadd(csetemp37,ToReal(-8),kmul(csetemp12,ToReal(7))))))),csetemp28))),kmadd(csetemp12,kmul(csetemp47,ksub(kmadd(csetemp9,kadd(csetemp10,kmadd(kmadd(csetemp37,ToReal(-8),csetemp12),ToReal(3),csetemp11)),kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,ToReal(18)))),csetemp28)),kmadd(csetemp49,kmul(csetemp57,kmul(rXYZ,kmul(ToReal(-5),ToReal(M)))),kmadd(csetemp36,kmul(ToReal(-4),ToReal(M)),kmadd(csetemp43,kmul(kadd(csetemp10,kadd(csetemp11,kmadd(csetemp9,ToReal(-3),kmul(csetemp12,ToReal(-2))))),kmul(ToReal(2),ToReal(M))),kmadd(csetemp13,kmul(csetemp33,kmul(csetemp39,kmul(kadd(csetemp9,kmadd(csetemp12,ToReal(-3),kmadd(csetemp10,ToReal(-2),kmul(csetemp11,ToReal(-2))))),kmul(ToReal(4),ToReal(M))))),kmadd(csetemp27,kmul(csetemp28,kmul(csetemp41,kmul(ToReal(2),kmul(kmadd(csetemp9,ToReal(-2),kmadd(csetemp10,ToReal(4),kmadd(csetemp11,ToReal(4),kmul(csetemp12,ToReal(7))))),ToReal(M))))),kmadd(csetemp12,kmul(csetemp32,kmul(csetemp9,kmul(kmadd(csetemp28,ToReal(-4),kmadd(ToReal(-4),kadd(csetemp54,kadd(csetemp55,kmsub(csetemp10,kmul(csetemp11,ToReal(2)),csetemp41))),kmul(csetemp9,kmadd(csetemp10,ToReal(8),kmadd(csetemp11,ToReal(8),kmul(csetemp12,ToReal(9))))))),ToReal(M)))),kmadd(csetemp12,kmul(csetemp51,kmul(csetemp9,kmul(ToReal(2),kmul(kmadd(csetemp9,ToReal(-7),kmadd(csetemp10,ToReal(9),kmadd(csetemp11,ToReal(9),kmul(csetemp12,ToReal(11))))),ToReal(M))))),kmadd(csetemp46,kmul(kmadd(csetemp28,ToReal(-2),kmadd(csetemp9,kmul(kadd(csetemp10,kmadd(csetemp12,ToReal(-7),csetemp11)),ToReal(2)),kmul(csetemp12,kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),ToReal(11))))),ToReal(M)),kmul(csetemp41,kmsub(csetemp25,kmul(csetemp9,knmsub(csetemp9,kadd(csetemp10,kadd(csetemp11,kmadd(csetemp37,ToReal(2),kmul(csetemp12,ToReal(3))))),kmadd(csetemp28,ToReal(5),kmul(kadd(csetemp10,kadd(csetemp11,csetemp12)),kmul(csetemp37,ToReal(6)))))),kmul(csetemp26,kmul(csetemp28,kmul(kmadd(csetemp28,ToReal(2),kmadd(csetemp54,ToReal(4),kmadd(csetemp55,ToReal(4),kmadd(csetemp41,ToReal(7),kmadd(csetemp9,kmul(ToReal(-2),kmadd(csetemp10,ToReal(3),kmadd(csetemp11,ToReal(3),kmul(csetemp12,ToReal(7))))),kmadd(csetemp11,kmul(csetemp12,ToReal(11)),kmul(csetemp10,kmadd(csetemp11,ToReal(8),kmul(csetemp12,ToReal(11)))))))))),ToReal(M)))))))))))))))))))))))));
    
    CCTK_REAL_VEC betap1 = 
      kmul(csetemp13,kmul(kadd(csetemp16,csetemp17),kmul(INV(kmul(kadd(csetemp15,csetemp9),kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M))))))),kmul(ToReal(2),ToReal(M)))));
    
    CCTK_REAL_VEC betap2 = 
      kmul(csetemp13,kmul(kadd(csetemp19,csetemp20),kmul(INV(kmul(kadd(csetemp15,csetemp9),kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M))))))),kmul(ToReal(2),ToReal(M)))));
    
    CCTK_REAL_VEC betap3 = 
      kmul(csetemp15,kmul(Z,kmul(INV(kadd(csetemp14,kmadd(csetemp12,csetemp9,kmul(csetemp13,kmul(ToReal(2),ToReal(M)))))),kmul(ToReal(2),ToReal(M)))));
    
    CCTK_REAL_VEC dtbetap1 = ToReal(0);
    
    CCTK_REAL_VEC dtbetap2 = ToReal(0);
    
    CCTK_REAL_VEC dtbetap3 = ToReal(0);
    
    CCTK_REAL_VEC csetemp58 = SQR(Jac11);
    
    CCTK_REAL_VEC csetemp59 = SQR(Jac21);
    
    CCTK_REAL_VEC csetemp60 = SQR(Jac31);
    
    CCTK_REAL_VEC gxxL = 
      kmadd(csetemp58,G11,kmadd(csetemp59,G22,kmadd(csetemp60,G33,kmul(kmadd(G32,kmul(Jac21,Jac31),kmul(Jac11,kmadd(G21,Jac21,kmul(G31,Jac31)))),ToReal(2)))));
    
    CCTK_REAL_VEC gxyL = 
      kmadd(Jac12,kmadd(G11,Jac11,kmadd(G21,Jac21,kmul(G31,Jac31))),kmadd(Jac22,kmadd(G21,Jac11,kmadd(G22,Jac21,kmul(G32,Jac31))),kmul(kmadd(G31,Jac11,kmadd(G32,Jac21,kmul(G33,Jac31))),Jac32)));
    
    CCTK_REAL_VEC gxzL = 
      kmadd(Jac13,kmadd(G11,Jac11,kmadd(G21,Jac21,kmul(G31,Jac31))),kmadd(Jac23,kmadd(G21,Jac11,kmadd(G22,Jac21,kmul(G32,Jac31))),kmul(kmadd(G31,Jac11,kmadd(G32,Jac21,kmul(G33,Jac31))),Jac33)));
    
    CCTK_REAL_VEC csetemp61 = SQR(Jac12);
    
    CCTK_REAL_VEC csetemp62 = SQR(Jac22);
    
    CCTK_REAL_VEC csetemp63 = SQR(Jac32);
    
    CCTK_REAL_VEC gyyL = 
      kmadd(csetemp61,G11,kmadd(csetemp62,G22,kmadd(csetemp63,G33,kmul(kmadd(G32,kmul(Jac22,Jac32),kmul(Jac12,kmadd(G21,Jac22,kmul(G31,Jac32)))),ToReal(2)))));
    
    CCTK_REAL_VEC gyzL = 
      kmadd(Jac13,kmadd(G11,Jac12,kmadd(G21,Jac22,kmul(G31,Jac32))),kmadd(Jac23,kmadd(G21,Jac12,kmadd(G22,Jac22,kmul(G32,Jac32))),kmul(kmadd(G31,Jac12,kmadd(G32,Jac22,kmul(G33,Jac32))),Jac33)));
    
    CCTK_REAL_VEC csetemp64 = SQR(Jac13);
    
    CCTK_REAL_VEC csetemp65 = SQR(Jac23);
    
    CCTK_REAL_VEC csetemp66 = SQR(Jac33);
    
    CCTK_REAL_VEC gzzL = 
      kmadd(csetemp64,G11,kmadd(csetemp65,G22,kmadd(csetemp66,G33,kmul(kmadd(G32,kmul(Jac23,Jac33),kmul(Jac13,kmadd(G21,Jac23,kmul(G31,Jac33)))),ToReal(2)))));
    
    CCTK_REAL_VEC kxxL = 
      kmadd(csetemp58,K11,kmadd(csetemp59,K22,kmadd(csetemp60,K33,kmul(kmadd(Jac11,kmadd(Jac21,K21,kmul(Jac31,K31)),kmul(Jac21,kmul(Jac31,K32))),ToReal(2)))));
    
    CCTK_REAL_VEC kxyL = 
      kmadd(Jac11,kmadd(Jac12,K11,kmadd(Jac22,K21,kmul(Jac32,K31))),kmadd(Jac21,kmadd(Jac12,K21,kmadd(Jac22,K22,kmul(Jac32,K32))),kmul(Jac31,kmadd(Jac12,K31,kmadd(Jac22,K32,kmul(Jac32,K33))))));
    
    CCTK_REAL_VEC kxzL = 
      kmadd(Jac11,kmadd(Jac13,K11,kmadd(Jac23,K21,kmul(Jac33,K31))),kmadd(Jac21,kmadd(Jac13,K21,kmadd(Jac23,K22,kmul(Jac33,K32))),kmul(Jac31,kmadd(Jac13,K31,kmadd(Jac23,K32,kmul(Jac33,K33))))));
    
    CCTK_REAL_VEC kyyL = 
      kmadd(csetemp61,K11,kmadd(csetemp62,K22,kmadd(csetemp63,K33,kmul(kmadd(Jac12,kmadd(Jac22,K21,kmul(Jac32,K31)),kmul(Jac22,kmul(Jac32,K32))),ToReal(2)))));
    
    CCTK_REAL_VEC kyzL = 
      kmadd(Jac12,kmadd(Jac13,K11,kmadd(Jac23,K21,kmul(Jac33,K31))),kmadd(Jac22,kmadd(Jac13,K21,kmadd(Jac23,K22,kmul(Jac33,K32))),kmul(Jac32,kmadd(Jac13,K31,kmadd(Jac23,K32,kmul(Jac33,K33))))));
    
    CCTK_REAL_VEC kzzL = 
      kmadd(csetemp64,K11,kmadd(csetemp65,K22,kmadd(csetemp66,K33,kmul(kmadd(Jac13,kmadd(Jac23,K21,kmul(Jac33,K31)),kmul(Jac23,kmul(Jac33,K32))),ToReal(2)))));
    
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
  LC_ENDLOOP3VEC(KerrSchild_initial);
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
