
/*

T01 = ( C1   -S1   0   0 )
      ( S1    C1   0   0 )
      (  0     0   1   0 )

T12 = ( C2   -S2   0   0 )
      (  0     0  -1   0 )
      ( S2    C2   0   0 )

T23 = ( C3   -S3   0   0 )
      (  0     0   1  r3 )
      (-S3   -C3   0   0 )

T34 = ( C4   -S4   0   0 )
      (  0     0   1   0 )
      (-S4   -C4   0   0 )

T45 = ( C5   -S5   0   0 )
      (  0     0  -1  r5 )
      ( S5    C5   0   0 )

T56 = ( C6   -S6   0   0 )
      (  0     0  -1   0 )
      ( S6    C6   0   0 )

T67 = ( C7   -S7   0   0 )
      (  0     0   1   0 )
      ( -S7  -C7   0   0 )

                  (  C7 -S7   0    0 )
                  (   0   0   1    0 )
                  ( -S7 -C7   0    0 )
		    	        
( C6 -S6   0   0 )(  d1  *  -S6    0 )
(  0   0  -1   0 )(  S7  *    0    0 )
( S6  C6   0   0 )(  d2  *   C6    0 )
		    	        
( C5 -S5   0   0 )(  d3  *   d5    0 )
(  0   0  -1 -r5 )( -d2  *  -C6  -r5 )
( S5  C5   0   0 )(  d4  *   d6    0 )
		    	        
( C4 -S4   0   0 )(  d7  *   d9  d11 )
(  0   0   1   0 )(  d4  *   d6    0 )
(-S4 -C4   0   0 )(  d8  *  d10  d12 )
		    	        
( C3 -S3   0   0 )( d13  *  d15  d17 )
(  0   0   1  r3 )(  d8  *  d10  d12 + r3)
(-S3 -C3   0   0 )( d14  *  d16  d18 )
		    	        
( C2 -S2   0   0 )( d19  *  d21  d23 )
(  0   0  -1   0 )(-d14  * -d16 -d18 )
( S2  C2   0   0 )( d20  *  d22  d24 )

( C1 -S1   0   0 )( d25  *  d27  d29 )
( S1  C1   0   0 )( d26  *  d28  d30 )
(  0   0   1   0 )( d20  *  d22  d24 )
*/

#include "math.h"
#include <stdio.h>
#include "gb.h"

// to eventually put in gb.c
void Gb_th_print(Gb_th *th, char* s)
{
  printf("%s= %g %g %g %g\n", s, th->vx.x, th->vy.x, th->vz.x, th->vp.x);
  printf("%s= %g %g %g %g\n", s, th->vz.y, th->vx.y, th->vy.y, th->vp.y);
  printf("%s= %g %g %g %g\n", s, th->vx.z, th->vy.z, th->vz.z, th->vp.z);
}

// end to eventually put in gb.c

/*
 * compute the configuration e1 e2 e3
 * Inputs : 
 *   r3, r5 : constants of the arm
 *   eq : articulation coordinates
 * output :
 *   e1 : 1 ou -1
 *   e2 : 1 ou -1
 *   e3 : 1 ou -1
 */
void kukaLBR_gete1e2e3(double r3, double r5, Gb_q7* eq,
		       int* e1, int* e2, int* e3)
{
  double s2 = sin(eq->q2);
  double c2 = cos(eq->q2);
  double c3 = cos(eq->q3);
  double s4 = sin(eq->q4);
  double c4 = cos(eq->q4);
  double s6 = sin(eq->q6);
  double d23 = c2 *c3 * s4 * r5 - s2 * (c4 * r5 + r3);
  *e1 = (d23 >= 0) ? 1. : -1.;
  *e2 = (s4  >= 0) ? 1. : -1.;
  *e3 = (s6  >= 0) ? 1. : -1.;
}

void kukaLBR_mgd(Gb_q7* Q, double r3, double r5, Gb_th* th07)
{
  double m;
  double C1 = cos (Q->q1);
  double C2 = cos (Q->q2);
  double C3 = cos (Q->q3);
  double C4 = cos (Q->q4);
  double C5 = cos (Q->q5);
  double C6 = cos (Q->q6);
  double C7 = cos (Q->q7);
  double S1 = sin (Q->q1);
  double S2 = sin (Q->q2);
  double S3 = sin (Q->q3);
  double S4 = sin (Q->q4);
  double S5 = sin (Q->q5);
  double S6 = sin (Q->q6);
  double S7 = sin (Q->q7);

  double d1  = C6 * C7;
  double d2  = S6 * C7;
  double d3  = C5 * d1 - S5 * S7;
  double d4  = S5 * d1 + C5 * S7;
  double d5  =-C5 * S6;
  double d6  =-S5 * S6;
  double d7  = C4 * d3 + S4 * d2;
  double d8  =-S4 * d3 + C4 * d2;
  double d9  = C4 * d5 + S4 * C6;
  double d10 =-S4 * d5 + C4 * C6;
  double d11 = S4 * r5;
  double d12 = C4 * r5;
  double d13 = C3 * d7 - S3 * d4;
  double d14 =-S3 * d7 - C3 * d4; 
  double d15 = C3 * d9 - S3 * d6;
  double d16 =-S3 * d9 - C3 * d6; 
  double d17 = C3 * d11;
  double d18 =-S3 * d11;
  double d19 = C2 * d13 - S2 * d8; 
  double d20 = S2 * d13 + C2 * d8;
  double d21 = C2 * d15 - S2 * d10;
  double d22 = S2 * d15 + C2 * d10;
  double d23 = C2 * d17 - S2 * (d12 + r3);
  double d24 = S2 * d17 + C2 * (d12 + r3);
  double d25 = C1 * d19 + S1 * d14;
  double d26 = S1 * d19 - C1 * d14;
  double d27 = C1 * d21 + S1 * d16;
  double d28 = S1 * d21 - C1 * d16;
  double d29 = C1 * d23 + S1 * d18;
  double d30 = S1 * d23 - C1 * d18;
  
//  printf("d1-10= %g %g %g %g %g  %g %g %g %g %g\n", 
//	 d1, d2, d3, d4, d5, d6, d7, d8, d9, d10);
//  printf("d11-20= %g %g %g %g %g  %g %g %g %g %g\n", 
//	 d11, d12, d13, d14, d15, d16, d17, d18, d19, d20);
//  printf("d21-30= %g %g %g %g %g  %g %g %g %g %g\n", 
//	 d21, d22, d23, d24, d25, d26, d27, d28, d29, d30);
  th07->vx.x = d25;
  th07->vx.y = d26;
  th07->vx.z = d20;
  th07->vz.x = d27;
  th07->vz.y = d28;
  th07->vz.z = d22;
  th07->vp.x = d29;
  th07->vp.y = d30;
  th07->vp.z = d24;
  m = th07->vx.x * th07->vx.x + th07->vx.y * th07->vx.y + th07->vx.z * th07->vx.z;
  th07->vx.x /= m;
  th07->vx.y /= m;
  th07->vx.z /= m;
  Gb_v3_cross_product (&(th07->vz), &(th07->vx), &(th07->vy));
  m = th07->vy.x * th07->vy.x + th07->vy.y * th07->vy.y + th07->vy.z * th07->vy.z;
  th07->vy.x /= m;
  th07->vy.y /= m;
  th07->vy.z /= m;
  Gb_v3_cross_product(&(th07->vx), &(th07->vy), &(th07->vz));
}


/* kukaLBR_direct computes geometry and differential direct models
//  T01 = ( C1   -S1   0   0 )
//        ( S1    C1   0   0 )
//        (  0     0   1   0 )
//  
//  T12 = ( C2   -S2   0   0 )
//        (  0     0  -1   0 )
//        ( S2    C2   0   0 )
//  
//  T23 = ( C3   -S3   0   0 )
//        (  0     0   1  r3 )
//        (-S3   -C3   0   0 )
//  
//  T34 = ( C4   -S4   0   0 )
//        (  0     0   1   0 )
//        (-S4   -C4   0   0 )
//  
//  T45 = ( C5   -S5   0   0 )
//        (  0     0  -1 -r5 )
//        ( S5    C5   0   0 )
//  
//  T56 = ( C6   -S6   0   0 )
//        (  0     0  -1   0 )
//        ( S6    C6   0   0 )
//  
//  T67 = ( C7   -S7   0   0 )
//        (  0     0   1   0 )
//        ( -S7  -C7   0   0 )
//  
//  
//                      T12 = ( C2   -S2   0   0 )
//                            (  0     0  -1   0 )
//                            ( S2    C2   0   0 )   f1= C1*C2
//  T01 = ( C1   -S1   0   0 )( f1    f3  S1   0 )   f2= S1*C2
//        ( S1    C1   0   0 )( f2    f4 -C1   0 )   f3=-C1*S2
//        (  0     0   1   0 )( S2    C2   0   0 )   f4=-S1*S2
//  
//                      T23 = ( C3   -S3   0   0 )   f5= f1*C3-S1*S3 
//                            (  0     0   1  r3 )   f6= f2*C3+C1*S3
//                            (-S3   -C3   0   0 )   f7= S2*C3
//  T02 = ( f1    f3  S1   0 )( f5    f8  f3 f11 )   f8=-f1*S3-S1*C3
//        ( f2    f4 -C1   0 )( f6    f9  f4 f12 )   f9=-f2*S3+C1*C3
//        ( S2    C2   0   0 )( f7   f10  C2 f13 )   f10=-S2*S3
//  						 f11= f3*r3
//                                                   f12= f4*r3
//  						 f13= C2*r3
//  
//                      T34 = ( C4   -S4   0   0 )   f14= f5*C4-f3*S4
//                            (  0     0   1   0 )   f15= f6*C4-f4*S4
//                            (-S4   -C4   0   0 )   f16= f7*C4-C2*S4
//  T03 = ( f5    f8  f3 f11 )( f14  f17  f8 f11 )   f17=-f5*S4-f3*C4 
//        ( f6    f9  f4 f12 )( f15  f18  f9 f12 )   f18=-f6*S4-f4*C4
//        ( f7   f10  C2 f13 )( f16  f19 f10 f13 )   f19=-f7*S4-C2*C4
//  
//                      T45 = ( C5   -S5   0   0 )   f20= f14*C5+f8*S5
//                            (  0     0  -1 -r5 )   f21= f15*C5+f9*S5
//                            ( S5    C5   0   0 )   f22= f16*C5+F10*S5
//  T04 = ( f14  f17  f8 f11 )( f20  f23-f17 f26 )   f23=-f14*S5+f8*C5
//        ( f15  f18  f9 f12 )( f21  f24-f18 f27 )   f24=-f15*S5+f9*C5
//        ( f16  f19 f10 f13 )( f22  f25-f19 f28 )   f25=-f16*S5+f10*C5
//                                                   f26=-f17*r5+f11
//                                                   f27=-f18*r5+f12
//                                                   f28=-f19*r5+f13
//  
//  
//                      T56 = ( C6   -S6   0   0 )   f29= f20*C6-f17*S6
//                            (  0     0  -1   0 )   f30= f21*C6-f18*S6
//                            ( S6    C6   0   0 )   f31= f22*C6-f19*S6
//  T05 = ( f20  f23-f17 f26 )( f29  f32-f23 f26 )   f32=-f20*S6-f17*C6 
//        ( f21  f24-f18 f27 )( f30  f33-f24 f27 )   f33=-f21*S6-f18*C6
//        ( f22  f25-f19 f28 )( f31  f34-f25 f28 )   f34=-f22*S6-f19*C6
//  
//  
//  
//                      T67 = ( C7   -S7   0   0 )   f35= f29*C7+f23*S7
//                            (  0     0   1   0 )   f36= f30*C7+f24*S7
//                            ( -S7  -C7   0   0 )   f37= f31*C7+f25*S7
//  T06 = ( f29  f32-f23 f26 )( f35  f38 f32 f26 )   f38=-f29*S7+f23*C7
//        ( f30  f33-f24 f27 )( f36  f39 f33 f27 )   f39=-f30*S7+f24*C7
//        ( f31  f34-f25 f28 )( f37  f40 f34 f28 )   f40=-f31*S7+f25*C7
//  
//  
//  P7= O_0O_7 (position of the wrist center relatively to frame R_0)
//  Pu= 0_3O_7 (position of the wrist center relatively to frame R_3)
//  
//  Jac=( z1^P7  z2^P7  z3^P7  z4^Pu    0    0   0  )
//      (   z1     z2     z3     z4    z5   z6  z7  )
//  
//  z1^P7 = | 0   |f26   |-f27
//          | 0 ^ |f27 = |f26
//          | 1   |f28   |0
//  
//  z2^P7 = |S1     |f26   |-C1*f28
//          |-C1  ^ |f27 = |-S1*f28
//          |0      |f28   |S1*f27+C1*f26
//  
//  z3^P7 = |f3    |f26   |f4*f28-C2f27
//          |f4  ^ |f27 = |C2*f26-f3*f28
//          |C2    |f28   |f3*f27-f4*f26
//  
//  z4^Pu = z4 ^(-r5 y4) = r5 y4 ^ z4 = r5 x4 = r5 * [f14  f15  f16]
//  
 */
void kukaLBR_direct(Gb_q7* Q, double r3, double r5, Gb_th* th07, Gb_jac7 jac7)
{
  double m;
  double C1 = cos (Q->q1);
  double C2 = cos (Q->q2);
  double C3 = cos (Q->q3);
  double C4 = cos (Q->q4);
  double C5 = cos (Q->q5);
  double C6 = cos (Q->q6);
  double C7 = cos (Q->q7);
  double S1 = sin (Q->q1);
  double S2 = sin (Q->q2);
  double S3 = sin (Q->q3);
  double S4 = sin (Q->q4);
  double S5 = sin (Q->q5);
  double S6 = sin (Q->q6);
  double S7 = sin (Q->q7);

  double  f1 = C1*C2;
  double  f2 = S1*C2;
  double  f3 =-C1*S2;
  double  f4 =-S1*S2;
  double  f5 = f1*C3-S1*S3;
  double  f6 = f2*C3+C1*S3;
  double  f7 = S2*C3;
  double  f8 =-f1*S3-S1*C3;
  double  f9 =-f2*S3+C1*C3;
  double f10 =-S2*S3;
  double f11 = f3*r3;
  double f12 = f4*r3;
  double f13 = C2*r3;
  double f14 = f5*C4-f3*S4;
  double f15 = f6*C4-f4*S4;
  double f16 = f7*C4-C2*S4;
  double f17 =-f5*S4-f3*C4;
  double f18 =-f6*S4-f4*C4;
  double f19 =-f7*S4-C2*C4;
  double f20 = f14*C5+f8*S5;
  double f21 = f15*C5+f9*S5;
  double f22 = f16*C5+f10*S5;
  double f23 =-f14*S5+f8*C5;
  double f24 =-f15*S5+f9*C5;
  double f25 =-f16*S5+f10*C5;
  double f26 =-f17*r5+f11;
  double f27 =-f18*r5+f12;
  double f28 =-f19*r5+f13;
  double f29 = f20*C6-f17*S6;
  double f30 = f21*C6-f18*S6;
  double f31 = f22*C6-f19*S6;
  double f32 =-f20*S6-f17*C6;
  double f33 =-f21*S6-f18*C6;
  double f34 =-f22*S6-f19*C6;
  double f35 = f29*C7+f23*S7;
  double f36 = f30*C7+f24*S7;
  double f37 = f31*C7+f25*S7;
  double f38 =-f29*S7+f23*C7;
  double f39 =-f30*S7+f24*C7;
  double f40 =-f31*S7+f25*C7;

  jac7.c1.x  = -f27;
  jac7.c1.y  = f26;
  jac7.c1.z  = 0;
  jac7.c1.rx = 0;
  jac7.c1.ry = 0;
  jac7.c1.rz = 1;
  jac7.c2.x  = -C1*f28;
  jac7.c2.y  = -S1*f28;
  jac7.c2.z  = S1*f27+C1*f26;
  jac7.c2.rx = S1;
  jac7.c2.ry = -C1;
  jac7.c2.rz = 0;
  jac7.c3.x  = f4*f28-C2*f27;
  jac7.c3.y  = C2*f26-f3*f28;
  jac7.c3.z  = f3*f27-f4*f26;
  jac7.c3.rx = f3;
  jac7.c3.ry = f4;
  jac7.c3.rz = C2;
  jac7.c4.x  = r5 * f14;
  jac7.c4.y  = r5 * f15;
  jac7.c4.z  = r5 * f16;
  jac7.c4.rx = f8;
  jac7.c4.ry = f9;
  jac7.c4.rz = f10;
  jac7.c5.x  = 0;
  jac7.c5.y  = 0;
  jac7.c5.z  = 0;
  jac7.c5.rx = -f17;
  jac7.c5.ry = -f18;
  jac7.c5.rz = -f19;
  jac7.c6.x  = 0;
  jac7.c6.y  = 0;
  jac7.c6.z  = 0;
  jac7.c6.rx = -f23;
  jac7.c6.ry = -f24;
  jac7.c6.rz = -f25;
  jac7.c7.x  = 0;
  jac7.c7.y  = 0;
  jac7.c7.z  = 0;
  jac7.c7.rx = f32;
  jac7.c7.ry = f33;
  jac7.c7.rz = f34;

  th07->vx.x = f35;
  th07->vx.y = f36;
  th07->vx.z = f37;
  th07->vy.x = f38;
  th07->vy.y = f39;
  th07->vy.z = f40;
  th07->vz.x = f32;
  th07->vz.y = f33;
  th07->vz.z = f34;
  th07->vp.x = f26;
  th07->vp.y = f27;
  th07->vp.z = f28;
}

/*
 X*X + Y*Y + Z*Z = d29*d29 + d30*d30 + d24 *d24
    =  C1*C1* d23*d23 + S1*S1* d18*d18 + 2* C1*d23*S1*d18
      +S1*S1* d23*d23 + C1*C1* d18*d18 - 2* S1*d23*C1*d18
      + d24*d24
    = d23*d23 + d18*d18 + d24*d24
    = d18*d18 + C2*C2 * d17*d17 + S2*S2 * (d12+r3)(d12+r3) -2* C2*d17*S2*(d12+r3)
              + S2*S2 * d17*d17 + C2*C2 * (d12+r3)(d12+r3) +2* S2*d17*C2*(d12+r3)
    = d18*d18 + d17*d17 + (d12+r3)(d12+r3)
    = S3*S3 * d11*d11 + C3*C3 * d11*d11 + (d12+r3)(d12+r3)
    = d11*d11 + (d12+r3)(d12+r3)
    = S4*S4 * r5*r5 + (C4 * r5 + r3)(C4 * r5 + r3)
    = S4*S4 * r5*r5 + C4*C4 * r5*r5 + r3*r3 + 2 * C4 * r5 * r3
    = r5*r5 + r3*r3 + 2 * C4 * r5 * r3
C4 = d24*d24 + d29*d29 + d30*d30 - r5*r5 - r3*r3 / 2. / r5 / r3
===> S4 et Q4 avec Deux solutions

  d1  = C6 * C7;
  d2  = S6 * C7;
  d3  = C5 * d1 - S5 * S7;
  d4  = S5 * d1 + C5 * S7;
  d5  =-C5 * S6;
  d6  =-S5 * S6;
  d7  = C4 * d3 + S4 * d2;
  d8  =-S4 * d3 + C4 * d2;
  d9  = C4 * d5 + S4 * C6;
  d10 =-S4 * d5 + C4 * C6;
  d11 = S4 * r5;
  d12 = C4 * r5;
  d13 = C3 * d7 - S3 * d4;
  d14 =-S3 * d7 - C3 * d4; 
  d15 = C3 * d9 - S3 * d6;
  d16 =-S3 * d9 - C3 * d6; 
  d17 = C3 * d11;
  d18 =-S3 * d11;
  d19 = C2 * d13 - S2 * d8; 
  d20 = S2 * d13 + C2 * d8;
  d21 = C2 * d15 - S2 * d10;
  d22 = S2 * d15 + C2 * d10;
  d23 = C2 * d17 - S2 * (d12 + r3);
  d24 = S2 * d17 + C2 * (d12 + r3);
  d25 = C1 * d19 + S1 * d14;
  d26 = S1 * d19 - C1 * d14;
  d27 = C1 * d21 + S1 * d16;
  d28 = S1 * d21 - C1 * d16;
  d29 = C1 * d23 + S1 * d18;
  d30 = S1 * d23 - C1 * d18;
 */

int kukaLBR_mgi_q3(Gb_th* th07, Gb_q7* Qp, double r3, double r5, double epsilon,
		   Gb_q7* q) 
{
  double C1, C2, C3, C4, C5, C6, C7;
  double S1, S2, S3, S4, S5, S6, S7;
  double d25 = th07->vx.x;
  double d26 = th07->vx.y;
  double d20 = th07->vx.z;
  double d27 = th07->vz.x;
  double d28 = th07->vz.y;
  double d22 = th07->vz.z;
  double d29 = th07->vp.x;
  double d30 = th07->vp.y;
  double d24 = th07->vp.z;
  double d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15;
  double d16, d17, d18, d19, d21, d23;
  q->q3 = Qp->q3;
  C3 = cos(q->q3);
  S3 = sin(q->q3);

  C4 = ( d24*d24 + d29*d29 + d30*d30 - r5*r5 - r3*r3 ) / 2. / r5 / r3;
  if ( C4 < -1.-epsilon ) { 
    return -1;
  } else if ( C4 < -1 )  { 
    C4 = -1;
    S4 = 0;
    q->q4 = M_PI;
  } else if ( C4 < 1 )  {
    S4 = sqrt(1 - C4*C4);
    if (sin(Qp->q4) < 0.) S4 = -S4;
    q->q4 = atan2(S4, C4);
    //if ?? autre solution
  } else if ( C4 < 1 + epsilon ) {
    C4 = 1;
    S4 = 0;
    q->q4 = 0;
  } else { 
    return -1;
  }
  d11 = S4 * r5;
  d12 = C4 * r5;
  d17 = C3 * d11;
  d18 = -S3 * d11;
  if ( (d29*d29 + d30*d30 -d18*d18) < -epsilon) {
    return -2;
  } else if ( (d29*d29 + d30*d30 -d18*d18) < 0) {
    d23 = 0;
    // approximated ??
  } else {
    d23 = sqrt(d29*d29 + d30*d30 -d18*d18);
  }
  //if ?? autre solution
  if ((d23*d23 + d18*d18) == 0) {
    // return value  singular configuration
    q->q1 = Qp->q1;
    C1 = cos(Qp->q1);
    S1 = sin(Qp->q1);
  } else {
    C1 = (d23*d29 - d18*d30) / (d23*d23 + d18*d18);
    S1 = (d18*d29 + d23*d30) / (d23*d23 + d18*d18);
    q->q1 = atan2(S1, C1);
  }
  d21 = C1 * d27 + S1 * d28;
  d16 = S1 * d27 - C1 * d28;
  d19 = C1 * d25 + S1 * d26;
  d14 = S1 * d25 - C1 * d26;
  if (( d17*d17+(d12+r3)*(d12+r3) ) == 0) {
    // return value  singular configuration
    q->q2 = Qp->q2;
    C2 = cos(Qp->q2);
    S2 = sin(Qp->q2);
  } else {
    C2 = (d17*d23 + (d12+r3)*d24) / (d17*d17+(d12+r3)*(d12+r3));
    S2 = (-(d12+r3)*d23 + d17*d24) / (d17*d17+(d12+r3)*(d12+r3));
    q->q2 = atan2(S2, C2);
  }
  d15 = C2*d21+S2*d22;
  d10 =-S2*d21+C2*d22;
  d13 = C2*d19 + S2*d20;
  d8 = -S2*d19 + C2*d20;
  d9 = C3*d15 - S3*d16;
  d6 = -S3*d15 - C3*d16;
  d7 = C3*d13 - S3*d14;
  d4 = -S3*d13 - C3*d14;
  d5 = C4*d9 - S4*d10;
  C6 = S4*d9 + C4*d10;
  d3 = C4*d7 - S4*d8;
  d2 = S4*d7 + C4*d8;
  S6= sqrt(d5*d5 + d6*d6);
  if (sin(Qp->q6) < 0) S6 = -S6;
  q->q6 = atan2(S6, C6);
  if (S6 == 0) {
    // return value  singular configuration
    q->q6 = Qp->q6;
    C6 = cos(Qp->q6);
    S6 = sin(Qp->q6);
  } else {
    C5 = -d5/S6;
    S5 = -d6/S6;
    q->q5 = atan2(S5, C5);
  }
  d1 = C5*d3 + S5*d4;
  S7 = -S5*d3 + C5*d4;
  C7 = C6*d1 + S6*d2;
  q->q7 = atan2(S7, C7);
}

Gb_statusMGI kukaLBR_mgi_q_e(Gb_th* th07, Gb_q7* Qp, double r3, double r5,
			     double epsilon, int e1, int e2, int e3,
			     Gb_q7* q) 
{
  Gb_statusMGI ret = MGI_OK;
  double C1, C2, C3, C4, C5, C6, C7;
  double S1, S2, S3, S4, S5, S6, S7;
  double d25 = th07->vx.x;
  double d26 = th07->vx.y;
  double d20 = th07->vx.z;
  double d27 = th07->vz.x;
  double d28 = th07->vz.y;
  double d22 = th07->vz.z;
  double d29 = th07->vp.x;
  double d30 = th07->vp.y;
  double d24 = th07->vp.z;
  double d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15;
  double d16, d17, d18, d19, d21, d23;
  q->q3 = Qp->q3;
  C3 = cos(q->q3);
  S3 = sin(q->q3);

  C4 = ( d24*d24 + d29*d29 + d30*d30 - r5*r5 - r3*r3 ) / 2. / r5 / r3;
  if ( C4 < -1.-epsilon ) { 
    return MGI_ERROR;
  } else if ( C4 < -1 )  { 
    C4 = -1;
    S4 = 0;
    q->q4 = M_PI;
  } else if ( C4 < 1 )  {
    S4 = sqrt(1 - C4*C4);
    S4 *= e2;
    q->q4 = atan2(S4, C4);
  } else if ( C4 < 1 + epsilon ) {
    C4 = 1;
    S4 = 0;
    q->q4 = 0;
  } else { 
    return MGI_ERROR;
    // here it is possible to reteurn an approximated q...
  }
  d11 = S4 * r5;
  d12 = C4 * r5;
  d17 = C3 * d11;
  d18 = -S3 * d11;
  if ( (d29*d29 + d30*d30 -d18*d18) < -epsilon) {
    return MGI_ERROR;
  } else if ( (d29*d29 + d30*d30 -d18*d18) < 0) {
    d23 = 0;
  } else {
    d23 = sqrt(d29*d29 + d30*d30 -d18*d18);
    d23 *= e1;
  }
  if ((d23*d23 + d18*d18) == 0) {
    ret = MGI_SINGULAR;
    q->q1 = Qp->q1;
    C1 = cos(Qp->q1);
    S1 = sin(Qp->q1);
  } else {
    C1 = (d23*d29 - d18*d30) / (d23*d23 + d18*d18);
    S1 = (d18*d29 + d23*d30) / (d23*d23 + d18*d18);
    q->q1 = atan2(S1, C1);
  }
  d21 = C1 * d27 + S1 * d28;
  d16 = S1 * d27 - C1 * d28;
  d19 = C1 * d25 + S1 * d26;
  d14 = S1 * d25 - C1 * d26;
  if (( d17*d17+(d12+r3)*(d12+r3) ) == 0) {
    ret = MGI_SINGULAR;
    q->q2 = Qp->q2;
    C2 = cos(Qp->q2);
    S2 = sin(Qp->q2);
  } else {
    C2 = (d17*d23 + (d12+r3)*d24) / (d17*d17+(d12+r3)*(d12+r3));
    S2 = (-(d12+r3)*d23 + d17*d24) / (d17*d17+(d12+r3)*(d12+r3));
    q->q2 = atan2(S2, C2);
  }
  d15 = C2*d21+S2*d22;
  d10 =-S2*d21+C2*d22;
  d13 = C2*d19 + S2*d20;
  d8 = -S2*d19 + C2*d20;
  d9 = C3*d15 - S3*d16;
  d6 = -S3*d15 - C3*d16;
  d7 = C3*d13 - S3*d14;
  d4 = -S3*d13 - C3*d14;
  d5 = C4*d9 - S4*d10;
  C6 = S4*d9 + C4*d10;
  d3 = C4*d7 - S4*d8;
  d2 = S4*d7 + C4*d8;
  S6= sqrt(d5*d5 + d6*d6);
  S6 *= e3;
  q->q6 = atan2(S6, C6);
  if (S6 == 0) {
    ret = MGI_SINGULAR;
    q->q6 = Qp->q6;
    C6 = cos(Qp->q6);
    S6 = sin(Qp->q6);
  } else {
    C5 = -d5/S6;
    S5 = -d6/S6;
    q->q5 = atan2(S5, C5);
  }
  d1 = C5*d3 + S5*d4;
  S7 = -S5*d3 + C5*d4;
  C7 = C6*d1 + S6*d2;
  q->q7 = atan2(S7, C7);
  return ret;
}

/*
int main(int argc, char** argv) 
{
  Gb_q7 q, qs;
  double r3 = 0.4;
  double r5 = 0.39;
  Gb_th th07;
  Gb_th thp;
  double epsilon = 1e-7;

//  q.q1 = 0;
//  q.q2 = 0;
//  q.q3 = 0;
//  q.q4 = 0;
//  q.q5 = 0;
//  q.q6 = 0;
//  q.q7 = 0;
  q.q1 = M_PI / 7.;
  q.q2 =-M_PI / 5.;
  q.q3 = M_PI / 7.;
  q.q4 =-M_PI / 15.;
  q.q5 = M_PI / 7.;
  q.q6 =-M_PI / 8.;
  q.q7 = M_PI / 7.;
//  q.q1 = M_PI / 7.;
//  q.q2 =-M_PI / 5.;
//  q.q3 = M_PI / 7.;
//  q.q4 =-M_PI / 5.;
//  q.q5 = M_PI / 7.;
//  q.q6 =-M_PI / 8.;
//  q.q7 = M_PI / 7.;
//
  printf("q= %g %g %g  %g %g %g  %g\n", q.q1, q.q2, q.q3, q.q4, q.q5, q.q6, q.q7);

  kukaLBR_mgd(&q, r3, r5, &th07);
  Gb_th_print(&th07, "th07");
  kukaLBR_mgi_q3(&th07, &q, r3, r5, epsilon, &qs);

  printf("qs= %g %g %g  %g %g %g  %g\n", 
	 qs.q1, qs.q2, qs.q3, qs.q4, qs.q5, qs.q6, qs.q7);

  int e1, e2, e3;
  Gb_statusMGI status;
  kukaLBR_gete1e2e3(r3, r5, &q, &e1, &e2, &e3);
  status = kukaLBR_mgi_q_e(&th07, &q, r3, r5, epsilon, e1, e2, e3, &qs);
  printf("%s: %d%d%d qs= %g %g %g  %g %g %g  %g\n", 
	 Gb_statusMGI_s(status), e1, e2, e3,
	 qs.q1, qs.q2, qs.q3, qs.q4, qs.q5, qs.q6, qs.q7);
  kukaLBR_mgd(&qs, r3, r5, &thp);   Gb_th_print(&thp, "th07");
  printf("\n");
  e1 = 1; e2 = 1; e3 = 1;
  status = kukaLBR_mgi_q_e(&th07, &q, r3, r5, epsilon, e1, e2, e3, &qs);
  printf("%s: %d%d%d qs= %g %g %g  %g %g %g  %g\n", 
	 Gb_statusMGI_s(status), e1, e2, e3,
	 qs.q1, qs.q2, qs.q3, qs.q4, qs.q5, qs.q6, qs.q7);
  kukaLBR_mgd(&qs, r3, r5, &thp);   Gb_th_print(&thp, "th07");
  e1 = 1; e2 = 1; e3 =-1;
  status = kukaLBR_mgi_q_e(&th07, &q, r3, r5, epsilon, e1, e2, e3, &qs);
  printf("%s: %d%d%d qs= %g %g %g  %g %g %g  %g\n", 
	 Gb_statusMGI_s(status), e1, e2, e3,
	 qs.q1, qs.q2, qs.q3, qs.q4, qs.q5, qs.q6, qs.q7);
  kukaLBR_mgd(&qs, r3, r5, &thp);   Gb_th_print(&thp, "th07");
  e1 = 1; e2 =-1; e3 = 1;
  status = kukaLBR_mgi_q_e(&th07, &q, r3, r5, epsilon, e1, e2, e3, &qs);
  printf("%s: %d%d%d qs= %g %g %g  %g %g %g  %g\n", 
	 Gb_statusMGI_s(status), e1, e2, e3,
	 qs.q1, qs.q2, qs.q3, qs.q4, qs.q5, qs.q6, qs.q7);
  kukaLBR_mgd(&qs, r3, r5, &thp);   Gb_th_print(&thp, "th07");
  e1 = 1; e2 =-1; e3 =-1;
  status = kukaLBR_mgi_q_e(&th07, &q, r3, r5, epsilon, e1, e2, e3, &qs);
  printf("%s: %d%d%d qs= %g %g %g  %g %g %g  %g\n", 
	 Gb_statusMGI_s(status), e1, e2, e3,
	 qs.q1, qs.q2, qs.q3, qs.q4, qs.q5, qs.q6, qs.q7);
  kukaLBR_mgd(&qs, r3, r5, &thp);   Gb_th_print(&thp, "th07");
  e1 =-1; e2 = 1; e3 = 1;
  status = kukaLBR_mgi_q_e(&th07, &q, r3, r5, epsilon, e1, e2, e3, &qs);
  printf("%s: %d%d%d qs= %g %g %g  %g %g %g  %g\n", 
	 Gb_statusMGI_s(status), e1, e2, e3,
	 qs.q1, qs.q2, qs.q3, qs.q4, qs.q5, qs.q6, qs.q7);
  kukaLBR_mgd(&qs, r3, r5, &thp);   Gb_th_print(&thp, "th07");
  e1 =-1; e2 = 1; e3 =-1;
  status = kukaLBR_mgi_q_e(&th07, &q, r3, r5, epsilon, e1, e2, e3, &qs);
  printf("%s: %d%d%d qs= %g %g %g  %g %g %g  %g\n", 
	 Gb_statusMGI_s(status), e1, e2, e3,
	 qs.q1, qs.q2, qs.q3, qs.q4, qs.q5, qs.q6, qs.q7);
  kukaLBR_mgd(&qs, r3, r5, &thp);   Gb_th_print(&thp, "th07");
  e1 =-1; e2 =-1; e3 = 1;
  status = kukaLBR_mgi_q_e(&th07, &q, r3, r5, epsilon, e1, e2, e3, &qs);
  printf("%s: %d%d%d qs= %g %g %g  %g %g %g  %g\n", 
	 Gb_statusMGI_s(status), e1, e2, e3,
	 qs.q1, qs.q2, qs.q3, qs.q4, qs.q5, qs.q6, qs.q7);
  kukaLBR_mgd(&qs, r3, r5, &thp);   Gb_th_print(&thp, "th07");
  e1 =-1; e2 =-1; e3 =-1;
  status = kukaLBR_mgi_q_e(&th07, &q, r3, r5, epsilon, e1, e2, e3, &qs);
  printf("%s: %d%d%d qs= %g %g %g  %g %g %g  %g\n", 
	 Gb_statusMGI_s(status), e1, e2, e3,
	 qs.q1, qs.q2, qs.q3, qs.q4, qs.q5, qs.q6, qs.q7);
  kukaLBR_mgd(&qs, r3, r5, &thp);   Gb_th_print(&thp, "th07");
  return 0;
}*/

