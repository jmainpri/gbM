
/*

T01 = ( C1   -S1   0   0 )
      ( S1    C1   0   0 )
      (  0     0   1   0 )

T12 = ( C2   -S2   0   0 )
      (  0     0  -1   0 )
      ( S2    C2   0   0 )

T23 = ( C3   -S3   0   0 )
      (  0     0   1   0 )
      (-S3   -C3   0   0 )

T34 = ( C4   -S4   0   0 )
      (  0     0   1   0 )
      (-S4   -C4   0   0 )

T45 = ( C5   -S5   0   0 )
      (  0     0  -1   0 )
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
typedef struct Gb_q7 {
  double q1;
  double q2;
  double q3;
  double q4;
  double q5;
  double q6;
  double q7;
} Gb_q7;
typedef struct Gb_v3 {
  double x;
  double y;
  double z;
} Gb_v3;
typedef struct Gb_th {
  Gb_v3 vx;
  Gb_v3 vy;
  Gb_v3 vz;
  Gb_v3 vp;
} Gb_th;

#include "math.h"

void kukaLBR_mgd(Gb_q7* Q, Gb_q7* Q0, double r3, double r5, Gb_th* th07)
{
  double C1 = cos (Q->q1 - Q0->q1);
  double C2 = cos (Q->q2 - Q0->q2);
  double C3 = cos (Q->q3 - Q0->q3);
  double C4 = cos (Q->q4 - Q0->q4);
  double C5 = cos (Q->q5 - Q0->q5);
  double C6 = cos (Q->q6 - Q0->q6);
  double C7 = cos (Q->q7 - Q0->q7);
  double S1 = sin (Q->q1 - Q0->q1);
  double S2 = sin (Q->q2 - Q0->q2);
  double S3 = sin (Q->q3 - Q0->q3);
  double S4 = sin (Q->q4 - Q0->q4);
  double S5 = sin (Q->q5 - Q0->q5);
  double S6 = sin (Q->q6 - Q0->q6);
  double S7 = sin (Q->q7 - Q0->q7);

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
  
  th07->vx.x = d25;
  th07->vx.y = d26;
  th07->vx.z = d20;
  th07->vy.x = d28 * d20 - d22 * d26;
  th07->vy.y = d22 * d25 - d27 * d20;
  th07->vy.z = d27 * d26 - d28 * d25;
  th07->vz.x = d27;
  th07->vz.y = d28;
  th07->vz.z = d22;
  th07->vp.x = d29;
  th07->vp.y = d30;
  th07->vp.z = d24;
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
    if (Qp->q4 < 0.) S4 = -S4;
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
  if ( (d29*d29 + d30*d30 -d18*d18) < -epsilon) {
    return -2;
  } else if ( (d29*d29 + d30*d30 -d18*d18) < 0) {
    d23 = 0;
    // approximated ??
  } else {
    d23 = sqrt(d29*d29 + d30*d30 -d18*d18);
   pas top     +d18 puis -d18
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
  if (( d17*d17+(d12+r3)*(d12+r3) ) == 0) {
    // return value  singular configuration
    q->q2 = Qp->q2;
    C2 = cos(Qp->q2);
    S2 = sin(Qp->q2);
  } else {
    C2 = (d17*d23 + (d12+r3)*d24) / (d17*d17+(d12+r3)*(d12+r3));
    S2 = ((d12+r3)*d23 + d17*d24) / (d17*d17+(d12+r3)*(d12+r3));
    q->q2 = atan2(S2, C2);
  }
  d15 = C2*d21+S2*d22;
  d10 =-S2*d21+C2*d22;
  d13 = C2*d19 + S2*d20;
  d8 = -S2*d19 + C2*d20;
  d9 = C3*d15 - S3*d16;
  d6 = -S3*d15 - C3*d6;
  d7 = C3*d13 - S3*d14;
  d4 = -S3*d13 - C3*d14;
  d5 = C4*d9 - S4*d10;
  C6 = S4*d9 + C4*d10;
  d3 = C4*d7 - S4*d8;
  d2 = S4*d7 + C4*d8;
  S6= sqrt(d5*d5 + d6*d6);
  // if ?? autre solution
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

#include <stdio.h>

int main(int argc, char** argv) 
{
  Gb_q7 q, q0, qs;
  double r3 = 0.4;
  double r5 = 0.39;
  Gb_th th07;
  double epsilon = 1e-7;

  q0.q1 = 0.0;
  q0.q2 = 0.0;
  q0.q3 = 0.0;
  q0.q4 = 0.0;
  q0.q5 = 0.0;
  q0.q6 = 0.0;
  q0.q7 = 0.0;

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
  q.q4 =-M_PI / 5.;
  q.q5 = M_PI / 7.;
  q.q6 =-M_PI / 8.;
  q.q7 = M_PI / 7.;

  printf("q= %g %g %g  %g %g %g  %g\n", q.q1, q.q2, q.q3, q.q4, q.q5, q.q6, q.q7);

  kukaLBR_mgd(&q, &q0, r3, r5, &th07);
  kukaLBR_mgi_q3(&th07, &q, r3, r5, epsilon, &qs);

  printf("qs= %g %g %g  %g %g %g  %g\n", 
	 qs.q1, qs.q2, qs.q3, qs.q4, qs.q5, qs.q6, qs.q7);
  return 0;
}
