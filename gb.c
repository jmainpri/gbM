/*
 * Copyright (c) 2002 LAAS/CNRS -- RIA --
 * Daniel SIDOBRE -- mai 2002
 */

#include <math.h>
#include <stdio.h>

#include "gb.h"

#ifdef VXWORKS 
#define	M_PI		3.14159265358979323846
#define	M_PI_2		1.57079632679489661923
#define	M_PI_4		0.78539816339744830962
#include <stdlib.h> 
#endif

#include "SProto_gb.h"

/* 
 * Gb_v3_norme : Calcule le vecteur normé s defini par e, 
 *   si e=0 retourne (0,0,1)
 * accepte s=e
 */
double Gb_v3_norme(const Gb_v3* e, Gb_v3* s)
{
  double w = sqrt(e->x*e->x + e->y*e->y + e->z*e->z);
  if(w == 0.) {
    s->x = 0;
    s->y = 0;
    s->z = 1.;
  } else {
    s->x = e->x / w;
    s->y = e->y / w;
    s->z = e->z / w;
  }
  return w;
}

double Gb_prs(const double w1[3], const double w2[3])
{
  return(w1[0]*w2[0]+w1[1]*w2[1]+w1[2]*w2[2]);
}

double Gb_v3_prs(const Gb_v3* v1, const Gb_v3* v2)
{
  return(v1->x * v2->x + v1->y * v2->y + v1->z * v2->z);
}

/*
 * Gb_v3_cross_product : compute the cros product w of u and v
 * accepte output=u ou v
 */
void Gb_v3_cross_product(const Gb_v3* u, const Gb_v3* v, Gb_v3* output)
{
  double x, y, z; /* Pour output = u ou v */
  x = u->y * v->z - u->z * v->y;
  y = u->z * v->x - u->x * v->z;
  z = u->x * v->y - u->y * v->x;
  output->x = x;
  output->y = y;
  output->z = z;
}

double Gb_v3_product(const Gb_v3* u, const Gb_v3* v)
{
  return u->x * v->x + u->y * v->y + u->z * v->z;
}

double Gb_v3_module(const Gb_v3* u)
{
  return sqrt(u->x * u->x + u->y * u->y + u->z * u->z);
}

void Gb_v3_plus(const Gb_v3* u, const Gb_v3* v, Gb_v3* output)
{
  output->x = u->x + v->x;
  output->y = u->y + v->y;
  output->z = u->z + v->z;
}

void Gb_v3_moins(const Gb_v3* u, const Gb_v3* v, Gb_v3* output)
{
  output->x = u->x - v->x;
  output->y = u->y - v->y;
  output->z = u->z - v->z;
}

void Gb_v3_product_r(const Gb_v3* u, double r, Gb_v3* output)
{
  output->x = u->x * r;
  output->y = u->y * r;
  output->z = u->z * r;
} 

void Gb_dep_set(Gb_dep *dep, double x, double y, double z, 
		double rx, double ry, double rz, double a)
{
  dep->x  = x;
  dep->y  = y;
  dep->z  = z;
  dep->rx = rx;
  dep->ry = ry;
  dep->rz = rz;
  dep->a  = a;
}

void Gb_dep_get(Gb_dep *dep, double *x, double *y, double *z, 
		double *rx, double *ry, double *rz, double *a)
{
  *x  = dep->x ;
  *y  = dep->y ;
  *z  = dep->z ;
  *rx = dep->rx;
  *ry = dep->ry;
  *rz = dep->rz;
  *a  = dep->a ;
}


/* Calcul d'une matrice de préproduit vectoriel */
void Gb_m33_ppv(const Gb_v3 *ve, Gb_m33 *ms)
{
  ms->vx.x = 0;
  ms->vx.y = ve->z;
  ms->vx.z =-ve->y;
  ms->vy.x =-ve->z;
  ms->vy.y = 0;
  ms->vy.z = ve->x;
  ms->vz.x = ve->y;
  ms->vz.y =-ve->x;
  ms->vz.z = 0;
}

void Gb_dep_quat(const Gb_dep* dep, Gb_quat* q)
{
  double m, sa2;
  q->x = dep->x;
  q->y = dep->y;
  q->z = dep->z;
  m = sqrt(dep->rx*dep->rx + dep->ry*dep->ry + dep->rz*dep->rz);
  if (m == 0) {
    q->vx = 0;
    q->vy = 0;
    q->vz = 0;
    q->w  = 1;
    return;
  }
  q->w = cos(dep->a / 2.);
  sa2  = sin(dep->a / 2.);
  q->vx = dep->rx / m * sa2;
  q->vy = dep->ry / m * sa2;
  q->vz = dep->rz / m * sa2;
}

void Gb_quat_dep(const Gb_quat* q, Gb_dep* dep)
{
  double sinq;
  dep->x = q->x;
  dep->y = q->y;
  dep->z = q->z;
  dep->a = acos(q->w) * 2.;
  sinq = sqrt(1. - q->w * q->w);
  if (sinq == 0) {
    dep->rx = 0;
    dep->ry = 0;
    dep->rz = 0;
  } else {
    dep->rx = q->vx / sinq;
    dep->ry = q->vy / sinq;
    dep->rz = q->vz / sinq;
  }
}


void Gb_quat_th(const Gb_quat* q, Gb_th* th)
{
  /* peu teste */
  /* la version de la faq matrice et quaternion comporte une 
   *   multiplication de moins (ww est inutile) mais au moins
   *   3 multiplication par deux de plus.
   */
  double ww, wx, wy, wz, xy, xz, yz;
  ww = q->w * q->w;
  wx =-q->w * q->vx;
  wy =-q->w * q->vy;
  wz =-q->w * q->vz;
  xy = q->vx * q->vy;
  xz = q->vx * q->vz;
  yz = q->vy * q->vz;
  th->vx.x = 2 * (ww + q->vx * q->vx) - 1.;
  th->vx.y = 2 * (xy - wz);
  th->vx.z = 2 * (xz + wy);
  th->vy.x = 2 * (xy + wz);
  th->vy.y = 2 * (ww + q->vy * q->vy) - 1.;
  th->vy.z = 2 * (yz - wx);
  th->vz.x = 2 * (xz - wy);
  th->vz.y = 2 * (yz + wx);
  th->vz.z = 2 * (ww + q->vz * q->vz) -1.;
  th->vp.x = q->x;
  th->vp.y = q->y;
  th->vp.z = q->z;
}

void Gb_dep_th(const Gb_dep* dep, Gb_th* th)
{
  Gb_quat q;
  Gb_dep_quat(dep, &q);
  Gb_quat_th(&q, th);
}

void Gb_th_dep(const Gb_th* th, Gb_dep* dep)
{
  Gb_quat q;
  Gb_th_quat(th, &q);
  Gb_quat_dep(&q, dep);
}

void Gb_th_quat(const Gb_th* th, Gb_quat* q)
{
  double s0, s1, s2, s3;
  q->x = th->vp.x;
  q->y = th->vp.y;
  q->z = th->vp.z;
  s0 = th->vx.x + th->vy.y + th->vz.z;
  if (s0 > -0.1) {
    q->w = sqrt(1 + s0) / 2.;
    q->vx = (th->vy.z - th->vz.y) / 4. / q->w;
    q->vy = (th->vz.x - th->vx.z) / 4. / q->w;
    q->vz = (th->vx.y - th->vy.x) / 4. / q->w;
  } else {
    s1 = th->vx.x - th->vy.y - th->vz.z;
    if (s1 > -0.1) {
      q->vx = sqrt(1 + s1) / 2.;
      q->w  = (th->vy.z - th->vz.y) / 4. / q->vx;
      q->vy = (th->vx.y + th->vy.x) / 4. / q->vx;
      q->vz = (th->vx.z + th->vz.x) / 4. / q->vx;
    } else {
      s2 = - th->vx.x + th->vy.y - th->vz.z;
      if (s2 > -0.1) {
	q->vy = sqrt(1 + s2) / 2.;
	q->w  = (th->vz.x - th->vx.z) / 4. / q->vy;
	q->vx = (th->vx.y + th->vy.x) / 4. / q->vy;
	q->vz = (th->vz.y + th->vy.z) / 4. / q->vy;
      } else {
	s3 = - th->vx.x - th->vy.y + th->vz.z;
	if (s3 < -0.1) {
	  fprintf(stderr, "Pb Gb_th_quat s3 < -0.1 !!!\n");
	}
	q->vz = sqrt(1 + s3) / 2.;
	q->w  = (th->vx.y - th->vy.x) / 4. / q->vz;
	q->vx = (th->vx.z + th->vz.x) / 4. / q->vz;
	q->vy = (th->vz.y + th->vy.z) / 4. / q->vz;
      }
    }
  }
}


void Gb_v6_plus(const Gb_v6* a, const Gb_v6* b, Gb_v6* s)
{
  s->x  = a->x  + b->x ;
  s->y  = a->y  + b->y ;
  s->z  = a->z  + b->z ;
  s->rx = a->rx + b->rx;
  s->ry = a->ry + b->ry;
  s->rz = a->rz + b->rz;
}

void Gb_v6_moins(const Gb_v6* a, const Gb_v6* b, Gb_v6* s)
{
  s->x  = a->x  - b->x ;
  s->y  = a->y  - b->y ;
  s->z  = a->z  - b->z ;
  s->rx = a->rx - b->rx;
  s->ry = a->ry - b->ry;
  s->rz = a->rz - b->rz;
}

double Gb_v6_module (const Gb_v6* a)
{
  return sqrt(a->x  * a->x  +
	      a->y  * a->y  +
	      a->z  * a->z  +
	      a->rx * a->rx +
	      a->ry * a->ry +
	      a->rz * a->rz);
}

Gb_v3* Gb_v6_get_t(Gb_v6* v)
{
  return (Gb_v3*) &v->x;
}

/*const*/
Gb_v3* Gb_v6c_get_t(const Gb_v6* v)
{
  return (Gb_v3*) &v->x;
}

Gb_v3* Gb_v6_get_r(Gb_v6* v)
{
  return (Gb_v3*) &v->rx;
}

/*const*/
Gb_v3* Gb_v6c_get_r(const Gb_v6* v)
{
  return (Gb_v3*) &v->rx;
}



void Gb_th_x_v3(const Gb_th* th, const Gb_v3* v, Gb_v3* vs)
{
  Gb_v3 v3;
  Gb_v3* vu;
  if (v == vs)
    vu = &v3;
  else
    vu = vs;
  vu->x = th->vx.x * v->x + th->vy.x * v->y + th->vz.x * v->z;
  vu->y = th->vx.y * v->x + th->vy.y * v->y + th->vz.y * v->z;
  vu->z = th->vx.z * v->x + th->vy.z * v->y + th->vz.z * v->z;
  if (vu != vs) {
    vs->x = vu->x;
    vs->y = vu->y;
    vs->z = vu->z;
  }
}

void Gb_thInv_x_v3(const Gb_th* th, const Gb_v3* v, Gb_v3* vs)
{
  Gb_v3 v3;
  Gb_v3* vu;
  if (v == vs)
    vu = &v3;
  else
    vu = vs;
  vs->x = th->vx.x * v->x + th->vx.y * v->y + th->vx.z * v->z;
  vs->y = th->vy.x * v->x + th->vy.y * v->y + th->vy.z * v->z;
  vs->z = th->vz.x * v->x + th->vz.y * v->y + th->vz.z * v->z;
  if (vu != vs) {
    vs->x = vu->x;
    vs->y = vu->y;
    vs->z = vu->z;
  }
}

void Gb_th_x_v6(const Gb_th* th, const Gb_v6* v, Gb_v6* vs)
{
  Gb_th_x_v3(th, Gb_v6c_get_t(v), Gb_v6_get_t(vs));
  Gb_th_x_v3(th, Gb_v6c_get_r(v), Gb_v6_get_r(vs));
}

void Gb_thInv_x_v6(const Gb_th* th, const Gb_v6* v, Gb_v6* vs)
{
  Gb_thInv_x_v3(th, Gb_v6c_get_t(v), Gb_v6_get_t(vs));
  Gb_thInv_x_v3(th, Gb_v6c_get_r(v), Gb_v6_get_r(vs));
}

void Gb_th_produit(const Gb_th* a, const Gb_th* b, Gb_th* s)
{
  Gb_th temp;
  Gb_th *u;
  if (a == s || b == s) {
    u = &temp;
  } else {
    u = s;
  }
  u->vx.x = a->vx.x * b->vx.x + a->vy.x * b->vx.y + a->vz.x * b->vx.z;
  u->vx.y = a->vx.y * b->vx.x + a->vy.y * b->vx.y + a->vz.y * b->vx.z;
  u->vx.z = a->vx.z * b->vx.x + a->vy.z * b->vx.y + a->vz.z * b->vx.z;
  u->vy.x = a->vx.x * b->vy.x + a->vy.x * b->vy.y + a->vz.x * b->vy.z;
  u->vy.y = a->vx.y * b->vy.x + a->vy.y * b->vy.y + a->vz.y * b->vy.z;
  u->vy.z = a->vx.z * b->vy.x + a->vy.z * b->vy.y + a->vz.z * b->vy.z;
  u->vz.x = a->vx.x * b->vz.x + a->vy.x * b->vz.y + a->vz.x * b->vz.z;
  u->vz.y = a->vx.y * b->vz.x + a->vy.y * b->vz.y + a->vz.y * b->vz.z;
  u->vz.z = a->vx.z * b->vz.x + a->vy.z * b->vz.y + a->vz.z * b->vz.z;
  u->vp.x = a->vx.x * b->vp.x + a->vy.x * b->vp.y + a->vz.x * b->vp.z +a->vp.x;
  u->vp.y = a->vx.y * b->vp.x + a->vy.y * b->vp.y + a->vz.y * b->vp.z +a->vp.y;
  u->vp.z = a->vx.z * b->vp.x + a->vy.z * b->vp.y + a->vz.z * b->vp.z +a->vp.z;
  if (u != s) {
    s->vx.x = u->vx.x;
    s->vx.y = u->vx.y;
    s->vx.z = u->vx.z;
    s->vy.x = u->vy.x;
    s->vy.y = u->vy.y;
    s->vy.z = u->vy.z;
    s->vz.x = u->vz.x;
    s->vz.y = u->vz.y;
    s->vz.z = u->vz.z;
    s->vp.x = u->vp.x;
    s->vp.y = u->vp.y;
    s->vp.z = u->vp.z;
  }
}

/*
 * th représente la matrice de passage de th_12
 * f représente une force exprimée dans le repère 2  f_2
 * fs est la force calculée dans le repère 1  fs_1
 *  ( fs3_1  ) = th_12 * f3_2
 *  ( fsr3_1 ) = th_12 * fr3_2  + th_12->vp ^ fs3_1
 */
void Gb_th_x_force(const Gb_th* th, const Gb_force* f, Gb_force* fs)
{
  Gb_v3* f3 = (Gb_v3*) f;
  Gb_v3* fs3 = (Gb_v3*) fs;
  Gb_v3* fr3 = (Gb_v3*) &f->rx;
  Gb_v3* fsr3 = (Gb_v3*) &fs->rx;
  Gb_th_x_v3(th, f3, fs3);
  Gb_th_x_v3(th, fr3, fsr3);
  fs->rx = fsr3->x  + th->vp.y * fs3->z - th->vp.z * fs3->y;
  fs->ry = fsr3->y  + th->vp.z * fs3->x - th->vp.x * fs3->z;
  fs->rz = fsr3->z  + th->vp.x * fs3->y - th->vp.y * fs3->x;
}

void Gb_th_inverse(const Gb_th* th, Gb_th* ths)
{
  double tmp;
  Gb_v3 v;
  if (th != ths) {
    ths->vx.x = th->vx.x;
    ths->vx.y = th->vy.x;
    ths->vx.z = th->vz.x;
    ths->vy.x = th->vx.y;
    ths->vy.y = th->vy.y;
    ths->vy.z = th->vz.y;
    ths->vz.x = th->vx.z;
    ths->vz.y = th->vy.z;
    ths->vz.z = th->vz.z;
  } else {
    tmp = ths->vx.y;
    ths->vx.y = ths->vy.x;
    ths->vy.x = tmp;
    tmp = ths->vx.z;
    ths->vx.z = ths->vz.x;
    ths->vz.x = tmp;
    tmp = ths->vy.z;
    ths->vy.z = ths->vz.y;
    ths->vz.y = tmp;
  }
  Gb_th_x_v3(ths, &th->vp, &v);
  ths->vp.x = -v.x;
  ths->vp.y = -v.y;
  ths->vp.z = -v.z;
}

void Gb_thInv_x_force(const Gb_th* th, const Gb_force* f, Gb_force* fs)
{
  Gb_th inverse;
  /* a ame'liorer ? */
  Gb_th_inverse(th, &inverse);
  Gb_th_x_force(&inverse, f, fs);
}

/*
 *  ( vs3_1  ) = ( th_12 * vs3_2 + th_12->vp ^ ( th_12 * vsr3_2 ) )
 *  ( vsr3_1 )   ( th_12 * vsr3_2                                 )
 */
void Gb_th_x_vitesse(const Gb_th* th, const Gb_vitesse* v, Gb_vitesse* vs)
{
  Gb_v3* v3 = (Gb_v3*) v;
  Gb_v3* vs3 = (Gb_v3*) vs;
  Gb_v3* vr3 = (Gb_v3*) &v->rx;
  Gb_v3* vsr3 = (Gb_v3*) &vs->rx;
  Gb_th_x_v3(th, v3, vs3);
  Gb_th_x_v3(th, vr3, vsr3);
  vs->x = vs3->x  + th->vp.y * vsr3->z - th->vp.z * vsr3->y;
  vs->y = vs3->y  + th->vp.z * vsr3->x - th->vp.x * vsr3->z;
  vs->z = vs3->z  + th->vp.x * vsr3->y - th->vp.y * vsr3->x;
}

void Gb_quat_x_v3(const Gb_quat* q, const Gb_v3* v, Gb_v3* vs)
{
  /*  ?????*/
  /*
   * soit v_q le quaternion ( 0 v)
   * soit vs_q le quaternion (0 vs)
   * vs_q = q * v_q * q_inv
   *
   */
}
void Gb_quat_x_quat(const Gb_quat* q1, const Gb_quat* q2, Gb_quat* qs)
{
  /*
   * qs.w = q1.w * q2.w - q1.v . q2.v
   * qs.v = q1.w * q2.v + q2.w * q1.v + q1.v /\ q2.v
   */
}

