/*
 * Copyright (c) 2002 LAAS/CNRS -- RIA --
 * Daniel SIDOBRE -- mai 2002
 */
#ifndef gb_h
#define gb_h

/*
 *  Pour Genom il faut séparer la déclaration des structures de celle
 *    des fonctions.
 */
#include "gbStruct.h"

#ifndef GB_GENOM
#include "Proto_gb.h"
#include "Proto_gbModeles.h"
#endif

/* double gb_v3_norme(const Gb_v3* e, Gb_v3* s); */

/* extern void gb_v3_cross_product ( Gb_v3* u, Gb_v3* v, Gb_v3* output ); */
/* extern double gb_v3_product ( Gb_v3* u, Gb_v3* v ); */
/* extern double gb_v3_module ( Gb_v3* u ); */
/* extern void gb_v3_plus ( Gb_v3* u, Gb_v3* v, Gb_v3* output ); */
/* extern void gb_v3_moins ( Gb_v3* u, Gb_v3* v, Gb_v3* output ); */
/* extern void gb_v3_product_r ( Gb_v3* u, double r, Gb_v3* output ); */


/* void gb_v6_plus(Gb_v6* a, Gb_v6* b, Gb_v6* s); */
/* void gb_v6_moins(Gb_v6* a, Gb_v6* b, Gb_v6* s); */
/* Gb_v3* Gb_v6_get_t(Gb_v6* v); */
/* Gb_v3* Gb_v6_get_r(Gb_v6* v); */

/* void Gb_th_vitesse(Gb_th* th, Gb_vitesse* v, Gb_vitesse* vs); */
/* void Gb_th_force(Gb_th* th, Gb_force* f, Gb_force* fs); */
/* void Gb_th_produit(Gb_th* a, Gb_th* b, Gb_th* s); */
/* void Gb_th_v3(Gb_th* th, Gb_v3* v, Gb_v3* vs); */
/* void Gb_th_v3_inverse(Gb_th* th, Gb_v3* v, Gb_v3* vs); */
/* void Gb_th_v6(Gb_th* th, Gb_v6* v, Gb_v6* vs); */
/* void Gb_th_v6_inverse(Gb_th* th, Gb_v6* v, Gb_v6* vs); */
/* void Gb_th_inverse(Gb_th* th, Gb_th* ths); */
/* void Gb_th_force_inverse(Gb_th* th, Gb_force* f, Gb_force* fs); */


/* void gb_dep_set(Gb_dep *dep, double x, double y, double z,  */
/* 		double rx, double ry, double rz, double a); */

/* void gb_dep_get(Gb_dep *dep, double *x, double *y, double *z,  */
/* 		double *rx, double *ry, double *rz, double *a); */



/* void gb_m33_ppv(Gb_m33 *ms, Gb_v3 *ve); */

/* void gb_dep_th(Gb_dep* dep, Gb_th* th); */


/* void gb_dep_quat(Gb_dep* dep, Gb_quat* q); */
/* void gb_quat_dep(Gb_quat* q, Gb_dep* dep); */
/* void gb_th_quat(Gb_th* th, Gb_quat* q); */
/* void gb_quat_th(Gb_quat* q, Gb_th* th); */


/* extern Arm_robot* new_Arm_robot(int nb); */
/* extern void delete_Arm_robot(Arm_robot* r); */

/* extern Arm_liaison* Arm_robot_define_axe(Arm_robot* r, int axe, */
/* 					 Arm_liaison_type type, double a_1, */
/* 					 double alpha_1, double theta_ou_r); */


/* extern int Arm_robot_masse_set(Arm_robot* arm, int i, double masse); */
/* extern int Arm_robot_masse_get(Arm_robot* arm, int i, double* masse); */

#endif
