/*
 *   This file was automatically generated by version 1.7 of cextract.
 *   Manual editing not recommended.
 *
 *   Created: Wed Apr 27 10:57:03 2005
 */
#ifndef __CEXTRACT__
#if __STDC__

extern double Gb_v3_norme ( const Gb_v3* e, Gb_v3* s );
extern double Gb_prs ( const double w1[3], const double w2[3] );
extern double Gb_v3_prs ( const Gb_v3* v1, const Gb_v3* v2 );
extern void Gb_v3_cross_product ( const Gb_v3* u, const Gb_v3* v, Gb_v3* output );
extern double Gb_v3_product ( const Gb_v3* u, const Gb_v3* v );
extern double Gb_v3_module ( const Gb_v3* u );
extern void Gb_v3_plus ( const Gb_v3* u, const Gb_v3* v, Gb_v3* output );
extern void Gb_v3_moins ( const Gb_v3* u, const Gb_v3* v, Gb_v3* output );
extern void Gb_v3_product_r ( const Gb_v3* u, double r, Gb_v3* output );
extern void Gb_v3_div_r(const Gb_v3* u, double r, Gb_v3* output);
extern void Gb_v3_dist(const Gb_v3* u, const Gb_v3* v, double* dist);
extern double Gb_v3_dist_droite(const Gb_v3* u, const Gb_v3* v, const Gb_v3* a);
extern void Gb_v3_set( Gb_v3* u, double x, double y, double z);
extern void Gb_v3_copy_into(const Gb_v3* e, Gb_v3* s);
extern double Gb_v3_angle(const Gb_v3* u, const Gb_v3* v);
extern void Gb_v3_print(const Gb_v3* u);
extern void Gb_dep_set ( Gb_dep* dep, double x, double y, double z, double rx, double ry, double rz, double a );
extern void Gb_dep_get ( Gb_dep* dep, double* x, double* y, double* z, double* rx, double* ry, double* rz, double* a );
extern void Gb_dep_print(const Gb_dep* dep);
extern void Gb_m33_ppv ( const Gb_v3* ve, Gb_m33* ms );
extern void Gb_dep_quat ( const Gb_dep* dep, Gb_quat* q );
extern void Gb_quat_dep ( const Gb_quat* q, Gb_dep* dep );
extern void Gb_quat_th ( const Gb_quat* q, Gb_th* th );
extern void Gb_dep_th ( const Gb_dep* dep, Gb_th* th );
extern void Gb_th_dep ( const Gb_th* th, Gb_dep* dep );
extern void Gb_th_quat ( const Gb_th* th, Gb_quat* q );
extern void Gb_v6_plus ( const Gb_v6* a, const Gb_v6* b, Gb_v6* s );
extern void Gb_v6_moins ( const Gb_v6* a, const Gb_v6* b, Gb_v6* s );
extern double Gb_v6_module ( const Gb_v6* a );
extern Gb_v3* Gb_v6_get_t ( Gb_v6* v );
extern Gb_v3* Gb_v6c_get_t ( const Gb_v6* v );
extern Gb_v3* Gb_v6_get_r ( Gb_v6* v );
extern Gb_v3* Gb_v6c_get_r ( const Gb_v6* v );
extern void Gb_th_x_v3 ( const Gb_th* th, const Gb_v3* v, Gb_v3* vs );
extern void Gb_thInv_x_v3 ( const Gb_th* th, const Gb_v3* v, Gb_v3* vs );
extern void Gb_th_x_v6 ( const Gb_th* th, const Gb_v6* v, Gb_v6* vs );
extern void Gb_thInv_x_v6 ( const Gb_th* th, const Gb_v6* v, Gb_v6* vs );
extern void Gb_th_produit ( const Gb_th* a, const Gb_th* b, Gb_th* s );
extern void Gb_th_x_force ( const Gb_th* th, const Gb_force* f, Gb_force* fs );
extern void Gb_th_inverse ( const Gb_th* th, Gb_th* ths );
extern void Gb_thInv_x_force ( const Gb_th* th, const Gb_force* f, Gb_force* fs );
extern void Gb_th_x_vitesse ( const Gb_th* th, const Gb_vitesse* v, Gb_vitesse* vs );
extern void Gb_th_print(const Gb_th* th, const char* label);
extern void Gb_quat_x_v3 ( const Gb_quat* q, const Gb_v3* u, Gb_v3* vs );
extern void Gb_quat_x_quat ( const Gb_quat* q1, const Gb_quat* q2, Gb_quat* qs );
extern void Gb_quat_inverse ( const Gb_quat* qi, Gb_quat* qo );
extern void Gb_quat_conjugue ( const Gb_quat* qi, Gb_quat* qo );
extern void Gb_quat_interpole ( const Gb_quat* q1, const Gb_quat* q2, double s, Gb_quat* qo );
extern int Gb_quat_interpole_dep ( const Gb_dep* d1, const Gb_dep* d2, double s, Gb_dep* d_o );
extern int Gb_quat_interpole_dep2 ( const Gb_quat* q1, const Gb_quat* q2, double s, Gb_quat* qo );
extern void Gb_quat_interpole_diff ( const Gb_quat* q1, const Gb_quat* q2, double s, Gb_quat* qq, Gb_quat* qd, Gb_quat* qdiff );
extern void Gb_q6_set(Gb_q6* s, double q1, double q2, double q3, double q4, double q5, double q6);
extern void Gb_q6_get(const Gb_q6* e, double* q1, double* q2, double* q3, double* q4, double* q5, double* q6);
extern void Gb_q6_print(const Gb_q6* e);
extern void Gb_quat_compute_relativeDep_to_interpole(const Gb_quat* q1, const Gb_quat* q2, Gb_dep* relDep);
extern void Gb_quat_interpole_depRel(const Gb_quat* q1, const Gb_quat* q2, double s, Gb_quat* qo, Gb_dep *relDep);

#else /* __STDC__ */

extern double Gb_v3_norme (/* const Gb_v3* e, Gb_v3* s */);
extern double Gb_prs (/* const double w1[3], const double w2[3] */);
extern double Gb_v3_prs (/* const Gb_v3* v1, const Gb_v3* v2 */);
extern void Gb_v3_cross_product (/* const Gb_v3* u, const Gb_v3* v, Gb_v3* output */);
extern double Gb_v3_product (/* const Gb_v3* u, const Gb_v3* v */);
extern double Gb_v3_module (/* const Gb_v3* u */);
extern void Gb_v3_plus (/* const Gb_v3* u, const Gb_v3* v, Gb_v3* output */);
extern void Gb_v3_moins (/* const Gb_v3* u, const Gb_v3* v, Gb_v3* output */);
extern void Gb_v3_product_r (/* const Gb_v3* u, double r, Gb_v3* output */);
extern void Gb_v3_dist(/*const Gb_v3* u, const Gb_v3* v, double* dist*/);
extern void Gb_v3_print(/*const Gb_v3* u*/);
extern void Gb_dep_set (/* Gb_dep* dep, double x, double y, double z, double rx, double ry, double rz, double a */);
extern void Gb_dep_get (/* Gb_dep* dep, double* x, double* y, double* z, double* rx, double* ry, double* rz, double* a */);
extern void Gb_dep_print(/*const Gb_dep* dep*/);
extern void Gb_m33_ppv (/* const Gb_v3* ve, Gb_m33* ms */);
extern void Gb_dep_quat (/* const Gb_dep* dep, Gb_quat* q */);
extern void Gb_quat_dep (/* const Gb_quat* q, Gb_dep* dep */);
extern void Gb_quat_th (/* const Gb_quat* q, Gb_th* th */);
extern void Gb_dep_th (/* const Gb_dep* dep, Gb_th* th */);
extern void Gb_th_dep (/* const Gb_th* th, Gb_dep* dep */);
extern void Gb_th_quat (/* const Gb_th* th, Gb_quat* q */);
extern void Gb_v6_plus (/* const Gb_v6* a, const Gb_v6* b, Gb_v6* s */);
extern void Gb_v6_moins (/* const Gb_v6* a, const Gb_v6* b, Gb_v6* s */);
extern double Gb_v6_module (/* const Gb_v6* a */);
extern Gb_v3* Gb_v6_get_t (/* Gb_v6* v */);
extern Gb_v3* Gb_v6c_get_t (/* const Gb_v6* v */);
extern Gb_v3* Gb_v6_get_r (/* Gb_v6* v */);
extern Gb_v3* Gb_v6c_get_r (/* const Gb_v6* v */);
extern void Gb_th_x_v3 (/* const Gb_th* th, const Gb_v3* v, Gb_v3* vs */);
extern void Gb_thInv_x_v3 (/* const Gb_th* th, const Gb_v3* v, Gb_v3* vs */);
extern void Gb_th_x_v6 (/* const Gb_th* th, const Gb_v6* v, Gb_v6* vs */);
extern void Gb_thInv_x_v6 (/* const Gb_th* th, const Gb_v6* v, Gb_v6* vs */);
extern void Gb_th_produit (/* const Gb_th* a, const Gb_th* b, Gb_th* s */);
extern void Gb_th_x_force (/* const Gb_th* th, const Gb_force* f, Gb_force* fs */);
extern void Gb_th_inverse (/* const Gb_th* th, Gb_th* ths */);
extern void Gb_thInv_x_force (/* const Gb_th* th, const Gb_force* f, Gb_force* fs */);
extern void Gb_th_x_vitesse (/* const Gb_th* th, const Gb_vitesse* v, Gb_vitesse* vs */);
extern void Gb_quat_x_v3 (/* const Gb_quat* q, const Gb_v3* u, Gb_v3* vs */);
extern void Gb_quat_x_quat (/* const Gb_quat* q1, const Gb_quat* q2, Gb_quat* qs */);
extern void Gb_quat_inverse (/* const Gb_quat* qi, Gb_quat* qo */);
extern void Gb_quat_conjugue (/* const Gb_quat* qi, Gb_quat* qo */);
extern void Gb_quat_interpole (/* const Gb_quat* q1, const Gb_quat* q2, double s, Gb_quat* qo */);
extern int Gb_quat_interpole_dep (/* const Gb_dep* d1, const Gb_dep* d2, double s, Gb_dep* d_o */);
extern int Gb_quat_interpole_dep2 (/* const Gb_quat* q1, const Gb_quat* q2, double s, Gb_quat* qo */);
extern void Gb_quat_interpole_diff (/* const Gb_quat* q1, const Gb_quat* q2, double s, Gb_quat* qq, Gb_quat* qd, Gb_quat* qdiff */);
extern void Gb_q6_set(/*Gb_q6* s, double q1, double q2, double q3, double q4, double q5, double q6*/);
extern void Gb_q6_get(/*const Gb_q6* e, double* q1, double* q2, double* q3, double* q4, double* q5, double* q6*/);
extern void Gb_q6_print(/*const Gb_q6* e*/);
extern void Gb_quat_compute_relativeDep_to_interpole(/*const Gb_quat* q1, const Gb_quat* q2, Gb_dep* relDep*/);
extern void Gb_quat_interpole_depRel();
#endif /* __STDC__ */
#endif /* __CEXTRACT__ */
