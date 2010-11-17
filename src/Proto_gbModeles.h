/*
 *   This file was automatically generated by version 1.7 of cextract.
 *   Manual editing not recommended.
 *
 *   Created: Wed Jun 29 16:41:44 2005
 */
#ifndef __CEXTRACT__
#if __STDC__

extern char* Gb_statusMGI_s ( Gb_statusMGI u );
extern double Gb_atan2 ( double y, double x );
extern void Gb_dataMGD_print ( FILE* output, Gb_dataMGD *d );
extern void Gb_MGD6rTh ( Gb_6rParameters* bras, Gb_q6* eq, Gb_dataMGD* d, Gb_th* th );
extern void Gb_MGD6r_gete1e2e3(Gb_6rParameters* bras, Gb_q6* eq, int* e1, int* e2, int* e3);
extern void Gb_MGD6r_6Th ( Gb_6rParameters* bras, Gb_q6* eq, Gb_dataMGD* d, Gb_th* th01, Gb_th* th02, Gb_th* th03, Gb_th* th04, Gb_th* th05, Gb_th* th06 );
extern Gb_statusMGI Gb_MGI6rTh ( Gb_6rParameters* bras, Gb_th* eth, int e1, int e2, int e3, Gb_q6* old_q, Gb_dataMGD* d, Gb_q6* sq );
Gb_statusMGI Gb_MGI6rTh_O(Gb_6rParameters* bras, Gb_th* eth, Gb_q6* old_q, Gb_dataMGD* d, Gb_q6* sq);
extern void Gb_MDD6r ( Gb_6rParameters* bras, Gb_dataMGD* d, Gb_th* t06, Gb_jac* jac );
extern void kukaLBR_mgd(Gb_q7* Q, double r3, double r5, Gb_th* th07);
extern Gb_statusMGI kukaLBR_mgi_q_e(Gb_th* th07, Gb_q7* Qp, double r3, double r5, double epsilon, int e1, int e2, int e3, Gb_q7* q) ;
extern void kukaLBR_gete1e2e3(double r3, double r5, Gb_q7* eq, int* e1, int* e2, int* e3);
extern Gb_statusMGI pr2_mgi_q3_8(Gb_th* th07, Gb_q7* Qp,
                                 double a1, double r3, double r5, Gb_q7 *qMin, Gb_q7 *qMax, 
                                 double epsilon, Gb_q7 qsol[32], int* nbsolution);

#else /* __STDC__ */

extern char* Gb_statusMGI_s (/* Gb_statusMGI u */);
extern double Gb_atan2 (/* double y, double x */);
extern void Gb_dataMGD_print (/* FILE* output, Gb_dataMGD *d */);
extern void Gb_MGD6rTh (/* Gb_6rParameters* bras, Gb_q6* eq, Gb_dataMGD* d, Gb_th* th */);
extern void Gb_MGD6r_gete1e2e3(/* Gb_6rParameters* bras, Gb_q6* eq, int* e1, int* e2, int* e3 */);
extern void Gb_MGD6r_6Th (/* Gb_6rParameters* bras, Gb_q6* eq, Gb_dataMGD* d, Gb_th* th01, Gb_th* th02, Gb_th* th03, Gb_th* th04, Gb_th* th05, Gb_th* th06 */);
extern Gb_statusMGI Gb_MGI6rTh (/* Gb_6rParameters* bras, Gb_th* eth, int e1, int e2, int e3, Gb_q6* old_q, Gb_dataMGD* d, Gb_q6* sq */);
extern void Gb_MDD6r (/* Gb_6rParameters* bras, Gb_dataMGD* d, Gb_th* t06, Gb_jac* jac */);
extern void kukaLBR_mgd(Gb_q7* Q, double r3, double r5, Gb_th* th07);
extern Gb_statusMGI kukaLBR_mgi_q_e(Gb_th* th07, Gb_q7* Qp, double r3, double r5, double epsilon, int e1, int e2, int e3, Gb_q7* q) ;
extern void kukaLBR_gete1e2e3(double r3, double r5, Gb_q7* eq, int* e1, int* e2, int* e3);
extern Gb_statusMGI pr2_mgi_q3_8(Gb_th* th07, Gb_q7* Qp,
                                 double a1, double r3, double r5, Gb_q7 *qMin, Gb_q7 *qMax, 
                                 double epsilon, Gb_q7 qsol[32], int* nbsolution) ;

#endif /* __STDC__ */
#endif /* __CEXTRACT__ */
