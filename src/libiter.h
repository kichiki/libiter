/* header file for library 'iter' -- mygmres.c, bi-cgstab.c, and orthomin.c.
 * Copyright (C) 2001 Kengo Ichiki <ichiki@haloumi.tn.utwente.nl>
 * $Id: libiter.h,v 2.3 2001/10/19 14:34:59 ichiki Exp $
 */

/* from blas.h */
void
daxpy (int n, double da, double *dx, int incx, double *dy, int incy);
void
daxpyz (int n, double da, double *dx, int incx, double *dy, int incy,
	double *dz, int incz);
void
dcopy (int n, double *dx, int incx,
       double * dy, int incy);
double
dnrm2_ (int n, double *x, int incx);
double
dnrm2 (int n, double *x, int incx);
double
ddot (int n, double *dx, int incx, double *dy, int incy);
void
dscal (int n, double da, double *dx, int incx);
void
dscalz (int n, double da, double *dx, int incx,
	double *dz, int incz);


/** from mygmres.h **/
/* header file of mygmres.c
 * Copyright (C) 1998-2001 Kengo Ichiki <ichiki@haloumi.tn.utwente.nl>
 * Id: mygmres.h,v 2.3 2001/10/19 14:29:41 ichiki Exp
 */

/** global variables **/
extern int ITER_mygmres_debug; /* [0|1]: [not print/print] iter and res */


void
solve_iter_gmres (int n,
		  double *b, double *x,
		  void (*atimes) (int, double *, double *, void *),
		  void * user_data,
		  int it_max, int it_restart, double eps);
void
mygmres_m (int n, double *f, double *x,
	   int m, double tol, int itmax,
	   int *iter, double *res,
	   void (*myatimes) (int, double *, double *, void *),
	   void * user_data);
void
mygmres (int n, double *f, double *x,
	 double tol, int itmax,
	 int *iter, double *res,
	 void (*myatimes) (int, double *, double *, void *),
	 void * user_data);


/** from bi-cgstab.h **/
/* header file of bi-cgstab.c (wrapper for iterative solver routines)
 * Copyright (C) 1999-2001 Kengo Ichiki <ichiki@haloumi.tn.utwente.nl>
 * Id: bi-cgstab.h,v 2.2 2001/10/19 14:31:22 ichiki Exp
 */

/** global variables **/
extern int ITER_stab_debug; /* [0|1]: [not print/print] iter and res */

void
solve_iter_stab (int n,
		 double *b, double *x,
		 void (*atimes) (int, double *, double *, void *),
		 void * user_data,
		 void (*solver) (int, double *, double *, int,
				 double, double, int *, double *,
				 void (*)
				 (int, double *, double *, void *),
				 void *),
		 int it_max, double log10_eps);
void
sta (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *, void *),
     void * user_data);
void
sta2 (int m, double *b, double *x, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, double *, double *, void *),
      void * user_data);
void
st2_chk (int m, double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, double *, double *, void *),
	 void * user_data);
void
gpb (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *, void *),
     void * user_data);
void
gpb_chk (int m, double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, double *, double *, void *),
	 void * user_data);


/** from orthomin.h **/
/* header file of orthomin.c
 * Copyright (C) 1999-2001 Kengo Ichiki <ichiki@haloumi.tn.utwente.nl>
 * Id: orthomin.h,v 2.2 2001/10/19 14:32:55 ichiki Exp
 */

/** global variables **/
extern int ITER_otmk_debug; /* [0|1]: [not print/print] iter and res */

void
solve_iter_otmk (int n,
		 double *b, double *x,
		 void (*atimes) (int, double *, double *, void *),
		 void * user_data,
		 void (*solver) (int, double *, double *,
				 int, int, double,
				 double, int *, double *,
				 void (*)
				 (int, double *, double *, void *),
				 void *),
		 int it_max, double log10_eps,
		 int it_restart);
void
otmk (int m, double *b, double *x,
      int kres, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, double *, double *, void *),
      void * user_data);
