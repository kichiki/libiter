/* header file for library 'iter' -- mygmres.c, bi-cgstab.c, and orthomin.c.
 * Copyright (C) 2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: libiter.h,v 1.1 2001/01/07 05:50:15 ichiki Exp $
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

/* from mygmres.h */
void
mygmres_m (int n, double *f, double *x,
	   int m, double tol, int itmax,
	   int *iter, double *err,
	   void (*myatimes) (int, double *, double *));
void
mygmres (int n, double *f, double *x,
	 double tol, int itmax,
	 int *iter, double *err,
	 void (*myatimes) (int, double *, double *));
void
back_sub (int m, int nn,
	  double *r, double *g, double *y);

/* from bi-cgstab.h */
void
sta (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *));
void
st2 (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *));
void
gpb (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *));

/* from orthomin.h */
void
otmk (int m, double *b, double *x,
      int kres, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, double *, double *));
