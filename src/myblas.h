/* BLAS
 * Copyright (C) 1999,2003 Kengo Ichiki <ichiki@jhu.edu>
 * $Id: myblas.h,v 1.1 2003/04/26 02:49:28 ichiki Exp $
 *
 * Excerpted from : BLAS package
 */

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
