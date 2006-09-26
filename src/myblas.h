/* header file for blas.c --
 * BLAS
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: myblas.h,v 1.2 2006/09/26 05:17:41 ichiki Exp $
 *
 * Excerpted from : BLAS package
 */

/*    constant times a vector plus a vector.
 *    uses unrolled loops for increments equal to one.
 *    jack dongarra, linpack, 3/11/78.
 *    modified 12/3/93, array(1) declarations changed to array(*)
 * calc dy [i] = dy [i] + da * dx [i]
 * INPUT
 *    n : dimension
 *    da : constant
 *    dx [n * incx] : vector
 *    dy [n * incy] : vector
 * OUTPUT
 *    dy [n * incy] : vector
 */
void
daxpy (int n, double da, const double *dx, int incx, double *dy, int incy);

/* calc dz [i] = dy [i] + da * dx [i]
 * INPUT
 *    n : dimension
 *    da : constant
 *    dx [n * incx] : vector
 *    dy [n * incy] : vector
 * OUTPUT
 *    dz [n * incy] : vector
 */
void
daxpyz (int n, double da, const double *dx, int incx,
	const double *dy, int incy,
	double *dz, int incz);

/*     copies a vector, x, to a vector, y.
 *     uses unrolled loops for increments equal to one.
 *     jack dongarra, linpack, 3/11/78.
 *     modified 12/3/93, array(1) declarations changed to array(*)
 * calc dy [i] = dx [i]
 * INPUT
 *    n : dimension
 *    dx [n * incx] : vector
 * OUTPUT
 *    dy [n * incy] : vector
 */
void
dcopy (int n, const double *dx, int incx,
       double * dy, int incy);

/*  DNRM2 returns the euclidean norm of a vector via the function
 *  name, so that
 *
 *     DNRM2 := sqrt( x'*x )
 *
 *  -- This version written on 25-October-1982.
 *     Modified on 14-October-1993 to inline the call to DLASSQ.
 *     Sven Hammarling, Nag Ltd.
 */
double
dnrm2_ (int n, const double *x, int incx);

double
dnrm2 (int n, const double *x, int incx);

/*     forms the dot product of two vectors.
 *     uses unrolled loops for increments equal to one.
 *     jack dongarra, linpack, 3/11/78.
 *     modified 12/3/93, array(1) declarations changed to array(*)
 * calc (dy, dx)
 * INPUT
 *    n : dimension
 *    dx [n * incx] : vector
 *    dy [n * incy] : vector
 * OUTPUT (return value)
 */
double
ddot (int n, const double *dx, int incx, const double *dy, int incy);

/*     scales a vector by a constant.
 *     uses unrolled loops for increment equal to one.
 *     jack dongarra, linpack, 3/11/78.
 *     modified 3/93 to return if incx .le. 0.
 *     modified 12/3/93, array(1) declarations changed to array(*)
 * calc dx [i] = da * dx [i]
 */
void
dscal (int n, double da, double *dx, int incx);
/*     scales a vector by a constant.
 *     uses unrolled loops for increment equal to one.
 *     jack dongarra, linpack, 3/11/78.
 *     modified 3/93 to return if incx .le. 0.
 *     modified 12/3/93, array(1) declarations changed to array(*)
 * calc dz [i] = da * dx [i]
 */
void
dscalz (int n, double da, const double *dx, int incx,
	double *dz, int incz);
