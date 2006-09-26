/* header file for library 'iter' -- gmres.c, bi-cgstab.c, and orthomin.c.
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: libiter.h,v 2.5 2006/09/26 17:14:32 ichiki Exp $
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

/** from myblas.h **/
/* header file for blas.c --
 * Excerpted from BLAS package
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


/** from gmres.h **/
/* header file of mygmres.c --
 * generalized minimum residual method
 *
 * Reference :
 *   GMRES(m) : Y.Saad & M.H.Schultz, SIAM J.Sci.Stat.Comput.
 *   vol7 (1986) pp.856-869
 */

/** global variables **/
extern int ITER_gmres_debug; /* [0|1]: [not print/print] iter and res */


/* wrapper routine for gmres_m ()
 * INPUT
 *   n : size of vectors v[] and f[] -- expected to be np * nelm for red-sym
 *   b [n] : given vector
 *   atimes (n, x, b) : routine to calc A.x and return b[]
 *   it_max : max # iterations
 *   it_restart : # iterations to restart
 *   eps : the accuracy
 * OUTPUT
 *   x [n] : solution
 */
void
solve_iter_gmres (int n,
		  const double *b, double *x,
		  void (*atimes) (int, const double *, double *, void *),
		  void * user_data,
		  int it_max, int it_restart, double eps);
void
gmres_m (int n, const double *f, double *x,
	 int m, double tol, int itmax,
	 int *iter, double *res,
	 void (*myatimes) (int, const double *, double *, void *),
	 void * user_data);
void
gmres (int n, const double *f, double *x,
       double tol, int itmax,
       int *iter, double *res,
       void (*myatimes) (int, const double *, double *, void *),
       void * user_data);

/** from bi-cgstab.h **/
/* header file of bi-cgstab.c --
 * wrapper for iterative solver routines
 *
 * solver routines are translated into C by K.I. from fortran code
 * originally written by martin h. gutknecht
 *             ===============rog.f==============
 *             problem given by martin h. gutknecht
 *             numerical method: bi-cgsta        method -- sta ()
 *             numerical method: bi-cgsta2       method -- st2 ()
 *             numerical method: gpbi-cg         method -- gpb ()
 *             ver. 1.0 aug. 08 1995  s. l. zhang
 *    at urus.slzhang.fort.iter.real.gpbcg.gutknecht-p(final.f)
 */

/** global variables **/
extern int ITER_stab_debug; /* [0|1]: [not print/print] iter and res */

void
solve_iter_stab (int n, const double *b,
		 double *x,
		 void (*atimes) (int, const double *, double *, void *),
		 void * user_data,
		 void (*solver) (int, const double *, double *, int,
				 double, double, int *, double *,
				 void (*)
				 (int, const double *, double *, void *),
				 void *),
		 int it_max, double log10_eps);
void
sta (int m, const double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, const double *, double *, void *),
     void * user_data);
void
sta2 (int m, const double *b, double *x, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, const double *, double *, void *),
      void * user_data);
void
st2_chk (int m, const double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, const double *, double *, void *),
	 void * user_data);
void
gpb (int m, const double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, const double *, double *, void *),
     void * user_data);
void
gpb_chk (int m, const double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, const double *, double *, void *),
	 void * user_data);

/** from orthomin.h **/
/* header file of orthomin.c --
 * orthomin scheme
 *
 * solver routines are translated into C by K.I. from fortran code
 * originally written by martin h. gutknecht
 *             ===============rog.f==============
 *             problem given by martin h. gutknecht
 *             numerical method: orthomin(k) method
 *             ver. 1.0 jul. 28 1995 by s. l. zhang
 *             ver. 1.1 aug. 31 1995 by s. l. zhang
 *    at slzhang.fort.iter.complex.orthomin.gutknecht-problem(gut.f)
 */

/** global variables **/
extern int ITER_otmk_debug; /* [0|1]: [not print/print] iter and res */

void
solve_iter_otmk (int n, const double *b,
		 double *x,
		 void (*atimes) (int, const double *, double *, void *),
		 void * user_data,
		 void (*solver) (int, const double *, double *,
				 int, int, double,
				 double, int *, double *,
				 void (*)
				 (int, const double *, double *, void *),
				 void *),
		 int it_max, double log10_eps,
		 int it_restart);
void
otmk (int m, const double *b, double *x,
      int kres, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, const double *, double *, void *),
      void * user_data);
