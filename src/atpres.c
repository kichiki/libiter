/* ATPRES -- Weiss, Algorithm 5
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: atpres.c,v 1.6 2007/11/25 18:47:36 kichiki Exp $
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
#include "../config.h"

#include <stdlib.h>
#include <math.h>
#include "libiter.h" // struct iter
#include "myblas.h"
#include "memory-check.h" // CHECK_MALLOC


#ifdef HAVE_CBLAS_H
/* use ATLAS' CBLAS routines */

#include <cblas.h>

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
/* use Fortran BLAS routines */

double
ddot_(int* N, 
      double* X, int* incX, 
      double* Y, int* incY);
double
dnrm2_(int* N, 
       double* X, int* incX);
int
dcopy_(int* N,
       double* X, int* incX,
       double* Y, int* incY);
int
daxpy_(int* N,
       double* alpha,
       double* X, int* incX,
       double* Y, int* incY);
int
dscal_(int* N,
       double* alpha,
       double* X, int* incX);

# else // !HAVE_BLAS_H
/* use local BLAS routines */

#include "myblas.h"

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


/*  for details to the algorithm see
 *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
 *  Algorithm 5 : ATPRES
 * INPUT
 *   n : dimension of the problem
 *   b[n] : r-h-s vector
 *   atimes (int n, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_t (int n, static double *x, double *b, void *param) :
 *        calc matrix-vector product A^T.x = b.
 *   atimes_param : parameters for atimes() and atimes_t().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps  : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x[n] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
atpres (int n, const double *b, double *x,
	void (*atimes) (int, const double *, double *, void *),
	void (*atimes_t) (int, const double *, double *, void *),
	void *atimes_param,
	struct iter *it)
{
#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  double eps = it->eps;
  int maxiter = it->max;

  double *alpha = (double *)malloc (sizeof (double) * 2);
  double *r     = (double *)malloc (sizeof (double) * n);
  double *xx    = (double *)malloc (sizeof (double) * 2 * n);
  double *rr    = (double *)malloc (sizeof (double) * 2 * n);
  double *tmp   = (double *)malloc (sizeof (double) * n);
  double *tmp0  = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (alpha, "atpres");
  CHECK_MALLOC (r,     "atpres");
  CHECK_MALLOC (xx,    "atpres");
  CHECK_MALLOC (rr,    "atpres");
  CHECK_MALLOC (tmp,   "atpres");
  CHECK_MALLOC (tmp0,  "atpres");

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  //double b_norm = cblas_dnrm2 (n, b, 1);
  double b_norm = sqrt (cblas_ddot (n, b, 1, b, 1));

  /* r = A * x(i) - b */
  double tiny = 1.0e-15;
  //if (fabs (cblas_dnrm2 (n, x, 1) / (double) n) > tiny)
  double x_norm = sqrt (cblas_ddot (n, x, 1, x, 1));
  if (fabs (x_norm / (double) n) > tiny)
    {
      atimes (n, x, tmp, atimes_param);
      //my_daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
      cblas_daxpy (n, -1.0, b, 1, tmp, 1);
      cblas_dcopy (n, tmp, 1, r, 1);
    }
  else
    {
      //my_dscalz (n, - 1.0, b, 1, r, 1);
      cblas_dcopy (n, b, 1, r, 1);
      cblas_dscal (n, -1.0, r, 1);
    }

  /* initial condition */
  //double res = cblas_dnrm2 (n, r, 1);
  double res = sqrt (cblas_ddot (n, r, 1, r, 1));
  cblas_dcopy (n, r, 1, rr, 1);
  cblas_dcopy (n, x, 1, xx, 1);

  int iter = 0;
  while (res > eps * b_norm
	 && iter < maxiter)
    {
      if (it->debug == 2)
	{
	  fprintf (it->out, "ATPRES %d %e\n", iter, res);
	}
      iter ++;

      /* set min (k, 2) */
      int mink2;
      if (iter < 2) mink2 = iter;
      else          mink2 = 2;

      /* set alpha */
      int i;
      for (i=0; i<mink2; i++)
	{
	  atimes_t (n, & rr [i * n], tmp, atimes_param);
	  atimes_t (n, & rr [0], tmp0, atimes_param);
	  alpha [i] =
	    - cblas_ddot (n, tmp, 1, tmp0, 1)
	    / cblas_ddot (n, & rr [i * n], 1, & rr [i * n], 1);
	}

      /* set phi */
      double sum = 0.0;
      for (i=0; i<mink2; i++) sum += alpha [i];
      double phi = 1.0 / sum;

      /* update xx */
      cblas_dcopy (n, tmp0, 1, tmp, 1);
      for (i=0; i<mink2; i++)
	{
	  cblas_daxpy (n, alpha [i], & xx [i * n], 1, tmp, 1);
	}
      cblas_dcopy (n, & xx [0], 1, & xx [n], 1);
      //my_dscalz (n, phi, tmp, 1, xx, 1);
      cblas_dcopy (n, tmp, 1, xx, 1);
      cblas_dscal (n, phi, xx, 1);

      /* update rr */
      atimes (n, tmp0, tmp, atimes_param);
      for (i=0; i<mink2; i++)
	{
	  cblas_daxpy (n, alpha [i], & rr [i * n], 1, tmp, 1);
	}
      cblas_dcopy (n, & rr [0], 1, & rr [n], 1);
      //my_dscalz (n, phi, tmp, 1, rr, 1);
      cblas_dcopy (n, tmp, 1, rr, 1);
      cblas_dscal (n, phi, rr, 1);

      /* set gamma */
      //my_daxpyz (n, -1.0, r, 1, rr, 1, tmp, 1);
      cblas_dcopy (n, rr, 1, tmp, 1);
      cblas_daxpy (n, -1.0, r, 1, tmp, 1);
      double gamma
	= - cblas_ddot (n, r, 1, tmp, 1) / cblas_ddot (n, tmp, 1, tmp, 1);

      /* update r */
      cblas_daxpy (n, gamma, tmp, 1, r, 1);
      /* update x */
      //my_daxpyz (n, -1.0, x, 1, xx, 1, tmp, 1);
      cblas_dcopy (n, xx, 1, tmp, 1);
      cblas_daxpy (n, -1.0, x, 1, tmp, 1);

      cblas_daxpy (n, gamma, tmp, 1, x, 1);

      //res = cblas_dnrm2 (n, r, 1);
      res = sqrt (cblas_ddot (n, r, 1, r, 1));
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b_norm = dnrm2_ (&n, b, &i_1);

  /* r = A * x(i) - b */
  double tiny = 1.0e-15;
  if (fabs (dnrm2_ (&n, x, &i_1) / (double) n) > tiny)
    {
      atimes (n, x, tmp, atimes_param);
      //my_daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
      daxpy_ (&n, &d_m1, b, &i_1, tmp, &i_1);
      dcopy_ (&n, tmp, &i_i, r, &i_1);
    }
  else
    {
      //my_dscalz (n, - 1.0, b, 1, r, 1);
      dcopy_ (&n, b, &i_1, r, &i_1);
      dscal_ (&n, &d_m1, r, &i_i);
    }

  /* initial condition */
  double res = dnrm2_ (&n, r, &i_1);
  dcopy_ (&n, r, &i_1, rr, &i_1);
  dcopy_ (&n, x, &i_1, xx, &i_1);

  int iter = 0;
  while (res > eps * b_norm
	 && iter < maxiter)
    {
      if (it->debug == 2)
	{
	  fprintf (it->out, "ATPRES %d %e\n", iter, res);
	}
      iter ++;

      /* set min (k, 2) */
      int mink2;
      if (iter < 2) mink2 = iter;
      else          mink2 = 2;

      /* set alpha */
      int i;
      for (i=0; i<mink2; i++)
	{
	  atimes_t (n, & rr [i * n], tmp, atimes_param);
	  atimes_t (n, & rr [0], tmp0, atimes_param);
	  alpha [i] =
	    - ddot_ (&n, tmp, &i_1, tmp0, &i_1)
	    / ddot_ (&n, &rr[i * n], &i_1, &rr[i * n], &i_1);
	}

      /* set phi */
      double sum = 0.0;
      for (i=0; i<mink2; i++) sum += alpha [i];
      double phi = 1.0 / sum;

      /* update xx */
      dcopy_ (&n, tmp0, &i_1, tmp, &i_1);
      for (i=0; i<mink2; i++)
	{
	  daxpy_ (&n, &alpha[i], &xx[i * n], &i_1, tmp, &i_1);
	}
      dcopy_ (&n, &xx[0], &i_1, &xx[n], &i_1);
      //my_dscalz (n, phi, tmp, 1, xx, 1);
      dcopy_ (&n, tmp, &i_1, xx, &i_1);
      dscal_ (&n, &phi, xx, &i_1);

      /* update rr */
      atimes (n, tmp0, tmp, atimes_param);
      for (i=0; i<mink2; i++)
	{
	  daxpy_ (&n, &alpha[i], &rr[i * n], &i_1, tmp, &i_1);
	}
      dcopy_ (&n, &rr[0], &i_1, &rr[n], &i_1);
      //my_dscalz (n, phi, tmp, 1, rr, 1);
      dcopy_ (&n, tmp, &i_1, rr, &i_1);
      dscal_ (&n, &phi, rr, &i_1);

      /* set gamma */
      //my_daxpyz (n, -1.0, r, 1, rr, 1, tmp, 1);
      dcopy_ (&n, rr, &i_1, tmp, &i_1);
      daxpy_ (&n, &d_m1, r, &i_1, tmp, &i_1);
      double gamma
	= - ddot_ (&n, r, &i_1, tmp, &i_1)
	/ ddot_ (&n, tmp, &i_1, tmp, &i_1);

      /* update r */
      daxpy_ (&n, &gamma, tmp, &i_1, r, &i_1);
      /* update x */
      //my_daxpyz (n, -1.0, x, 1, xx, 1, tmp, 1);
      dcopy_ (&n, xx, &i_1, tmp, &i_1);
      daxpy_ (&n, &d_m1, x, &i_1, tmp, &i_1);

      daxpy_ (&n, &gamma, tmp, &i_1, x, &i_1);

      res = dnrm2_ (&n, r, &i_1);
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b_norm = my_dnrm2 (n, b, 1);

  /* r = A * x(i) - b */
  double tiny = 1.0e-15;
  if (fabs (my_dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp, atimes_param);
      my_daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    {
      my_dscalz (n, - 1.0, b, 1, r, 1);
    }

  /* initial condition */
  double res = my_dnrm2 (n, r, 1);
  my_dcopy (n, r, 1, rr, 1);
  my_dcopy (n, x, 1, xx, 1);

  int iter = 0;
  while (res > eps * b_norm
	 && iter < maxiter)
    {
      if (it->debug == 2)
	{
	  fprintf (it->out, "ATPRES %d %e\n", iter, res);
	}
      iter ++;

      /* set min (k, 2) */
      int mink2;
      if (iter < 2) mink2 = iter;
      else          mink2 = 2;

      /* set alpha */
      int i;
      for (i=0; i<mink2; i++)
	{
	  atimes_t (n, & rr [i * n], tmp, atimes_param);
	  atimes_t (n, & rr [0], tmp0, atimes_param);
	  alpha [i] =
	    - my_ddot (n, tmp, 1, tmp0, 1)
	    / my_ddot (n, & rr [i * n], 1, & rr [i * n], 1);
	}

      /* set phi */
      double sum = 0.0;
      for (i=0; i<mink2; i++) sum += alpha [i];
      double phi = 1.0 / sum;

      /* update xx */
      my_dcopy (n, tmp0, 1, tmp, 1);
      for (i=0; i<mink2; i++)
	{
	  my_daxpy (n, alpha [i], & xx [i * n], 1, tmp, 1);
	}
      my_dcopy (n, & xx [0], 1, & xx [n], 1);
      my_dscalz (n, phi, tmp, 1, xx, 1);

      /* update rr */
      atimes (n, tmp0, tmp, atimes_param);
      for (i=0; i<mink2; i++)
	{
	  my_daxpy (n, alpha [i], & rr [i * n], 1, tmp, 1);
	}
      my_dcopy (n, & rr [0], 1, & rr [n], 1);
      my_dscalz (n, phi, tmp, 1, rr, 1);

      /* set gamma */
      my_daxpyz (n, -1.0, r, 1, rr, 1, tmp, 1);
      double gamma
	= - my_ddot (n, r, 1, tmp, 1) / my_ddot (n, tmp, 1, tmp, 1);

      /* update r */
      my_daxpy (n, gamma, tmp, 1, r, 1);
      /* update x */
      my_daxpyz (n, -1.0, x, 1, xx, 1, tmp, 1);
      my_daxpy (n, gamma, tmp, 1, x, 1);

      res = my_dnrm2 (n, r, 1);
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  int ret = -1;
  if (res <= eps * b_norm) ret = 0; // success

  free (alpha);
  free (r);
  free (xx);
  free (rr);
  free (tmp);
  free (tmp0);

  it->niter = iter;
  it->res2  = res * res / (b_norm * b_norm);
  return (ret);
}
