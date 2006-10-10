/* CGS -- Weiss, Algorithm 11
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: cgs.c,v 2.2 2006/10/10 18:08:06 ichiki Exp $
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

#include <stdio.h>
#include <stdlib.h>
#include "libiter.h"


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


/* Ref: Weiss, Algorithm 11 CGS
 */
void
cgs (int n, const double *b, double *x,
     void (*atimes) (int, const double *, double *, void *),
     void * user_data,
     struct iter * it_param)
{
  double tol, tol2;
  int itmax;

  int i;

  double * r;
  double * r0;
  double * p;
  double * u;
  double * ap; // A.p
  double * q;
  double * t;

  //double * bqu;
  double *qu;

  double r0ap; // (r0*, A.p)

  double rho, rho1;
  double delta;
  double beta;

  double res2 = 0.0;

#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;
  double d_2 = 2.0;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H



  tol = it_param->eps;
  tol2 = tol * tol;
  itmax = it_param->max;

  r  = (double *) malloc (sizeof (double) * n);
  r0 = (double *) malloc (sizeof (double) * n);
  p  = (double *) malloc (sizeof (double) * n);
  u  = (double *) malloc (sizeof (double) * n);
  ap = (double *) malloc (sizeof (double) * n);
  q  = (double *) malloc (sizeof (double) * n);
  t  = (double *) malloc (sizeof (double) * n);

#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */

  // initial guess
  cblas_dcopy (n, b, 1, x, 1);

  // initial residue
  atimes (n, x, r, user_data); // r = A.x
  cblas_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b

  cblas_dcopy (n, r, 1, r0, 1); // r0* = r
  cblas_dcopy (n, r, 1, p, 1); // p = r
  cblas_dcopy (n, r, 1, u, 1); // u = r

  rho = cblas_ddot (n, r0, 1, r, 1); // rho = (r0*, r)

  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, user_data); // ap = A.p
      r0ap = cblas_ddot (n, r0, 1, ap, 1); // r0ap = (r0*, A.p)
      delta = - rho / r0ap;

      cblas_dcopy (n, u, 1, q, 1); // q = u
      cblas_dscal (n, 2.0, q, 1); // q = 2 u
      cblas_daxpy (n, delta, ap, 1, q, 1); // q = 2 u + delta A.p

      atimes (n, q, t, user_data); // t = A.q

      cblas_daxpy (n, delta, t, 1, r, 1); // r = r + delta t
      cblas_daxpy (n, delta, q, 1, x, 1); // x = x + delta q

      res2 = cblas_ddot (n, r, 1, r, 1);
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-cgs %d %e\n", i, res2);
	}
      if (res2 < tol2) break;

      rho1 = cblas_ddot (n, r0, 1, r, 1); // rho = (r0*, r)
      beta = rho1 / rho;
      rho = rho1;

      qu = t; // here t is not used so that this is used for working area.
      cblas_dcopy (n, q, 1, qu, 1); // qu = q
      cblas_daxpy (n, -1.0, u, 1, qu, 1); // qu = q - u
      cblas_dcopy (n, r, 1, u, 1); // u = r
      cblas_daxpy (n, beta, qu, 1, u, 1); // u = r + beta (q - u)

      cblas_daxpy (n, beta, p, 1, qu, 1); // qu = q - u + beta * p
      cblas_dcopy (n, u, 1, p, 1); // p = u
      cblas_daxpy (n, beta, qu, 1, p, 1); // p = u + beta (q - u + b * p)
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  // initial guess
  dcopy_ (&n, b, &i_1, x, &i_1);

  // initial residue
  atimes (n, x, r, user_data); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b

  dcopy_ (&n, r, &i_1, r0, &i_1); // r0* = r
  dcopy_ (&n, r, &i_1, p, &i_1); // p = r
  dcopy_ (&n, r, &i_1, u, &i_1); // u = r

  rho = ddot_ (&n, r0, &i_1, r, &i_1); // rho = (r0*, r)

  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, user_data); // ap = A.p
      r0ap = ddot_ (&n, r0, &i_1, ap, &i_1); // r0ap = (r0*, A.p)
      delta = - rho / r0ap;

      dcopy_ (&n, u, &i_1, q, &i_1); // q = u
      dscal_ (&n, &d_2, q, &i_1); // q = 2 u
      daxpy_ (&n, &delta, ap, &i_1, q, &i_1); // q = 2 u + delta A.p

      atimes (n, q, t, user_data); // t = A.q

      daxpy_ (&n, &delta, t, &i_1, r, &i_1); // r = r + delta t
      daxpy_ (&n, &delta, q, &i_1, x, &i_1); // x = x + delta q

      res2 = ddot_ (&n, r, &i_1, r, &i_1);
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-cgs %d %e\n", i, res2);
	}
      if (res2 < tol2) break;

      rho1 = ddot_ (&n, r0, &i_1, r, &i_1); // rho = (r0*, r)
      beta = rho1 / rho;
      rho = rho1;

      qu = t; // here t is not used so that this is used for working area.
      dcopy_ (&n, q, &i_1, qu, &i_1); // qu = q
      daxpy_ (&n, &d_m1, u, &i_1, qu, &i_1); // qu = q - u
      dcopy_ (&n, r, &i_1, u, &i_1); // u = r
      daxpy_ (&n, &beta, qu, &i_1, u, &i_1); // u = r + beta (q - u)

      daxpy_ (&n, &beta, p, &i_1, qu, &i_1); // qu = q - u + beta * p
      dcopy_ (&n, u, &i_1, p, &i_1); // p = u
      daxpy_ (&n, &beta, qu, &i_1, p, &i_1); // p = u + beta (q - u + b * p)
    }

# else // !HAVE_BLAS_H
  /* use local BLAS routines */

  // initial guess
  my_dcopy (n, b, 1, x, 1);

  // initial residue
  atimes (n, x, r, user_data); // r = A.x
  my_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b

  my_dcopy (n, r, 1, r0, 1); // r0* = r
  my_dcopy (n, r, 1, p, 1); // p = r
  my_dcopy (n, r, 1, u, 1); // u = r

  rho = my_ddot (n, r0, 1, r, 1); // rho = (r0*, r)

  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, user_data); // ap = A.p
      r0ap = my_ddot (n, r0, 1, ap, 1); // r0ap = (r0*, A.p)
      delta = - rho / r0ap;

      my_dcopy (n, u, 1, q, 1); // q = u
      my_dscal (n, 2.0, q, 1); // q = 2 u
      my_daxpy (n, delta, ap, 1, q, 1); // q = 2 u + delta A.p

      atimes (n, q, t, user_data); // t = A.q

      my_daxpy (n, delta, t, 1, r, 1); // r = r + delta t
      my_daxpy (n, delta, q, 1, x, 1); // x = x + delta q

      res2 = my_ddot (n, r, 1, r, 1);
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-cgs %d %e\n", i, res2);
	}
      if (res2 < tol2) break;

      rho1 = my_ddot (n, r0, 1, r, 1); // rho = (r0*, r)
      beta = rho1 / rho;
      rho = rho1;

      qu = t; // here t is not used so that this is used for working area.
      my_dcopy (n, q, 1, qu, 1); // qu = q
      my_daxpy (n, -1.0, u, 1, qu, 1); // qu = q - u
      my_dcopy (n, r, 1, u, 1); // u = r
      my_daxpy (n, beta, qu, 1, u, 1); // u = r + beta (q - u)

      my_daxpy (n, beta, p, 1, qu, 1); // qu = q - u + beta * p
      my_dcopy (n, u, 1, p, 1); // p = u
      my_daxpy (n, beta, qu, 1, p, 1); // p = u + beta (q - u + b * p)
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (r);
  free (r0);
  free (p);
  free (u);
  free (ap);
  free (q);
  free (t);

  if (it_param->debug == 1)
    {
      fprintf (it_param->out, "libiter-cgs %d %e\n", i, res2);
    }
}
