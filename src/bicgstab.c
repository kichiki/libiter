/* BiCGSTAB - Weiss, Algorithm 12
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bicgstab.c,v 2.2 2006/10/10 18:08:59 ichiki Exp $
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


/* Ref: Weiss, Algorithm 12 BiCGSTAB
 */
void
bicgstab (int n, const double *b, double *x,
	  void (*atimes) (int, const double *, double *, void *),
	  void * user_data,
	  struct iter * it_param)
{
  double tol, tol2;
  int itmax;

  int i;

  double * r;
  double * rs;
  double * p;
  double * ap; // A.p
  double * s;
  double * t;

  double rsap; // (r*, A.p)
  double st;
  double t2;

  double rho, rho1;
  double delta;
  double gamma;
  double beta;

  double res2 = 0.0;

#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


  tol = it_param->eps;
  tol2 = tol * tol;
  itmax = it_param->max;

  r  = (double *) malloc (sizeof (double) * n);
  rs = (double *) malloc (sizeof (double) * n);
  p  = (double *) malloc (sizeof (double) * n);
  ap = (double *) malloc (sizeof (double) * n);
  s  = (double *) malloc (sizeof (double) * n);
  t  = (double *) malloc (sizeof (double) * n);

#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */

  // initial guess
  cblas_dcopy (n, b, 1, x, 1);
  
  atimes (n, x, r, user_data); // r = A.x
  cblas_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b

  cblas_dcopy (n, r, 1, rs, 1); // r* = r
  cblas_dcopy (n, r, 1, p, 1); // p = r

  rho = cblas_ddot (n, rs, 1, r, 1); // rho = (r*, r)

  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, user_data); // ap = A.p
      rsap = cblas_ddot (n, rs, 1, ap, 1); // rsap = (r*, A.p)
      delta = - rho / rsap;

      cblas_dcopy (n, r, 1, s, 1); // s = r
      cblas_daxpy (n, delta, ap, 1, s, 1); // s = r + delta A.p
      atimes (n, s, t, user_data); // t = A.s

      st = cblas_ddot (n, s, 1, t, 1); // st = (s, t)
      t2 = cblas_ddot (n, t, 1, t, 1); // t2 = (t, t)
      gamma = - st / t2;

      cblas_dcopy (n, s, 1, r, 1); // r = s
      cblas_daxpy (n, gamma, t, 1, r, 1); // r = s + gamma t

      cblas_daxpy (n, delta, p, 1, x, 1); // x = x + delta p...
      cblas_daxpy (n, gamma, s, 1, x, 1); //       + gamma s

      res2 = cblas_ddot (n, r, 1, r, 1);
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-bicgstab %d %e\n", i, res2);
	}
      if (res2 < tol2) break;

      rho1 = cblas_ddot (n, rs, 1, r, 1); // rho = (r*, r)
      beta = rho1 / rho * delta / gamma;
      rho = rho1;

      cblas_daxpy (n, gamma, ap, 1, p, 1); // p = p + gamma A.p
      cblas_dscal (n, beta, p, 1); // p = beta (p + gamma A.p)
      cblas_daxpy (n, 1.0, r, 1, p, 1); // p = r + beta(p + gamma A.p)
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  // initial guess
  dcopy_ (&n, b, &i_1, x, &i_1);
  
  atimes (n, x, r, user_data); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b

  dcopy_ (&n, r, &i_1, rs, &i_1); // r* = r
  dcopy_ (&n, r, &i_1, p, &i_1); // p = r

  rho = ddot_ (&n, rs, &i_1, r, &i_1); // rho = (r*, r)

  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, user_data); // ap = A.p
      rsap = ddot_ (&n, rs, &i_1, ap, &i_1); // rsap = (r*, A.p)
      delta = - rho / rsap;

      dcopy_ (&n, r, &i_1, s, &i_1); // s = r
      daxpy_ (&n, &delta, ap, &i_1, s, &i_1); // s = r + delta A.p
      atimes (n, s, t, user_data); // t = A.s

      st = ddot_ (&n, s, &i_1, t, &i_1); // st = (s, t)
      t2 = ddot_ (&n, t, &i_1, t, &i_1); // t2 = (t, t)
      gamma = - st / t2;

      dcopy_ (&n, s, &i_1, r, &i_1); // r = s
      daxpy_ (&n, &gamma, t, &i_1, r, &i_1); // r = s + gamma t

      daxpy_ (&n, &delta, p, &i_1, x, &i_1); // x = x + delta p...
      daxpy_ (&n, &gamma, s, &i_1, x, &i_1); //       + gamma s

      res2 = ddot_ (&n, r, &i_1, r, &i_1);
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-bicgstab %d %e\n", i, res2);
	}
      if (res2 < tol2) break;

      rho1 = ddot_ (&n, rs, &i_1, r, &i_1); // rho = (r*, r)
      beta = rho1 / rho * delta / gamma;
      rho = rho1;

      daxpy_ (&n, &gamma, ap, &i_1, p, &i_1); // p = p + gamma A.p
      dscal_ (&n, &beta, p, &i_1); // p = beta (p + gamma A.p)
      daxpy_ (&n, &d_1, r, &i_1, p, &i_1); // p = r + beta(p + gamma A.p)
    }

# else // !HAVE_BLAS_H
  /* use local BLAS routines */

  // initial guess
  my_dcopy (n, b, 1, x, 1);
  
  atimes (n, x, r, user_data); // r = A.x
  my_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b

  my_dcopy (n, r, 1, rs, 1); // r* = r
  my_dcopy (n, r, 1, p, 1); // p = r

  rho = my_ddot (n, rs, 1, r, 1); // rho = (r*, r)

  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, user_data); // ap = A.p
      rsap = my_ddot (n, rs, 1, ap, 1); // rsap = (r*, A.p)
      delta = - rho / rsap;

      my_dcopy (n, r, 1, s, 1); // s = r
      my_daxpy (n, delta, ap, 1, s, 1); // s = r + delta A.p
      atimes (n, s, t, user_data); // t = A.s

      st = my_ddot (n, s, 1, t, 1); // st = (s, t)
      t2 = my_ddot (n, t, 1, t, 1); // t2 = (t, t)
      gamma = - st / t2;

      my_dcopy (n, s, 1, r, 1); // r = s
      my_daxpy (n, gamma, t, 1, r, 1); // r = s + gamma t

      my_daxpy (n, delta, p, 1, x, 1); // x = x + delta p...
      my_daxpy (n, gamma, s, 1, x, 1); //       + gamma s

      res2 = my_ddot (n, r, 1, r, 1);
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-bicgstab %d %e\n", i, res2);
	}
      if (res2 < tol2) break;

      rho1 = my_ddot (n, rs, 1, r, 1); // rho = (r*, r)
      beta = rho1 / rho * delta / gamma;
      rho = rho1;

      my_daxpy (n, gamma, ap, 1, p, 1); // p = p + gamma A.p
      my_dscal (n, beta, p, 1); // p = beta (p + gamma A.p)
      my_daxpy (n, 1.0, r, 1, p, 1); // p = r + beta(p + gamma A.p)
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (r);
  free (rs);
  free (p);
  free (ap);
  free (s);
  free (t);

  if (it_param->debug == 1)
    {
      fprintf (it_param->out, "libiter-bicgstab %d %e\n", i, res2);
    }
}
