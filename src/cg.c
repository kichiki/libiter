/* Classical CG method -- Weiss' Algorithm 2
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: cg.c,v 2.2 2006/10/10 18:07:10 ichiki Exp $
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



/* Classical CG method -- Weiss' Algorithm 2
 */
void
cg (int n, const double *b, double *x,
    void (*atimes) (int, const double *, double *, void *),
    void * user_data,
    struct iter * it_param)
{
  double tol, tol2;
  int itmax;

  int i;

  double * p;
  double * r;
  double * ap;

  double r2;
  double rr2 = 0.0;
  double pap;

  double gamma;
  double beta;

#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;
  double scale;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


  tol = it_param->eps;
  tol2 = tol * tol;
  itmax = it_param->max;

  p  = (double *) malloc (sizeof (double) * n);
  r  = (double *) malloc (sizeof (double) * n);
  ap = (double *) malloc (sizeof (double) * n);

#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */

  // initial guess
  cblas_dcopy (n, b, 1, x, 1);

  atimes (n, x, r, user_data); // r = A.x
  cblas_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  
  cblas_dcopy (n, r, 1, p, 1); // p = r

  for (i = 0; i < itmax; i ++)
    {
      r2 = cblas_ddot (n, r, 1, r, 1); // r2 = (r, r)
      
      atimes (n, p, ap, user_data); // ap = A.p
      pap = cblas_ddot (n, p, 1, ap, 1); // pap = (p, A.p)

      gamma = - r2 / pap; // gamma = - (r, r) / (p, A.p)
      
      cblas_daxpy (n, gamma, p, 1, x, 1); // x += gamma p
      cblas_daxpy (n, gamma, ap, 1, r, 1); // r += gamma Ap

      // new norm of r
      rr2 = cblas_ddot (n, r, 1, r, 1); // (r, r)
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-cg %d %e\n", i, rr2);
	}
      if (rr2 <= tol2) break;

      beta = rr2 / r2; // beta = (r, r) / (r0, r0)
      
      cblas_dscal (n, beta, p, 1); // p *= beta
      cblas_daxpy (n, 1.0, r, 1, p, 1); // p = r + beta p

      r2 = rr2;
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  // initial guess
  dcopy_ (&n, b, &i_1, x, &i_1);

  atimes (n, x, r, user_data); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b
  
  dcopy_ (&n, r, &i_1, p, &i_1); // p = r

  for (i = 0; i < itmax; i ++)
    {
      r2 = ddot_ (&n, r, &i_1, r, &i_1); // r2 = (r, r)
      
      atimes (n, p, ap, user_data); // ap = A.p
      pap = ddot_ (&n, p, &i_1, ap, &i_1); // pap = (p, A.p)

      gamma = - r2 / pap; // gamma = - (r, r) / (p, A.p)
      
      daxpy_ (&n, &gamma, p, &i_1, x, &i_1); // x += gamma p
      daxpy_ (&n, &gamma, ap, &i_1, r, &i_1); // r += gamma Ap

      // new norm of r
      rr2 = ddot_ (&n, r, &i_1, r, &i_1); // (r, r)
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-cg %d %e\n", i, rr2);
	}
      if (rr2 <= tol2) break;

      beta = rr2 / r2; // beta = (r, r) / (r0, r0)
      
      dscal_ (&n, &beta, p, &i_1); // p *= beta
      daxpy_ (&n, &d_1, r, &i_1, p, &i_1); // p = r + beta p

      r2 = rr2;
    }

# else // !HAVE_BLAS_H
  /* use local BLAS routines */

  // initial guess
  my_dcopy (n, b, 1, x, 1);

  atimes (n, x, r, user_data); // r = A.x
  my_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  
  my_dcopy (n, r, 1, p, 1); // p = r

  for (i = 0; i < itmax; i ++)
    {
      r2 = my_ddot (n, r, 1, r, 1); // r2 = (r, r)
      
      atimes (n, p, ap, user_data); // ap = A.p
      pap = my_ddot (n, p, 1, ap, 1); // pap = (p, A.p)

      gamma = - r2 / pap; // gamma = - (r, r) / (p, A.p)
      
      my_daxpy (n, gamma, p, 1, x, 1); // x += gamma p
      my_daxpy (n, gamma, ap, 1, r, 1); // r += gamma Ap

      // new norm of r
      rr2 = my_ddot (n, r, 1, r, 1); // (r, r)
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-cg %d %e\n", i, rr2);
	}
      if (rr2 <= tol2) break;

      beta = rr2 / r2; // beta = (r, r) / (r0, r0)
      
      my_dscal (n, beta, p, 1); // p *= beta
      my_daxpy (n, 1.0, r, 1, p, 1); // p = r + beta p

      r2 = rr2;
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (p);
  free (r);
  free (ap);

  if (it_param->debug == 1)
    {
      fprintf (it_param->out, "libiter-cg %d %e\n", i, rr2);
    }
}
