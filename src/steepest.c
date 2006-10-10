/* Steepest Descent -- Weiss' Algorithm 1
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: steepest.c,v 2.2 2006/10/10 18:06:16 ichiki Exp $
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



/* Steepest Descent -- Weiss' Algorithm 1
 */
void
steepest (int n, const double *b, double *x,
	  void (*atimes) (int, const double *, double *, void *),
	  void * user_data,
	  struct iter * it_param)
{
  double tol, tol2;
  int itmax;

  int i;

  double * r;
  double * ar;

  double r2 = 0.0;
  double rar;

  double delta;

#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_m1 = -1.0;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


  tol = it_param->eps;
  tol2 = tol * tol;
  itmax = it_param->max;

  r  = (double *) malloc (sizeof (double) * n);
  ar = (double *) malloc (sizeof (double) * n);

#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */

  // initial guess
  cblas_dcopy (n, b, 1, x, 1);

  atimes (n, x, r, user_data); // r = A.x
  cblas_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  
  for (i = 0; i < itmax; i ++)
    {
      r2 = cblas_ddot (n, r, 1, r, 1); // r2 = (r, r)
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-steepest %d %e\n", i, r2);
	}
      if (r2 <= tol2) break;
      
      atimes (n, r, ar, user_data); // ar = A.r
      rar = cblas_ddot (n, r, 1, ar, 1); // rar = (r, A.r)

      delta = - r2 / rar; // delta = - (r, r) / (r, A.r)
      
      cblas_daxpy (n, delta, r, 1, x, 1); // x += delta r
      cblas_daxpy (n, delta, ar, 1, r, 1); // r += delta A.r

    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  // initial guess
  dcopy_ (&n, b, &i_1, x, &i_1);

  atimes (n, x, r, user_data); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b
  
  for (i = 0; i < itmax; i ++)
    {
      r2 = ddot_ (&n, r, &i_1, r, &i_1); // r2 = (r, r)
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-steepest %d %e\n", i, r2);
	}
      if (r2 <= tol2) break;
      
      atimes (n, r, ar, user_data); // ar = A.r
      rar = ddot_ (&n, r, &i_1, ar, &i_1); // rar = (r, A.r)

      delta = - r2 / rar; // delta = - (r, r) / (r, A.r)
      
      daxpy_ (&n, &delta, r, &i_1, x, &i_1); // x += delta r
      daxpy_ (&n, &delta, ar, &i_1, r, &i_1); // r += delta A.r

    }

# else // !HAVE_BLAS_H
  /* use local BLAS routines */

  // initial guess
  my_dcopy (n, b, 1, x, 1);

  atimes (n, x, r, user_data); // r = A.x
  my_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  
  for (i = 0; i < itmax; i ++)
    {
      r2 = my_ddot (n, r, 1, r, 1); // r2 = (r, r)
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-steepest %d %e\n", i, r2);
	}
      if (r2 <= tol2) break;
      
      atimes (n, r, ar, user_data); // ar = A.r
      rar = my_ddot (n, r, 1, ar, 1); // rar = (r, A.r)

      delta = - r2 / rar; // delta = - (r, r) / (r, A.r)
      
      my_daxpy (n, delta, r, 1, x, 1); // x += delta r
      my_daxpy (n, delta, ar, 1, r, 1); // r += delta A.r

    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


  free (r);
  free (ar);

  if (it_param->debug == 1)
    {
      fprintf (it_param->out, "libiter-steepest %d %e\n", i, r2);
    }
}
