/* Steepest Descent -- Weiss' Algorithm 1
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: steepest.c,v 2.1 2006/10/09 20:09:24 ichiki Exp $
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* BLAS functions */
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

/* Steepest Descent -- Weiss' Algorithm 1
 */
void
steepest (int n, const double *b, double *x,
	  double tol, int itmax,
	  int *iter, double *res,
	  void (*atimes) (int, const double *, double *, void *),
	  void * user_data)
{
  int i;

  int i_1 = 1;
  double d_m1 = -1.0;

  double * r;
  double * ar;

  double r2 = 0.0;
  double rar;

  double delta;

  r  = (double *) malloc (sizeof (double) * n);
  ar = (double *) malloc (sizeof (double) * n);

  // initial guess
  dcopy_ (&n, b, &i_1, x, &i_1);

  atimes (n, x, r, user_data); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b
  
  for (i = 0; i < itmax; i ++)
    {
      r2 = ddot_ (&n, r, &i_1, r, &i_1); // r2 = (r, r)
      if (r2 <= tol) break;
      
      atimes (n, r, ar, user_data); // ar = A.r
      rar = ddot_ (&n, r, &i_1, ar, &i_1); // rar = (r, A.r)

      delta = - r2 / rar; // delta = - (r, r) / (r, A.r)
      
      daxpy_ (&n, &delta, r, &i_1, x, &i_1); // x += delta r
      daxpy_ (&n, &delta, ar, &i_1, r, &i_1); // r += delta A.r

    }
  free (r);
  free (ar);

  *iter = i;
  *res = sqrt (r2);
}
