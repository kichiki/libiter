/* Classical CG method -- Weiss' Algorithm 2
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: cg.c,v 2.1 2006/10/09 20:09:24 ichiki Exp $
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


/* Classical CG method -- Weiss' Algorithm 2
 */
void
cg (int n, const double *b, double *x,
    double tol, int itmax,
    int *iter, double *res,
    void (*atimes) (int, const double *, double *, void *),
    void * user_data)
{
  int i;

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;
  double * p;
  double * r;
  double * ap;

  double r2, rr2;
  double pap;

  double gamma;
  double beta;

  p  = (double *) malloc (sizeof (double) * n);
  r  = (double *) malloc (sizeof (double) * n);
  ap = (double *) malloc (sizeof (double) * n);

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
      //rr2 = dnrm2_ (&n, r, &i_1);
      rr2 = ddot_ (&n, r, &i_1, r, &i_1);
      if (rr2 <= tol) break;

      beta = rr2 / r2; // beta = (r, r) / (r0, r0)
      
      dscal_ (&n, &beta, p, &i_1); // p *= beta
      daxpy_ (&n, &d_1, r, &i_1, p, &i_1); // p = r + beta p

      r2 = rr2;
    }

  free (p);
  free (r);
  free (ap);

  *iter = i;
  *res = rr2;
}
