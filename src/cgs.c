/* CGS -- Weiss, Algorithm 11
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: cgs.c,v 2.1 2006/10/09 20:09:24 ichiki Exp $
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


/* Ref: Weiss, Algorithm 11 CGS
 */
void
cgs (int n, const double *b, double *x,
     double tol, int itmax,
     int *iter, double *res,
     void (*atimes) (int, const double *, double *, void *),
     void * user_data)
{
  int i;

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;
  double d_2 = 2.0;

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

  double tol2;
  double res2 = 0.0;


  tol2 = tol * tol;

  r  = (double *) malloc (sizeof (double) * n);
  r0 = (double *) malloc (sizeof (double) * n);
  p  = (double *) malloc (sizeof (double) * n);
  u  = (double *) malloc (sizeof (double) * n);
  ap = (double *) malloc (sizeof (double) * n);
  q  = (double *) malloc (sizeof (double) * n);
  t  = (double *) malloc (sizeof (double) * n);

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

  free (r);
  free (r0);
  free (p);
  free (u);
  free (ap);
  free (q);
  free (t);

  *iter = i;
  *res = sqrt (res2);
}
