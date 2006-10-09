/* BiCGSTAB - Weiss, Algorithm 12
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bicgstab.c,v 2.1 2006/10/09 20:09:24 ichiki Exp $
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


/* Ref: Weiss, Algorithm 12 BiCGSTAB
 */
void
bicgstab (int n, const double *b, double *x,
	  double tol, int itmax,
	  int *iter, double *res,
	  void (*atimes) (int, const double *, double *, void *),
	  void * user_data)
{
  int i;

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;

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

  double tol2;
  double res2 = 0.0;


  tol2 = tol * tol;

  r  = (double *) malloc (sizeof (double) * n);
  rs = (double *) malloc (sizeof (double) * n);
  p  = (double *) malloc (sizeof (double) * n);
  ap = (double *) malloc (sizeof (double) * n);
  s  = (double *) malloc (sizeof (double) * n);
  t  = (double *) malloc (sizeof (double) * n);

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
      if (res2 < tol2) break;

      rho1 = ddot_ (&n, rs, &i_1, r, &i_1); // rho = (r*, r)
      beta = rho1 / rho * delta / gamma;
      rho = rho1;

      daxpy_ (&n, &gamma, ap, &i_1, p, &i_1); // p = p + gamma A.p
      dscal_ (&n, &beta, p, &i_1); // p = beta (p + gamma A.p)
      daxpy_ (&n, &d_1, r, &i_1, p, &i_1); // p = r + beta(p + gamma A.p)
    }

  free (r);
  free (rs);
  free (p);
  free (ap);
  free (s);
  free (t);

  *iter = i;
  *res = sqrt (res2);
}
