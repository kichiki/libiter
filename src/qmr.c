/* QMR -- Weiss' Algorithm 10
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: qmr.c,v 2.1 2006/10/09 20:09:24 ichiki Exp $
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


/* QMR -- Weiss' Algorithm 10
 */
void
qmr (int n, const double *b, double *x,
     double tol, int itmax,
     int *iter, double *res,
     void (*atimes) (int, const double *, double *, void *),
     void (*atimes_trans) (int, const double *, double *, void *),
     void * user_data)
{
  int i;

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;

  double * xt;
  double * r;
  double * rt;
  double * rs;
  double * p;
  double * ps;
  double * atps;
  double * ap;

  double t2;
  double rt2;
  double patps;
  double rtrs, rtrs1;

  double delta;
  double beta;
  double gamma;
  double d_1_gamma;

  double tol2;
  double res2 = 0.0;


  tol2 = tol * tol;

  xt = (double *) malloc (sizeof (double) * n);
  r  = (double *) malloc (sizeof (double) * n);
  rt = (double *) malloc (sizeof (double) * n);
  rs = (double *) malloc (sizeof (double) * n);
  p  = (double *) malloc (sizeof (double) * n);
  ps = (double *) malloc (sizeof (double) * n);
  atps = (double *) malloc (sizeof (double) * n);
  ap = (double *) malloc (sizeof (double) * n);

  // initial guess
  dcopy_ (&n, b, &i_1, x, &i_1);

  dcopy_ (&n, x, &i_1, xt, &i_1); // xt = x
  
  atimes (n, x, r, user_data); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b

  dcopy_ (&n, r, &i_1, rt, &i_1); // rt = r
  dcopy_ (&n, r, &i_1, rs, &i_1); // r* = r
  dcopy_ (&n, r, &i_1, p,  &i_1); // p  = r
  dcopy_ (&n, r, &i_1, ps, &i_1); // p* = r

  t2 = ddot_ (&n, rt, &i_1, rt, &i_1); // t2 = (rt, rt)

  for (i = 0; i < itmax; i ++)
    {
      atimes_trans (n, ps, atps, user_data); // atps = At.p*
      patps = ddot_ (&n, p, &i_1, atps, &i_1); // patps = (p, At.p*)
      rtrs = ddot_ (&n, rt, &i_1, rs, &i_1); // rtrs = (rt, r*)
      delta = - rtrs / patps;

      atimes (n, p, ap, user_data); // ap = A.p
      daxpy_ (&n, &delta, ap, &i_1, rt, &i_1); // rt += delta A.ps
      daxpy_ (&n, &delta, atps, &i_1, rs, &i_1); // r* += delta At.ps

      rtrs1 = ddot_ (&n, rt, &i_1, rs, &i_1); // rtrs = (rt, r*) for news
      beta = rtrs1 / rtrs;
      rtrs = rtrs1;

      daxpy_ (&n, &delta, p, &i_1, xt, &i_1); // xt += delta p(old)

      dscal_ (&n, &beta, p, &i_1); // p = beta p(old)
      daxpy_ (&n, &d_1, rt, &i_1, p, &i_1); // p = rt + beta p(old)

      dscal_ (&n, &beta, ps, &i_1); // p* = beta p*(old)
      daxpy_ (&n, &d_1, rs, &i_1, ps, &i_1); // p* = r* + beta p*(old)

      rt2 = ddot_ (&n, rt, &i_1, rt, &i_1); // rt2 = (rt, rt)
      t2 = 1.0 / ((1.0 / t2) + (1.0 / rt2));
      gamma = t2 / rt2;

      d_1_gamma = 1.0 - gamma;
      dscal_ (&n, &d_1_gamma, x, &i_1); // x = (1-gamma) x(old)
      daxpy_ (&n, &gamma, xt, &i_1, x, &i_1); // x += gamma(xt - x(old))

      dscal_ (&n, &d_1_gamma, r, &i_1); // r = (1-gamma) r(old)
      daxpy_ (&n, &gamma, rt, &i_1, r, &i_1); // r += gamma(rt - r(old))

      res2 = ddot_ (&n, r, &i_1, r, &i_1);
      if (res2 < tol2) break;
    }

  free (xt);
  free (r);
  free (rt);
  free (rs);
  free (p);
  free (ps);
  free (atps);
  free (ap);

  *iter = i;
  *res = sqrt (res2);
}
