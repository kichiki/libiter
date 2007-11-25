/* BICO -- Weiss' Algorithm 9
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bico.c,v 2.3 2007/11/25 18:44:42 kichiki Exp $
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
#include "libiter.h" // struct iter
#include "memory-check.h" // CHECK_MALLOC


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


/* BICO -- Weiss' Algorithm 9
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
bico (int n, const double *b, double *x,
      void (*atimes) (int, const double *, double *, void *),
      void (*atimes_t) (int, const double *, double *, void *),
      void *atimes_param,
      struct iter *it)
{
  int ret = -1;
  int itmax = it->max;
  double eps2 = it->eps * it->eps;

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;

  double *xt   = (double *)malloc (sizeof (double) * n);
  double *r    = (double *)malloc (sizeof (double) * n);
  double *rt   = (double *)malloc (sizeof (double) * n);
  double *rs   = (double *)malloc (sizeof (double) * n);
  double *p    = (double *)malloc (sizeof (double) * n);
  double *ps   = (double *)malloc (sizeof (double) * n);
  double *atps = (double *)malloc (sizeof (double) * n);
  double *ap   = (double *)malloc (sizeof (double) * n);
  double *rtmr = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (xt,   "bico");
  CHECK_MALLOC (r,    "bico");
  CHECK_MALLOC (rt,   "bico");
  CHECK_MALLOC (rs,   "bico");
  CHECK_MALLOC (p,    "bico");
  CHECK_MALLOC (ps,   "bico");
  CHECK_MALLOC (atps, "bico");
  CHECK_MALLOC (ap,   "bico");
  CHECK_MALLOC (rtmr, "bico");

  double b2 = ddot_ (&n, b, &i_1, b, &i_1);
  eps2 *= b2;
  double res2 = 0.0;

  dcopy_ (&n, x, &i_1, xt, &i_1); // xt = x
  
  atimes (n, x, r, atimes_param); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b

  dcopy_ (&n, r, &i_1, rt, &i_1); // rt = r
  dcopy_ (&n, r, &i_1, rs, &i_1); // r* = r
  dcopy_ (&n, r, &i_1, p,  &i_1); // p  = r
  dcopy_ (&n, r, &i_1, ps, &i_1); // p* = r


  int i;
  for (i = 0; i < itmax; i ++)
    {
      atimes_t (n, ps, atps, atimes_param); // atps = At.p*
      double patps = ddot_ (&n, p, &i_1, atps, &i_1); // patps = (p, At.p*)
      double rtrs = ddot_ (&n, rt, &i_1, rs, &i_1); // rtrs = (rt, r*)
      double delta = - rtrs / patps;

      atimes (n, p, ap, atimes_param); // ap = A.p
      daxpy_ (&n, &delta, ap, &i_1, rt, &i_1); // rt += delta A.ps
      daxpy_ (&n, &delta, atps, &i_1, rs, &i_1); // r* += delta At.ps

      double rtrs1 = ddot_ (&n, rt, &i_1, rs, &i_1); // rtrs = (rt, r*) for news
      double beta = rtrs1 / rtrs;
      rtrs = rtrs1;

      daxpy_ (&n, &delta, p, &i_1, xt, &i_1); // xt += delta p(old)

      dscal_ (&n, &beta, p, &i_1); // p = beta p(old)
      daxpy_ (&n, &d_1, rt, &i_1, p, &i_1); // p = rt + beta p(old)

      dscal_ (&n, &beta, ps, &i_1); // p* = beta p*(old)
      daxpy_ (&n, &d_1, rs, &i_1, ps, &i_1); // p* = r* + beta p*(old)

      dcopy_ (&n, rt, &i_1, rtmr, &i_1); // rtmr = rt
      daxpy_ (&n, &d_m1, r, &i_1, rtmr, &i_1); // rtmr = rt - r
      // rtmr2 = (rt - r, rt - r)
      double rtmr2 = ddot_ (&n, rtmr, &i_1, rtmr, &i_1);
      double rrtmr = ddot_ (&n, r, &i_1, rtmr, &i_1); // rrtmr = (r, rt - r)
      double gamma = - rrtmr / rtmr2;

      double d_1_gamma = 1.0 - gamma;
      dscal_ (&n, &d_1_gamma, x, &i_1); // x = (1-gamma) x(old)
      daxpy_ (&n, &gamma, xt, &i_1, x, &i_1); // x += gamma(xt - x(old))

      dscal_ (&n, &d_1_gamma, r, &i_1); // r = (1-gamma) r(old)
      daxpy_ (&n, &gamma, rt, &i_1, r, &i_1); // r += gamma(rt - r(old))

      res2 = ddot_ (&n, r, &i_1, r, &i_1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-bico %d %e\n", i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}
    }

  free (xt);
  free (r);
  free (rt);
  free (rs);
  free (p);
  free (ps);
  free (atps);
  free (ap);
  free (rtmr);

  if (it->debug == 1)
    {
      fprintf (it->out, "libiter-bico %d %e\n", i, res2 / b2);
    }

  it->niter = i;
  it->res2  = res2 / b2;
  return (ret);
}
