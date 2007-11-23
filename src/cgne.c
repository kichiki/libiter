/* CGNE -- Weiss, Algorithm 6
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: cgne.c,v 1.5 2007/11/23 05:01:38 kichiki Exp $
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
#include <stdlib.h>
#include <math.h>
#include "libiter.h" // struct iter
#include "myblas.h"
#include "memory-check.h" // CHECK_MALLOC


/*  for details to the algorithm see
 *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
 *  Algorithm 6 : CGNE
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
 *   x[n] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
void
cgne (int n, const double *b, double *x,
      void (*atimes) (int, const double *, double *, void *),
      void (*atimes_t) (int, const double *, double *, void *),
      void *atimes_param,
      struct iter *it)
{
  double eps = it->eps;
  int maxiter = it->max;

  double *alpha = (double *)malloc (sizeof (double) * 2);
  double *r     = (double *)malloc (sizeof (double) * 2 * n);
  double *xx    = (double *)malloc (sizeof (double) * 2 * n);
  double *tmp   = (double *)malloc (sizeof (double) * n);
  double *tmp0  = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (alpha, "atpres");
  CHECK_MALLOC (r,     "atpres");
  CHECK_MALLOC (xx,    "atpres");
  CHECK_MALLOC (tmp,   "atpres");
  CHECK_MALLOC (tmp0,  "atpres");

  double b_norm = my_dnrm2 (n, b, 1);

  /* r = A * x(i) - b */
  double tiny = 1.0e-15;
  if (fabs (my_dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp, atimes_param);
      my_daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    {
      my_dscalz (n, - 1.0, b, 1, r, 1);
    }

  /* initial condition */
  double res = my_dnrm2 (n, r, 1);
  my_dcopy (n, x, 1, xx, 1);

  int iter = 0;
  while (res > eps * b_norm
	 && iter < maxiter)
    {
      if (it->debug == 2)
	{
	  fprintf (it->out, "CGNE %d %e\n", iter, res);
	}
      iter ++;

      /* set min (k, 2) */
      int mink2;
      if (iter < 2) mink2 = iter;
      else          mink2 = 2;

      /* set alpha */
      int i;
      for (i=0; i<mink2; i++)
	{
	  atimes_t (n, & r [i * n], tmp, atimes_param);
	  atimes_t (n, & r [0], tmp0, atimes_param);
	  alpha [i] =
	    - my_ddot (n, tmp, 1, tmp0, 1)
	    / my_ddot (n, & r [i * n], 1, & r [i * n], 1);
	}


      /* set phi */
      double sum = 0.0;
      for (i=0; i<mink2; i++) sum += alpha [i];
      double phi = 1.0 / sum;

      /* update xx */
      my_dcopy (n, tmp0, 1, tmp, 1);
      for (i=0; i<mink2; i++)
	{
	  my_daxpy (n, alpha [i], & xx [i * n], 1, tmp, 1);
	}
      my_dcopy (n, & xx [0], 1, & xx [n], 1);
      my_dscalz (n, phi, tmp, 1, xx, 1);

      /* update r */
      atimes (n, tmp0, tmp, atimes_param);
      for (i=0; i<mink2; i++)
	{
	  my_daxpy (n, alpha [i], & r [i * n], 1, tmp, 1);
	}
      my_dcopy (n, & r [0], 1, & r [n], 1);
      my_dscalz (n, phi, tmp, 1, r, 1);

      res = my_dnrm2 (n, r, 1);
    }

  my_dcopy (n, xx, 1, x, 1);

  free (alpha);
  free (r);
  free (xx);
  free (tmp);
  free (tmp0);

  it->niter = iter;
  it->res2  = res * res / (b_norm * b_norm);
}
