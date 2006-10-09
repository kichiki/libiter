/* CGNE -- Weiss, Algorithm 6
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: cgne.c,v 1.3 2006/10/09 22:02:11 ichiki Exp $
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
#include "myblas.h"


/*  for details to the algorithm see
 *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
 *  Algorithm 6 : CGNE
 */
void
cgne (int n, const double *b, double *x,
      double eps, int maxiter,
      int *iter, double *res,
      void (*atimes) (int, const double *, double *, void *),
      void (*atimes_t) (int, const double *, double *, void *),
      void * user_data)
{
  double phi;
  double sum;
  double *alpha;
  double b_norm;
  double *r;
  double *xx;
  double *tmp;
  double *tmp0;
  double tiny = 1.0e-15;
  int i;
  int mink2;


  alpha = (double *) malloc (sizeof (double) * 2);
  r     = (double *) malloc (sizeof (double) * 2 * n);
  xx    = (double *) malloc (sizeof (double) * 2 * n);
  tmp   = (double *) malloc (sizeof (double) * n);
  tmp0  = (double *) malloc (sizeof (double) * n);


  b_norm = my_dnrm2 (n, b, 1);

  (*iter) = 0;
  /* r = A * x(i) - b */
  if (fabs (my_dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp, user_data);
      my_daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    {
      my_dscalz (n, - 1.0, b, 1, r, 1);
    }


  /* initial condition */
  (*res) = my_dnrm2 (n, r, 1);
  my_dcopy (n, x, 1, xx, 1);

  while ((*res) > eps * b_norm
	 && (*iter) < maxiter)
    {
      /*fprintf (stderr, "#CGNE %d %e\n", (*iter), (*res));*/
      (*iter) ++;

      /* set min (k, 2) */
      if ((*iter) < 2) mink2 = (*iter);
      else             mink2 = 2;

      /* set alpha */
      for (i=0; i<mink2; i++)
	{
	  atimes_t (n, & r [i * n], tmp, user_data);
	  atimes_t (n, & r [0], tmp0, user_data);
	  alpha [i] =
	    - my_ddot (n, tmp, 1, tmp0, 1)
	    / my_ddot (n, & r [i * n], 1, & r [i * n], 1);
	}


      /* set phi */
      sum = 0.0;
      for (i=0; i<mink2; i++) sum += alpha [i];
      phi = 1.0 / sum;

      /* update xx */
      my_dcopy (n, tmp0, 1, tmp, 1);
      for (i=0; i<mink2; i++)
	{
	  my_daxpy (n, alpha [i], & xx [i * n], 1, tmp, 1);
	}
      my_dcopy (n, & xx [0], 1, & xx [n], 1);
      my_dscalz (n, phi, tmp, 1, xx, 1);

      /* update r */
      atimes (n, tmp0, tmp, user_data);
      for (i=0; i<mink2; i++)
	{
	  my_daxpy (n, alpha [i], & r [i * n], 1, tmp, 1);
	}
      my_dcopy (n, & r [0], 1, & r [n], 1);
      my_dscalz (n, phi, tmp, 1, r, 1);

      (*res) = my_dnrm2 (n, r, 1);
    }

  my_dcopy (n, xx, 1, x, 1);

  free (alpha);
  free (r);
  free (xx);
  free (tmp);
  free (tmp0);
}
