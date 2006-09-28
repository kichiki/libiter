/* generalized minimum residual method
 * Copyright (C) 1998-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: gmres.c,v 2.7 2006/09/28 04:25:11 kichiki Exp $
 *
 * Reference :
 *   GMRES(m) : Y.Saad & M.H.Schultz, SIAM J.Sci.Stat.Comput.
 *   vol7 (1986) pp.856-869
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
#include <stdlib.h> /* malloc (), free() */
#include <stdio.h>
#include <math.h>
#include "myblas.h"

#include "gmres.h"


/* m  : number of iteration */
/* nn : dimension of matrix r [] (nnxnn) */
static void
back_sub (int m, int nn,
	  const double *r, const double *g,
	  double *y)
{
  int i, j, jj;

  /*for (j = m - 1;j >= 0;j --){*/
  /* above for-loop fail, because j is unsigned!! */
  for (j = m - 1, jj = 0; jj < m; j --, jj ++)
    {
      y [j] = g [j];
      for (i = j + 1; i < m; i ++)
	{
	  y [j] -= r [j * nn + i] * y [i];
	}
      y [j] /= r [j * nn + j];
    }
}
void
gmres_m (int n, const double *f, double *x,
	 int m, double tol, int itmax,
	 int *iter, double *res,
	 void (*myatimes) (int, const double *, double *, void *),
	 void * user_data)
{
  /* solve linear system A.x = f */
  /* n: dimension of this sysmtem */
  /* m: # of iteration at once */
  int i, j, k;
  double hv;
  double rr, hh;
  double r1, r2;
  double g0;
  double * v, * h, * g, * c, * s;


  v   = (double *) malloc (sizeof (double) * (m + 1) * n);
  h   = (double *) malloc (sizeof (double) * m * m);
  g   = (double *) malloc (sizeof (double) * m + 1);
  c   = (double *) malloc (sizeof (double) * m);
  s   = (double *) malloc (sizeof (double) * m);


  (*iter) = 0;
  /* 1. start: */
  /* compute r0 */
  myatimes (n, x, v + 0, user_data); /* use v [0] temporaliry */
  daxpyz (n, -1.0, v + 0, 1, f, 1, v + 0, 1);
  /* compute v1 */
  /* beta */
  g [0] = dnrm2 (n, v + 0, 1);
  dscal (n, 1.0 / g [0], v + 0, 1);

  /* main loop */
  while ((*iter) <= itmax)
    {
      ++(*iter);
      /* 2. iterate: */
      for (j = 0; j < m; j ++)
	{
	  /* tmp = A.vj : use v [(j + 1) * n] directly */
	  myatimes (n, v + j * n, v + (j + 1) * n, user_data);
	  /* h_i,j (i=1,...,j) */
	  for (i = 0; i <= j; i ++)
	    {
	      h [i * m + j] = ddot (n, v + (j + 1) * n, 1,
				    v + i * n, 1);
	    }
	  /* vv_j+1 */
	  for (k = 0; k < n; k ++)
	    {
	      hv = 0.0;
	      for (i = 0; i <= j; i ++)
		{
		  hv += h [i * m + j] * v [i * n + k];
		}
	      v [(j + 1) * n + k] -= hv;
	    }
	  /* h_j+1,j */
	  hh = dnrm2 (n, v + (j + 1) * n, 1);
	  /* v_j+1 */
	  dscal (n, 1.0 / hh, v + (j + 1) * n, 1);
	  /* rotate */
	  for (i = 0; i < j; i ++)
	    {
	      r1 = h [ i      * m + j];
	      r2 = h [(i + 1) * m + j];
	      h [ i      * m + j] = c [i] * r1 - s [i] * r2;
	      h [(i + 1) * m + j] = s [i] * r1 + c [i] * r2;
	    }
	  rr = h [j * m + j];
	  hv = sqrt (rr * rr + hh * hh); /* temporary variable */
	  c [j] =  rr / hv;
	  s [j] = -hh / hv;
	  h [j * m + j] = hv; /* resultant (after rotated) element */

	  g0 = g [j];
	  g [j    ] = c [j] * g0;
	  g [j + 1] = s [j] * g0;
	}
      /* 3. form the approximate solution */
      /* solve y_k */
      back_sub (j/*m*/, m, h, g, c); /* use c [] as y_k */
      /* x_m */
      for (i = 0; i < n; i ++)
	{
	  for (k = 0; k < j/*m*/; k ++)
	    {
	      x [i] += v [k * n + i] * c [k];
	    }
	}

      /* 4. restart */
      (*res) = fabs (g [j/*m*/]); /* residual */
      /*fprintf (stderr, "# iter %d res %e\n", *iter, *res);*/
      /* if satisfied, */
      if ((*res) <= tol) break;
      /* else */
      /* compute r_m */
      myatimes (n, x, v + 0, user_data);
      /* r_m */
      daxpyz (n, -1.0, v + 0, 1, f, 1, v + 0, 1);
      /* compute v1 */
      g [0] = dnrm2 (n, v + 0, 1);
      dscal (n, 1.0 / g [0], v + 0, 1);
    }

  /* adjust iter */
  (*iter) *= m;

  free (v);
  free (h);
  free (g);
  free (c);
  free (s);
}

void
gmres (int n, const double *f, double *x,
       double tol, int itmax,
       int *iter, double *res,
       void (*myatimes) (int, const double *, double *, void *),
       void * user_data)
{
  /* solve linear system A.x = f */
  /* n: dimension of this sysmtem */
  int i, j, k, m;
  double hv;
  double rr, hh;
  double r1, r2;
  double g0;
  double * v, * h, * g, * c, * s;


  m = itmax;

  v   = (double *) malloc (sizeof (double) * (m + 1) * n);
  h   = (double *) malloc (sizeof (double) * m * m);
  g   = (double *) malloc (sizeof (double) * m + 1);
  c   = (double *) malloc (sizeof (double) * m);
  s   = (double *) malloc (sizeof (double) * m);


  /* 1. start: */
  /* compute r0 */
  myatimes (n, x, v + 0, user_data); /* use v [0] temporaliry */
  daxpyz (n, -1.0, v + 0, 1, f, 1, v + 0, 1);
  /* compute v1 */
  /* beta */
  g [0] = dnrm2 (n, v + 0, 1);
  dscal (n, 1.0 / g [0], v + 0, 1);

  /* main loop */
  /* 2. iterate: */
  for (j = 0; j < m; j ++)
    {
      /* tmp = A.vj : use v [(j + 1) * n] directly */
      myatimes (n, v + j * n, v + (j + 1) * n, user_data);
      /* h_i,j (i=1,...,j) */
      for (i = 0; i <= j; i ++)
	{
	  h [i * m + j] = ddot (n, v + (j + 1) * n, 1,
				v + i * n, 1);
	}
      /* vv_j+1 */
      for (k = 0; k < n; k ++)
	{
	  hv = 0.0;
	  for (i = 0; i <= j; i ++)
	    {
	      hv += h [i * m + j] * v [i * n + k];
	    }
	  v [(j + 1) * n + k] -= hv;
	}
      /* h_j+1,j */
      hh = dnrm2 (n, v + (j + 1) * n, 1);
      /* v_j+1 */
      dscal (n, 1.0 / hh, v + (j + 1) * n, 1);
      /* rotate */
      for (i = 0; i < j; i ++)
	{
	  r1 = h [ i      * m + j];
	  r2 = h [(i + 1) * m + j];
	  h [ i      * m + j] = c [i] * r1 - s [i] * r2;
	  h [(i + 1) * m + j] = s [i] * r1 + c [i] * r2;
	}
      rr = h [j * m + j];
      hv = sqrt (rr * rr + hh * hh); /* temporary variable */
      c [j] =  rr / hv;
      s [j] = -hh / hv;
      h [j * m + j] = hv; /* resultant (after rotated) element */

      g0 = g [j];
      g [j    ] = c [j] * g0;
      g [j + 1] = s [j] * g0;

      (*res) = fabs (g [j + 1]); /* residual */
      fprintf (stderr, "# iter %d res %e\n", j, *res);
      /* if satisfied, */
      if ((*res) <= tol)
	{
	  j ++; /* this is because ++(*iter) in gmres(m) */
	  break;
	}
    }
  (*iter) = j;

  /* 3. form the approximate solution */
  /* solve y_k */
  back_sub (j, m, h, g, c); /* use c [] as y_k */
  /* x_m */
  for (i = 0; i < n; i ++)
    {
      for (k = 0; k < j; k ++)
	{
	  x [i] += v [k * n + i] * c [k];
	}
    }

  free (v);
  free (h);
  free (g);
  free (c);
  free (s);
}

