/* generalized minimum residual method
 * Copyright (C) 1998-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: gmres.c,v 2.8 2006/10/09 21:56:59 ichiki Exp $
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
#include "../config.h"

#include <stdlib.h> /* malloc (), free() */
#include <stdio.h>
#include <math.h>


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

#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;
  double scale;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


  v   = (double *) malloc (sizeof (double) * (m + 1) * n);
  h   = (double *) malloc (sizeof (double) * m * m);
  g   = (double *) malloc (sizeof (double) * m + 1);
  c   = (double *) malloc (sizeof (double) * m);
  s   = (double *) malloc (sizeof (double) * m);


  (*iter) = 0;
  /* 1. start: */
  /* compute r0 */
  /* compute v1 */
  /* beta */
  myatimes (n, x, v + 0, user_data); /* use v [0] temporaliry */

#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */

  // v = f - v
  cblas_dscal (n, -1.0, v, 1); // v = - v
  cblas_daxpy (n, 1.0, f, 1, v, 1); // v = f - v

  //g[0] = cblas_dnrm2 (n, v, 1);
  g[0] = sqrt (cblas_ddot (n, v, 1, v, 1));
  cblas_dscal (n, 1.0 / g[0], v, 1);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  // v = f - v
  dscal_ (&n, &d_m1, v, &i_1); // v = - v
  daxpy_ (&n, &d_1, f, &i_1, v, &i_1); // v = f - v

  g[0] = sqrt (ddot_ (&n, v, &i_1, v, &i_1));
  scale = 1.0 / g[0];
  dscal_ (&n, &scale, v, &i_1);

# else // !HAVE_BLAS_H
  /* use local BLAS routines */

  my_daxpyz (n, -1.0, v + 0, 1, f, 1, v + 0, 1);

  g [0] = my_dnrm2 (n, v + 0, 1);
  my_dscal (n, 1.0 / g [0], v + 0, 1);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


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
#ifdef HAVE_CBLAS_H
	      /* use ATLAS' CBLAS routines */

	      h [i * m + j] =
		cblas_ddot (n, v + (j + 1) * n, 1,
			    v + i * n, 1);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
	      /* use Fortran BLAS routines */

	      h [i * m + j] =
		ddot_ (&n, v + (j + 1) * n, &i_1,
		       v + i * n, &i_1);

# else // !HAVE_BLAS_H
	      /* use local BLAS routines */

	      h [i * m + j] =
		my_ddot (n, v + (j + 1) * n, 1,
			 v + i * n, 1);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
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
	  /* v_j+1 */
#ifdef HAVE_CBLAS_H
	  /* use ATLAS' CBLAS routines */

	  //hh = cblas_dnrm2 (n, v + (j + 1) * n, 1);
	  hh = sqrt (cblas_ddot (n, v + (j + 1) * n, 1, v + (j + 1) * n, 1));
	  cblas_dscal (n, 1.0 / hh, v + (j + 1) * n, 1);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
	  /* use Fortran BLAS routines */

	  //hh = dnrm2_ (&n, v + (j + 1) * n, &i_1);
	  hh = sqrt (ddot_ (&n, v + (j + 1) * n, &i_1, v + (j + 1) * n, &i_1));
	  scale = 1.0 / hh;
	  dscal_ (&n, &scale, v + (j + 1) * n, &i_1);

# else // !HAVE_BLAS_H
	  /* use local BLAS routines */

	  hh = my_dnrm2 (n, v + (j + 1) * n, 1);
	  my_dscal (n, 1.0 / hh, v + (j + 1) * n, 1);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

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
      /* compute v1 */

#ifdef HAVE_CBLAS_H
      /* use ATLAS' CBLAS routines */

      // v = f - v
      cblas_dscal (n, -1.0, v, 1); // v = - v
      cblas_daxpy (n, 1.0, f, 1, v, 1); // v = f - v

      //g [0] = cblas_dnrm2 (n, v, 1);
      g[0] = sqrt (cblas_ddot (n, v, 1, v, 1));
      cblas_dscal (n, 1.0 / g[0], v, 1);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
      /* use Fortran BLAS routines */

      // v = f - v
      dscal_ (&n, &d_m1, v, &i_1); // v = - v
      daxpy_ (&n, &d_1, f, &i_1, v, &i_1); // v = f - v

      g[0] = sqrt (ddot_ (&n, v, &i_1, v, &i_1));
      scale = 1.0 / g[0];
      dscal_ (&n, &scale, v, &i_1);

# else // !HAVE_BLAS_H
      /* use local BLAS routines */

      my_daxpyz (n, -1.0, v + 0, 1, f, 1, v + 0, 1);

      g [0] = my_dnrm2 (n, v + 0, 1);
      my_dscal (n, 1.0 / g [0], v + 0, 1);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
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

#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;
  double scale;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


  m = itmax;

  v   = (double *) malloc (sizeof (double) * (m + 1) * n);
  h   = (double *) malloc (sizeof (double) * m * m);
  g   = (double *) malloc (sizeof (double) * m + 1);
  c   = (double *) malloc (sizeof (double) * m);
  s   = (double *) malloc (sizeof (double) * m);


  /* 1. start: */
  /* compute r0 */
  /* compute v1 */
  /* beta */
  myatimes (n, x, v + 0, user_data); /* use v [0] temporaliry */

#ifdef HAVE_CBLAS_H
  /* use ATLAS' CBLAS routines */

  // v = f - v
  cblas_dscal (n, -1.0, v, 1); // v = - v
  cblas_daxpy (n, 1.0, f, 1, v, 1); // v = f - v

  g[0] = sqrt (cblas_ddot (n, v, 1, v, 1));
  cblas_dscal (n, 1.0 / g[0], v, 1);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  // v = f - v
  dscal_ (&n, &d_m1, v, &i_1); // v = - v
  daxpy_ (&n, &d_1, f, &i_1, v, &i_1); // v = f - v

  g[0] = sqrt (ddot_ (&n, v, &i_1, v, &i_1));
  scale = 1.0 / g[0];
  dscal_ (&n, &scale, v, &i_1);

# else // !HAVE_BLAS_H
  /* use local BLAS routines */

  my_daxpyz (n, -1.0, v + 0, 1, f, 1, v + 0, 1);

  g [0] = my_dnrm2 (n, v + 0, 1);
  my_dscal (n, 1.0 / g [0], v + 0, 1);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  /* main loop */
  /* 2. iterate: */
  for (j = 0; j < m; j ++)
    {
      /* tmp = A.vj : use v [(j + 1) * n] directly */
      myatimes (n, v + j * n, v + (j + 1) * n, user_data);
      /* h_i,j (i=1,...,j) */
      for (i = 0; i <= j; i ++)
	{
#ifdef HAVE_CBLAS_H
	  /* use ATLAS' CBLAS routines */

	  h [i * m + j] =
	    cblas_ddot (n, v + (j + 1) * n, 1,
			v + i * n, 1);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
	  /* use Fortran BLAS routines */

	  h [i * m + j] =
	    ddot_ (&n, v + (j + 1) * n, &i_1,
		   v + i * n, &i_1);

# else // !HAVE_BLAS_H
	  /* use local BLAS routines */

	  h [i * m + j] =
	    my_ddot (n, v + (j + 1) * n, 1,
		     v + i * n, 1);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
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
      /* v_j+1 */
#ifdef HAVE_CBLAS_H
      /* use ATLAS' CBLAS routines */

      //hh = cblas_dnrm2 (n, v + (j + 1) * n, 1);
      hh = sqrt (cblas_ddot (n, v + (j + 1) * n, 1, v + (j + 1) * n, 1));
      cblas_dscal (n, 1.0 / hh, v + (j + 1) * n, 1);

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
      /* use Fortran BLAS routines */

      //hh = dnrm2_ (&n, v + (j + 1) * n, &i_1);
      hh = sqrt (ddot_ (&n, v + (j + 1) * n, &i_1, v + (j + 1) * n, &i_1));
      scale = 1.0 / hh;
      dscal_ (&n, &scale, v + (j + 1) * n, &i_1);

# else // !HAVE_BLAS_H
      /* use local BLAS routines */

      hh = my_dnrm2 (n, v + (j + 1) * n, 1);
      my_dscal (n, 1.0 / hh, v + (j + 1) * n, 1);

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

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

