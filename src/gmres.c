/* generalized minimum residual method
 * Copyright (C) 1998-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: gmres.c,v 2.11 2007/10/27 03:23:01 kichiki Exp $
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
#include "libiter.h"
#include "memory-check.h" // CHECK_MALLOC


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
	 void (*myatimes) (int, const double *, double *, void *),
	 void * user_data,
	 struct iter * it_param)
{
  double res = 0.0;

  /* solve linear system A.x = f */
  /* n: dimension of this sysmtem */
  /* m: # of iteration at once */
  int i, j, k;
  double hv;
  double rr, hh;
  double r1, r2;
  double g0;

#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;
  double scale;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


  int m = it_param->restart;
  int itmax = it_param->max;
  double tol = it_param->eps;

  double *v = (double *)malloc (sizeof (double) * (m + 1) * n);
  double *h = (double *)malloc (sizeof (double) * m * m);
  double *g = (double *)malloc (sizeof (double) * m + 1);
  double *c = (double *)malloc (sizeof (double) * m);
  double *s = (double *)malloc (sizeof (double) * m);
  CHECK_MALLOC (v, "gmres_m");
  CHECK_MALLOC (h, "gmres_m");
  CHECK_MALLOC (g, "gmres_m");
  CHECK_MALLOC (c, "gmres_m");
  CHECK_MALLOC (s, "gmres_m");


  int iter = 0;
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
  while (iter <= itmax)
    {
      ++iter;
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
      res = fabs (g [j/*m*/]); /* residual */
      /*fprintf (stderr, "# iter %d res %e\n", iter, *res);*/
      /* if satisfied, */
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-gmres(%d) %d %d %e\n",
		   m, iter, j, res*res);
	}
      if (res <= tol) break;
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

  free (v);
  free (h);
  free (g);
  free (c);
  free (s);

  /* adjust iter */
  iter *= m;

  if (it_param->debug == 1)
    {
      fprintf (it_param->out, "libiter-gmres(%d) it= %d res^2= %e\n",
	       m, iter, res*res);
    }
}

void
gmres (int n, const double *f, double *x,
       void (*myatimes) (int, const double *, double *, void *),
       void * user_data,
       struct iter * it_param)
{
  double res = 0.0;

  /* solve linear system A.x = f */
  /* n: dimension of this sysmtem */
  int i, j, k;
  double hv;
  double rr, hh;
  double r1, r2;
  double g0;

#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;
  double scale;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


  int m = it_param->max;
  double tol = it_param->eps;

  double *v = (double *) malloc (sizeof (double) * (m + 1) * n);
  double *h = (double *) malloc (sizeof (double) * m * m);
  double *g = (double *) malloc (sizeof (double) * m + 1);
  double *c = (double *) malloc (sizeof (double) * m);
  double *s = (double *) malloc (sizeof (double) * m);
  CHECK_MALLOC (v, "gmres");
  CHECK_MALLOC (h, "gmres");
  CHECK_MALLOC (g, "gmres");
  CHECK_MALLOC (c, "gmres");
  CHECK_MALLOC (s, "gmres");


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

      res = fabs (g [j + 1]); /* residual */
      if (it_param->debug == 2)
	{
	  fprintf (it_param->out, "libiter-gmres %d %e\n",
		   j, res*res);
	}
      /* if satisfied, */
      if (res <= tol)
	{
	  j ++; /* this is because ++(*iter) in gmres(m) */
	  break;
	}
    }

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

  if (it_param->debug == 1)
    {
      fprintf (it_param->out, "libiter-gmres it= %d res^2= %e\n", j, res*res);
    }
}

