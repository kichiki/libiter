/* generalized minimum residual method
 * Copyright (C) 1998,1999 Kengo ICHIKI <kengo@caltech.edu>
 * $Id: gmres.c,v 1.10 2000/06/27 07:23:54 ichiki Exp $
 *
 * Reference :
 *   GMRES(m) : Y.Saad & M.H.Schultz, SIAM J.Sci.Stat.Comput.
 *   vol7 (1986) pp.856-869
 */

#include <stdlib.h> /* malloc (), free() */
#include <stdio.h>
#include <math.h>
#include "../myroutines.h" /* my_d_malloc() */
#include "blas.h"

#include "mygmres.h"

void
mygmres_m (int n, double *f, double *x,
	   int m, double tol, int itmax,
	   int *iter, double *err,
	   void (*myatimes) (int, double *, double *))
{
  /* solve linear system A.x = f */
  /* n: dimension of this sysmtem */
  /* m: # of iteration at once */
  int i, j, k;
  double hv;
  double rr, hh;
  double r1, r2;
  double g0;
  double *tmp, *v, *h, *g, *c, *s;


  tmp = my_d_malloc (n, "tmp");
  v   = my_d_malloc ((m + 1) * n, "v");
  h   = my_d_malloc (m * m, "h");
  g   = my_d_malloc (m + 1, "g");
  c   = my_d_malloc (m, "c");
  s   = my_d_malloc (m, "s");

  *iter = 0;
  /* 1. start: */
  /* compute r0 */
  myatimes (n, x, tmp);
  daxpyz (n, -1.0, tmp, 1, f, 1, tmp, 1);
  /* compute v1 */
  /* beta */
  g [0] = dnrm2 (n, tmp, 1);
  dscalz (n, 1.0 / g [0], tmp, 1, & v [0], 1);

  /* main loop */
  while (*iter <= itmax)
    {
      ++(*iter);
      /* 2. iterate: */
      for (j=0; j<m; j++)
	{
	  /* tmp = A.vj */
	  myatimes (n, &v [j * n], tmp);
	  /* h_i,j (i=1,...,j) */
	  for (i=0; i<=j; i++)
	    h [i * m + j] = ddot (n, tmp, 1, & v [i * n], 1);
	  /* vv_j+1 */
	  for (k=0; k<n; k++)
	    {
	      hv = 0.0;
	      for (i=0; i<=j; i++)
		hv += h [i * m + j] * v [i * n + k];
	      tmp [k] -= hv;
	    }
	  /* h_j+1,j */
	  hh = dnrm2 (n, tmp, 1);
	  /* v_j+1 */
	  dscalz (n, 1.0 / hh, tmp, 1, & v [(j + 1) * n], 1);
	  /* rotate */
	  for (i=0; i<j; i++)
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
      back_sub (m, m, h, g, c); /* use c [] as y_k */
      /* x_m */
      for (i=0; i<n; i++)
	for (k=0; k<m; k++)
	  x [i] += v [k * n + i] * c [k];

      /* 4. restart */
      *err = fabs (g [m]); /* residual */
      /*fprintf (stderr, "# iter %d res %e\n", *iter, *err);*/
      /* if satisfied, */
      if (*err <= tol) break;
      /* else */
      /* compute r_m */
      /* tmp = A.x_m */
      myatimes (n, x, tmp);
      /* r_m */
      daxpyz (n, -1.0, tmp, 1, f, 1, tmp, 1);
      /* compute v1 */
      g [0] = dnrm2 (n, tmp, 1);
      dscalz (n, 1.0 / g [0], tmp, 1, & v [0], 1);
    }

  /* adjust iter */
  (*iter) *= m;

  free (tmp);
  free (v);
  free (h);
  free (g);
  free (c);
  free (s);
}

void
mygmres (int n, double *f, double *x,
	 double tol, int itmax,
	 int *iter, double *err,
	 void (*myatimes) (int, double *, double *))
{
  /* solve linear system A.x = f */
  /* n: dimension of this sysmtem */
  int i, j, k, m;
  double hv;
  double rr, hh;
  double r1, r2;
  double g0;
  double *tmp, *v, *h, *g, *c, *s;


  m = itmax;

  tmp = my_d_malloc (n, "tmp");
  v   = my_d_malloc ((m + 1) * n, "v");
  h   = my_d_malloc (m * m, "h");
  g   = my_d_malloc (m + 1, "g");
  c   = my_d_malloc (m, "c");
  s   = my_d_malloc (m, "s");


  /* 1. start: */
  /* compute r0 */
  myatimes (n, x, tmp);
  daxpyz (n, -1.0, tmp, 1, f, 1, tmp, 1);
  /* compute v1 */
  g [0] = dnrm2 (n, tmp, 1); /* beta */
  dscalz (n, 1.0 / g [0], tmp, 1, & v [0], 1);

  /* main loop */
  /* 2. iterate: */
  for (j=0; j<m; j++)
    {
      /* tmp = A.vj */
      myatimes (n, &v [j * n], tmp);
      /* h_i,j (i=1,...,j) */
      for (i=0; i<=j; i++)
	h [i * m + j] = ddot (n, tmp, 1, &v [i * n], 1);
      /* vv_j+1 */
      for (k=0; k<n; k++)
	{
	  hv = 0.0;
	  for (i=0; i<=j; i++)
	    hv += h [i * m + j] * v [i * n + k];
	  tmp [k] -= hv;
	}
      /* h_j+1,j */
      hh = dnrm2 (n, tmp, 1);
      /* v_j+1 */
      dscalz (n, 1.0 / hh, tmp, 1, & v [(j + 1) * n], 1);
      /* rotate */
      for (i=0; i<j; i++)
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

      *err = fabs (g [j + 1]); /* residual */
      /* if satisfied, */
      if (*err <= tol)
	{
	  j ++; /* this is because ++(*iter) in gmres(m) */
	  break;
	}
    }
  *iter = j;

  /* 3. form the approximate solution */
  /* solve y_k */
  back_sub (j, m, h, g, c); /* use c [] as y_k */
  /* x_m */
  for (i=0; i<n; i++)
    for (k=0; k<j; k++)
      x [i] += v [k * n + i] * c [k];

  free (tmp);
  free (v);
  free (h);
  free (g);
  free (c);
  free (s);
}

void
back_sub (int m, int nn,
	  double *r, double *g, double *y)
/* m  : number of iteration */
/* nn : dimension of matrix r [] (nnxnn) */
{
  int i, j, jj;

  /*for (j=m-1;j>=0;j--){*/
  /* above for-loop fail, because j is unsigned!! */
  for (jj=0; jj<m; jj++)
    {
      j = m - 1 - jj;
      y [j] = 0.0;
      for (i= j + 1; i<m; i++)
	{
	  y [j] -= r [j * nn + i] * y [i];
	}
      y [j] += g [j];
      y [j] = y [j] / r [j * nn + j];
    }
}
