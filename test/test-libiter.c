/* test code for libiter solvers
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: test-libiter.c,v 1.1 2006/10/10 19:50:51 ichiki Exp $
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

#include <libiter.h>
#include <dgetri_c.h>
#include <bench.h>

#include "../src/steepest.h"
#include "../src/cg.h"
#include "../src/bicgstab.h"
#include "../src/cgs.h"
#include "../src/qmr.h"
#include "../src/bico.h"
#include "../src/bicg.h"

#include "../src/atpres.h"
#include "../src/cgne.h"


/* BLAS functions */
void dgemv_(char *trans, int *m, int *n, double *alpha, double *a,
	    int *lda, double *x, int *incx,
	    double *beta, double *y, int *incy);
int
dcopy_(int* N,
       double* X, int* incX,
       double* Y, int* incY);


/*
 * INPUT
 *  gamma
 */
void
make_matrix_Toeplitz (int n, double gamma, double * mat)
{
  int i;

  /* Toeplitz matrix */
  for (i = 0; i < n * n; i ++)
    {
      mat[i] = 0.0;
    }

  for (i = 0; i < n; i ++)
    {
      mat[i * n + i] = 2.0;
    }
  for (i = 0; i < n-1; i ++)
    {
      mat[i * n + (i+1)] = 1.0;
    }
  for (i = 2; i < n; i ++)
    {
      mat[i * n + (i-2)] = gamma;
    }
}
/*
 * INPUT
 *  *user_data = (double) gamma
 */
void
atimes_Toeplitz (int n, const double * x, double * b, void * user_data)
{
  int i;
  double *gamma;

  /* Toeplitz matrix */
  gamma = (double *)user_data;
  for (i = 0; i < n; i ++)
    {
      b[i] = 2.0 * x[i];
    }
  for (i = 0; i < n-1; i ++)
    {
      b[i] += x[i+1];
    }
  for (i = 2; i < n; i ++)
    {
      b[i] += (*gamma) * x[i-2];
    }
}
void
atimes_t_Toeplitz (int n, const double * x, double * b, void * user_data)
{
  int i;
  double *gamma;

  /* Toeplitz matrix */
  gamma = (double *)user_data;
  for (i = 0; i < n; i ++)
    {
      b[i] = 2.0 * x[i];
    }
  for (i = 0; i < n-2; i ++)
    {
      b[i] += x[i+2];
    }
  for (i = 1; i < n; i ++)
    {
      b[i] += (*gamma) * x[i-1];
    }
}


/*
 * INPUT
 *  *user_data = (double *) mat
 */
void
atimes_by_matrix (int n, const double * x, double * b, void * user_data)
{
  double * mat;

  char trans = 'T'; /* fortran's memory allocation is transposed */
  int i_1 = 1;
  double d_1 = 1.0;
  double d_0 = 0.0;

  mat = (double *)user_data;
  dgemv_ (&trans, &n, &n, &d_1, mat, &n,
	  x, &i_1,
	  &d_0, b, &i_1);
}
void
atimes_t_by_matrix (int n, const double * x, double * b, void * user_data)
{
  double * mat;

  char trans = 'N'; /* fortran's memory allocation is transposed */
  int i_1 = 1;
  double d_1 = 1.0;
  double d_0 = 0.0;

  mat = (double *)user_data;
  dgemv_ (&trans, &n, &n, &d_1, mat, &n,
	  x, &i_1,
	  &d_0, b, &i_1);
}


/*
 * INPUT
 *   x0[n]: initial guess
 *   x[n] : exact solution.
 *          if NULL is given, no checking
 */
void
test_all (int n,
	  const double * b, const double * x0,
	  void (* atimes)(int, const double *, double *, void *),
	  void (* atimes_trans)(int, const double *, double *, void *),
	  void * user_data,
	  int itmax, double eps,
	  double * x)
{
  struct iter * iter; // for libiter

  int i;

  double * x1 = NULL;

  int it;
  double res;

  double t0, t;
  double tiny = 1.0e-6;


  x1 = (double *) malloc (sizeof (double) * n);


  // reset initial guess
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }
  iter = iter_init ("steepest", itmax, 20, eps, 1, stdout);
  t0 = ptime_ms_d();
  steepest (n, b, x1,
	    atimes, user_data,
	    iter);
  t = ptime_ms_d();
  fprintf (stdout, "SteepestDescent: cpu time %f\n", t-t0);
  if (x != NULL)
    {
      for (i = 0; i < n; i ++)
	{
	  if (fabs (x1[i] - x[i]) > tiny)
	    {
	      fprintf (stdout, "x[%d] : %e, %e %e\n",
		       i,
		       fabs (x1[i] - x[i]),
		       x[i], x1[i]);
	    }
	}
    }
  fprintf (stdout, "\n");

  // reset initial guess
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }
  iter = iter_init ("cg", itmax, 20, eps, 1, stdout);
  t0 = ptime_ms_d();
  cg (n, b, x1,
      atimes, user_data,
      iter);
  t = ptime_ms_d();
  fprintf (stdout, "CG: cpu time %f\n", t-t0);
  if (x != NULL)
    {
      for (i = 0; i < n; i ++)
	{
	  if (fabs (x1[i] - x[i]) > tiny)
	    {
	      fprintf (stdout, "x[%d] : %e, %e %e\n",
		       i,
		       fabs (x1[i] - x[i]),
		       x[i], x1[i]);
	    }
	}
    }
  fprintf (stdout, "\n");

  // reset initial guess
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }
  t0 = ptime_ms_d();
  atpres (n, b, x1, eps, itmax, &it, &res, atimes, atimes_trans, user_data);
  t = ptime_ms_d();
  fprintf (stdout, "ATPRES (cpu time %f): it = %d, res = %e\n",
	   t-t0, it, res);
  if (x != NULL)
    {
      for (i = 0; i < n; i ++)
	{
	  if (fabs (x1[i] - x[i]) > tiny)
	    {
	      fprintf (stdout, "x[%d] : %e, %e %e\n",
		       i,
		       fabs (x1[i] - x[i]),
		       x[i], x1[i]);
	    }
	}
    }
  fprintf (stdout, "\n");

  // reset initial guess
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }
  t0 = ptime_ms_d();
  cgne (n, b, x1, eps, itmax, &it, &res, atimes, atimes_trans, user_data);
  t = ptime_ms_d();
  fprintf (stdout, "CGNE (cpu time %f): it = %d, res = %e\n",
	   t-t0, it, res);
  if (x != NULL)
    {
      for (i = 0; i < n; i ++)
	{
	  if (fabs (x1[i] - x[i]) > tiny)
	    {
	      fprintf (stdout, "x[%d] : %e, %e %e\n",
		       i,
		       fabs (x1[i] - x[i]),
		       x[i], x1[i]);
	    }
	}
    }
  fprintf (stdout, "\n");

  // reset initial guess
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }
  t0 = ptime_ms_d();
  bicg (n, b, x1, eps, itmax, &it, &res,
	atimes,
	atimes_trans,
	user_data);
  t = ptime_ms_d();
  fprintf (stdout, "BiCG (cpu time %f): it = %d, res = %e\n",
	   t-t0, it, res);
  if (x != NULL)
    {
      for (i = 0; i < n; i ++)
	{
	  if (fabs (x1[i] - x[i]) > tiny)
	    {
	      fprintf (stdout, "x[%d] : %e, %e %e\n",
		       i,
		       fabs (x1[i] - x[i]),
		       x[i], x1[i]);
	    }
	}
    }
  fprintf (stdout, "\n");

  // reset initial guess
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }
  t0 = ptime_ms_d();
  bico (n, b, x1, eps, itmax, &it, &res,
	atimes,
	atimes_trans,
	user_data);
  t = ptime_ms_d();
  fprintf (stdout, "BICO (cpu time %f): it = %d, res = %e\n",
	   t-t0, it, res);
  if (x != NULL)
    {
      for (i = 0; i < n; i ++)
	{
	  if (fabs (x1[i] - x[i]) > tiny)
	    {
	      fprintf (stdout, "x[%d] : %e, %e %e\n",
		       i,
		       fabs (x1[i] - x[i]),
		       x[i], x1[i]);
	    }
	}
    }
  fprintf (stdout, "\n");

  // reset initial guess
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }
  t0 = ptime_ms_d();
  qmr (n, b, x1, eps, itmax, &it, &res,
       atimes,
       atimes_trans,
       user_data);
  t = ptime_ms_d();
  fprintf (stdout, "QMR (cpu time %f): it = %d, res = %e\n",
	   t-t0, it, res);
  if (x != NULL)
    {
      for (i = 0; i < n; i ++)
	{
	  if (fabs (x1[i] - x[i]) > tiny)
	    {
	      fprintf (stdout, "x[%d] : %e, %e %e\n",
		       i,
		       fabs (x1[i] - x[i]),
		       x[i], x1[i]);
	    }
	}
    }
  fprintf (stdout, "\n");

  // reset initial guess
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }
  iter = iter_init ("cgs", itmax, 20, eps, 1, stdout);
  t0 = ptime_ms_d();
  cgs (n, b, x1,
       atimes, user_data,
       iter);
  t = ptime_ms_d();
  fprintf (stdout, "CGS: cpu time %f\n", t-t0);
  if (x != NULL)
    {
      for (i = 0; i < n; i ++)
	{
	  if (fabs (x1[i] - x[i]) > tiny)
	    {
	      fprintf (stdout, "x[%d] : %e, %e %e\n",
		       i,
		       fabs (x1[i] - x[i]),
		       x[i], x1[i]);
	    }
	}
    }
  fprintf (stdout, "\n");

  // reset initial guess
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }
  iter = iter_init ("bicgstab", itmax, 20, eps, 1, stdout);
  t0 = ptime_ms_d();
  bicgstab (n, b, x1,
	    atimes, user_data,
	    iter);
  t = ptime_ms_d();
  fprintf (stdout, "BiCGSTAB: cpu time %f\n", t-t0);
  if (x != NULL)
    {
      for (i = 0; i < n; i ++)
	{
	  if (fabs (x1[i] - x[i]) > tiny)
	    {
	      fprintf (stdout, "x[%d] : %e, %e %e\n",
		       i,
		       fabs (x1[i] - x[i]),
		       x[i], x1[i]);
	    }
	}
    }
  fprintf (stdout, "\n");

  // reset initial guess
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }
  iter = iter_init ("sta",   itmax, 20, eps, 1, stdout);
  t0 = ptime_ms_d();
  solve_iter (n, b, x1,
	      atimes, user_data,
	      iter);
  t = ptime_ms_d();
  fprintf (stdout, "sta: cpu time %f\n", t-t0);
  iter_free (iter);
  if (x != NULL)
    {
      for (i = 0; i < n; i ++)
	{
	  if (fabs (x1[i] - x[i]) > tiny)
	    {
	      fprintf (stdout, "x[%d] : %e, %e %e\n",
		       i,
		       fabs (x1[i] - x[i]),
		       x[i], x1[i]);
	    }
	}
    }
  fprintf (stdout, "\n");

  // reset initial guess
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }
  iter = iter_init ("gmres", itmax, 20, eps, 1, stdout);
   t0 = ptime_ms_d();
  solve_iter (n, b, x1,
	      atimes, user_data,
	      iter);
  t = ptime_ms_d();
  fprintf (stdout, "GMRES: cpu time %f\n", t-t0);
  iter_free (iter);
  if (x != NULL)
    {
      for (i = 0; i < n; i ++)
	{
	  if (fabs (x1[i] - x[i]) > tiny)
	    {
	      fprintf (stdout, "x[%d] : %e, %e %e\n",
		       i,
		       fabs (x1[i] - x[i]),
		       x[i], x1[i]);
	    }
	}
    }
  fprintf (stdout, "\n");

  free (x1);
}



/* main program */
int
main (int argc, char** argv)
{
  double * b;  // given vector
  double * x;  // solution
  double * mat_Toeplitz = NULL;
  int n;
  int nn;
  double gamma;

  double eps;
  int itmax;

  int i;
  int j;

  int i_1 = 1;

  double * mat = NULL;
  double * inv = NULL;



  //itmax = 2000;
  itmax = 10000;
  eps = 1.0e-12;


  /* some symmetric matrix */
  fprintf (stdout, "*** some symmetric matrix ***\n");
  n = 10;
  nn = n * n;

  b  = (double *) malloc (sizeof (double) * n);
  x  = (double *) malloc (sizeof (double) * n);

  mat = (double *) malloc (sizeof (double) * nn);
  inv = (double *) malloc (sizeof (double) * nn);

  for (i = 0; i < n; i ++)
    {
      mat[i * n + i] = 10.0;
    }
  for (i = 0; i < n; i ++)
    {
      for (j = i+1; j < n; j ++)
	{
	  mat[i * n + j] = (double) (i+j);
	  mat[j * n + i] = (double) (i+j);
	}
    }
  dcopy_(&nn, mat, &i_1, inv, &i_1);
  lapack_inv_ (n, inv);

  /* given vector */
  for (i = 0; i < n; i ++)
    {
      b[i] = 1.0;
    }

  atimes_by_matrix (n, b, x, (void *) inv);
  free (inv);
  // x is the solution

  test_all (n, b, b/* initial guess*/,
	    atimes_by_matrix,
	    atimes_t_by_matrix,
	    (void *) mat,
	    itmax, eps,
	    NULL/* x */);

  free (b);
  free (x);
  free (mat);


  /* Toeplitz matrix */
  fprintf (stdout, "*** Toeplitz matrix ***\n");
  //n = 200;
  n = 1000;
  nn = n * n;

  gamma = 1.5;

  b  = (double *) malloc (sizeof (double) * n);
  x  = (double *) malloc (sizeof (double) * n);
  mat_Toeplitz = (double *) malloc (sizeof (double) * n * n);

  /* given vector */
  for (i = 0; i < n; i ++)
    {
      b[i] = 1.0;
    }

  make_matrix_Toeplitz (n, gamma, mat_Toeplitz);
  lapack_inv_ (n, mat_Toeplitz);
  atimes_by_matrix (n, b, x, (void *) mat_Toeplitz);
  free (mat_Toeplitz);
  // x is the solution

  test_all (n, b, b/* initial guess*/,
	    atimes_Toeplitz,
	    atimes_t_Toeplitz,
	    (void *) &gamma,
	    itmax, eps,
	    NULL/*x*/);

  free (b);
  free (x);

  return 0;
}
