/* Toeplitz matrix related code to check libiter solvers
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: toeplitz.c,v 1.2 2007/11/25 19:11:23 kichiki Exp $
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
#include <dgetri_c.h>
#include "memory-check.h" // macro CHECK_MALLOC

#include "check-iter-gen.h" // check_iter_gen()


/* BLAS functions */
void dgemv_(char *trans, int *m, int *n, double *alpha, double *a,
	    int *lda, double *x, int *incx,
	    double *beta, double *y, int *incy);


/* make Toeplitz matrix
 * INPUT
 *  gamma
 * OUTPUT
 *  mat[n*n]
 */
void
Toeplitz_make_matrix (int n, double gamma, double *mat)
{
  /* Toeplitz matrix */
  int i;
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

/* calc atimes b = A.x, where A is the Toeplitz matrix
 * INPUT
 *  x[n] : given vector
 *  *user_data = (double) gamma
 * OUTPUT
 *  b[n] := A.x
 */
void
Toeplitz_atimes (int n, const double *x, double *b, void *user_data)
{
  double *gamma = (double *)user_data;

  /* Toeplitz matrix */
  int i;
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

/* calc atimes b = A^T.x, where A is the Toeplitz matrix
 * INPUT
 *  x[n] : given vector
 *  *user_data = (double) gamma
 * OUTPUT
 *  b[n] := A^T.x
 */
void
Toeplitz_atimes_t (int n, const double *x, double *b, void *user_data)
{
  double *gamma = (double *)user_data;

  /* Toeplitz matrix */
  int i;
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
static void
atimes_by_matrix (int n, const double *x, double *b, void *user_data)
{
  char trans = 'T'; /* fortran's memory allocation is transposed */
  int i_1 = 1;
  double d_1 = 1.0;
  double d_0 = 0.0;

  double *mat = (double *)user_data;
  dgemv_ (&trans, &n, &n, &d_1, mat, &n,
	  x, &i_1,
	  &d_0, b, &i_1);
}

int
Toeplitz_check_all (int n, double gamma,
		    int it_max, int it_restart, double it_eps,
		    int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "Toeplitz_check_all n = %d gamma = %e (eps = %e): start\n",
	       n, gamma, tiny);
    }

  int check = 0;


  double *b = (double *)malloc (sizeof (double) * n);
  double *x = (double *)malloc (sizeof (double) * n);
  double *mat = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (b,   "Toeplitz_check_all");
  CHECK_MALLOC (x,   "Toeplitz_check_all");
  CHECK_MALLOC (mat, "Toeplitz_check_all");

  /* given vector */
  int i;
  for (i = 0; i < n; i ++)
    {
      b[i] = 1.0;
    }

  Toeplitz_make_matrix (n, gamma, mat);
  lapack_inv_ (n, mat);
  atimes_by_matrix (n, b, x, (void *)mat);
  // x[] is the solution
  free (mat);

  /* the following schemes failed...
   * Toeplitz problem is hard for them.
  // Steepest
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, NULL, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "steepest", it_max, it_restart, it_eps,
		    verbose, tiny);
  // CG
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, NULL, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "cg", it_max, it_restart, it_eps,
		    verbose, tiny);
  // CG, another implementation
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, NULL, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "cg_", it_max, it_restart, it_eps,
		    verbose, tiny);
  // CGS
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, NULL, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "cgs", it_max, it_restart, it_eps,
		    verbose, tiny);
  */

  /* Bi-CGSTAB may failed in high tolerance -- 
   * first converging up to a certain point,
   * then start diverging,
   * and after that, start converging again,
   * but the converged value has enormous errors
   * although the residual looks small....

  // Bi-CGSTAB (in Weiss)
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, NULL, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "bicgstab", it_max, it_restart, it_eps,
		    verbose, tiny);
  // Bi-CGSTAB
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, NULL, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "sta", it_max, it_restart, it_eps,
		    verbose, tiny);
  */
  // Bi-CGSTAB2
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, NULL, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "sta2", it_max, it_restart, it_eps,
		    verbose, tiny);

  // GPBi-CG
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, NULL, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "gpb", it_max, it_restart, it_eps,
		    verbose, tiny);

  // ORTHOMIN
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, NULL, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "otmk", it_max, it_restart, it_eps,
		    verbose, tiny);

  // GMRES
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, NULL, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "gmres", it_max, it_restart, it_eps,
		    verbose, tiny);

  /**
   * schemes with atimes_t()
   */
  // ATPRES
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, Toeplitz_atimes_t, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "atpres", it_max, it_restart, it_eps,
		    verbose, tiny);
  // CGNE
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, Toeplitz_atimes_t, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "cgne", it_max, it_restart, it_eps,
		    verbose, tiny);
  /* the following schemes failed...
   * maybe, Toeplitz problem is hard for them.
  // BICG
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, Toeplitz_atimes_t, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "bicg", it_max, it_restart, it_eps,
		    verbose, tiny);
  // BICO
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, Toeplitz_atimes_t, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "bico", it_max, it_restart, it_eps,
		    verbose, tiny);
  // QMR
  check +=
    check_iter_gen (n, b, b, // initial guess
		    Toeplitz_atimes, Toeplitz_atimes_t, (void *) &gamma,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "QMR", it_max, it_restart, it_eps,
		    verbose, tiny);
  */

  free (b);
  free (x);

  if (verbose != 0)
    {
      fprintf (stdout,
	       "--------------------------------------------------\n"
	       "Toeplitz_check_all n = %d gamma = %e (eps = %e): finished\n",
	       n, gamma, tiny);

      if (check == 0)
	fprintf (stdout, " => PASSED\n");
      else
	fprintf (stdout, " => FAILED\n");

      fprintf (stdout,
	       "==================================================\n\n");
    }

  return (check);
}
