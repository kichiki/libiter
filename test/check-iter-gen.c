/* test code for iterative schemes in general
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-iter-gen.c,v 1.1 2007/11/23 04:43:15 kichiki Exp $
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
#include <string.h> // strcmp()
#include <math.h>   // sqrt()
#include <bench.h> // ptime_ms_d()
#include "check.h" // compare()
#include "memory-check.h" // macro CHECK_MALLOC

#include <dgetri_c.h> // lapack_inv_()

#include "libiter.h"

#include "steepest.h" // steepest()
#include "cg.h"       // cg()
#include "cgs.h"      // cgs()

#include "gmres.h"    // 

#include "bi-cgstab.h"// sta(), sta2(), gpb()
#include "orthomin.h" // otmk()

// schemes require atimes_t()
// the following is implemented by my BLAS.
#include "atpres.h"
#include "cgne.h"
// the following is implemented by fortran BLAS
#include "qmr.h"
#include "bico.h"
#include "bicg.h"


/*
#include "bicgstab.h"
*/

/* BLAS functions */
void dgemv_(char *trans, int *m, int *n, double *alpha, double *a,
	    int *lda, double *x, int *incx,
	    double *beta, double *y, int *incy);
void dcopy_ (int *, double *, int *, double *, int *);


/*
 * INPUT
 *   x0[n]: initial guess
 *   x[n] : exact solution.
 *          if NULL is given, no checking
 *   solver : currently the following solver can handle:
 *            the following schemes need only atimes()
 *              "steepest"
 *              "cg"
 *              "cgs"
 *              "bicgstab"
 *              "sta"
 *              "sta2"
 *              "gpb"
 *              "otmk"
 *              "gmres"
 *            the following schemes also need atimes_t()
 *              "atpres"
 *              "cgne"
 *              "bicg"
 *              "bico"
 *              "qmr"
 *            the following schemes need the preconditioner inv()
 *              "sta_pc"
 *              "sta2_pc"
 *              "gpb_pc"
 *              "otmk_pc"
 *              "gmres_pc"
 *   itmax    :
 *   nrestart :
 *   eps      : tolerance for the iteration scheme
 */
int
check_iter_gen (int n, const double *b, const double *x0,
		void (*atimes)(int, const double *, double *, void *),
		void (*atimes_t)(int, const double *, double *, void *),
		void *atimes_param,
		void (*inv) (int, const double *, double *, void *),
		void *inv_param,
		double *x,
		const char *solver, int itmax, int nrestart, double eps,
		int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "--------------------------------------------------\n"
	       "check_iter_gen with %s (eps = %e): start\n",
	       solver, eps);
    }

  int check = 0;


  // set struct iter
  struct iter *it = iter_init (solver, itmax, nrestart, eps,
			       n, NULL, 0, // guess
			       0, NULL);
  CHECK_MALLOC (it, "check_iter_gen");

  double *x1 = (double *) malloc (sizeof (double) * n);
  CHECK_MALLOC (x1, "check_iter_gen");

  // initial guess
  int i;
  for (i = 0; i < n; i ++)
    {
      x1[i] = x0[i];
    }

  double t0 = 0.0;
  double t1 = 0.0;

  if (strcmp (solver, "steepest") == 0)
    {
      t0 = ptime_ms_d();
      steepest (n, b, x1, atimes, atimes_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "cg") == 0)
    {
      t0 = ptime_ms_d();
      cg (n, b, x1, atimes, atimes_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "cgs") == 0)
    {
      t0 = ptime_ms_d();
      cgs (n, b, x1, atimes, atimes_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "bicgstab") == 0)
    {
      fprintf (stderr, "not yet implemented\n");
    }
  else if (strcmp (solver, "sta") == 0)
    {
      t0 = ptime_ms_d();
      sta (n, b, x1, atimes, atimes_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "sta2") == 0)
    {
      t0 = ptime_ms_d();
      sta2 (n, b, x1, atimes, atimes_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "gpb") == 0)
    {
      t0 = ptime_ms_d();
      gpb (n, b, x1, atimes, atimes_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "otmk") == 0)
    {
      t0 = ptime_ms_d();
      otmk (n, b, x1, atimes, atimes_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "gmres") == 0)
    {
      t0 = ptime_ms_d();
      gmres_m (n, b, x1, atimes, atimes_param, it);
      //gmres_m_ (n, b, x1, atimes, atimes_param, it);
      t1 = ptime_ms_d();
    }
  /**
   * the following scheme needs atimes_t()
   */
  else if (strcmp (solver, "atpres") == 0)
    {
      t0 = ptime_ms_d();
      atpres (n, b, x1, atimes, atimes_t, atimes_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "cgne") == 0)
    {
      t0 = ptime_ms_d();
      cgne (n, b, x1, atimes, atimes_t, atimes_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "bicg") == 0)
    {
      t0 = ptime_ms_d();
      bicg (n, b, x1, atimes, atimes_t, atimes_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "bico") == 0)
    {
      t0 = ptime_ms_d();
      bico (n, b, x1, atimes, atimes_t, atimes_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "qmr") == 0)
    {
      t0 = ptime_ms_d();
      qmr (n, b, x1, atimes, atimes_t, atimes_param, it);
      t1 = ptime_ms_d();
    }
  /**
   * with preconditioner
   */
  else if (strcmp (solver, "sta_pc") == 0)
    {
      t0 = ptime_ms_d();
      sta_pc (n, b, x1, atimes, atimes_param, inv, inv_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "sta2_pc") == 0)
    {
      t0 = ptime_ms_d();
      sta2_pc (n, b, x1, atimes, atimes_param, inv, inv_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "gpb_pc") == 0)
    {
      t0 = ptime_ms_d();
      gpb_pc (n, b, x1, atimes, atimes_param, inv, inv_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "otmk_pc") == 0)
    {
      t0 = ptime_ms_d();
      otmk_pc (n, b, x1, atimes, atimes_param, inv, inv_param, it);
      t1 = ptime_ms_d();
    }
  else if (strcmp (solver, "gmres_pc") == 0)
    {
      t0 = ptime_ms_d();
      gmres_m_pc (n, b, x1, atimes, atimes_param, inv, inv_param, it);
      t1 = ptime_ms_d();
    }
  else
    {
      fprintf (stderr, "invalid solver %s\n", solver);
    }

  iter_free (it);

  double res = sqrt (it->res2);
  if (res >= eps)
    {
      if (verbose != 0)
	{
	  fprintf (stdout, "couldn't converge : res = %e, eps = %e\n",
		   res, eps);
	}
      check ++;
    }
  else
    {
      if (verbose != 0)
	{
	  fprintf (stdout, "converged : res = %e, eps = %e\n",
		   res, eps);
	}
      char label[80];
      for (i = 0; i < n; i ++)
	{
	  sprintf (label, "x[%d]", i);
	  check += compare (x1[i], x[i], label, verbose, tiny);
	}
    }

  free (x1);

  if (verbose != 0)
    {
      fprintf (stdout, "niter : %d\n", it->niter);
      fprintf (stdout, "(|r^2| / |b^2|)^{1/2} = %e\n", res);
      fprintf (stdout, "CPU time : %f\n\n", t1 - t0);
      if (check == 0)
	fprintf (stdout, " => PASSED\n\n");
      else
	fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
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
static void
atimes_t_by_matrix (int n, const double *x, double *b, void *user_data)
{
  char trans = 'N'; /* fortran's memory allocation is transposed */
  int i_1 = 1;
  double d_1 = 1.0;
  double d_0 = 0.0;

  double *mat = (double *)user_data;
  dgemv_ (&trans, &n, &n, &d_1, mat, &n,
	  x, &i_1,
	  &d_0, b, &i_1);
}

/* check small problem for all schemes
 * the problem is
 *  [1 2 3]^{-1}  [-3/4  1/2   1/4 ]
 *  [2 1 4]     = [ 1/2 -2/5   1/10]
 *  [3 4 1]       [ 1/4  1/10 -3/20]
 * thereofre,
 *  the solution of A.x = b for b = [4 10 20] is
 *  x = [-3+5+5, 2-4+2, 1+1-3] = [7, 0, -1]
 */
int
check_3_all (int it_max, int it_restart, double it_eps,
	     int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_3_gen (tiny = %e): start\n\n", tiny);
    }

  int check = 0;


  int n = 3;
  double *a = (double *)malloc (sizeof (double) * n * n);
  double *b = (double *)malloc (sizeof (double) * n);
  double *x = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (a, "check_3_all");
  CHECK_MALLOC (b, "check_3_all");
  CHECK_MALLOC (x, "check_3_all");

  a[0] = 1.0;
  a[1] = 2.0;
  a[2] = 3.0;

  a[3] = 2.0;
  a[4] = 1.0;
  a[5] = 4.0;

  a[6] = 3.0;
  a[7] = 4.0;
  a[8] = 1.0;

  b[0] = 4.0;
  b[1] = 10.0;
  b[2] = 20.0;

  x[0] = 7.0;
  x[1] = 0.0;
  x[2] = -1.0;

  /* steepest descent method has some problem...
  // Steepest
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "steepest", it_max, it_restart, it_eps,
		    verbose, tiny);
  */

  // CG
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "cg", it_max, it_restart, it_eps,
		    verbose, tiny);
  // CGS
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "cgs", it_max, it_restart, it_eps,
		    verbose, tiny);

  // Bi-CGSTAB
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "sta", it_max, it_restart, it_eps,
		    verbose, tiny);
  // Bi-CGSTAB2
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "sta2", it_max, it_restart, it_eps,
		    verbose, tiny);

  // GPBi-CG
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "gpb", it_max, it_restart, it_eps,
		    verbose, tiny);

  // ORTHOMIN
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "otmk", it_max, it_restart, it_eps,
		    verbose, tiny);

  // GMRES
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL,
		    (void *)a,
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
		    atimes_by_matrix, atimes_t_by_matrix,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "atpres", it_max, it_restart, it_eps,
		    verbose, tiny);
  // CGNE
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, atimes_t_by_matrix,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "cgne", it_max, it_restart, it_eps,
		    verbose, tiny);

  // BICG
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, atimes_t_by_matrix,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "bicg", it_max, it_restart, it_eps,
		    verbose, tiny);
  // BICO
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, atimes_t_by_matrix,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "bico", it_max, it_restart, it_eps,
		    verbose, tiny);
  // QMR
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, atimes_t_by_matrix,
		    (void *)a,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "qmr", it_max, it_restart, it_eps,
		    verbose, tiny);

  free (a);
  free (b);
  free (x);

  if (verbose != 0)
    {
      fprintf (stdout,
	       "--------------------------------------------------\n"
	       "check_3_gen (tiny = %e): finished\n", tiny);

      if (check == 0)
	fprintf (stdout, " => PASSED\n");
      else
	fprintf (stdout, " => FAILED\n");

      fprintf (stdout,
	       "==================================================\n\n");
    }

  return (check);
}

int
check_symmetric_all (int n, int it_max, int it_restart, double it_eps,
		     int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_symmetric_all n = %d (eps = %e): start\n",
	       n, tiny);
    }

  int check = 0;


  double *b = (double *)malloc (sizeof (double) * n);
  double *x = (double *)malloc (sizeof (double) * n);
  double *mat = (double *)malloc (sizeof (double) * n * n);
  double *inv = (double *)malloc (sizeof (double) * n * n);

  int i, j;
  for (i = 0; i < n; i ++)
    {
      mat[i * n + i] = 1.0;
      for (j = i+1; j < n; j ++)
	{
	  mat[i * n + j] = 1.0/(double)(i+j);
	  mat[j * n + i] = 1.0/(double)(i+j);
	}
    }
  int i_1 = 1;
  int nn = n * n;
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


  /*
  // Steepest
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "steepest", it_max, it_restart, it_eps,
		    verbose, tiny);
  */

  // CG
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "cg", it_max, it_restart, it_eps,
		    verbose, tiny);

  // CGS
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "cgs", it_max, it_restart, it_eps,
		    verbose, tiny);

  // Bi-CGSTAB
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "sta", it_max, it_restart, it_eps,
		    verbose, tiny);

  // Bi-CGSTAB2
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "sta2", it_max, it_restart, it_eps,
		    verbose, tiny);

  // GPBi-CG
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "gpb", it_max, it_restart, it_eps,
		    verbose, tiny);

  // ORTHOMIN
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "otmk", it_max, it_restart, it_eps,
		    verbose, tiny);

  // GMRES
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, NULL, (void *)mat,
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
		    atimes_by_matrix, atimes_t_by_matrix, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "atpres", it_max, it_restart, it_eps,
		    verbose, tiny);
  // CGNE
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, atimes_t_by_matrix, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "cgne", it_max, it_restart, it_eps,
		    verbose, tiny);
  // BICG
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, atimes_t_by_matrix, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "bicg", it_max, it_restart, it_eps,
		    verbose, tiny);
  // BICO
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, atimes_t_by_matrix, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "bico", it_max, it_restart, it_eps,
		    verbose, tiny);
  // QMR
  check +=
    check_iter_gen (n, b, b, // initial guess
		    atimes_by_matrix, atimes_t_by_matrix, (void *)mat,
		    NULL, NULL, // for preconditioner
		    x, // exact solution
		    "qmr", it_max, it_restart, it_eps,
		    verbose, tiny);

  free (b);
  free (x);
  free (mat);

  if (verbose != 0)
    {
      fprintf (stdout,
	       "--------------------------------------------------\n"
	       "check_symmetric_all n = %d (eps = %e): finished\n",
	       n, tiny);

      if (check == 0)
	fprintf (stdout, " => PASSED\n");
      else
	fprintf (stdout, " => FAILED\n");

      fprintf (stdout,
	       "==================================================\n\n");
    }

  return (check);
}
