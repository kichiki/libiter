/* test code for otmk_pc() in orthomin.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-otmk-pc.c,v 1.1 2007/12/01 18:07:01 kichiki Exp $
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
#include <math.h>  // sqrt()
#include <bench.h> // ptime_ms_d()
#include "check.h" // compare()
#include "memory-check.h" // macro CHECK_MALLOC

#include "libiter.h"
#include "toeplitz.h" // Toeplitz_make_matrix()

#include "ilu.h" // inv_diag(), inv_ILU()

#include "orthomin.h"


static void
local_atimes (int n, const double *x,
	      double *b, void *param)
{
  double *a = (double *)param;

  int i, j;
  for (i = 0; i < n; i ++)
    {
      b[i] = 0.0;
      for (j = 0; j < n; j ++)
	{
	  b[i] += a[i*n+j] * x[j];
	}
    }
}

int
check_otmk_pc (int n,
	       int it_max, int it_restart, double it_eps,
	       int verbose, double tiny)
{
  if (verbose != 0)
    {
      fprintf (stdout,
	       "==================================================\n"
	       "check_otmk_pc(n=%d): start\n", n);
    }

  int check = 0;
  double max = 0.0;

  double *a  = (double *)malloc (sizeof (double) * n * n);
  double *b  = (double *)malloc (sizeof (double) * n);
  double *x  = (double *)malloc (sizeof (double) * n);
  double *x_ = (double *)malloc (sizeof (double) * n);
  double *lu = (double *)malloc (sizeof (double) * n * n);
  CHECK_MALLOC (a,  "check_otmk_pc");
  CHECK_MALLOC (b,  "check_otmk_pc");
  CHECK_MALLOC (x,  "check_otmk_pc");
  CHECK_MALLOC (x_, "check_otmk_pc");
  CHECK_MALLOC (lu, "check_otmk_pc");


  Toeplitz_make_matrix (n, 1.5, a);
  int i;
  srand48(0);
  for (i = 0; i < n; i ++)
    {
      b[i] = drand48();
    }

  struct iter *it
    = iter_init ("otmk", it_max, it_restart, it_eps, // solver param
		 n, NULL, 0, // guess
		 0, NULL); // debug info
  CHECK_MALLOC (it, "check_otmk_pc");

  char label[80];


  // by solve_iter
  for (i = 0; i < n; i ++)
    {
      x[i] = 0.0;
    }
  double t0 = ptime_ms_d ();
  solve_iter (n, b, x, local_atimes, (void *)a, it);
  double t1 = ptime_ms_d ();
  if (verbose != 0)
    {
      fprintf (stdout, "solve_iter(otmk) "
	       "iter = %d, res = %e, CPU time = %f\n",
	       it->niter, sqrt(it->res2), t1 - t0);
    }

  // otmk
  for (i = 0; i < n; i ++)
    {
      x_[i] = 0.0;
    }
  t0 = ptime_ms_d ();
  otmk (n, b, x_,
	local_atimes, (void *)a,
	it);
  t1 = ptime_ms_d ();
  if (verbose != 0)
    {
      fprintf (stdout, "otmk             "
	       "iter = %d, res = %e, CPU time = %f\n",
	       it->niter, sqrt(it->res2), t1 - t0);
    }
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "x[%d]", i);
      check += compare_max (x[i], x_[i], label, verbose, tiny, &max);
    }

  // otmk_pc with inv_diag()
  for (i = 0; i < n; i ++)
    {
      x_[i] = 0.0;
    }
  t0 = ptime_ms_d ();
  otmk_pc (n, b, x_,
	   local_atimes, (void *)a,
	   inv_diag, (void *)a,
	   it);
  t1 = ptime_ms_d ();
  if (verbose != 0)
    {
      fprintf (stdout, "otmk_pc(diag)    "
	       "iter = %d, res = %e, CPU time = %f\n",
	       it->niter, sqrt(it->res2), t1 - t0);
    }
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "x[%d]", i);
      check += compare_max (x[i], x_[i], label, verbose, tiny, &max);
    }

  // otmk_pc with inv_ILU()
  for (i = 0; i < n; i ++)
    {
      x_[i] = 0.0;
    }
  LU_Gauss (n, a, lu);
  t0 = ptime_ms_d ();
  otmk_pc (n, b, x_,
	   local_atimes, (void *)a,
	   inv_ILU, (void *)lu,
	   it);
  t1 = ptime_ms_d ();
  if (verbose != 0)
    {
      fprintf (stdout, "otmk_pc(LU)      "
	       "iter = %d, res = %e, CPU time = %f\n",
	       it->niter, sqrt(it->res2), t1 - t0);
    }
  for (i = 0; i < n; i ++)
    {
      sprintf (label, "x[%d]", i);
      check += compare_max (x[i], x_[i], label, verbose, tiny, &max);
    }

  iter_free (it);

  free (a);
  free (b);
  free (x);
  free (x_);
  free (lu);

  if (verbose != 0)
    {
      fprintf (stdout, " max error = %e vs tiny = %e\n", max, tiny);
      if (check == 0) fprintf (stdout, " => PASSED\n\n");
      else            fprintf (stdout, " => FAILED\n\n");
    }

  return (check);
}
