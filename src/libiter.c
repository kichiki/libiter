/* overall wrapper for iterative solver routines
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: libiter.c,v 1.6 2007/11/23 04:56:55 kichiki Exp $
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
#include <stdio.h> /* fprintf() */
#include <stdlib.h> /* malloc(), free() */
#include <string.h> /* strcmp() */

#include "libiter.h"

#include "steepest.h"
#include "cg.h"
#include "cgs.h"
#include "bicgstab.h"

#include "gmres.h"
#include "bi-cgstab.h"
#include "orthomin.h"

#include "memory-check.h" // CHECK_MALLOC


/* initialize parameters
 * INPUT
 *   solver : string indicating the solver
 *            sta, sta2, gpb, otmk, or gmres (default)
 *   max, restart, eps : iteration parameters
 *   n          : dimension of the problem
 *   guess[n]   : initial guess
 *                if NULL is given, set zero for the guess
 *   flag_guess = 0 : don't use the guess[] in solve_iter()
 *              = 1 : use guess[] for the initial guess in solve_iter()
 *   debug = 0 : no debug info
 *         = 1 : iteration numbs and residue
 *         = 2 : residue for each iteration step
 *   out   : FILE * to output debug info.
 * OUTPUT
 *   (struct iter *) returned value :
 */
struct iter *
iter_init (const char *solver,
	   int max, int restart, double eps,
	   int n, const double *guess, int flag_guess,
	   int debug, FILE *out)
{
  struct iter *param = (struct iter *)malloc (sizeof (struct iter));
  CHECK_MALLOC (param, "iter_init");

  param->solver = (char *)malloc (sizeof (char) * (strlen (solver) + 1));
  CHECK_MALLOC (param->solver, "iter_init");
  strcpy (param->solver, solver);

  param->max = max;
  param->restart = restart;
  param->eps = eps;

  /* to keep good initial guess for the next step */
  param->n = n;
  param->guess = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (param->solver, "iter_init");
  int i;
  for (i = 0; i < n; i ++)
    {
      if (guess == NULL)
	{
	  param->guess[i] = 0.0;
	}
      else
	{
	  param->guess[i] = guess[i];
	}
    }
  param->flag_guess = flag_guess;

  /* results of the iteration */
  param->niter = 0;
  param->res2  = 0.0;

  /* to report the debug informations */
  param->out   = out;
  param->debug = debug;

  return (param);
}

void
iter_free (struct iter *param)
{
  if (param == NULL) return;

  if (param->solver != NULL) free (param->solver);
  if (param->guess  != NULL) free (param->guess);
  free (param);
}

/* wrapper routine for iterative solvers
 * INPUT
 *   n : size of vectors v[] and f[] -- expected to be np * nelm for red-sym
 *   b [n] : given vector
 *   atimes (n, x, b, user_data) : routine to calc A.x and return b[]
 *   user_data : pointer to be passed to solver and atimes routines
 *   it_param : parameters for iterative solvers
 *              solver : string indicating the solver
 *                "steepest" : steepest descent method
 *                "cg"       : conjugate gradient
 *                "cgs"      : conjugate gradient squared
 *                "bicgstab" : bi-conjugate gradient stabilized
 *                "sta", "sta2", "gpb", "otmk" :
 *                "gmres"    : generalized minimum residual method  (default)
 *              max, restart, eps
 *              n, guess[n] : the result at the last process
 *              flag_guess : 0 == don't keep the results,
 *                           1 == keep the results for the next.
 * OUTPUT
 *   x [n] : solution
 */
void
solve_iter (int n, const double *b,
	    double *x,
	    void (*atimes) (int, const double *, double *, void *),
	    void *atimes_param,
	    struct iter *it_param)
{
  int i;
  if (it_param->n != n)
    {
      fprintf (stderr, "libiter solve_iter : n is different %d != %d\n",
	       it_param->n, n);
      exit (1);
    }

  // set the initial guess
  if (it_param->flag_guess == 0)
    {
      for (i = 0; i < n; ++i)
	{
	  x[i] = 0.0;
	}
    }
  else
    {
      for (i = 0; i < n; ++i)
	{
	  x[i] = it_param->guess[i];
	}
    }


  if (strcmp (it_param->solver, "steepest") == 0)
    {
      steepest
	(n, b, x,
	 atimes, atimes_param,
	 it_param);
    }
  else if (strcmp (it_param->solver, "cg") == 0)
    {
      cg (n, b, x,
	  atimes, atimes_param,
	  it_param);
    }
  else if (strcmp (it_param->solver, "cgs") == 0)
    {
      cgs (n, b, x,
	   atimes, atimes_param,
	   it_param);
    }
  else if (strcmp (it_param->solver, "bicgstab") == 0)
    {
      bicgstab (n, b, x,
		atimes, atimes_param,
		it_param);
    }
  else if (strcmp (it_param->solver, "sta") == 0)
    {
      sta (n, b, x,
	   atimes, atimes_param,
	   it_param);
    }
  else if (strcmp (it_param->solver, "sta2") == 0)
    {
      sta2 (n, b, x,
	    atimes, atimes_param,
	    it_param);
    }
  else if (strcmp (it_param->solver, "gpb") == 0)
    {
      gpb (n, b, x,
	   atimes, atimes_param,
	   it_param);
    }
  else if (strcmp (it_param->solver, "otmk") == 0)
    {
      otmk (n, b, x,
	    atimes, atimes_param,
	    it_param);
    }
  else //if (strcmp (it_param->solver, "gmres") == 0)
    {
      gmres_m (n, b, x,
	       atimes, atimes_param,
	       it_param);
    }

  /* keep the results for the next first guess */
  if (it_param->flag_guess != 0)
    {
      for (i = 0; i < n; ++i)
	{
	  it_param->guess[i] = x[i];
	}
    }
}
