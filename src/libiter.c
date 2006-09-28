/* overall wrapper for iterative solver routines
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: libiter.c,v 1.1 2006/09/28 04:26:50 kichiki Exp $
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
#include <math.h> /* log10() */
#include <stdlib.h> /* malloc(), free() */
#include <string.h> /* strcmp() */

#include "gmres.h"
#include "bi-cgstab.h"
#include "orthomin.h"

#include "libiter.h"


/* initialize parameters
 *   solver : string indicating the solver
 *            sta, sta2, gpb, otmk, or gmres (default)
 *   eps and log10_eps
 *   max (and restart)
 */
struct iter *
iter_init (const char * solver,
	   int max,
	   int restart,
	   double eps,
	   int debug)
{
  struct iter * param = NULL;

  param = (struct iter *) malloc (sizeof (struct iter));

  param->solver = (char *) malloc (sizeof (char) * (strlen (solver) + 1));
  strcpy (param->solver, solver);

  param->max = max;
  param->restart = restart;
  param->eps = eps;
  param->log10_eps = log10 (eps);
  param->debug = debug;

  return (param);
}

void
iter_free (struct iter * param)
{
  if (param != NULL)
    {
      if (param->solver != NULL)
	{
	  free (param->solver);
	}
      free (param);
    }
}

/* wrapper routine for iterative solvers
 * INPUT
 *   n : size of vectors v[] and f[] -- expected to be np * nelm for red-sym
 *   b [n] : given vector
 *   atimes (n, x, b, user_data) : routine to calc A.x and return b[]
 *   user_data : pointer to be passed to solver and atimes routines
 *   it_param : parameters for iterative solvers
 *              solver : string indicating the solver
 *                       sta, sta2, gpb, otmk, or gmres (default)
 *              eps and log10_eps
 *              max (and restart)
 * OUTPUT
 *   x [n] : solution
 */
void
solve_iter (int n, const double *b,
	    double *x,
	    void (*atimes) (int, const double *, double *, void *),
	    void * user_data,
	    struct iter * it_param)
{
  int i;

  double hnor;
  double residual;
  int iter;


  /* some preparation */
  hnor = 0.0;
  for (i = 0; i < n; ++i)
    {
      hnor += b [i] * b [i];
    }
  if (hnor != 0.0)
    {
      hnor = log10 (hnor) / 2.0;
    }


  if (strcmp (it_param->solver, "sta") == 0)
    {
      sta (n, b, x,
	   it_param->max,
	   it_param->log10_eps,
	   hnor, &iter, &residual,
	   atimes, user_data);
    }
  else if (strcmp (it_param->solver, "sta2") == 0)
    {
      sta2 (n, b, x,
	    it_param->max,
	    it_param->log10_eps,
	    hnor, &iter, &residual,
	    atimes, user_data);
    }
  else if (strcmp (it_param->solver, "gpb") == 0)
    {
      gpb (n, b, x,
	   it_param->max,
	   it_param->log10_eps,
	   hnor, &iter, &residual,
	   atimes, user_data);
    }
  else if (strcmp (it_param->solver, "otmk") == 0)
    {
      otmk (n, b, x,
	    it_param->restart,
	    it_param->max,
	    it_param->log10_eps,
	    hnor, &iter, &residual,
	    atimes, user_data);
    }
  else //if (strcmp (it_param->solver, "gmres") == 0)
    {
      gmres_m (n, b, x,
	       it_param->restart,
	       it_param->eps,
	       it_param->max,
	       &iter, &residual,
	       atimes, user_data);
    }

  if (it_param->debug != 0)
    {
      fprintf (stderr, "libiter: iter=%d res=%e\n", iter, residual);
    }
}
