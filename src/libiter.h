/* header file for library 'iter' -- gmres.c, bi-cgstab.c, and orthomin.c.
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: libiter.h,v 2.8 2006/10/10 18:14:01 ichiki Exp $
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
#ifndef	_LIBITER_H_
#define	_LIBITER_H_

struct iter {
  char * solver;
  int max;
  int restart;
  double eps;
  double log10_eps;

  int debug; /* = 0 : no debug info
	      * = 1 : iteration numbs and residue
	      * = 2 : residue for each iteration step
	      */
  FILE * out;
};

/* initialize parameters
 * INPUT
 *   solver : string indicating the solver
 *            sta, sta2, gpb, otmk, or gmres (default)
 *   eps and log10_eps
 *   max (and restart)
 *   debug = 0 : no debug info
 *         = 1 : iteration numbs and residue
 *   out   : FILE * to output debug info.
 * OUTPUT
 *   (struct iter *) returned value :
 */
struct iter *
iter_init (const char * solver,
	   int max,
	   int restart,
	   double eps,
	   int debug,
	   FILE * out);

void
iter_free (struct iter * param);

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
	    struct iter * it_param);

#endif /* !_LIBITER_H_ */
