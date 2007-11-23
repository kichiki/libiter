/* header file for library 'iter' -- gmres.c, bi-cgstab.c, and orthomin.c.
 * Copyright (C) 1999-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: libiter.h,v 2.11 2007/11/23 04:57:22 kichiki Exp $
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

#include <stdio.h> // FILE

struct iter {
  char *solver; /* for solve_iter() routine.
		 * currently following solvers are implemented:
		 * "steepest" : steepest descent method
		 * "cg"       : conjugate gradient
		 * "cgs"      : conjugate gradient squared
		 * "bicgstab" : bi-conjugate gradient stabilized
		 * "sta", "sta2", "gpb", "otmk" :
		 * "gmres"    : generalized minimum residual method  (default)
		 */
  int max;
  int restart;
  double eps;
  //double log10_eps;

  /* to keep good initial guess for the next step */
  int n;
  double *guess;
  int flag_guess; /* 0 : don't keep the results for the guess
		   * 1 : keep the results for the next guess
		   */
  /* results of the iteration */
  int niter;   // number of iteration
  double res2; // |r^2| / |b^2|

  /* to report the debug informations */
  FILE *out;
  int debug; /* = 0 : no debug info
	      * = 1 : iteration numbs and residue
	      * = 2 : residue for each iteration step
	      */
};

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
	   int debug, FILE *out);

void
iter_free (struct iter *param);

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
	    struct iter *it_param);


#endif /* !_LIBITER_H_ */
