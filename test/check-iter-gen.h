/* header file for check-gen.c --
 * test code for iterative schemes in general
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-iter-gen.h,v 1.1 2007/11/23 04:43:41 kichiki Exp $
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
#ifndef	_CHECK_GEN_H_
#define	_CHECK_GEN_H_


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
		int verbose, double tiny);

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
	     int verbose, double tiny);

int
check_symmetric_all (int n, int it_max, int it_restart, double it_eps,
		     int verbose, double tiny);


#endif /* !_CHECK_GEN_H_ */
