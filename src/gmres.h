/* header file of mygmres.c --
 * generalized minimum residual method
 * Copyright (C) 1998-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: gmres.h,v 2.5 2006/09/26 17:08:52 ichiki Exp $
 *
 * Reference :
 *   GMRES(m) : Y.Saad & M.H.Schultz, SIAM J.Sci.Stat.Comput.
 *   vol7 (1986) pp.856-869
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

/** global variables **/
extern int ITER_gmres_debug; /* [0|1]: [not print/print] iter and res */


/* wrapper routine for gmres_m ()
 * INPUT
 *   n : size of vectors v[] and f[] -- expected to be np * nelm for red-sym
 *   b [n] : given vector
 *   atimes (n, x, b) : routine to calc A.x and return b[]
 *   it_max : max # iterations
 *   it_restart : # iterations to restart
 *   eps : the accuracy
 * OUTPUT
 *   x [n] : solution
 */
void
solve_iter_gmres (int n,
		  const double *b, double *x,
		  void (*atimes) (int, const double *, double *, void *),
		  void * user_data,
		  int it_max, int it_restart, double eps);
void
gmres_m (int n, const double *f, double *x,
	 int m, double tol, int itmax,
	 int *iter, double *res,
	 void (*myatimes) (int, const double *, double *, void *),
	 void * user_data);
void
gmres (int n, const double *f, double *x,
       double tol, int itmax,
       int *iter, double *res,
       void (*myatimes) (int, const double *, double *, void *),
       void * user_data);
