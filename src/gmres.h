/* header file of mygmres.c --
 * generalized minimum residual method
 * Copyright (C) 1998-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: gmres.h,v 2.4 2006/09/26 05:22:20 ichiki Exp $
 *
 * Reference :
 *   GMRES(m) : Y.Saad & M.H.Schultz, SIAM J.Sci.Stat.Comput.
 *   vol7 (1986) pp.856-869
 */

/** global variables **/
extern int ITER_mygmres_debug; /* [0|1]: [not print/print] iter and res */


void
solve_iter_gmres (int n,
		  const double *b, double *x,
		  void (*atimes) (int, const double *, double *, void *),
		  void * user_data,
		  int it_max, int it_restart, double eps);
void
mygmres_m (int n, const double *f, double *x,
	   int m, double tol, int itmax,
	   int *iter, double *res,
	   void (*myatimes) (int, const double *, double *, void *),
	   void * user_data);
void
mygmres (int n, const double *f, double *x,
	 double tol, int itmax,
	 int *iter, double *res,
	 void (*myatimes) (int, const double *, double *, void *),
	 void * user_data);
