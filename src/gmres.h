/* header file of mygmres.c
 * Copyright (C) 1998-2001 Kengo Ichiki <ichiki@haloumi.tn.utwente.nl>
 * $Id: gmres.h,v 1.1 2001/10/13 12:10:31 ichiki Exp $
 */

extern int GMRES_it_max;
extern int GMRES_it_restart;
extern double GMRES_eps;

void
solve_iter_gmres (int n,
		  double *b, double *x,
		  void (*atimes) (int, double *, double *));
void
mygmres_m (int n, double *f, double *x,
	   int m, double tol, int itmax,
	   int *iter, double *err,
	   void (*myatimes) (int, double *, double *));
void
mygmres (int n, double *f, double *x,
	 double tol, int itmax,
	 int *iter, double *err,
	 void (*myatimes) (int, double *, double *));
