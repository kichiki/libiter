/* header file of mygmres.c
 * Copyright (C) 1998-2001 Kengo Ichiki <ichiki@haloumi.tn.utwente.nl>
 * $Id: gmres.h,v 2.1 2001/10/13 12:12:18 ichiki Exp $
 */

void
solve_iter_gmres (int n,
		  double *b, double *x,
		  void (*atimes) (int, double *, double *, void *),
		  void * user_data,
		  int it_max, int it_restart, double eps);
void
mygmres_m (int n, double *f, double *x,
	   int m, double tol, int itmax,
	   int *iter, double *err,
	   void (*myatimes) (int, double *, double *, void *),
	   void * user_data);
void
mygmres (int n, double *f, double *x,
	 double tol, int itmax,
	 int *iter, double *err,
	 void (*myatimes) (int, double *, double *, void *),
	 void * user_data);
