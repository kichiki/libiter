/* header file of orthomin.c --
 * orthomin scheme
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: orthomin.h,v 2.3 2006/09/26 05:33:45 ichiki Exp $
 *
 * translated from fortran code
 *             ===============rog.f==============
 *             problem given by martin h. gutknecht
 *             numerical method: orthomin(k) method
 *             ver. 1.0 jul. 28 1995 by s. l. zhang
 *             ver. 1.1 aug. 31 1995 by s. l. zhang
 *    at slzhang.fort.iter.complex.orthomin.gutknecht-problem(gut.f)
 */

/** global variables **/
extern int ITER_otmk_debug; /* [0|1]: [not print/print] iter and res */

void
solve_iter_otmk (int n, const double *b,
		 double *x,
		 void (*atimes) (int, const double *, double *, void *),
		 void * user_data,
		 void (*solver) (int, const double *, double *,
				 int, int, double,
				 double, int *, double *,
				 void (*)
				 (int, const double *, double *, void *),
				 void *),
		 int it_max, double log10_eps,
		 int it_restart);
void
otmk (int m, const double *b, double *x,
      int kres, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, const double *, double *, void *),
      void * user_data);
