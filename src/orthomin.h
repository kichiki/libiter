/* header file of orthomin.c
 * Copyright (C) 1999-2001 Kengo Ichiki <ichiki@haloumi.tn.utwente.nl>
 * $Id: orthomin.h,v 2.2 2001/10/19 14:32:55 ichiki Exp $
 */

/** global variables **/
extern int ITER_otmk_debug; /* [0|1]: [not print/print] iter and res */

void
solve_iter_otmk (int n,
		 double *b, double *x,
		 void (*atimes) (int, double *, double *, void *),
		 void * user_data,
		 void (*solver) (int, double *, double *,
				 int, int, double,
				 double, int *, double *,
				 void (*)
				 (int, double *, double *, void *),
				 void *),
		 int it_max, double log10_eps,
		 int it_restart);
void
otmk (int m, double *b, double *x,
      int kres, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, double *, double *, void *),
      void * user_data);
