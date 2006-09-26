/* header file of bi-cgstab.c --
 * wrapper for iterative solver routines
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bi-cgstab.h,v 2.3 2006/09/26 05:35:12 ichiki Exp $
 *
 * (solver routines themselves are originally written by martin h. gutknecht)
 */

/** global variables **/
extern int ITER_stab_debug; /* [0|1]: [not print/print] iter and res */

void
solve_iter_stab (int n, const double *b,
		 double *x,
		 void (*atimes) (int, const double *, double *, void *),
		 void * user_data,
		 void (*solver) (int, const double *, double *, int,
				 double, double, int *, double *,
				 void (*)
				 (int, const double *, double *, void *),
				 void *),
		 int it_max, double log10_eps);
void
sta (int m, const double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, const double *, double *, void *),
     void * user_data);
void
sta2 (int m, const double *b, double *x, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, const double *, double *, void *),
      void * user_data);
void
st2_chk (int m, const double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, const double *, double *, void *),
	 void * user_data);
void
gpb (int m, const double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, const double *, double *, void *),
     void * user_data);
void
gpb_chk (int m, const double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, const double *, double *, void *),
	 void * user_data);
