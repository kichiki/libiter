/* header file of bi-cgstab.c (wrapper for iterative solver routines)
 * Copyright (C) 1999-2001 Kengo Ichiki <ichiki@haloumi.tn.utwente.nl>
 * $Id: bi-cgstab.h,v 2.2 2001/10/19 14:31:22 ichiki Exp $
 */

/** global variables **/
extern int ITER_stab_debug; /* [0|1]: [not print/print] iter and res */

void
solve_iter_stab (int n,
		 double *b, double *x,
		 void (*atimes) (int, double *, double *, void *),
		 void * user_data,
		 void (*solver) (int, double *, double *, int,
				 double, double, int *, double *,
				 void (*)
				 (int, double *, double *, void *),
				 void *),
		 int it_max, double log10_eps);
void
sta (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *, void *),
     void * user_data);
void
sta2 (int m, double *b, double *x, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, double *, double *, void *),
      void * user_data);
void
st2_chk (int m, double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, double *, double *, void *),
	 void * user_data);
void
gpb (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *, void *),
     void * user_data);
void
gpb_chk (int m, double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, double *, double *, void *),
	 void * user_data);
