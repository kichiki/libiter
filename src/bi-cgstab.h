/* header file of bi-cgstab.c (wrapper for iterative solver routines)
 * Copyright (C) 1999-2001 Kengo Ichiki <ichiki@haloumi.tn.utwente.nl>
 * $Id: bi-cgstab.h,v 1.1 2001/10/13 11:44:42 ichiki Exp $
 */

/** global variables **/
extern int STAB_it_max;
extern double STAB_log10_exp;


void
solve_iter_stab (int n,
		 double *b, double *x,
		 void (*atimes) (int, double *, double *),
		 void (*solver) (int, double *, double *, int,
				 double, double, int *, double *,
				 void (*) (int, double *, double *)));

void
sta (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *));
void
st2 (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *));
void
st2_chk (int m, double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, double *, double *));
void
gpb (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *));
void
gpb_chk (int m, double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, double *, double *));
