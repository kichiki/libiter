/*             ===============rog.f==============
 *             problem given by martin h. gutknecht
 *             numerical method: orthomin(k) method
 *             ver. 1.0 jul. 28 1995 by s. l. zhang
 *             ver. 1.1 aug. 31 1995 by s. l. zhang
 *    at slzhang.fort.iter.complex.orthomin.gutknecht-problem(gut.f)
 *
 * translated from fortran to C by Kengo Ichiki <ichiki@haloumi.tn.utwente.nl>
 * $Id: orthomin.c,v 2.1 2001/10/13 12:01:27 ichiki Exp $
 */
#include <stdio.h> /* fprintf() */
#include <math.h> /* log10() */
#include <stdlib.h> /* malloc(), free() */

#include "orthomin.h"


/* wrapper routine for solvers below
 * INPUT
 *   n : size of vectors v[] and f[] -- expected to be np * nelm for red-sym
 *   b [n] : given vector
 *   atimes (n, x, b, user_data) : routine to calc A.x and return b[]
 *   user_data : pointer to be passed to solver and atimes routines
 *   solver : solver routine to call
 *   it_max : max # iterations
 *   log10_eps : log10(eps), where eps is the accuracy
 * OUTPUT
 *   x [n] : solution
 */
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
		 int it_restart)
{
  int i;

  double hnor;
  double residual;
  int iter;


  hnor = 0.0;
  for (i = 0; i < n; ++i)
    hnor += b [i] * b [i];
  if (hnor != 0.0)
    hnor = log10 (hnor) / 2.0;

  solver (n, b, x,
	  it_restart, it_max, log10_eps,
	  hnor, &iter, &residual, atimes,
	  user_data);

  fprintf (stderr, "# iter=%d res=%e\n", iter, residual);
}


/* orthomin(k) method
 * INPUT
 *   m : dimension of the problem
 *   kres : steps for restart
 *   kend : 
 *   b[m] : r-h-s vector
 *   eps : log10 of cutoff
 *   hnor : log10 of norm of b[]
 *   myatimes (int m, double *x, double *b) : calc matrix-vector product
 *   user_data : pointer to be passed to atimes routines
 * OUTPUT
 *   x[m] : solution
 *   *iter : # of iteration
 *   *hg : log10(residual)
 */
void
otmk (int m, double *b, double *x,
      int kres, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, double *, double *, void *),
      void * user_data)
{
  double *r; /* r[m] */
  double *p; /* p[(kres+1) * m] */
  double *ap; /* ap[(kres+1) * m] */
  double *beta; /* beta[kres+1] */
  double *pap; /* pap[kres+1] */


  int i, j;
  int k1, k2, k3;

  double tbs;
  double res;
  double rap;
  double alpha;

  /* for min0() */
  int jj;

  /* myatimes */
  double *tmp; /* tmp[m] */


  /* allocation of matrices
  r    = my_d_malloc (m, "r");
  p    = my_d_malloc (m * (kres + 1), "p");
  ap   = my_d_malloc (m * (kres + 1), "ap");
  beta = my_d_malloc ((kres + 1), "beta");
  pap  = my_d_malloc ((kres + 1), "pap");
  tmp  = my_d_malloc (m, "tmp"); */
  r    = (double *) malloc (sizeof (double) * m);
  p    = (double *) malloc (sizeof (double) * m * (kres + 1));
  ap   = (double *) malloc (sizeof (double) * m * (kres + 1));
  beta = (double *) malloc (sizeof (double) * (kres + 1));
  pap  = (double *) malloc (sizeof (double) * (kres + 1));
  tmp  = (double *) malloc (sizeof (double) * m);
  if (r == NULL
      || p == NULL
      || ap == NULL
      || beta == NULL
      || pap == NULL
      || tmp == NULL)
    {
      fprintf (stderr, "malloc in otmk ()");
      exit (1);
    }


  myatimes (m, x, tmp, user_data);
  for (i=0; i<m; i++) /* 110 */
    {
      tbs = b[i] - tmp[i];
      r[i] = tbs;
      p[0 * m + i] = tbs;
    }/* 110 */
  myatimes (m, r, tmp, user_data);
  for (i=0; i<m; i++) /* 120 */
    {
      ap[0 * m + i] = tmp[i];
    }/* 120 */

  for ((*iter) = 0; (*iter) <= kend; (*iter) ++) /* 10 */
    {
      k1 = (*iter) % (kres + 1);
      res = 0.0;
      rap = 0.0;
      pap[k1] = 0.0;
      for (i=0; i<m; i++) /* 210 */
	{
          res += r[i] * r[i];
	  rap += r[i] * ap[k1 * m + i];
	  pap[k1] += ap[k1 * m + i] * ap[k1 * m + i];
	}/* 210 */

      /* call resd(62,k) */
      (*hg) = log10 (res) / 2.0 - hnor;
      /*fprintf (stdout, "32 %d %e\n", k, hg);*/

      if((*hg) <= eps) goto end_otmk;
      alpha = rap / pap[k1];
      for (i=0; i<m; i++) /* 220 */
	{
	  r[i] -= alpha * ap[k1 * m + i];
	  x[i] += alpha * p[k1 * m + i];
	}/* 220 */
      k2 = ((*iter) + 1) % (kres + 1);
      myatimes (m, r, tmp, user_data);
      for (i=0; i<m; i++) /* 230 */
	{
	  p[k2 * m + i] = r[i];
	  ap[k2 * m + i] = tmp[i];
	}/* 230 */
      /*for (j = 0; j <= min0(kres-1,k); j ++)*/ /* 240 */
      jj = (*iter);
      if ((*iter) > (kres - 1))
	jj = kres - 1;
      for (j = 0; j <= jj; j ++) /* 240 */
	{
	  k3 = ((*iter) - j) % (kres + 1);
	  beta[k3] = 0.0;
	  for (i=0; i<m; i++) /* 310 */
	    {
	      beta[k3] += ap[k2 * m + i] * ap[k3 * m + i];
	    }/* 310 */
	  beta[k3] = - beta[k3] / pap[k3];
	  for (i=0; i<m; i++) /* 320 */
	    {
	      p[k2 * m + i] += beta[k3] * p[k3 * m + i];
	      ap[k2 * m + i] += beta[k3] * ap[k3 * m + i];
	    }/* 320 */
	}/* 240 */
    }/* 10 */

end_otmk:
  free (r);
  free (p);
  free (ap);
  free (beta);
  free (pap);
  free (tmp);
}
