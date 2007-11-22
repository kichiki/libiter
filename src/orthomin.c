/* orthomin scheme
 * Copyright (C) 1999-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: orthomin.c,v 2.7 2007/11/22 05:48:31 kichiki Exp $
 *
 * solver routines are translated into C by K.I. from fortran code
 * originally written by martin h. gutknecht
 *             ===============rog.f==============
 *             problem given by martin h. gutknecht
 *             numerical method: orthomin(k) method
 *             ver. 1.0 jul. 28 1995 by s. l. zhang
 *             ver. 1.1 aug. 31 1995 by s. l. zhang
 *    at slzhang.fort.iter.complex.orthomin.gutknecht-problem(gut.f)
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
#include "../config.h"

#include <stdio.h> /* fprintf() */
#include <math.h> /* log10() */
#include <stdlib.h> /* malloc(), free() */
#include "memory-check.h" // CHECK_MALLOC
#include "libiter.h" // struct iter

#ifdef HAVE_CBLAS_H
/* use ATLAS' CBLAS routines */

#include <cblas.h>

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
/* use Fortran BLAS routines */

double
ddot_(int* N, 
      double* X, int* incX, 
      double* Y, int* incY);
int
dcopy_(int* N,
       double* X, int* incX,
       double* Y, int* incY);
int
daxpy_(int* N,
       double* alpha,
       double* X, int* incX,
       double* Y, int* incY);
int
dscal_(int* N,
       double* alpha,
       double* X, int* incX);

# else // !HAVE_BLAS_H
/* use local BLAS routines */

#include "myblas.h"

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


#include "orthomin.h"


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
otmk (int m, const double *b, double *x,
      int kres, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, const double *, double *, void *),
      void * user_data)
{
  int i, j;
  int k1, k2, k3;

  double tbs;
  double res;
  double rap;
  double alpha;

  /* for min0() */
  int jj;

  /**
   * allocation of matrices
   * r   [m]
   * p   [(kres+1) * m]
   * ap  [(kres+1) * m]
   * beta[kres+1]
   * pap [kres+1]
   * tmp [m] for myatimes
   */
  double *r    = (double *)malloc (sizeof (double) * m);
  double *p    = (double *)malloc (sizeof (double) * m * (kres + 1));
  double *ap   = (double *)malloc (sizeof (double) * m * (kres + 1));
  double *beta = (double *)malloc (sizeof (double) * (kres + 1));
  double *pap  = (double *)malloc (sizeof (double) * (kres + 1));
  double *tmp  = (double *)malloc (sizeof (double) * m);
  CHECK_MALLOC (r   , "otmk");
  CHECK_MALLOC (p   , "otmk");
  CHECK_MALLOC (ap  , "otmk");
  CHECK_MALLOC (beta, "otmk");
  CHECK_MALLOC (pap , "otmk");
  CHECK_MALLOC (tmp , "otmk");


  myatimes (m, x, tmp, user_data);
  for (i = 0; i < m; i++) /* 110 */
    {
      tbs = b[i] - tmp[i];
      r[i] = tbs;
      p[0 * m + i] = tbs;
    }/* 110 */
  myatimes (m, r, tmp, user_data);
  for (i = 0; i < m; i++) /* 120 */
    {
      ap[0 * m + i] = tmp[i];
    }/* 120 */

  for ((*iter) = 0; (*iter) <= kend; (*iter) ++) /* 10 */
    {
      k1 = (*iter) % (kres + 1);
      res = 0.0;
      rap = 0.0;
      pap[k1] = 0.0;
      for (i = 0; i < m; i++) /* 210 */
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
      for (i = 0; i < m; i++) /* 220 */
	{
	  r[i] -= alpha * ap[k1 * m + i];
	  x[i] += alpha * p[k1 * m + i];
	}/* 220 */
      k2 = ((*iter) + 1) % (kres + 1);
      myatimes (m, r, tmp, user_data);
      for (i = 0; i < m; i++) /* 230 */
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
	  for (i = 0; i < m; i++) /* 310 */
	    {
	      beta[k3] += ap[k2 * m + i] * ap[k3 * m + i];
	    }/* 310 */
	  beta[k3] = - beta[k3] / pap[k3];
	  for (i = 0; i < m; i++) /* 320 */
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

/* orthomin(k) method with BLAS/ATLAS
 * INPUT
 *   m : dimension of the problem
 *   b[m] : r-h-s vector
 *   atimes (int m, double *x, double *b) : calc matrix-vector product
 *   atimes_param : pointer to be passed to atimes routines
 *   it : struct iter. max, restart, log10_eps are used.
 * OUTPUT
 *   x[m] : solution
 *   *iter : # of iteration
 *   *hg : log10(residual)
 */
void
otmk_ (int m, const double *b, double *x,
       int *iter, double *hg,
       void (*atimes) (int, const double *, double *, void *),
       void *atimes_param,
       struct iter *it)
{
  int kend = it->max;
  int kres = it->restart;
  double eps2 = it->eps * it->eps;

#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_m1 = -1.0;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  /**
   * allocation of matrices
   * r   [m]
   * p   [(kres+1) * m]
   * ap  [(kres+1) * m]
   * beta[kres+1]
   * pap [kres+1]
   */
  double *r    = (double *)malloc (sizeof (double) * m);
  double *p    = (double *)malloc (sizeof (double) * m * (kres + 1));
  double *ap   = (double *)malloc (sizeof (double) * m * (kres + 1));
  double *beta = (double *)malloc (sizeof (double) * (kres + 1));
  double *pap  = (double *)malloc (sizeof (double) * (kres + 1));
  CHECK_MALLOC (r   , "otmk_");
  CHECK_MALLOC (p   , "otmk_");
  CHECK_MALLOC (ap  , "otmk_");
  CHECK_MALLOC (beta, "otmk_");
  CHECK_MALLOC (pap , "otmk_");


#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double res2;
  double b2 = cblas_ddot (m, b, 1, b, 1); // (b,b)

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  cblas_daxpy (m, -1.0, b, 1, r, 1);
  cblas_dscal (m, -1.0, r, 1);

  // p(0) = r(0)
  cblas_dcopy (m, r, 1, p, 1);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  for ((*iter) = 0; (*iter) <= kend; (*iter) ++)
    {
      res2 = cblas_ddot (m, r, 1, r, 1); // (r, r)
      res2 /= b2;
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk_ %d %e\n", (*iter), res2);
	}
      if(res2 <= eps2) goto end_otmk_;

      int k1 = (*iter) % (kres + 1);
      double rap = cblas_ddot (m, r, 1, ap + k1 * m, 1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1] = cblas_ddot (m, ap + k1 * m, 1, ap + k1 * m, 1);

      double alpha = rap / pap[k1];
      cblas_daxpy (m, +alpha, p  + k1 * m, 1, x, 1); // x = x + alpha * p(k1)
      cblas_daxpy (m, -alpha, ap + k1 * m, 1, r, 1); // r = r - alpha * ap(k1)

      int k2 = ((*iter) + 1) % (kres + 1);
      cblas_dcopy (m, r, 1, p + k2 * m, 1);     // p(k2) = r
      atimes (m, r, ap + k2 * m, atimes_param); // ap(k2) = A.r

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = (*iter);
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = ((*iter) - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.r
	  beta[k3] = cblas_ddot (m, ap + k2 * m, 1, ap + k3 * m, 1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  cblas_daxpy (m, beta [k3], p + k3 * m, 1, p + k2 * m, 1);
	  // ap(k2) += beta(k3) ap(k3)
	  cblas_daxpy (m, beta [k3], ap + k3 * m, 1, ap + k2 * m, 1);
	}
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double res2;
  double b2 = ddot_ (&m, b, &i_1, b, &i_1); // (b,b)

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  daxpy_ (&m, &d_m1, b, &i_1, r, &i_1);
  dscal_ (&m, &d_m1, r, &i_1);

  // p(0) = r(0)
  dcopy_ (&m, r, &i_1, p, &i_1);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  for ((*iter) = 0; (*iter) <= kend; (*iter) ++)
    {
      res2 = ddot_ (&m, r, &i_1, r, &i_1); // (r, r)
      res2 /= b2;
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk_ %d %e\n", (*iter), res2);
	}
      if(res2 <= eps2) goto end_otmk_;

      int k1 = (*iter) % (kres + 1);
      double rap = ddot_ (&m, r, &i_1, ap + k1 * m, &i_1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1] = ddot_ (&m, ap + k1 * m, &i_1, ap + k1 * m, &i_1);

      double alpha = rap / pap[k1];
      double malpha = -alpha;
      daxpy_ (&m, &alpha,  p  + k1 * m, &i_1, x, &i_1); // x = x + alpha*p(k1)
      daxpy_ (&m, &malpha, ap + k1 * m, &i_1, r, &i_1); // r = r - alpha*ap(k1)

      int k2 = ((*iter) + 1) % (kres + 1);
      dcopy_ (&m, r, &i_1, p + k2 * m, &i_1);   // p(k2) = r
      atimes (m, r, ap + k2 * m, atimes_param); // ap(k2) = A.r

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = (*iter);
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = ((*iter) - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.r
	  beta[k3] = ddot_ (&m, ap + k2 * m, &i_1, ap + k3 * m, &i_1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  daxpy_ (&m, beta + k3, p + k3 * m,  &i_1, p + k2 * m,  &i_1);
	  // ap(k2) += beta(k3) ap(k3)
	  daxpy_ (&m, beta + k3, ap + k3 * m, &i_1, ap + k2 * m, &i_1);
	}
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double res2;
  double b2 = my_ddot (m, b, 1, b, 1); // (b,b)

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  my_daxpy (m, -1.0, b, 1, r, 1);
  my_dscal (m, -1.0, r, 1);

  // p(0) = r(0)
  my_dcopy (m, r, 1, p, 1);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  for ((*iter) = 0; (*iter) <= kend; (*iter) ++)
    {
      res2 = my_ddot (m, r, 1, r, 1); // (r, r)
      res2 /= b2;
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk_ %d %e\n", (*iter), res2);
	}
      if(res2 <= eps2) goto end_otmk_;

      int k1 = (*iter) % (kres + 1);
      double rap = my_ddot (m, r, 1, ap + k1 * m, 1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1] = my_ddot (m, ap + k1 * m, 1, ap + k1 * m, 1);

      double alpha = rap / pap[k1];
      my_daxpy (m, +alpha, p  + k1 * m, 1, x, 1); // x = x + alpha * p(k1)
      my_daxpy (m, -alpha, ap + k1 * m, 1, r, 1); // r = r - alpha * ap(k1)

      int k2 = ((*iter) + 1) % (kres + 1);
      my_dcopy (m, r, 1, p + k2 * m, 1);     // p(k2) = r
      atimes (m, r, ap + k2 * m, atimes_param); // ap(k2) = A.r

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = (*iter);
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = ((*iter) - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.r
	  beta[k3] = my_ddot (m, ap + k2 * m, 1, ap + k3 * m, 1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  my_daxpy (m, beta [k3], p + k3 * m, 1, p + k2 * m, 1);
	  // ap(k2) += beta(k3) ap(k3)
	  my_daxpy (m, beta [k3], ap + k3 * m, 1, ap + k2 * m, 1);
	}
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

end_otmk_:
  (*hg) = log10 (res2) / 2.0;
  free (r);
  free (p);
  free (ap);
  free (beta);
  free (pap);

  if (it->debug == 1)
    {
      fprintf (it->out, "otmk_ %d %e\n", (*iter), res2);
    }
}

/* orthomin(k) method with preconditioning
 * INPUT
 *   m : dimension of the problem
 *   b[m] : r-h-s vector
 *   atimes (int m, double *x, double *b) : calc matrix-vector product
 *   atimes_param : pointer to be passed to atimes routines
 *   it : struct iter. max, restart, log10_eps are used.
 * OUTPUT
 *   x[m] : solution
 *   *iter : # of iteration
 *   *hg : log10(residual)
 */
void
otmk_pc (int m, const double *b, double *x,
	 int *iter, double *hg,
	 void (*atimes) (int, const double *, double *, void *),
	 void *atimes_param,
	 void (*inv) (int, const double *, double *, void *),
	 void *inv_param,
	 struct iter *it)
{
  int kend = it->max;
  int kres = it->restart;
  double eps2 = it->eps * it->eps;

#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_m1 = -1.0;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  /**
   * allocation of matrices
   * r   [m]
   * p   [(kres+1) * m]
   * ap  [(kres+1) * m]
   * beta[kres+1]
   * pap [kres+1]
   * Kr  [m] for K^{-1}.r, where K is the preconditioner.
   */
  double *r    = (double *)malloc (sizeof (double) * m);
  double *p    = (double *)malloc (sizeof (double) * m * (kres + 1));
  double *ap   = (double *)malloc (sizeof (double) * m * (kres + 1));
  double *beta = (double *)malloc (sizeof (double) * (kres + 1));
  double *pap  = (double *)malloc (sizeof (double) * (kres + 1));
  CHECK_MALLOC (r   , "otmk_pc");
  CHECK_MALLOC (p   , "otmk_pc");
  CHECK_MALLOC (ap  , "otmk_pc");
  CHECK_MALLOC (beta, "otmk_pc");
  CHECK_MALLOC (pap , "otmk_pc");


#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double res2;
  double b2 = cblas_ddot (m, b, 1, b, 1); // (b,b)

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  cblas_daxpy (m, -1.0, b, 1, r, 1);
  cblas_dscal (m, -1.0, r, 1);

  // p(0) = K^{-1}.r(0)
  inv (m, r, p, inv_param);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  for ((*iter) = 0; (*iter) <= kend; (*iter) ++)
    {
      res2 = cblas_ddot (m, r, 1, r, 1); // (r, r)
      res2 /= b2;
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk_pc %d %e\n", (*iter), res2);
	}
      if(res2 <= eps2) goto end_otmk_pc;

      int k1 = (*iter) % (kres + 1);
      double rap = cblas_ddot (m, r, 1, ap + k1 * m, 1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1]    = cblas_ddot (m, ap + k1 * m, 1, ap + k1 * m, 1);

      // alpha = (r, ap(k1)) / (ap(k1), ap(k1))
      double alpha = rap / pap[k1];
      cblas_daxpy (m, +alpha, p  + k1 * m, 1, x, 1); // x = x + alpha * p(k1)
      cblas_daxpy (m, -alpha, ap + k1 * m, 1, r, 1); // r = r - alpha * ap(k1)
      
      int k2 = ((*iter) + 1) % (kres + 1);
      inv (m, r, p + k2 * m, inv_param);                 // p(k2)  = K^{-1}.r
      atimes (m, p + k2 * m, ap + k2 * m, atimes_param); // ap(k2) = A.p(k2)

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = (*iter);
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = ((*iter) - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.K^{-1}.r
	  beta[k3] = cblas_ddot (m, ap + k2 * m, 1, ap + k3 * m, 1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  cblas_daxpy (m, beta [k3], p + k3 * m, 1, p + k2 * m, 1);
	  // ap(k2) += beta(k3) ap(k3)
	  cblas_daxpy (m, beta [k3], ap + k3 * m, 1, ap + k2 * m, 1);
	}
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double res2;
  double b2 = ddot_ (&m, b, &i_1, b, &i_1); // (b,b)

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  daxpy_ (&m, &d_m1, b, &i_1, r, &i_1);
  dscal_ (&m, &d_m1, r, &i_1);

  // p(0) = K^{-1}.r(0)
  inv (m, r, p, inv_param);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  for ((*iter) = 0; (*iter) <= kend; (*iter) ++)
    {
      res2 = ddot_ (&m, r, &i_1, r, &i_1); // (r, r)
      res2 /= b2;
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk_pc %d %e\n", (*iter), res2);
	}
      if(res2 <= eps2) goto end_otmk_pc;

      int k1 = (*iter) % (kres + 1);
      double rap = ddot_ (&m, r, &i_1, ap + k1 * m, &i_1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1] = ddot_ (&m, ap + k1 * m, &i_1, ap + k1 * m, &i_1);

      // alpha = (r, ap(k1)) / (ap(k1), ap(k1))
      double alpha = rap / pap[k1];
      double malpha = -alpha;
      daxpy_ (&m, &alpha,  p  + k1 * m, &i_1, x, &i_1); // x = x + alpha*p(k1)
      daxpy_ (&m, &malpha, ap + k1 * m, &i_1, r, &i_1); // r = r - alpha*ap(k1)
      
      int k2 = ((*iter) + 1) % (kres + 1);
      inv (m, r, p + k2 * m, inv_param);                 // p(k2)  = K^{-1}.r
      atimes (m, p + k2 * m, ap + k2 * m, atimes_param); // ap(k2) = A.p(k2)

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = (*iter);
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = ((*iter) - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.K^{-1}.r
	  beta[k3] = ddot_ (&m, ap + k2 * m, &i_1, ap + k3 * m, &i_1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  daxpy_ (&m, beta + k3, p + k3 * m,  &i_1, p + k2 * m,  &i_1);
	  // ap(k2) += beta(k3) ap(k3)
	  daxpy_ (&m, beta + k3, ap + k3 * m, &i_1, ap + k2 * m, &i_1);
	}
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double res2;
  double b2 = my_ddot (m, b, 1, b, 1); // (b,b)

  // r(0) = b - A.x
  atimes (m, x, r, atimes_param);
  my_daxpy (m, -1.0, b, 1, r, 1);
  my_dscal (m, -1.0, r, 1);

  // p(0) = K^{-1}.r(0)
  inv (m, r, p, inv_param);
  // ap(0) = A.p(0) (= A.r(0))
  atimes (m, p, ap, atimes_param);

  for ((*iter) = 0; (*iter) <= kend; (*iter) ++)
    {
      res2 = my_ddot (m, r, 1, r, 1); // (r, r)
      res2 /= b2;
      if (it->debug == 2)
	{
	  fprintf (it->out, "otmk_pc %d %e\n", (*iter), res2);
	}
      if(res2 <= eps2) goto end_otmk_pc;

      int k1 = (*iter) % (kres + 1);
      double rap = my_ddot (m, r, 1, ap + k1 * m, 1); // (r, ap(k1))
      // pap(k1) = (ap(k1), ap(k1))
      pap[k1]    = my_ddot (m, ap + k1 * m, 1, ap + k1 * m, 1);

      // alpha = (r, ap(k1)) / (ap(k1), ap(k1))
      double alpha = rap / pap[k1];
      my_daxpy (m, +alpha, p  + k1 * m, 1, x, 1); // x = x + alpha * p(k1)
      my_daxpy (m, -alpha, ap + k1 * m, 1, r, 1); // r = r - alpha * ap(k1)
      
      int k2 = ((*iter) + 1) % (kres + 1);
      inv (m, r, p + k2 * m, inv_param);                 // p(k2)  = K^{-1}.r
      atimes (m, p + k2 * m, ap + k2 * m, atimes_param); // ap(k2) = A.p(k2)

      /*for (j = 0; j <= min0(kres-1,k); j ++)*/
      int jj = (*iter);
      if (jj > (kres - 1))
	{
	  jj = kres - 1;
	}
      int j;
      for (j = 0; j <= jj; j ++)
	{
	  int k3 = ((*iter) - j) % (kres + 1);
	  // beta(k3) = - (ap(k2), ap(k3)) / (ap(k3), ap(k3))
	  // note that ap(k2) = A.K^{-1}.r
	  beta[k3] = my_ddot (m, ap + k2 * m, 1, ap + k3 * m, 1);
	  beta[k3] = - beta[k3] / pap[k3];

	  // p(k2) += beta(k3) p(k3)
	  my_daxpy (m, beta [k3], p + k3 * m, 1, p + k2 * m, 1);
	  // ap(k2) += beta(k3) ap(k3)
	  my_daxpy (m, beta [k3], ap + k3 * m, 1, ap + k2 * m, 1);
	}
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

end_otmk_pc:
  (*hg) = log10 (res2) / 2.0;
  free (r);
  free (p);
  free (ap);
  free (beta);
  free (pap);

  if (it->debug == 1)
    {
      fprintf (it->out, "otmk_pc %d %e\n", (*iter), res2);
    }
}
