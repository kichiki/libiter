/* iterative solver in Weiss (1996)
 * Copyright (C) 1999 Kengo ICHIKI <kengo@caltech.edu>
 * $Id: atpres.c,v 1.1 1999/09/05 21:01:24 ichiki Exp $
 *
 * Reference :
 *   R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h> /* free () */
#include "myroutines.h" /* my_d_malloc () */

#include "blas.h"

#include "weiss.h"


void
orthores (int n, double *b, double *x,
	  int nres,
	  int maxiter, double eps,
	  int *iter, double *res,
	  void (*atimes) (int, double *, double *))
{
  /*  for details to the algorithm see
   *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
   *  Algorithm 3 : ORTHORES
   */
  double phi;
  double sum;
  double *alpha;
  double b_norm;
  double *r;
  double *xx;
  double *rr;
  double *tmp;
  double tiny = 1.0e-15;
  int i;
  int ires;
  int ires1;

  alpha = my_d_malloc (nres, "alpha");
  r = my_d_malloc (n, "r");
  xx = my_d_malloc ((nres + 1) * n, "xx");
  rr = my_d_malloc ((nres + 1) * n, "rr");
  tmp = my_d_malloc (n, "tmp");


  b_norm = dnrm2 (n, b, 1);

  (*iter) = 0;
  /* r = A * x(i) - b */
  if (fabs (dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp);
      daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    dscalz (n, - 1.0, b, 1, r, 1);


  /* initial condition */
  (*res) = dnrm2 (n, r, 1);
  dcopy (n, r, 1, rr, 1);
  dcopy (n, x, 1, xx, 1);

  while ((*res) > eps * b_norm
	 && (*iter) < maxiter)
    {
      fprintf (stderr, "#ORTHORES %d %e\n", (*iter), (*res));
      ires = (*iter) % nres;
      ires1 = ires + 1;
      (*iter) ++;

      /* set alpha */
      atimes (n, & rr [ires * n], tmp);
      for (i=0; i<=ires; i++)
	alpha [i] =
	  - ddot (n, & rr [i * n], 1, tmp, 1)
	  / ddot (n, & rr [i * n], 1, & rr [i * n], 1);

      /* set phi */
      sum = 0.0;
      for (i=0; i<=ires; i++)
	sum += alpha [i];
      phi = 1.0 / sum;

      /* update xx */
      dcopy (n, & rr [ires * n], 1, & xx [ires1 * n], 1);
      for (i=0; i<=ires; i++)
	daxpy (n, alpha [i], & xx [i * n], 1, & xx [ires1 * n], 1);
      dscal (n, phi, & xx [ires1 * n], 1);

      /* update rr */
      dcopy (n, tmp, 1, & rr [ires1 * n], 1);
      for (i=0; i<=ires; i++)
	daxpy (n, alpha [i], & rr [i * n], 1, & rr [ires1 * n], 1);
      dscal (n, phi, & rr [ires1 * n], 1);

      dcopy (n, & xx [ires1 * n], 1, x, 1);
      dcopy (n, & rr [ires1 * n], 1, r, 1);

      (*res) = dnrm2 (n, r, 1);

      if (ires1 == nres)
	{
	  /* restart */
	  dcopy (n, r, 1, rr, 1);
	  dcopy (n, x, 1, xx, 1);
	}
    }

  free (alpha);
  free (r);
  free (xx);
  free (rr);
  free (tmp);
}

void
gmres (int n, double *b, double *x,
       int nres,
       int maxiter, double eps,
       int *iter, double *res,
       void (*atimes) (int, double *, double *))
{
  /*  for details to the algorithm see
   *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
   *  Algorithm 4 : GMRES
   */
  double gamma;
  double phi;
  double sum;
  double *alpha;
  double b_norm;
  double *r;
  double *xx;
  double *rr;
  double *tmp;
  double tiny = 1.0e-15;
  int i;
  int ires;
  int ires1;

  alpha = my_d_malloc (nres, "alpha");
  r = my_d_malloc (n, "r");
  xx = my_d_malloc ((nres + 1) * n, "xx");
  rr = my_d_malloc ((nres + 1) * n, "rr");
  tmp = my_d_malloc (n, "tmp");


  b_norm = dnrm2 (n, b, 1);

  (*iter) = 0;
  /* r = A * x(i) - b */
  if (fabs (dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp);
      daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    dscalz (n, - 1.0, b, 1, r, 1);


  /* initial condition */
  (*res) = dnrm2 (n, r, 1);
  dcopy (n, r, 1, rr, 1);
  dcopy (n, x, 1, xx, 1);

  while ((*res) > eps * b_norm
	 && (*iter) < maxiter)
    {
      fprintf (stderr, "#GMRES %d %e\n", (*iter), (*res));
      ires = (*iter) % nres;
      ires1 = ires + 1;
      (*iter) ++;

      /* set alpha */
      atimes (n, & rr [ires * n], tmp);
      for (i=0; i<=ires; i++)
	alpha [i] =
	  - ddot (n, & rr [i * n], 1, tmp, 1)
	  / ddot (n, & rr [i * n], 1, & rr [i * n], 1);

      /* set phi */
      sum = 0.0;
      for (i=0; i<=ires; i++)
	sum += alpha [i];
      phi = 1.0 / sum;

      /* update xx */
      dcopy (n, & rr [ires * n], 1, & xx [ires1 * n], 1);
      for (i=0; i<=ires; i++)
	daxpy (n, alpha [i], & xx [i * n], 1, & xx [ires1 * n], 1);
      dscal (n, phi, & xx [ires1 * n], 1);

      /* update rr */
      dcopy (n, tmp, 1, & rr [ires1 * n], 1);
      for (i=0; i<=ires; i++)
	daxpy (n, alpha [i], & rr [i * n], 1, & rr [ires1 * n], 1);
      dscal (n, phi, & rr [ires1 * n], 1);

      /* set gamma */
      daxpyz (n, -1.0, r, 1, & rr [ires1 * n], 1, tmp, 1);
      gamma = - ddot (n, r, 1, tmp, 1) / ddot (n, tmp, 1, tmp, 1);

      /* update r */
      daxpy (n, gamma, tmp, 1, r, 1);
      /* update x */
      daxpyz (n, -1.0, x, 1, & xx [ires1 * n], 1, tmp, 1);
      daxpy (n, gamma, tmp, 1, x, 1);

      (*res) = dnrm2 (n, r, 1);

      if (ires1 == nres)
	{
	  /* restart */
	  dcopy (n, r, 1, rr, 1);
	  dcopy (n, x, 1, xx, 1);
	}
    }

  free (alpha);
  free (r);
  free (xx);
  free (rr);
  free (tmp);
}

void
atpres (int n, double *b, double *x,
	int maxiter, double eps,
	int *iter, double *res,
	void (*atimes) (int, double *, double *),
	void (*atimes_t) (int, double *, double *))
{
  /*  for details to the algorithm see
   *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
   *  Algorithm 5 : ATPRES
   */
  double gamma;
  double phi;
  double sum;
  double *alpha;
  double b_norm;
  double *r;
  double *xx;
  double *rr;
  double *tmp;
  double *tmp0;
  double tiny = 1.0e-15;
  int i;
  int mink2;


  alpha = my_d_malloc (2, "alpha");
  r = my_d_malloc (n, "r");
  xx = my_d_malloc (2 * n, "xx");
  rr = my_d_malloc (2 * n, "rr");
  tmp = my_d_malloc (n, "tmp");
  tmp0 = my_d_malloc (n, "tmp0");


  b_norm = dnrm2 (n, b, 1);

  (*iter) = 0;
  /* r = A * x(i) - b */
  if (fabs (dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp);
      daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    dscalz (n, - 1.0, b, 1, r, 1);


  /* initial condition */
  (*res) = dnrm2 (n, r, 1);
  dcopy (n, r, 1, rr, 1);
  dcopy (n, x, 1, xx, 1);

  while ((*res) > eps * b_norm
	 && (*iter) < maxiter)
    {
      fprintf (stderr, "#ATPRES %d %e\n", (*iter), (*res));
      (*iter) ++;

      /* set min (k, 2) */
      if ((*iter) < 2)
	mink2 = (*iter);
      else
	mink2 = 2;

      /* set alpha */
      for (i=0; i<mink2; i++)
	{
	  atimes_t (n, & rr [i * n], tmp);
	  atimes_t (n, & rr [0], tmp0);
	  alpha [i] =
	    - ddot (n, tmp, 1, tmp0, 1)
	    / ddot (n, & rr [i * n], 1, & rr [i * n], 1);
	}

      /* set phi */
      sum = 0.0;
      for (i=0; i<mink2; i++)
	sum += alpha [i];
      phi = 1.0 / sum;

      /* update xx */
      dcopy (n, tmp0, 1, tmp, 1);
      for (i=0; i<mink2; i++)
	daxpy (n, alpha [i], & xx [i * n], 1, tmp, 1);
      dcopy (n, & xx [0], 1, & xx [n], 1);
      dscalz (n, phi, tmp, 1, xx, 1);

      /* update rr */
      atimes (n, tmp0, tmp);
      for (i=0; i<mink2; i++)
	daxpy (n, alpha [i], & rr [i * n], 1, tmp, 1);
      dcopy (n, & rr [0], 1, & rr [n], 1);
      dscalz (n, phi, tmp, 1, rr, 1);

      /* set gamma */
      daxpyz (n, -1.0, r, 1, rr, 1, tmp, 1);
      gamma = - ddot (n, r, 1, tmp, 1) / ddot (n, tmp, 1, tmp, 1);

      /* update r */
      daxpy (n, gamma, tmp, 1, r, 1);
      /* update x */
      daxpyz (n, -1.0, x, 1, xx, 1, tmp, 1);
      daxpy (n, gamma, tmp, 1, x, 1);

      (*res) = dnrm2 (n, r, 1);
    }

  free (alpha);
  free (r);
  free (xx);
  free (rr);
  free (tmp);
  free (tmp0);
}

void
cgne (int n, double *b, double *x,
      int maxiter, double eps,
      int *iter, double *res,
      void (*atimes) (int, double *, double *),
      void (*atimes_t) (int, double *, double *))
{
  /*  for details to the algorithm see
   *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
   *  Algorithm 6 : CGNE
   */
  double phi;
  double sum;
  double *alpha;
  double b_norm;
  double *r;
  double *xx;
  double *tmp;
  double *tmp0;
  double tiny = 1.0e-15;
  int i;
  int mink2;


  alpha = my_d_malloc (2, "alpha");
  r = my_d_malloc (2 * n, "r");
  xx = my_d_malloc (2 * n, "xx");
  tmp = my_d_malloc (n, "tmp");
  tmp0 = my_d_malloc (n, "tmp0");


  b_norm = dnrm2 (n, b, 1);

  (*iter) = 0;
  /* r = A * x(i) - b */
  if (fabs (dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp);
      daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    dscalz (n, - 1.0, b, 1, r, 1);


  /* initial condition */
  (*res) = dnrm2 (n, r, 1);
  dcopy (n, x, 1, xx, 1);

  while ((*res) > eps * b_norm
	 && (*iter) < maxiter)
    {
      fprintf (stderr, "#CGNE %d %e\n", (*iter), (*res));
      (*iter) ++;

      /* set min (k, 2) */
      if ((*iter) < 2)
	mink2 = (*iter);
      else
	mink2 = 2;

      /* set alpha */
      for (i=0; i<mink2; i++)
	{
	  atimes_t (n, & r [i * n], tmp);
	  atimes_t (n, & r [0], tmp0);
	  alpha [i] =
	    - ddot (n, tmp, 1, tmp0, 1)
	    / ddot (n, & r [i * n], 1, & r [i * n], 1);
	}

      /* set phi */
      sum = 0.0;
      for (i=0; i<mink2; i++)
	sum += alpha [i];
      phi = 1.0 / sum;

      /* update xx */
      dcopy (n, tmp0, 1, tmp, 1);
      for (i=0; i<mink2; i++)
	daxpy (n, alpha [i], & xx [i * n], 1, tmp, 1);
      dcopy (n, & xx [0], 1, & xx [n], 1);
      dscalz (n, phi, tmp, 1, xx, 1);

      /* update r */
      atimes (n, tmp0, tmp);
      for (i=0; i<mink2; i++)
	daxpy (n, alpha [i], & r [i * n], 1, tmp, 1);
      dcopy (n, & r [0], 1, & r [n], 1);
      dscalz (n, phi, tmp, 1, r, 1);

      (*res) = dnrm2 (n, r, 1);
    }

  dcopy (n, xx, 1, x, 1);

  free (alpha);
  free (r);
  free (xx);
  free (tmp);
  free (tmp0);
}

void
bcg (int n, double *b, double *x,
     int maxiter, double eps,
     int *iter, double *res,
     void (*atimes) (int, double *, double *),
     void (*atimes_t) (int, double *, double *))
{
  /*  for details to the algorithm see
   *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
   *  Algorithm 8 : BCG
   */

  double beta;
  double delta;
  double rho, rho_old = 1.0;
  double b_norm;
  double *r, *r_;
  double *p, *p_;
  double *tmp;
  double tiny = 1.0e-15;


  r = my_d_malloc (n, "r");
  r_ = my_d_malloc (n, "r_");
  p = my_d_malloc (n, "p");
  p_ = my_d_malloc (n, "p_");
  tmp = my_d_malloc (n, "tmp");


  b_norm = dnrm2 (n, b, 1);

  (*iter) = 0;
  /* r = A * x(i) - b */
  if (fabs (dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp);
      daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    dscalz (n, - 1.0, b, 1, r, 1);

  /* plain BiCG (z = r, z_ = r_) */
  (*res) = dnrm2 (n, r, 1);
  dcopy (n, r, 1, p, 1);
  dcopy (n, r, 1, r_, 1);
  dcopy (n, r, 1, p_, 1);
  rho = ddot (n, r, 1, r_, 1);

  while ((*res) > eps * b_norm
	 && (*iter) < maxiter)
    {
      fprintf (stderr, "#BCG %d %e\n", (*iter), (*res));
      (*iter) ++;

      /* set delta */
      atimes (n, p, tmp);
      delta = - rho / ddot (n, p_, 1, tmp, 1);

      /* update r, r_ */
      daxpy (n, delta, tmp, 1, r, 1);
      atimes_t (n, p_, tmp);
      daxpy (n, delta, tmp, 1, r_, 1);

      /* update x */
      daxpy (n, delta, p, 1, x, 1);

      /* set beta */
      rho_old = rho;
      rho = ddot (n, r, 1, r_, 1);
      beta = rho / rho_old;

      /* set p and p_ */
      daxpyz (n, beta, p, 1, r, 1, p, 1);
      daxpyz (n, beta, p_, 1, r_, 1, p_, 1);

      (*res) = dnrm2 (n, r, 1);
    }

  free (r);
  free (r_);
  free (p);
  free (p_);
  free (tmp);
}

void
bico (int n, double *b, double *x,
      int maxiter, double eps,
      int *iter, double *res,
      void (*atimes) (int, double *, double *),
      void (*atimes_t) (int, double *, double *))
{
  /*  for details to the algorithm see
   *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
   *  Algorithm 9 : BICO
   */
  double beta;
  double delta;
  double gamma;
  double rho, rho_old = 1.0;
  double b_norm;
  double *r, *r_;
  double *rr, *rr_;
  double *p, *p_;
  double *xx;
  double *q;
  double *tmp;
  double tiny = 1.0e-15;


  r = my_d_malloc (n, "r"); /* r */
  r_ = my_d_malloc (n, "r_"); /* r^* */
  rr = my_d_malloc (n, "rr"); /* ~r */
  rr_ = my_d_malloc (n, "rr_"); /* ~r^* */
  p = my_d_malloc (n, "p"); /* p */
  p_ = my_d_malloc (n, "p_"); /* p^* */
  xx = my_d_malloc (n, "xx"); /* ~x */
  q = my_d_malloc (n, "q"); /* A . p */
  tmp = my_d_malloc (n, "tmp");


  b_norm = dnrm2 (n, b, 1);

  (*iter) = 0;
  /* r = A * x(i) - b */
  if (fabs (dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp);
      daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    dscalz (n, - 1.0, b, 1, r, 1);


  (*res) = dnrm2 (n, r, 1);

  /* initial condition */
  dcopy (n, r, 1, rr, 1);
  dcopy (n, r, 1, p, 1);
  dcopy (n, r, 1, r_, 1);
  dcopy (n, r, 1, rr_, 1);
  dcopy (n, r, 1, p_, 1);
  dcopy (n, x, 1, xx, 1);
  rho = ddot (n, rr, 1, rr_, 1);

  while ((*res) > eps * b_norm
	 && (*iter) < maxiter)
    {
      fprintf (stderr, "#BICO %d %e\n", (*iter), (*res));
      (*iter) ++;

      /* set delta */
      atimes (n, p, q);
      atimes_t (n, p_, tmp);
      delta = - rho / ddot (n, p, 1, tmp, 1);

      /* update rr, rr_ */
      daxpy (n, delta, q, 1, rr, 1);
      daxpy (n, delta, tmp, 1, rr_, 1);

      /* update xx */
      daxpy (n, delta, p, 1, xx, 1);

      /* set beta */
      rho_old = rho;
      rho = ddot (n, rr, 1, rr_, 1);
      beta = rho / rho_old;

      /* update p, p_ */
      daxpyz (n, beta, p, 1, rr, 1, p, 1);
      daxpyz (n, beta, p_, 1, rr_, 1, p_, 1);

      /* set gamma */
      daxpyz (n, -1.0, r, 1, rr, 1, tmp, 1);
      gamma = - ddot (n, r, 1, tmp, 1) / ddot (n, tmp, 1, tmp, 1);

      /* update r */
      daxpy (n, gamma, tmp, 1, r, 1);
      /* update x */
      daxpyz (n, -1.0, x, 1, xx, 1, tmp, 1);
      daxpy (n, gamma, tmp, 1, x, 1);

      (*res) = dnrm2 (n, r, 1);
    }

  free (r);
  free (r_);
  free (rr);
  free (rr_);
  free (p);
  free (p_);
  free (xx);
  free (q);
  free (tmp);
}

void
qmr (int n, double *b, double *x,
     int maxiter, double eps,
     int *iter, double *res,
     void (*atimes) (int, double *, double *),
     void (*atimes_t) (int, double *, double *))
{
  /*  for details to the algorithm see
   *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
   *  Algorithm 10 : QMR
   */
  double beta;
  double delta;
  double gamma;
  double rho, rho_old = 1.0;
  double tau2, tau2_old, rr_norm2;
  double b_norm;
  double *r, *r_;
  double *rr, *rr_;
  double *p, *p_;
  double *xx;
  double *q;
  double *tmp;
  double tiny = 1.0e-15;


  r = my_d_malloc (n, "r"); /* r */
  r_ = my_d_malloc (n, "r_"); /* r^* */
  rr = my_d_malloc (n, "rr"); /* ~r */
  rr_ = my_d_malloc (n, "rr_"); /* ~r^* */
  p = my_d_malloc (n, "p"); /* p */
  p_ = my_d_malloc (n, "p_"); /* p^* */
  xx = my_d_malloc (n, "xx"); /* ~x */
  q = my_d_malloc (n, "q"); /* A . p */
  tmp = my_d_malloc (n, "tmp");


  b_norm = dnrm2 (n, b, 1);

  (*iter) = 0;
  /* r = A * x(i) - b */
  if (fabs (dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp);
      daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    dscalz (n, - 1.0, b, 1, r, 1);


  (*res) = dnrm2 (n, r, 1);

  /* initial condition */
  dcopy (n, r, 1, rr, 1);
  dcopy (n, r, 1, p, 1);
  dcopy (n, r, 1, r_, 1);
  dcopy (n, r, 1, rr_, 1);
  dcopy (n, r, 1, p_, 1);
  dcopy (n, x, 1, xx, 1);
  rho = ddot (n, rr, 1, rr_, 1);
  tau2 = ddot (n, rr, 1, rr, 1);

  while ((*res) > eps * b_norm
	 && (*iter) < maxiter)
    {
      fprintf (stderr, "#QMR %d %e\n", (*iter), (*res));
      (*iter) ++;

      /* set delta */
      atimes (n, p, q);
      atimes_t (n, p_, tmp);
      delta = - rho / ddot (n, p, 1, tmp, 1);

      /* update rr, rr_ */
      daxpy (n, delta, q, 1, rr, 1);
      daxpy (n, delta, tmp, 1, rr_, 1);

      /* update xx */
      daxpy (n, delta, p, 1, xx, 1);

      /* set beta */
      rho_old = rho;
      rho = ddot (n, rr, 1, rr_, 1);
      beta = rho / rho_old;

      /* update p, p_ */
      daxpyz (n, beta, p, 1, rr, 1, p, 1);
      daxpyz (n, beta, p_, 1, rr_, 1, p_, 1);

      /* set gamma */
      tau2_old = tau2;
      rr_norm2 = ddot (n, rr, 1, rr, 1);
      tau2 = 1.0 / (1.0 / tau2_old + 1.0 / rr_norm2);
      gamma = tau2 / rr_norm2;

      /* update x */
      daxpyz (n, -1.0, x, 1, xx, 1, tmp, 1);
      daxpy (n, gamma, tmp, 1, x, 1);
      /* update r */
      daxpyz (n, -1.0, r, 1, rr, 1, tmp, 1);
      daxpy (n, gamma, tmp, 1, r, 1);

      (*res) = dnrm2 (n, r, 1);
    }

  free (r);
  free (r_);
  free (rr);
  free (rr_);
  free (p);
  free (p_);
  free (xx);
  free (q);
  free (tmp);
}

void
cgs (int n, double *b, double *x,
     int maxiter, double eps,
     int *iter, double *res,
     void (*atimes) (int, double *, double *))
{
  /*  for details to the algorithm see
   *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
   *  Algorithm 11 : CGS
   */
  double beta;
  double delta;
  double rho, rho_old = 1.0;
  double b_norm;
  double *r, *r0;
  double *p;
  double *u;
  double *q;
  double *tmp;
  double tiny = 1.0e-15;


  r = my_d_malloc (n, "r"); /* r */
  r0 = my_d_malloc (n, "r0"); /* r0 */
  p = my_d_malloc (n, "p"); /* p */
  u = my_d_malloc (n, "u"); /* u */
  q = my_d_malloc (n, "q"); /* q */
  tmp = my_d_malloc (n, "tmp");


  b_norm = dnrm2 (n, b, 1);

  (*iter) = 0;
  /* r = A * x(i) - b */
  if (fabs (dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp);
      daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    dscalz (n, - 1.0, b, 1, r, 1);


  (*res) = dnrm2 (n, r, 1);

  /* initial condition */
  dcopy (n, r, 1, p, 1);
  dcopy (n, r, 1, u, 1);
  dcopy (n, r, 1, r0, 1);
  rho = ddot (n, r, 1, r0, 1);

  while ((*res) > eps * b_norm
	 && (*iter) < maxiter)
    {
      fprintf (stderr, "#CGS %d %e\n", (*iter), (*res));
      (*iter) ++;

      /* set delta */
      atimes (n, p, tmp);
      delta = - rho / ddot (n, r0, 1, tmp, 1);

      /* update q */
      dscal (n, delta, tmp, 1);
      daxpyz (n, 2.0, u, 1, tmp, 1, q, 1);

      /* update r, x */
      atimes (n, q, tmp);
      daxpy (n, delta, tmp, 1, r, 1);
      daxpy (n, delta, q, 1, x, 1);

      /* set beta */
      rho_old = rho;
      rho = ddot (n, r, 1, r0, 1);
      beta = rho / rho_old;

      /* update u, p */
      daxpyz (n, -1.0, u, 1, q, 1, tmp, 1);
      daxpyz (n, beta, tmp, 1, r, 1, u, 1);
      daxpy (n, beta, p, 1, tmp, 1);
      daxpyz (n, beta, tmp, 1, u, 1, p, 1);

      (*res) = dnrm2 (n, r, 1);
    }

  free (r);
  free (r0);
  free (p);
  free (u);
  free (q);
  free (tmp);
}

void
bicgstab (int n, double *b, double *x,
	  int maxiter, double eps,
	  int *iter, double *res,
	  void (*atimes) (int, double *, double *))
{
  /*  for details to the algorithm see
   *  R. Weiss : Parameter-Free Iterative Linear Solvers (1996)
   *  Algorithm 12 : BiCGSTAB
   */
  double beta;
  double delta;
  double gamma;
  double rho, rho_old = 1.0;
  double b_norm;
  double *r, *r0;
  double *p;
  double *v;
  double *s;
  double *t;
  double *tmp;
  double tiny = 1.0e-15;


  r = my_d_malloc (n, "r"); /* r */
  r0 = my_d_malloc (n, "r0"); /* r0^* */
  p = my_d_malloc (n, "p"); /* p */
  v = my_d_malloc (n, "v"); /* v */
  s = my_d_malloc (n, "s"); /* s */
  t = my_d_malloc (n, "t"); /* t */
  tmp = my_d_malloc (n, "tmp");


  b_norm = dnrm2 (n, b, 1);

  (*iter) = 0;
  /* r = A * x(i) - b */
  if (fabs (dnrm2 (n, x, 1) / (double) n) > tiny)
    {
      atimes (n, x, tmp);
      daxpyz (n, - 1.0, b, 1, tmp, 1, r, 1);
    }
  else
    dscalz (n, - 1.0, b, 1, r, 1);


  (*res) = dnrm2 (n, r, 1);

  /* initial condition */
  dcopy (n, r, 1, p, 1);
  dcopy (n, r, 1, r0, 1);
  rho = ddot (n, r, 1, r0, 1);

  while ((*res) > eps * b_norm
	 && (*iter) < maxiter)
    {
      fprintf (stderr, "#BiCGSTAB %d %e\n", (*iter), (*res));
      (*iter) ++;

      /* set delta */
      atimes (n, p, v);
      delta = - rho / ddot (n, r0, 1, v, 1);

      /* update s, t */
      daxpyz (n, delta, v, 1, r, 1, s, 1);
      atimes (n, s, t);

      /* set gamma */
      gamma = - ddot (n, s, 1, t, 1) / ddot (n, t, 1, t, 1);

      /* update r, x */
      daxpyz (n, gamma, t, 1, s, 1, r, 1);
      daxpy (n, gamma, s, 1, x, 1);
      daxpy (n, delta, p, 1, x, 1);

      /* set rho, beta */
      rho_old = rho;
      rho = ddot (n, r, 1, r0, 1);
      beta = rho / rho_old * delta / gamma;

      /* update p */
      daxpyz (n, gamma, v, 1, p, 1, tmp, 1);
      daxpyz (n, beta, tmp, 1, r, 1, p, 1);

      (*res) = dnrm2 (n, r, 1);
    }

  free (r);
  free (r0);
  free (p);
  free (v);
  free (s);
  free (t);
  free (tmp);
}
