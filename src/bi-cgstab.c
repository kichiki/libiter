/* wrapper for solver routines here
 * Copyright (C) 2001 Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 * $Id: bi-cgstab.c,v 1.6 2001/02/01 07:36:47 ichiki Exp $
 *
 * (solver routines themselves are originally written by martin h. gutknecht)
 */
#include <stdio.h> /* fprintf() */
#include <math.h> /* log10() */
#include <stdlib.h> /* malloc(), free() */

#include "bi-cgstab.h"


/** global variables **/
int STAB_it_max;
double STAB_log10_eps;

/* wrapper routine for solvers below
 * INPUT
 *   n : size of vectors v[] and f[] -- expected to be np * nelm for red-sym
 *   b [n] : given vector
 *   atimes (n, x, b) : routine to calc A.x and return b[]
 *   solver : solver routine to call
 *   (global) STAB_it_max : max # iterations
 *   (global) STAB_log10_eps : log10(eps), where eps is the accuracy
 * OUTPUT
 *   x [n] : solution
 */
void
solve_iter_stab (int n,
		 double *b, double *x,
		 void (*atimes) (int, double *, double *),
		 void (*solver) (int, double *, double *, int,
				 double, double, int *, double *,
				 void (*) (int, double *, double *)))
{
  extern int STAB_it_max;
  extern double STAB_log10_eps;

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
	  STAB_it_max, STAB_log10_eps,
	  hnor, &iter, &residual, atimes);

  fprintf (stderr, "# iter=%d res=%e\n", iter, residual);
}


/*             ===============rog.f==============
 *             problem given by martin h. gutknecht
 *             numerical method: bi-cgsta        method -- sta ()
 *             numerical method: bi-cgsta2       method -- st2 ()
 *             numerical method: gpbi-cg         method -- gpb ()
 *             ver. 1.0 aug. 08 1995  s. l. zhang
 *    at urus.slzhang.fort.iter.real.gpbcg.gutknecht-p(final.f)
 *
 * translated from fortran into C
 *   by Kengo Ichiki <ichiki@kona.jinkan.kyoto-u.ac.jp>
 */

/* bi-cgstab method
 *   m : dimension of the problem
 *   kend : max of iteration
 *   b [m] : r-h-s vector
 *   eps : log10 of cutoff
 *   hnor : log10 of norm of b []
 *   myatimes (int m, double *x, double *b) : calc matrix-vector product
 * OUTPUT
 *   x [m] : solution
 *   *iter : # of iteration
 *   *hg : log10(residual)
 */
void
sta (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *))
{
  int i;

  double tbs;

  double alpha;
  double beta;
  double omega;

  double rho0;
  double rho1;

  double *r0, *p, *q, *t, *r;
  /* dimension r0(m), p(-1:m+1), q(m), t(m), r(-1:m+1) */

  double hal;

  double pap;
  double rur;


  /* for myatimes () */
  double *tmp;


  /*r0 = my_d_malloc (m, "r0");
  p  = my_d_malloc (m, "p");
  q  = my_d_malloc (m, "q");
  t  = my_d_malloc (m, "t");
  r  = my_d_malloc (m, "r");
  tmp= my_d_malloc (m, "tmp");*/
  r0  = (double *) malloc (sizeof (double) * m);
  p   = (double *) malloc (sizeof (double) * m);
  q   = (double *) malloc (sizeof (double) * m);
  t   = (double *) malloc (sizeof (double) * m);
  r   = (double *) malloc (sizeof (double) * m);
  tmp = (double *) malloc (sizeof (double) * m);
  if (r0 == NULL
      || p == NULL
      || q == NULL
      || t == NULL
      || r == NULL
      || tmp == NULL)
    {
      fprintf (stderr, "malloc in sta ()");
      exit (1);
    }

  beta  = 0.0;
  omega = 0.0;
  rho0  = 0.0;
  myatimes (m, x, tmp);
  for (i=0; i<m; i++) /* 100 */
    {
      tbs = b [i] - tmp [i];
      r0 [i] = tbs;
      r [i] = tbs;
      p [i] = 0.0;
      q [i] = 0.0;
      rho0 += tbs * tbs;
    }/* 100 */
  for ((*iter) = 0; (*iter)<= kend; (*iter)++) /* 10 */
    {
      /* call resd(63,kk) */
      hal = 0.0;
      for (i=0; i<m; i++) /* 210 */
	{
          hal += r [i] * r [i];
          p [i] = r [i] + beta * (p [i] - omega * q [i]);
	} /* 210 */
      (*hg) = log10 (hal) / 2.0 - hnor;
      /* for check
      fprintf (stdout, "#STA %d %e\n", (* iter), (* hg)); */

      if((*hg) <= eps) goto end_sta;
      pap = 0.0;
      myatimes (m, p, tmp);
      for (i=0; i<m; i++) /* 220 */
	{
	  tbs = tmp [i];
          q [i] = tbs;
          pap  += tbs * r0 [i];
	}/* 220 */
      alpha = rho0 / pap;
      /* write(73,*) kk, ' alpha = ', alpha */
      for (i=0; i<m; i++) /* 230 */
	{
          r [i] -= alpha * q [i];
          x [i] += alpha * p [i];
	}/* 230 */
      pap = 0.0;
      rur = 0.0;
      myatimes (m, r, tmp);
      for (i=0; i<m; i++) /* 240 */
	{
	  tbs = tmp [i];
          t [i] = tbs;
          pap += tbs * r [i];
          rur += tbs * tbs;
	}/* 240 */
      omega = pap / rur;
      rho1 = 0.0;
      for (i=0; i<m; i++) /* 250 */
	{
          x [i] += omega * r [i];
	}/* 250 */
      for (i=0; i<m; i++) /* 260 */
	{
          tbs  = r [i] - omega * t [i];
          r [i] = tbs;
          rho1 += r0 [i] * tbs;
	}/* 260 */
      beta = (rho1 / rho0) * (alpha / omega);
      /* write(73,*) kk, ' beta = ', beta */
      rho0 = rho1;
    }/* 10 */

end_sta:
  free (r0);
  free (p);
  free (q);
  free (t);
  free (r);
  free (tmp);
}



/* bi-cgstab2 method
 *   m : dimension of the problem
 *   kend : max of iteration
 *   b [m] : r-h-s vector
 *   eps : log10 of cutoff
 *   hnor : log10 of norm of b []
 *   myatimes (int m, double *x, double *b) : calc matrix-vector product
 * OUTPUT
 *   x [m] : solution
 *   *iter : # of iteration
 *   *hg : log10(residual)
 */
void
st2 (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *))
{
  int i;

  double tbs;

  double alpha;
  double beta;
  double eta;
  double zeta;

  double rho0;
  double rho1;

  double *r0, *w, *q, *u, *z, *y;
  /*dimension r0(m), w(m), q(m), u(m), z(m), y(m) */
  double *r, *p;
  /*dimension r(-1:m+1), p(-1:m+1) */

  double taa, tbb, tcc, tdd, tee;

  double hal;

  double pap;


  /* for myatimes () */
  double *tmp;


  /*r0  = my_d_malloc (m, "r0");
  w   = my_d_malloc (m, "w");
  q   = my_d_malloc (m, "q");
  u   = my_d_malloc (m, "u");
  z   = my_d_malloc (m, "z");
  y   = my_d_malloc (m, "y");
  r   = my_d_malloc (m, "r");
  p   = my_d_malloc (m, "p");
  tmp = my_d_malloc (m, "tmp");*/
  r0  = (double *) malloc (sizeof (double) * m);
  w   = (double *) malloc (sizeof (double) * m);
  q   = (double *) malloc (sizeof (double) * m);
  u   = (double *) malloc (sizeof (double) * m);
  z   = (double *) malloc (sizeof (double) * m);
  y   = (double *) malloc (sizeof (double) * m);
  r   = (double *) malloc (sizeof (double) * m);
  p   = (double *) malloc (sizeof (double) * m);
  tmp = (double *) malloc (sizeof (double) * m);
  if (r0 == NULL
      || w == NULL
      || q == NULL
      || u == NULL
      || z == NULL
      || y == NULL
      || r == NULL
      || p == NULL
      || tmp == NULL)
    {
      fprintf (stderr, "malloc in st2 ()");
      exit (1);
    }

  beta = 0.0;
  rho0 = 0.0;
  myatimes (m, x, tmp);
  for (i=0; i<m; i++) /* 100 */
    {
      tbs = b [i] - tmp [i];
      r0 [i] =  tbs;
      r [i] =  tbs;
      y [i] = -tbs;
      p [i] =  0.0;
      u [i] =  0.0;
      /* my correction (maybe we should do this) */
      w [i] =  0.0;
      rho0  +=  tbs * tbs;
    }/* 100 */
  for ((*iter) = 0; (*iter)<= kend; (*iter)++) /* 10 */
    {
      /* call resd(64,kk)*/
      hal  = 0.0;
      for (i=0; i<m; i++) /* 210 */
	{
          hal += r [i] * r [i];
          p [i] = r [i] + beta * (p [i] - u [i]);
          u [i] = y [i] + beta * u [i];
	}/* 210 */
      (*hg) = log10 (hal) / 2.0 - hnor;
      /* for check
      fprintf (stdout, "#ST2 %d %e\n", (* iter), (* hg)); */

      if((*hg) <= eps) goto end_st2;
      pap = 0.0;
      myatimes (m, p, tmp);
      for (i=0; i<m; i++) /* 220 */
	{
	  tbs = tmp [i];
          pap += r0 [i] * tbs;
          q [i] = tbs;
	}/* 220 */
      alpha = rho0 / pap;
      /* write(74,*) kk, ' alpha = ', alpha */
      for (i=0; i<m; i++) /* 230 */
	{
          tbs  = alpha * q [i];
          y [i] += tbs - alpha * w [i];
          r [i] -= tbs;
          x [i] += alpha * p [i];
	}/* 230 */
      myatimes (m, r, w); /* w [] -> wn = A tn , tn <- r [] */
      /* 240 */
      taa = 0.0;
      tbb = 0.0;
      tcc = 0.0;
      tdd = 0.0;
      tee = 0.0;
      for (i=0; i<m; i++) /* 250 */
	{
          taa += w [i] * w [i];
          tbb += y [i] * y [i];
          tcc += w [i] * y [i];
          tdd += w [i] * r [i];
          tee += y [i] * r [i];
	}/* 250 */
      if (((*iter) / 2 * 2) == (*iter))
	{
	  zeta = tdd / taa;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (tbb * tdd - tcc * tee) / (taa * tbb - tcc * tcc);
	  eta  = (taa * tee - tdd * tcc) / (taa * tbb - tcc * tcc);
	}
      /* 500 */
      rho1 = 0.0;
      for (i=0; i<m; i++) /* 260 */
	{
          tbs  = zeta * r [i] + eta * (z [i] - alpha * u [i]);
          z [i] = tbs;
          x [i] += tbs;
	}/* 260 */
      for (i=0; i<m; i++) /* 270 */
	{
          tbs  = eta * y [i] + zeta * w [i];
          y [i] = tbs;
          r [i] -= tbs;
          u [i] = zeta * q [i] + eta * u [i];
          rho1 += r0 [i] * r [i];
	}/* 270 */
      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      /* write(74,*) kk, ' beta = ', beta */
      for (i=0; i<m; i++) /* 280 */
	{
          w [i] += beta * q [i];
	}/* 280 */
    }/* 10 */

end_st2:
  free (r0);
  free (w);
  free (q);
  free (u);
  free (z);
  free (y);
  free (r);
  free (p);
  free (tmp);
}

/* gpbi-cg method
 *   m : dimension of the problem
 *   kend : max of iteration
 *   b [m] : r-h-s vector
 *   eps : log10 of cutoff
 *   hnor : log10 of norm of b []
 *   myatimes (int m, double *x, double *b) : calc matrix-vector product
 * OUTPUT
 *   x [m] : solution
 *   *iter : # of iteration
 *   *hg : log10(residual)
 */
void
gpb (int m, double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, double *, double *))
{
  int i;

  double tbs;

  double alpha;
  double beta;
  double eta;
  double zeta;

  double rho0;
  double rho1;

  double *r0, *w, *q, *u, *z, *y;
  /*dimension r0(m), w(m), q(m), u(m), z(m), y(m) */
  double *r, *p;
  /*dimension r(-1:m+1), p(-1:m+1) */

  double taa, tbb, tcc, tdd, tee;

  double hal;

  double pap;


  /* for myatimes () */
  double *tmp;


  /*r0 = my_d_malloc (m, "r0");
  w  = my_d_malloc (m, "w");
  q  = my_d_malloc (m, "q");
  u  = my_d_malloc (m, "u");
  z  = my_d_malloc (m, "z");
  y  = my_d_malloc (m, "y");
  r  = my_d_malloc (m, "r");
  p  = my_d_malloc (m, "p");
  tmp= my_d_malloc (m, "tmp");*/
  r0  = (double *) malloc (sizeof (double) * m);
  w   = (double *) malloc (sizeof (double) * m);
  q   = (double *) malloc (sizeof (double) * m);
  u   = (double *) malloc (sizeof (double) * m);
  z   = (double *) malloc (sizeof (double) * m);
  y   = (double *) malloc (sizeof (double) * m);
  r   = (double *) malloc (sizeof (double) * m);
  p   = (double *) malloc (sizeof (double) * m);
  tmp = (double *) malloc (sizeof (double) * m);
  if (r0 == NULL
      || w == NULL
      || q == NULL
      || u == NULL
      || z == NULL
      || y == NULL
      || r == NULL
      || p == NULL
      || tmp == NULL)
    {
      fprintf (stderr, "malloc in gpb ()");
      exit (1);
    }

  beta = 0.0;
  rho0 = 0.0;
  myatimes (m, x, tmp);
  for (i=0; i<m; i++) /* 100 */
    {
      tbs = b [i] - tmp [i];
      r0 [i] =  tbs;
      r [i] =  tbs;
      y [i] = -tbs;
      p [i] =  0.0;
      u [i] =  0.0;
      /* my correction (maybe we should do this) */
      w [i] =  0.0;
      rho0 +=  tbs * tbs;
    }/* 100 */
  for ((*iter) = 0; (*iter)<= kend; (*iter)++) /* 10 */
    {
      /* call resd(66,kk) */
      hal  = 0.0;
      for (i=0; i<m; i++) /* 210 */
	{
          hal += r [i] * r [i];
          p [i] = r [i] + beta * (p [i] - u [i]);
          u [i] = y [i] + beta * u [i];
	}/* 210 */
      (*hg) = log10 (hal) / 2.0 - hnor;
      /* for check
      fprintf (stdout, "#GPB %d %e\n", (* iter), (* hg)); */

      if((*hg) <= eps) goto end_gpb;
      pap = 0.0;
      myatimes (m, p, tmp); /* p [] -> pn */
      for (i=0; i<m; i++) /* 220 */
	{
	  tbs = tmp [i];
          pap += r0 [i] * tbs;
          q [i] = tbs; /* q [] -> A pn */
	}/* 220 */
      alpha = rho0 / pap;
      /* write(76,*) kk, ' alpha = ', alpha */
      for (i=0; i<m; i++) /* 230 */
	{
          tbs  = alpha * q [i];
          y [i] += tbs - alpha * w [i];
          r [i] -= tbs;
          x [i] += alpha * p [i];
	}/* 230 */
      myatimes (m, r, w); /* w [] -> wn = A tn , tn <- r [] */
      /* 240 */
      taa = 0.0; /* taa -> (A tn, A tn) */
      tbb = 0.0; /* tbb -> (  yn,   yn) */
      tcc = 0.0; /* tcc -> (A tn,   yn) */
      tdd = 0.0; /* tdd -> (A tn,   tn) */
      tee = 0.0; /* tee -> (  yn,   tn) */
      for (i=0; i<m; i++) /* 250 */
	{
          taa += w [i] * w [i];
          tbb += y [i] * y [i];
          tcc += w [i] * y [i];
          tdd += w [i] * r [i];
          tee += y [i] * r [i];
	}/* 250 */
      if ((*iter) == 0)
	{
	  zeta = tdd / taa;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (tbb * tdd - tcc * tee) / (taa * tbb - tcc * tcc);
	  eta  = (taa * tee - tdd * tcc) / (taa * tbb - tcc * tcc);
	}
      /* 500 */
      rho1 = 0.0; /* rho1 -> (r*0, rn+1) */
      for (i=0; i<m; i++) /* 260 */
	{
          tbs  = zeta * r [i] + eta * (z [i] - alpha * u [i]);
          z [i] = tbs;
          x [i] += tbs;
	}/* 260 */
      for (i=0; i<m; i++) /* 270 */
	{
          tbs  = eta * y [i] + zeta * w [i];
          y [i] = tbs;
          r [i] -= tbs;
          u [i] = zeta * q [i] + eta * u [i];
          rho1 += r0 [i] * r [i];
	}/* 270 */
      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      /* write(76,*) kk, ' beta = ', beta */
      for (i=0; i<m; i++) /* 280 */
	{
          w [i] += beta * q [i]; /* wn += beta A pn */
	}/* 280 */
    }/* 10 */

end_gpb:
  free (r0);
  free (w);
  free (q);
  free (u);
  free (z);
  free (y);
  free (r);
  free (p);
  free (tmp);
}

/* gpbi-cg method with check (printing residual in iterations)
 */
void
gpb_chk (int m, double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, double *, double *))
{
  int i;

  double tbs;

  double alpha;
  double beta;
  double eta;
  double zeta;

  double rho0;
  double rho1;

  double *r0, *w, *q, *u, *z, *y;
  double *r, *p;

  double taa, tbb, tcc, tdd, tee;

  double hal;

  double pap;


  /* for myatimes () */
  double *tmp;


  r0  = (double *) malloc (sizeof (double) * m);
  w   = (double *) malloc (sizeof (double) * m);
  q   = (double *) malloc (sizeof (double) * m);
  u   = (double *) malloc (sizeof (double) * m);
  z   = (double *) malloc (sizeof (double) * m);
  y   = (double *) malloc (sizeof (double) * m);
  r   = (double *) malloc (sizeof (double) * m);
  p   = (double *) malloc (sizeof (double) * m);
  tmp = (double *) malloc (sizeof (double) * m);
  if (r0 == NULL
      || w == NULL
      || q == NULL
      || u == NULL
      || z == NULL
      || y == NULL
      || r == NULL
      || p == NULL
      || tmp == NULL)
    {
      fprintf (stderr, "malloc in gpb ()");
      exit (1);
    }

  beta = 0.0;
  rho0 = 0.0;
  myatimes (m, x, tmp);
  for (i=0; i<m; i++) /* 100 */
    {
      tbs = b [i] - tmp [i];
      r0 [i] =  tbs;
      r [i] =  tbs;
      y [i] = -tbs;
      p [i] =  0.0;
      u [i] =  0.0;
      /* my correction (maybe we should do this) */
      w [i] =  0.0;
      rho0 +=  tbs * tbs;
    }/* 100 */
  for ((*iter) = 0; (*iter)<= kend; (*iter)++) /* 10 */
    {
      /* call resd(66,kk) */
      hal  = 0.0;
      for (i=0; i<m; i++) /* 210 */
	{
          hal += r [i] * r [i];
          p [i] = r [i] + beta * (p [i] - u [i]);
          u [i] = y [i] + beta * u [i];
	}/* 210 */
      (*hg) = log10 (hal) / 2.0 - hnor;
      /* for check */
      fprintf (stdout, "#GPB %d %e\n", (* iter), (* hg));

      if((*hg) <= eps) goto end_gpb_chk;
      pap = 0.0;
      myatimes (m, p, tmp); /* p [] -> pn */
      for (i=0; i<m; i++) /* 220 */
	{
	  tbs = tmp [i];
          pap += r0 [i] * tbs;
          q [i] = tbs; /* q [] -> A pn */
	}/* 220 */
      alpha = rho0 / pap;
      /* write(76,*) kk, ' alpha = ', alpha */
      for (i=0; i<m; i++) /* 230 */
	{
          tbs  = alpha * q [i];
          y [i] += tbs - alpha * w [i];
          r [i] -= tbs;
          x [i] += alpha * p [i];
	}/* 230 */
      myatimes (m, r, w); /* w [] -> wn = A tn , tn <- r [] */
      /* 240 */
      taa = 0.0; /* taa -> (A tn, A tn) */
      tbb = 0.0; /* tbb -> (  yn,   yn) */
      tcc = 0.0; /* tcc -> (A tn,   yn) */
      tdd = 0.0; /* tdd -> (A tn,   tn) */
      tee = 0.0; /* tee -> (  yn,   tn) */
      for (i=0; i<m; i++) /* 250 */
	{
          taa += w [i] * w [i];
          tbb += y [i] * y [i];
          tcc += w [i] * y [i];
          tdd += w [i] * r [i];
          tee += y [i] * r [i];
	}/* 250 */
      if ((*iter) == 0)
	{
	  zeta = tdd / taa;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (tbb * tdd - tcc * tee) / (taa * tbb - tcc * tcc);
	  eta  = (taa * tee - tdd * tcc) / (taa * tbb - tcc * tcc);
	}
      /* 500 */
      rho1 = 0.0; /* rho1 -> (r*0, rn+1) */
      for (i=0; i<m; i++) /* 260 */
	{
          tbs  = zeta * r [i] + eta * (z [i] - alpha * u [i]);
          z [i] = tbs;
          x [i] += tbs;
	}/* 260 */
      for (i=0; i<m; i++) /* 270 */
	{
          tbs  = eta * y [i] + zeta * w [i];
          y [i] = tbs;
          r [i] -= tbs;
          u [i] = zeta * q [i] + eta * u [i];
          rho1 += r0 [i] * r [i];
	}/* 270 */
      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      /* write(76,*) kk, ' beta = ', beta */
      for (i=0; i<m; i++) /* 280 */
	{
          w [i] += beta * q [i]; /* wn += beta A pn */
	}/* 280 */
    }/* 10 */

end_gpb_chk:
  free (r0);
  free (w);
  free (q);
  free (u);
  free (z);
  free (y);
  free (r);
  free (p);
  free (tmp);
}
