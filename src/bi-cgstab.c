/* wrapper for iterative solver routines
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bi-cgstab.c,v 2.6 2006/09/28 04:21:36 kichiki Exp $
 *
 * solver routines are translated into C by K.I. from fortran code
 * originally written by martin h. gutknecht
 *             ===============rog.f==============
 *             problem given by martin h. gutknecht
 *             numerical method: bi-cgsta        method -- sta ()
 *             numerical method: bi-cgsta2       method -- st2 ()
 *             numerical method: gpbi-cg         method -- gpb ()
 *             ver. 1.0 aug. 08 1995  s. l. zhang
 *    at urus.slzhang.fort.iter.real.gpbcg.gutknecht-p(final.f)
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
#include <stdio.h> /* fprintf() */
#include <math.h> /* log10() */
#include <stdlib.h> /* malloc(), free() */

#include "bi-cgstab.h"


/* bi-cgstab method
 *   m : dimension of the problem
 *   kend : max of iteration
 *   b [m] : r-h-s vector
 *   eps : log10 of cutoff
 *   hnor : log10 of norm of b []
 *   myatimes (int m, double *x, double *b) : calc matrix-vector product
 *   user_data : pointer to be passed to atimes routines
 * OUTPUT
 *   x [m] : solution
 *   *iter : # of iteration
 *   *hg : log10(residual)
 */
void
sta (int m, const double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, const double *, double *, void *),
     void * user_data)
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
  myatimes (m, x, tmp, user_data);
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
      myatimes (m, p, tmp, user_data);
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
      myatimes (m, r, tmp, user_data);
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
 *   user_data : pointer to be passed to atimes routines
 * OUTPUT
 *   x [m] : solution
 *   *iter : # of iteration
 *   *hg : log10(residual)
 */
void
sta2 (int m, const double *b, double *x, int kend,
      double eps, double hnor,
      int *iter, double *hg,
      void (*myatimes) (int, const double *, double *, void *),
      void * user_data)
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
      fprintf (stderr, "malloc in st2 ()");
      exit (1);
    }

  beta = 0.0;
  rho0 = 0.0;
  myatimes (m, x, tmp, user_data);
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
      myatimes (m, p, tmp, user_data);
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
      myatimes (m, r, w, user_data); /* w [] -> wn = A tn , tn <- r [] */
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

/* bi-cgstab2 method with check (printing residual in iterations)
 */
void
st2_chk (int m, const double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, const double *, double *, void *),
	 void * user_data)
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
      fprintf (stderr, "malloc in st2 ()");
      exit (1);
    }

  beta = 0.0;
  rho0 = 0.0;
  myatimes (m, x, tmp, user_data);
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
      /* for check */
      fprintf (stdout, "#ST2 %d %e\n", (* iter), (* hg));

      if((*hg) <= eps) goto end_st2_chk;
      pap = 0.0;
      myatimes (m, p, tmp, user_data);
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
      myatimes (m, r, w, user_data); /* w [] -> wn = A tn , tn <- r [] */
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

end_st2_chk:
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
 *   user_data : pointer to be passed to atimes routines
 * OUTPUT
 *   x [m] : solution
 *   *iter : # of iteration
 *   *hg : log10(residual)
 */
void
gpb (int m, const double *b, double *x, int kend,
     double eps, double hnor,
     int *iter, double *hg,
     void (*myatimes) (int, const double *, double *, void *),
     void * user_data)
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
  myatimes (m, x, tmp, user_data);
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
      myatimes (m, p, tmp, user_data); /* p [] -> pn */
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
      myatimes (m, r, w, user_data); /* w [] -> wn = A tn , tn <- r [] */
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
gpb_chk (int m, const double *b, double *x, int kend,
	 double eps, double hnor,
	 int *iter, double *hg,
	 void (*myatimes) (int, const double *, double *, void *),
	 void * user_data)
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
  myatimes (m, x, tmp, user_data);
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
      myatimes (m, p, tmp, user_data); /* p [] -> pn */
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
      myatimes (m, r, w, user_data); /* w [] -> wn = A tn , tn <- r [] */
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
