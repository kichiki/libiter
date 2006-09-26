/* header file of bi-cgstab.c --
 * wrapper for iterative solver routines
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bi-cgstab.h,v 2.4 2006/09/26 17:10:06 ichiki Exp $
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
