/* header file of mygmres.c --
 * generalized minimum residual method
 * Copyright (C) 1998-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: gmres.h,v 2.8 2007/11/22 05:51:12 kichiki Exp $
 *
 * Reference :
 *   GMRES(m) : Y.Saad & M.H.Schultz, SIAM J.Sci.Stat.Comput.
 *   vol7 (1986) pp.856-869
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
#ifndef	_GMRES_H_
#define	_GMRES_H_

void
gmres_m (int n, const double *f, double *x,
	 void (*atimes) (int, const double *, double *, void *),
	 void *atimes_param,
	 struct iter *it_param);

void
gmres_m_ (int n, const double *f, double *x,
	  void (*atimes) (int, const double *, double *, void *),
	  void *atimes_param,
	  struct iter *it_param);

void
gmres_m_pc (int n, const double *f, double *x,
	    void (*atimes) (int, const double *, double *, void *),
	    void *atimes_param,
	    void (*inv) (int, const double *, double *, void *),
	    void *inv_param,
	    struct iter *it_param);

void
gmres (int n, const double *f, double *x,
       void (*atimes) (int, const double *, double *, void *),
       void *atimes_param,
       struct iter *it_param);


#endif /* !_GMRES_H_ */
