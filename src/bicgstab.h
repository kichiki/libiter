/* BiCGSTAB - Weiss, Algorithm 12
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bicgstab.h,v 2.3 2007/11/23 05:17:21 kichiki Exp $
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
#ifndef	_BICGSTAB_H_
#define	_BICGSTAB_H_


/* Ref: Weiss, Algorithm 12 BiCGSTAB
 * INPUT
 *   n : dimension of the problem
 *   b [n] : r-h-s vector
 *   atimes (int n, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps  : criteria for |r^2|/|b^2|
 * OUTPUT
 *   x [n] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
void
bicgstab (int n, const double *b, double *x,
	  void (*atimes) (int, const double *, double *, void *),
	  void *atimes_param,
	  struct iter *it);


#endif /* !_BICGSTAB_H_ */
