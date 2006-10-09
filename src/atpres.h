/* header file for atpres.c --
 * ATPRES -- Weiss, Algorithm 5
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: atpres.h,v 2.2 2006/10/09 22:01:37 ichiki Exp $
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
#ifndef	_ATPRES_H_
#define	_ATPRES_H_


void
atpres (int n, const double *b, double *x,
	double eps, int maxiter,
	int *iter, double *res,
	void (*atimes) (int, const double *, double *, void *),
	void (*atimes_t) (int, const double *, double *, void *),
	void * user_data);


#endif /* !_ATPRES_H_ */
