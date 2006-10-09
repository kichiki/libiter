/* header file for cg.c --
 * Classical CG method -- Weiss' Algorithm 2
 * my implementation of Classical CG method
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: cg.h,v 2.1 2006/10/09 20:09:24 ichiki Exp $
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
#ifndef	_CG_H_
#define	_CG_H_

void
cg (int n, const double *b, double *x,
    double tol, int itmax,
    int *iter, double *res,
    void (*atimes) (int, const double *, double *, void *),
    void * user_data);


#endif /* !_CG_H_ */