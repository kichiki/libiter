/* Steepest Descent -- Weiss' Algorithm 1
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: steepest.h,v 2.3 2007/11/23 05:03:03 kichiki Exp $
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
#ifndef	_STEEPEST_H_
#define	_STEEPEST_H_


/* Steepest Descent -- Weiss' Algorithm 1
 * INPUT
 *  it_param : eps, max, debug, out are used.
 */
void
steepest (int n, const double *b, double *x,
	  void (*atimes) (int, const double *, double *, void *),
	  void *atimes_param,
	  struct iter * it_param);


#endif /* !_STEEPEST_H_ */
