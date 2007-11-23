/* header file for toeplitz.c --
 * Toeplitz matrix related code to check libiter solvers
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: toeplitz.h,v 1.1 2007/11/23 04:42:38 kichiki Exp $
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
#ifndef	_TOEPLITZ_H_
#define	_TOEPLITZ_H_


/* make Toeplitz matrix
 * INPUT
 *  gamma
 * OUTPUT
 *  mat[n*n]
 */
void
Toeplitz_make_matrix (int n, double gamma, double *mat);

/* calc atimes b = A.x, where A is the Toeplitz matrix
 * INPUT
 *  x[n] : given vector
 *  *user_data = (double) gamma
 * OUTPUT
 *  b[n] := A.x
 */
void
Toeplitz_atimes (int n, const double *x, double *b, void *user_data);

/* calc atimes b = A^T.x, where A is the Toeplitz matrix
 * INPUT
 *  x[n] : given vector
 *  *user_data = (double) gamma
 * OUTPUT
 *  b[n] := A^T.x
 */
void
Toeplitz_atimes_t (int n, const double *x, double *b, void *user_data);

int
Toeplitz_check_all (int n, double gamma,
		    int it_max, int it_restart, double it_eps,
		    int verbose, double tiny);


#endif /* !_TOEPLITZ_H_ */
