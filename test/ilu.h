/* header file for ilu.c --
 * Incomplete LU decomposition
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ilu.h,v 1.1 2007/12/01 18:03:42 kichiki Exp $
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
#ifndef	_ILU_H_
#define	_ILU_H_


/*
 * INPUT
 *  a[n*n] : given matrix
 * OUTPUT
 *  lu[n*n]:
 */
void
ILU (int n, const double *a,
     double *lu);

/*
 * INPUT
 *  a[n*n] : given matrix
 * OUTPUT
 *  lu[n*n]:
 */
void
LU_Gauss (int n, const double *a,
	  double *lu);

/* approximation of x = A^{-1}.b for preconditioning
 * INPUT
 *  b[n] : given vector
 *  user_data : (double *)lu[n*n]
 * OUTPUT
 *  x[n] := A^{-1}.b = U^{-1}.L^{-1}.b
 */
void
inv_ILU (int n, const double *b,
	 double *x, void *user_data);

/* diagonal preconditioner
 * INPUT
 *  b[n] : given vector
 *  inv_param : (double *)a[n*n], the coefficient matrix
 * OUTPUT
 *  x[n] := D^{-1}.b, where D is the diagonal part of A.
 */
void
inv_diag (int n, const double *b,
	  double *x, void *inv_param);

void
mul_LU (int n, const double *lu,
	double *a);


#endif /* !_ILU_H_ */
