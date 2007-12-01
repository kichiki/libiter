/* Incomplete LU decomposition
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ilu.c,v 1.1 2007/12/01 18:02:42 kichiki Exp $
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
#include <math.h>
#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for exit() */
#include "memory-check.h" // macro CHECK_MALLOC


/*
 * INPUT
 *  a[n*n] : given matrix
 * OUTPUT
 *  lu[n*n]:
 */
void
ILU (int n, const double *a,
     double *lu)
{
  int k;
  for (k = 0; k < n*n; k ++)
    {
      lu[k] = a[k];
    }

  for (k = 0; k < (n-1); k ++)
    {
      if (lu[k*n+k] == 0.0)
	{
	  fprintf (stderr, "ILU : breakdown a[%d,%d] = 0\n", k, k);
	}
      double d = 1.0 / lu[k*n+k];
      int i;
      for (i = (k+1); i < n; i ++)
	{
	  double e = d * lu[i*n+k];
	  lu[i*n+k] = e;
	  int j;
	  for (j = (k+1); j < n; j ++)
	    {
	      lu[i*n+j] -= e * lu[k*n+j];
	    }
	}
    }
}

/*
 * INPUT
 *  a[n*n] : given matrix
 * OUTPUT
 *  lu[n*n]:
 */
void
LU_Gauss (int n, const double *a,
	  double *lu)
{
  int k;
  for (k = 0; k < n*n; k ++)
    {
      lu[k] = a[k];
    }

  for (k = 0; k < (n-1); k ++)
    {
      if (lu[k*n+k] == 0.0)
	{
	  fprintf (stderr, "LU_GAUSS : breakdown a[%d,%d] = 0\n", k, k);
	}
      int i, j;
      for (i = (k+1); i < n; i ++)
	{
	  lu[i*n+k] = lu[i*n+k] / lu[k*n+k];
	  for (j = (k+1); j < n; j ++)
	    {
	      lu[i*n+j] -= lu[i*n+k] * lu[k*n+j];
	    }
	}
    }
}

/* approximation of x = A^{-1}.b for preconditioning
 * INPUT
 *  b[n] : given vector
 *  user_data : (double *)lu[n*n]
 * OUTPUT
 *  x[n] := A^{-1}.b = U^{-1}.L^{-1}.b
 */
void
inv_ILU (int n, const double *b,
	 double *x, void *user_data)
{
  double *lu = (double *)user_data;

  double *z = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (z, "inv_ILU");

  int i, j;
  // solve z from L.z = b
  for (i = 0; i < n; i ++)
    {
      z[i] = b[i]; // L[i*n+i] = 1
      for (j = 0; j < i; j ++)
	{
	  // j runs from  0, ..., i-1.
	  z[i] -= lu[i*n+j] * z[j];
	}
    }

  // solve x from U.x = z
  for (i = n-1; i >= 0; i --)
    {
      // i runs from n-1, n-2, ..., 0.
      x[i] = z[i];
      for (j = i + 1; j < n; j ++)
	{
	  // j runs from  i+1, ..., n-1
	  x[i] -= lu[i*n+j] * x[j];
	}
      x[i] /= lu[i*n+i];
    }

  free (z);
}

/* diagonal preconditioner
 * INPUT
 *  b[n] : given vector
 *  inv_param : (double *)a[n*n], the coefficient matrix
 * OUTPUT
 *  x[n] := D^{-1}.b, where D is the diagonal part of A.
 */
void
inv_diag (int n, const double *b,
	  double *x, void *inv_param)
{
  double *a = (double *)inv_param;
  int i;
  for (i = 0; i < n; i ++)
    {
      if (a[i*n+i] != 0.0) x[i] = b[i] / a[i*n+i];
      else                 x[i] = b[i]; // just for filling
    }
}

void
mul_LU (int n, const double *lu,
	double *a)
{
  // calc aa = L . U
  int i;
  for (i = 0; i < n; i ++)
    {
      int j;
      for (j = 0; j < n; j ++)
	{
	  a[i*n+j] = 0.0;
	  int k;
	  for (k = 0; k < n; k ++)
	    {
	      if (i < k) continue; // l_{ik} = 0
	      if (k > j) continue; // u_{kj} = 0
	      if (i == k)
		{
		  // l_{ik} = 1
		  a[i*n+j] += lu[k*n+j];
		}
	      else
		{
		  a[i*n+j] += lu[i*n+k] * lu[k*n+j];
		}
	    }
	}
    }
}
