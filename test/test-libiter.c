/* test code for libiter solvers
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: test-libiter.c,v 1.3 2007/11/23 04:44:09 kichiki Exp $
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "check-iter-gen.h"
#include "toeplitz.h"


/* main program */
int
main (int argc, char** argv)
{
  int check = 0;

  check += check_3_all (2000, 20, 1.0e-7, // max, restart, eps
			1, 1.0e-6); // for check

  check += check_symmetric_all (500, // dimension
				10000, 20, 1.0e-12, // max, restart, eps
				1, 3.8e-9); // for check

  check += Toeplitz_check_all (1000, 1.5, // n, gamma
			       10000, 20, 1.0e-7, // max, restart, eps
			       1, 3.5e-6); // for check

  fprintf (stdout,
	   "==================================================\n"
	   "TOTAL ERRORS : %d\n", check);
  if (check == 0)
    {
      fprintf (stdout, "Conglaturation!! ALL TESTS PASSED\n");
    }
  else
    {
      fprintf (stdout, "sorry. some test(s) FAILED\n");
    }

  return 0;
}
