/* header file for check-gmres-m-pc.c --
 * test code for gmres_m_pc() in gmres.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-gmres-m-pc.h,v 1.1 2007/12/01 18:06:19 kichiki Exp $
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
#ifndef	_CHECK_GMRES_M_PC_H_
#define	_CHECK_GMRES_M_PC_H_


int
check_gmres_m_pc (int n,
		  int it_max, int it_restart, double it_eps,
		  int verbose, double tiny);


#endif /* !_CHECK_GMRES_M_PC_H_ */
