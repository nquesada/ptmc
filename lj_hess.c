/***********************************************************************
 *  lj_hess.c
 *
 *  (C) 2010 Nicolas Quesada, nquesada@pegasus.udea.edu.co   
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 **********************************************************************/

/*! \file lj_hess.c
   \brief See documentation of lj_hess.h - this file contains source code only, with function prototypes in lj_less.h
 */

#include"lj_hess.h"
#include "lj_params.h"

void
lj_hessian (int n, double *x, double *y, double *z, double *hess)
{

  int i, j;
  int tmp1, tmp2;
  int tn = 3 * n;
  int sn = 6 * n;
  int nn = 9 * n;
  double g, f, epssig;
  double xij, yij, zij;
  double rijs2, rijs4, rijs6, rijs12;

  for (j = 0; j < n; j++)
    {
      for (i = 0; i < j; i++)
	{
	  tmp1 = j * nn + 3 * i;
	  tmp2 = i * nn + 3 * j;
	  xij = x[i] - x[j];
	  yij = y[i] - y[j];
	  zij = z[i] - z[j];
	  rijs2 = lj_sigma2 (i, j) / (xij * xij + yij * yij + zij * zij);
	  rijs4 = rijs2 * rijs2;
	  rijs6 = rijs4 * rijs2;
	  rijs12 = rijs6 * rijs6;
	  epssig = lj_epsilon (i, j) / lj_sigma2 (i, j);
	  g = -epssig * 24.0 * rijs2 * (rijs6 - 2 * rijs12);
	  f =
	    -epssig * 4.0 * rijs4 * (168 * rijs12 -
				     48 * rijs6) / lj_sigma2 (i, j);

	  hess[tmp2] = hess[tmp1] = g + xij * xij * f;
	  hess[tmp1 + tn] = hess[tmp1 + 1] = hess[tmp2 + tn] =
	    hess[tmp2 + 1] = xij * yij * f;
	  hess[tmp1 + sn] = hess[tmp1 + 2] = hess[tmp2 + sn] =
	    hess[tmp2 + 2] = xij * zij * f;

	  hess[tmp1 + 1 + tn] = hess[tmp2 + 1 + tn] = g + yij * yij * f;
	  hess[tmp1 + 2 + tn] = hess[tmp1 + 1 + sn] = hess[tmp2 + 2 + tn] =
	    hess[tmp2 + 1 + sn] = yij * zij * f;

	  hess[tmp1 + 2 + sn] = hess[tmp2 + 2 + sn] = g + zij * zij * f;
	}
    }
  for (i = 0; i < n; i++)
    {
      tmp1 = i * nn + 3 * i;
      for (j = 0; j < i; j++)
	{
	  tmp2 = i * nn + 3 * j;
	  hess[tmp1] -= hess[tmp2];
	  hess[tmp1 + tn] -= hess[tmp2 + tn];
	  hess[tmp1 + sn] -= hess[tmp2 + sn];

	  hess[tmp1 + 1 + tn] -= hess[tmp2 + 1 + tn];
	  hess[tmp1 + 2 + tn] -= hess[tmp2 + 2 + tn];

	  hess[tmp1 + 2 + sn] -= hess[tmp2 + 2 + sn];
	}
      for (j = i + 1; j < n; j++)
	{
	  tmp2 = i * nn + 3 * j;
	  hess[tmp1] -= hess[tmp2];
	  hess[tmp1 + tn] -= hess[tmp2 + tn];
	  hess[tmp1 + sn] -= hess[tmp2 + sn];

	  hess[tmp1 + 1 + tn] -= hess[tmp2 + 1 + tn];
	  hess[tmp1 + 2 + tn] -= hess[tmp2 + 2 + tn];

	  hess[tmp1 + 2 + sn] -= hess[tmp2 + 2 + sn];
	}

      hess[tmp1 + 1] = hess[tmp1 + tn];
      hess[tmp1 + 2] = hess[tmp1 + sn];
      hess[tmp1 + 1 + sn] = hess[tmp1 + 2 + tn];
    }
}
