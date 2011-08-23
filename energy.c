/***********************************************************************
 *  energy.c
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

/*! \file energy.c
   \brief See documentation of energy.h - this file contains source code
    only, with function prototypes in energy.h
 */


#include"energy.h"
#include"lj_params.h"



double
interaction_energy (double *x, double *y, double *z, int i, int n)
{
  int j;
  double U = 0;
  double xij, yij, zij, rijm2, rijm6;

  for (j = 0; j < i; j++)
    {
      xij = x[j] - x[i];
      yij = y[j] - y[i];
      zij = z[j] - z[i];
      rijm2 = lj_sigma2 (i, j) / (xij * xij + yij * yij + zij * zij);
      rijm6 = rijm2 * rijm2 * rijm2;
      U += lj_epsilon (i, j) * rijm6 * (rijm6 - 1);
    }

  for (j = i + 1; j < n; j++)
    {
      xij = x[j] - x[i];
      yij = y[j] - y[i];
      zij = z[j] - z[i];
      rijm2 = lj_sigma2 (i, j) / (xij * xij + yij * yij + zij * zij);
      rijm6 = rijm2 * rijm2 * rijm2;
      U += lj_epsilon (i, j) * rijm6 * (rijm6 - 1);
    }
  U = 4 * U;
  return U;
}



double
cluster_energy (double *x, double *y, double *z, int n)
{
  int i, j;
  double U = 0;
  double xij, yij, zij, rijm2, rijm6;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < i; j++)
	{
	  xij = x[j] - x[i];
	  yij = y[j] - y[i];
	  zij = z[j] - z[i];
	  rijm2 = lj_sigma2 (i, j) / (xij * xij + yij * yij + zij * zij);
	  rijm6 = rijm2 * rijm2 * rijm2;
	  U += lj_epsilon (i, j) * rijm6 * (rijm6 - 1);
	}
    }
  U = 4 * U;
  return U;


}
