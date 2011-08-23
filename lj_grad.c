/***********************************************************************
 *  lj_grad.c
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

/*! \file lj_grad.c
   \brief See documentation of lj_grad.h - this file contains source code only, with function prototypes in lj_grad.h
 */


#include<gsl/gsl_vector.h>

#include"lj_params.h"
#include"lj_grad.h"

double
lj_f (const gsl_vector * r, void *params)
{
  int i, j;
  int n = ((int *) params)[0];
  double x, y, z, r2, r2s, rm6, uu = 0;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < i; j++)
	{
	  x = gsl_vector_get (r, 3 * i) - gsl_vector_get (r, 3 * j);
	  y = gsl_vector_get (r, 3 * i + 1) - gsl_vector_get (r, 3 * j + 1);
	  z = gsl_vector_get (r, 3 * i + 2) - gsl_vector_get (r, 3 * j + 2);
	  r2 = (x * x + y * y + z * z);
	  r2s = lj_sigma2 (i, j) / r2;
	  rm6 = r2s * r2s * r2s;
	  uu += lj_epsilon (i, j) * 4 * rm6 * (rm6 - 1);
	}
    }
  return uu;
}


void
lj_df (const gsl_vector * r, void *params, gsl_vector * g)
{
  int i, j;
  int n = ((int *) params)[0];
  double *dr;
  dr = g->data;
  for (i = 0; i < 3 * n; i++)
    {
      dr[i] = 0;
    }
  double x, y, z, r2, rm6, r2s, et;

  for (i = 0; i < n; i++)
    {
      for (j = 0; j < i; j++)
	{
	  x = gsl_vector_get (r, 3 * i) - gsl_vector_get (r, 3 * j);
	  y = gsl_vector_get (r, 3 * i + 1) - gsl_vector_get (r, 3 * j + 1);
	  z = gsl_vector_get (r, 3 * i + 2) - gsl_vector_get (r, 3 * j + 2);
	  r2 = 1 / (x * x + y * y + z * z);
	  r2s = lj_sigma2 (i, j) * r2;
	  rm6 = (r2s * r2s * r2s);
	  et = lj_epsilon (i, j) * 24 * rm6 * (1 - 2 * rm6) * r2;
	  dr[3 * i] += et * x;
	  dr[3 * i + 1] += et * y;
	  dr[3 * i + 2] += et * z;
	}
      for (j = i + 1; j < n; j++)
	{
	  x = gsl_vector_get (r, 3 * i) - gsl_vector_get (r, 3 * j);
	  y = gsl_vector_get (r, 3 * i + 1) - gsl_vector_get (r, 3 * j + 1);
	  z = gsl_vector_get (r, 3 * i + 2) - gsl_vector_get (r, 3 * j + 2);

	  r2 = 1 / (x * x + y * y + z * z);
	  r2s = lj_sigma2 (i, j) * r2;
	  rm6 = (r2s * r2s * r2s);
	  et = lj_epsilon (i, j) * 24 * rm6 * (1 - 2 * rm6) * r2;
	  dr[3 * i] += et * x;
	  dr[3 * i + 1] += et * y;
	  dr[3 * i + 2] += et * z;
	}
    }
}


void
lj_fdf (const gsl_vector * r, void *params, double *f, gsl_vector * g)
{
  int i, j;
  int n = ((int *) params)[0];
  double *dr;
  dr = g->data;
  for (i = 0; i < 3 * n; i++)
    {
      dr[i] = 0;
    }
  double uu = 0;
  double x, y, z, r2, rm6, r2s, et;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < i; j++)
	{
	  x = gsl_vector_get (r, 3 * i) - gsl_vector_get (r, 3 * j);
	  y = gsl_vector_get (r, 3 * i + 1) - gsl_vector_get (r, 3 * j + 1);
	  z = gsl_vector_get (r, 3 * i + 2) - gsl_vector_get (r, 3 * j + 2);
	  r2 = 1 / (x * x + y * y + z * z);
	  r2s = lj_sigma2 (i, j) * r2;
	  rm6 = (r2s * r2s * r2s);
	  et = lj_epsilon (i, j) * 24 * rm6 * (1 - 2 * rm6) * r2;
	  dr[3 * i] += et * x;
	  dr[3 * i + 1] += et * y;
	  dr[3 * i + 2] += et * z;
	  uu += lj_epsilon (i, j) * 4 * rm6 * (rm6 - 1);
	}
      for (j = i + 1; j < n; j++)
	{
	  x = gsl_vector_get (r, 3 * i) - gsl_vector_get (r, 3 * j);
	  y = gsl_vector_get (r, 3 * i + 1) - gsl_vector_get (r, 3 * j + 1);
	  z = gsl_vector_get (r, 3 * i + 2) - gsl_vector_get (r, 3 * j + 2);
	  r2 = 1 / (x * x + y * y + z * z);
	  r2s = lj_sigma2 (i, j) * r2;
	  rm6 = (r2s * r2s * r2s);
	  et = lj_epsilon (i, j) * 24 * rm6 * (1 - 2 * rm6) * r2;
	  dr[3 * i] += et * x;
	  dr[3 * i + 1] += et * y;
	  dr[3 * i + 2] += et * z;
	}
    }
  *f = uu;
}
