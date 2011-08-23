/*! \file example_hess.c
 * \brief This file provides and example of the usage of the functions
 * that calculate the Hessian of the Lennard Jones potential energy 
 * surface and the eigenvalue routines provided by GSL.
 * 
 * Input:		To execute the program the following arguments must be 
 * 				passed in the command line:
 * 				1. Number of atoms.
 * 				2. Name of the file with the geometry.
 * 
 * 				The geometry files are plain text files with 3
 * 				columns for the x-y-z coordinates of the n atoms.
 * 
 * 				example:
 * 				./example_grad.out 13 geometry.xyz 
 * 				
 * Output:		The program will print to the standard output the 
 * 				eigenvalues of the Hessian matrix and will also
 * 				print in the last line the value of the Geometric Mean
 * 				of such values.
 * 				It is assumed that all the eigenvalues of the Hessian
 * 				matrix are positive up to plus or minus 0.0001, if
 * 				there negative eigenvalues the Geometric Mean
 * 				might be complex and an error will occur in its 
 * 				calculation.
 * 
 * 				If the reader wants to know more about the 
 * 				eigenvalue routines should consult the GSL Reference
 * 				Manual.
 */ 


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_eigen.h>

#include"lj_hess.h"

int
main (int argc, char *argv[])
{
  double *x, *y, *z, *data;
  int i,gg;
  int n = atoi (argv[1]);
  FILE *pf;

  x = (double *) malloc (n * sizeof (double));
  y = (double *) malloc (n * sizeof (double));
  z = (double *) malloc (n * sizeof (double));
  data = (double *) malloc (9 * n * n * sizeof (double));

  pf = fopen (argv[2], "r");

  for (i = 0; i < n; i++)
    {
      gg = fscanf (pf, "%lf %lf %lf", &x[i], &y[i], &z[i]);
    }

  lj_hessian (n, x, y, z, data);

  gsl_matrix_view m = gsl_matrix_view_array (data, 3 * n, 3 * n);

  gsl_vector *eval = gsl_vector_alloc (3 * n);
  gsl_matrix *evec = gsl_matrix_alloc (3 * n, 3 * n);

  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (3 * n);

  gsl_eigen_symmv (&m.matrix, eval, evec, w);

  gsl_eigen_symmv_free (w);

  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

  {
    int i;
    double prod = 0;
    int count = 0;
    for (i = 0; i < 3 * n; i++)
      {
	double eval_i = gsl_vector_get (eval, i);
	gsl_vector_view evec_i = gsl_matrix_column (evec, i);

	printf ("%g\n", eval_i);
	if (eval_i > 0.0001)//Only eigenvalues gretaer than 0.0001
	  {
	    prod += log (eval_i);
	    count++;

	  }
      }
    fprintf (stdout, "\n GM=%lf \n", exp ( prod / count));
  }

  gsl_vector_free (eval);
  gsl_matrix_free (evec);

  free (x);
  free (y);
  free (z);

  return 0;
}
