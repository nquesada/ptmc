
/*! \file example_grad.c
 * \brief This file provides and example of the usage of the functions
 * that calculate the gradient of the Lennard Jones potential energy 
 * surface and the local optimization routines provided by GSL.
 * 
 * Input:		To execute the program the following arguments must be 
 * 				passed in the command line:
 * 				1. Number of atoms.
 * 				2. Name of the file with the initial geometry.
 * 				3. Name that will be given to the file with the final
 * 				geometry.
 * 
 * 				The geometry files are plain text files with 3
 * 				columns for the x-y-z coordinates of the n atoms.
 * 
 * 				example:
 * 				./example_grad.out 13 initial.xyz final.xyz
 * 				The program will locally optimize the geometry in file
 * 				initial.xyz by reading the first 13 lines (i.e.
 * 				the cluster is assumed to have 13 atoms) and will 
 * 				write in final.xyz the optimized geometry..
 * 
 * Output:		As it was mentioned before the program will produce
 * 				a file containing the optimized geometry. The geomtric
 * 				center of the geomtry will be translated to the
 * 				coordinate origin. Also it
 * 				will return to the standard output a message giving
 * 				the value of the initial energy and the norm of the 
 * 				gradient for that geometry and the energy and norm
 * 				of the gradient for the optimized geometry.
 * 				If the reader wants to know more about the local
 * 				optimization routines should consult the GSL Reference
 * 				Manual.
 */ 

//C Headers
#include<stdio.h>
#include<stdlib.h>
#include<math.h>


//GSL Minimization Headers
#include<gsl/gsl_multimin.h>	//Multidimensional minimization


#include"lj_grad.h"


int
main (int argc, char *argv[])
{
  int i, j;
  int n = atoi (argv[1]);
  double x, y, z;
  gsl_vector *r = gsl_vector_alloc (3 * n);
  gsl_vector *dr = gsl_vector_alloc (3 * n);
  double uu = 0;
  double est;
  double tmp1, tmp2;
  FILE *in, *out;

  //Reading input
  in = fopen (argv[2], "r");
  out = fopen (argv[3], "w");
  for (i = 0; i < n; i++)
    {
      j = fscanf (in, "%lf %lf %lf", &x, &y, &z);
      gsl_vector_set (r, 3 * i, x);
      gsl_vector_set (r, 3 * i + 1, y);
      gsl_vector_set (r, 3 * i + 2, z);
    }
  fclose (in);

  //Initial energies and forces before optimization
  lj_fdf (r, &n, &uu, dr);
  tmp1 = tmp2 = 0.0;
  for (i = 0; i < 3 * n; i++)
    {
      tmp1 = gsl_vector_get (dr, i);
      tmp2 += tmp1 * tmp1;
    }
  tmp2 = sqrt (tmp2);
  fprintf (stdout,
	   "\n \nInitial Energy in file %s was E=%lf and the value of the norm of the gradient was |df|=%lf\n",
	   argv[2], uu, tmp2);

  //Setting the optimizer
  size_t iter = 0;
  int status;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  gsl_multimin_function_fdf my_func;

  my_func.n = 3 * n;
  my_func.f = lj_f;
  my_func.df = lj_df;
  my_func.fdf = lj_fdf;
  my_func.params = &n;

  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc (T, 3 * n);
  gsl_multimin_fdfminimizer_set (s, &my_func, r, 0.01, 1e-4);

  //The actual optimization
  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);
      if (status)
	break;
      status = gsl_multimin_test_gradient (s->gradient, 1e-6);
    }
  while (status == GSL_CONTINUE && iter < 10000);


  est = gsl_multimin_fdfminimizer_minimum (s);
  gsl_vector_memcpy (r, gsl_multimin_fdfminimizer_x (s));

  //The values of the energies and forces after the optimization

  lj_fdf (r, &n, &uu, dr);
  tmp1 = tmp2 = 0.0;
  for (i = 0; i < 3 * n; i++)
    {
      tmp1 = gsl_vector_get (dr, i);
      tmp2 += tmp1 * tmp1;
    }
  tmp2 = sqrt (tmp2);

  fprintf (stdout,
	   "\n \nFinal Energy in file %s was E=%lf and the value of the norm of the gradient was |df|=%lf\n\n\n",
	   argv[3], uu, tmp2);
	   
  x = y = z = 0.0;

  for (i = 0; i < n; i++)
    {
      x+=gsl_vector_get (r, 3 * i);
      y+=gsl_vector_get (r, 3 * i + 1);
      y+=gsl_vector_get (r, 3 * i + 2);
    }

  for (i = 0; i < n; i++)
    {
      fprintf (out, "%.16e %.16e %.16e\n", gsl_vector_get (r, 3 * i)-x,
	       gsl_vector_get (r, 3 * i + 1)-y, gsl_vector_get (r, 3 * i + 2)-z);
    }
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (r);

  return 0;
}
