/*! \file lj_hess.h
   \brief Provides the headers and documentation for the function that 
 * computes values of the Hessian matrix of a Lennard-Jones potential 
 * energy surface for a given configuration using the structures 
 * required by the GSL Eigenvalue Routines 
 */

/*! \fn void lj_hessian(int n, double *x, double *y, double *z, 
 * 		double *hess);
 * \brief Calculates the Hessian matrix of the Lennard-Jones potential 
 * energy surface for the configuration in *x, *y, *z.
 * 
 * Input	:\n	int n, an integer that indicates the number of atoms. \n
 * 
 * 				double *x, *y, *z are arrays of doubles of dimension n 
 * 				that contain the x, y and z coordinates of the 
 * 				collection of n atoms that is being analized. \n
 * 
 * 				*hess is 9 n^2 double array that which the Hessian 
 * 				Hermitian matrix of the potential at the configuration 
 * 				stored in the arrays *x, *y, *z is returned by the 
 * 				function.
 * 				The columns (or quivalently the rows since the matrix is
 * 				hermitian) are stored one after the other in the array 
 * 				*hess.
*/

void lj_hessian(int n, double *x, double *y, double *z, double *hess);
