/*! \file lj_grad.h
   \brief Provides the headers and documentation for the functions that 
 * compute values of the Lennard-Jones potential energy surface and its
 * gradient for a given configuration using the structures required by 
 * the GSL Minimization Routines.
 * 
 * In lj_grad.c the Lennard-Jones potential and its gradient are 
 * implemented using the datatypes required by the GNU Scientific 
 * Library to perform numerical local minimization. \n
 * The Lennard-Jones potential for n atoms of different species is given
 *  by:\n \f$ V( R) = V(r_1, r_2,..., r_N)=\sum_{i=1}^n
 * \sum_{j>i}^n \epsilon_{i,j} V_{LJ}(|r_i-r_j|/\sigma_{i,j}) \f$\n
 * where \f$ V_{LJ}(x)=4 (1/r^{12}-1/r^6)\f$, \f$ \sigma_{i,j}\f$ 
 * is the distance for which the potential energy between particles 
 * i and j is zero and \f$ \epsilon_{i,j} \f$ is the depth of the 
 * potential well between particles i and j.
 * As it was mentioned before, the indices i and j are labels for atomic
 * species. The positions of the atoms are assumed to be in a 
 * gsl_vector *r of 3n components ordered as follows 
 * r={x_1, y_1, z_1, x_2, y_2, z_2, ..., x_n, y_n, z_n}. Notice that in 
 * C the first element of an array is the zeroth element, so x_1 is the 
 * 0 element ot the array whereas z_n is the 3n-1 element of such array. 
 * To perform the calculation the user needs to provide the functions
 * double lj_sigma2(int i, int j) and double lj_epsilon(int i, int j)
 * in the file lj_params.c.
*/

/*! \fn double lj_f(const gsl_vector *r, void *params);
 * \brief Returns the value of the Lennard-Jones potential for 
 * configuration r.
 *  
 * Input	:\n	gsl_vector * r configuration of the atoms. \n 
 * 
 *				void * params an integer should be passed by reference 
 * 				indicating the number of atoms, n. \n \n
 * Output	:\n	Returns the value of the potential for n atoms in 
 * 				the configuration specified by *r. \n
 */ 

/*! \fn void lj_df(const gsl_vector * r, void * params, gsl_vector * g);
 * \brief Calculates the value of the gradient of the Lennard-Jones
 * potential energy surface.
 * 
 * Input	:\n	gsl_vector * r configuration of the atoms. \n 
 * 
 *				void * params an integer should be passed by reference 
 * 				indicating the number of atoms, n. \n 
 * 
 * 				gsl_vector * g the value of the gradient will be 
 * 				returned in this vector of 3 n componets.  The ordering 
 * 				corresponds to the same ordering in the gsl_vector r, 
 * 				i.e. , g={dV/dx_1, dV/dy_1, dV/dz_1, dV/dx_2, dV/dy_2
 * 				, dV/dz_2, ..., dV/dx_n, dV/dy_n, dV/dz_n}.\n \n
*/

/*! \fn void lj_fdf(const gsl_vector *r, void *params, 
 * 		double *f, gsl_vector *g);
 * \brief Calculates the value of the Lennard-Jones potential energy 
 * surface and its gradient.
 * 
 * Input	:\n	gsl_vector * r configuration of the atoms. \n 
 * 
 *				void * params an integer should be passed by reference 
 * 				indicating the number of atoms, n. \n 
 * 
 * 				double *f the value of the potencial energy surface is 
 * 				returned in this variable.\n
 * 
 * 				gsl_vector * g the value of the gradient will be 
 * 				returned in this vector of 3 n components.  The ordering 
 * 				corresponds to the same ordering in the gsl_vector r, 
 * 				i.e. , g={dV/dx_1, dV/dy_1, dV/dz_1, dV/dx_2, dV/dy_2
 * 				, dV/dz_2, ..., dV/dx_n, dV/dy_n, dV/dz_n}.\n 
 */

double lj_f(const gsl_vector *r, void *params);
void lj_df(const gsl_vector * r, void * params, gsl_vector * g);
void lj_fdf(const gsl_vector *r, void *params, double *f, gsl_vector *g);
