/*! \file energy.h
 * \brief Provides the headers and documentation for the functions that 
 * computes the energy surface that is sampled using Parallel
 * Tempering Monte Carlo in the program ptmc.c
 * 
 * It is assumed that the energy surface can be written as a sum of two
 * body terms: 
 * * \f$ V(R) = V(r_1, r_2,....,r_N) =
 * \sum_{i=1}^n \sum_{j < i}^n V_{i,j}(r_i,r_j) \f$
 * with n being the total number of atoms in the cluster.
 * 
 * The code in ptmc.c stores the positions of the n particles in 3 
 * arrays of doubles of dimensions n, the *x for the x components, *y 
 * for the y components and *z for the z components.
 */

/*! \fn double interaction_energy(double *x, double *y, double *z, int i
 * 			, int n);
 * \brief Returns the constribution to the energy due to the interaction
 * of the ith atom and the other n-1 atoms, i.e.
 * \f$ \sum_{j != i}^n V_{i,j}(\vec r_i,\vec r_j) \f$.
 * 
 * Input	: 	double *x, double *y, double *z the xyz coordinates of
 * 				the particles.
 * 
 * 				int i, the particle whose interaction with the other 
 * 				particles is to be calculated.
 * 
 * 				int n, the number of particles in the cluster.
 * 
 * Output	:	The energy contribution due to the interaction between
 * 				the ith particle and the other n-1 particles.
 * 
 * The reason why this function is explicitely required is because 
 * ptmc.c does single particle moves and thus for binary potentials the 
 * change of the energy will not depend on all the possible pairs of 
 * particles but just on the pairs that involve the particle that is 
 * being moved.
 */

/*! \fn double cluster_energy(double *x, double *y, double *z, int n);
 * \brief Returns the value of the potential energy surface of the n
 * atoms at the coordinates specified by *x, *y and *z, i.e.:
 * \f$ V(\vec R) = V(\vec r_1, \vec r_2, \ldots, \vec r_N)
 * =\sum_{i=1}^n \sum_{j < i}^n V_{i,j}(\vec r_i,\vec r_j) \f$
 * 
 * Input	:	double *x, double *y, double *z the xyz coordinates of
 * 				the particles.
 * 
 * 				int n, the number of particles in the cluster.
 * 
 * Output	:	The value of the potential energy surface at the 
 * 				configuration given by *x, *y and *z.
 */

double interaction_energy (double *x, double *y, double *z, int i, int n);
double cluster_energy (double *x, double *y, double *z, int n);
