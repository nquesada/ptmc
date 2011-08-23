/*! \file lj_params.h
   \brief Provides the headers and documentation for the function that 
 * returns the Lennard-Jones parametes.
 */

/*! \fn double lj_sigma2(int i,int j);
 * \brief Returns the value of the squared of the Lennard-Jones
 * parameter sigma between atoms i and j. See the documentation of 
 * lj_grad.h
 * 
 * The function returns the square of sigma and not sigma itself simply
 * because the Lennard-Jones potential energy function depends only on
 * even powers of the separation between the atoms. Finally, notice that 
 * this funcion is defined as an inline function to save time, specially
 * in the case of homogeneous clusters in which it adopts the constant 
 * value sigma2(i,j)=1.
 */

/*! \fn double lj_epsilon(int i,int j);
 * \brief Returns the value of the Lennard-Jones
 * parameter epsilon between atoms i and j. See the documentation of 
 * lj_grad.h
 * 
 * Notice that this funcion is defined as an inline function to save 
 * time, specially in the case of a homogeneous cluster in which it 
 * adopts the constant value epsilon(i,j)=1.
 */

double lj_sigma2(int i,int j);
double lj_epsilon(int i,int j);
