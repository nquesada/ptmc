/*! \file ptmc.h
 * \brief Defines macros necessary for the compilation of ptmc.c 
 * the main routine to do Parallel Tempering Monte Carlo.
 * To process the data generated by this program see ptmc_data.c
 * 
 * Input	:	The program requieres 2 or 3 files depending on wether 
 * 				the radial distribution function will be calculated or
 * 				not. These files should be passes at the time of 
 * 				execution.\n\n\n
 * 
 * 				The first file (say, parameters.in) 
 * 				that should be passed to the program
 * 				should contain the following values in its firt line in
 * 				the following order (separated by spaces):\n
 * 
 * 				1. Number of atoms, this will be stored in int n.\n
 * 				2. Frequency (in monte carlo steps) at which swaps 
 * 				between replicas will be 
 * 				attempted, this will be stored in int swap_freq.\n
 * 				3. Number of Monte Carlo Steps, this will be stored 
 * 				in int mc_steps. [Added4Aug11] In monte carlo step
 * 				depending of the value of the macro NUMBER_OF_DOPANTS
 * 				two types of movements will be attempted. If 
 * 				NUMBER_OF_DOPANTS is zero single particle movements will
 * 				be attempted. In this case each of atom of the cluster
 * 				will be randomly translated. If NUMBER_OF_DOPANTS is not
 * 				zero then also swaps between atoms of different species
 * 				will be attempted. To satisfy detailed balance the swap
 * 				is done as follows: Two atoms of the cluster are 
 * 				selected at random. If the two atoms are of the same
 * 				type the swap will always be accepted if they are not
 * 				of the same type the swap will be attempted and the 
 * 				Metropolis comparison will be made.
 * 				\n
 * 				4. Frequency at which partial averages and cumulative 
 * 				sums will be written to disk, this will be stored in 
 * 				int save_freq.. \n
 * 				5. Number of Equilibration steps, this will be stored 
 * 				in int eq_steps. [Added4Aug11] 
 * 				Notice that this thermalization will 
 * 				not be done in parallel and the only type of Monte 
 * 				Carlo steps that will be done are single particle
 * 				movements.\n
 * 				6. Minimum temperature of the simulation, this will be 
 * 				stored in double t0.\n
 * 				7. Maximum temperature of the simulation, this will be 
 * 				stored in double tf. \n
 * 				8. Minimum value to be included in the energy 
 * 				histograms to be saved, this will be 
 * 				stored in double U0.\n
 * 				9. Maximum value to be included in the energy 
 * 				histograms to be saved, this will be 
 * 				stored in double Uf.\n
 * 				10. Size of the bins in the energy histograms to be 
 * 				saved, this will be stored in double dU.\n
 * 				11. Radius of the sphere in which the cluster will be 
 * 				constrained to move, this will be instored in 
 * 				double box_radius. Notice that the sphere in which the 
 * 				cluster will be confined is assumed to be centered at
 * 				the origin. [Added4Aug11] Also notice that the program 
 * 				will check that the initial configurations provided in
 * 				the initial configurations file (see next paragraph) fit
 * 				inside the simulation box it they do not the simulation
 * 				will be aborted. \n \n \n
 * 				
 * 				The second file (say, names.in)
 * 				should contain the names of the files
 * 				containing the initial coordinates (initial 
 * 				configurations) for each replica. Each line of the
 * 				file should contain the name of file at the name of 
 * 				such file should not contain more characters than
 * 				MAX_LENGTH_CHAR.
 * 				There should be at least as many names as replicas are
 * 				minus one, this is because the 0th node does not do any
 * 				monte carlo calculation only directs the other nodes.
 * 				The node 0 will send to ther other np-1 nodes the names
 * 				contained in names.in, and each of them will attempt
 * 				to read the first n lines (n being the number of atoms)
 * 				of such a file whih is assumed to have the coordinates 
 * 				in the following format:\n \n
 * 				x_1 y_1 z_1 \n
 * 				x_2 y_2 z_2 \n
 * 				....        \n
 * 				x_n y_n z_n \n
 * 				\n
 * 				i.e, each line containd the x, y and z components of the
 * 				coordinates of the atoms.
 * 				Once it read the initial configuration it will check
 * 				that the given coordinates are inside the sphere of
 * 				radius box_radius. \n \n \n
 * 
 * 
 * 				The third file (say rdfparameters.in)
 * 				should only be passed if the macro 
 * 				CALCULATE_RDF is set to a value different than zero 
 * 				which implies that the radial distribution functions
 * 				will be calculated. If this is the case the the file 
 * 				should contain in its first line the following 
 * 				parameters separated by spaces:
 * 
 * 				1. The frequency at which the RDF should be calculated.
 * 				this will be instored in int rdf_freq.\n
 * 				2. Minimum value to be included in the RDFs
 * 				histograms to be saved, this will be 
 * 				stored in double radialmin.\n
 * 				3. Maximum value to be included in the RDFs 
 * 				histograms to be saved, this will be 
 * 				stored in double radialmax.\n
 * 				4. Size of the bins in the RDF histograms to be 
 * 				saved, this will be stored in double dradial.\n 
 * 				
 * 				Notice that the RDFs are calculated with the origin of 
 * 				coordinates in the geometric center of the cluster.
 * 
 * \n\n\n
 * 
 * Output:		The program ptmc.c will produce a number of output files
 * 				\n
 * 				Node 0 will produce only two file.
 * 				The first file, exchange.c,  will contain
 * 				only two colums the first column will give the number
 * 				of times a swap between replicas i and i+1 was atempted.
 * 				The second column will contain the number of times
 * 				the exchange between replicas i and i+1 was succesfull.
 * 				The second file, nprocs.dat, will contain in the first
 * 				line the number processes that where used and then
 * 				each line will contain in increasing order the 
 * 				temperatures of each replica
 * 
 * 				Node i (i>0) will produce the following files which are 
 * 				added to the file every save_freq monte carlo steps 
 * 				(i.e. each [Added4Aug11] line of the files named below 
 * 				will contain an
 * 				 snapshop of the accumulated values):\n\n
 * 				- avs[i].dat contains the accumulated values of the 	
 * 				following quantities: \n
 * 				1. Number of monte carlo steps. \n
 * 				2. Number of monte carlo moves rejected.\n
 * 				3. Number of times the an atom tried to move outside the
 * 				 box. \n
 * 				4. Number of times the energy has been calculated normU. \n
 * 				5. Acumulated value of the energy, U. \n
 * 				6. Acumulated value of the energy squared U2 . \n
 * 				7. Number of exchanges attempted between two different
 * 				atoms.
 * 				8. Number of succesful exchanges between two different
 * 				atoms.\n\n
 * 
 * 				[Added4Aug11] The reason not only the last snapshot of 
 * 				the accumulated 
 * 				values of the observables is saved but rather the whole
 * 				dynamical evolution is because these snapshots can be used
 * 				to calculate partial averages that later can be used 
 * 				to calculate errors in the expected values of the 
 * 				observables.
 * 
 * 				Notice that expected values of the energy and it square
 * 				are given by <U>=U/normU and <U2>=U2/normU. \n
 * 				
 * 				The last two quantities in the list will only appear if
 * 				NUMBER_OF_DOPANTS is different from zero then:\n\n
 * 
 * 				- cos[i].dat in this file the coordinates that replica i
 * 				has will be saved every save_freq steps. \n\n
 * 
 * 				[Added4Aug11] Notice that a number of configurations 
 * 				will be saved not only the latest one which will be last
 * 				n lines of the file cos[i].dat (n being the number of atoms)
 * 				The purpose of saving configurations through the evolutions
 * 				is that later these can be used to do basin analysis,
 * 				for instance they can be locally optimized to know
 * 				the basins that the system is visiting at a given 
 * 				temperature.
 * 
 * 				- e[i].dat The statistical frequencies of the energy 
 * 				histogram will be saved in this file. The energies to
 * 				which these frequencies corresponds can be infered from
 * 				the input file.\n\n
 * 
 * 				Finally, if the radial distribution functions are being
 * 				calculated (i.e if CALCULATE_RDF is not set to zero) 
 * 				then the replica i will produce the following files:\n\n
 * 				
 * 				- ord[i].dat: contains the accumulated values of the 	
 * 				following quantities: \n
 * 				1. Number of times the distance of the atoms to the 
 * 				geometric center has been calculated normr
 * 				2. Acumulated value of the distance of the 
 * 				atoms to the geometric center, rA.\n
 * 				3. Acumulated value of the distance of the 
 * 				atoms to the geometric center squared rA2.\n
 * 
 * 				The mean distance of the atoms is given by 
 * 				<r>=rA/normr and the mean value of the square by:
 * 				<r^2>=rA2/normr.
 * 
 * 				- If also there are dopants then the first three items 
 * 				will correspond to the acumuted values for the atoms
 * 				of the matrix and 3 additional colums will appear giving
 * 				4. Number of times the distance of the dopant atoms to the 
 * 				geometric center has been calculated, normr \n
 * 				5. Acumulated value of the distance of the 
 * 				dopant atoms to the geometric center, rA.\n
 * 				6. Acumulated value of the distance of the 
 * 				dopant atoms to the geometric center squared rA2.\n
 * 
 * 				- rpfA[i].dat: will contain the statistical frequencies 
 * 				of the radial distribution function.
 * 
 * 				If NUMBER_OF_DOPANTS is different from zero then 
 * 				rpfA[i].dat will contain the frequency hisotgram of the
 * 				matrix atoms and:\n\n
 * 
 * 				-rpfB[i].dat: will contain the statistical frequencies 
 * 				of the radial distribution function of the dopant atom. 
 * \n\n\n
 * 
 * 
 * 
 * To execute the program type something like:
 * \n
 * mpirun -np nprocs -machinefile machines ptmc.out input.in names.in 
 * \n
 * or
 * \n
 * mpirun -np nprocs -machinefile machines ptmc.out input.in names.in order.in
 * \n
 * if the RDF will be calculated.



 */
 
/*! \def CALCULATE_RDF 
 * \brief Specifies wether or not the Radial Distribution Function (RDF) 
 * will be calculated, if it set to zero the RDFs will NOT be calculated
 * any positive value will instrunct ptmc.c to calculate them.
 */ 


/*! \def NUMBER_OF_DOPANTS
 * \brief Specifies the number of dopant atoms in the cluster. 
 * This implies that for the first M atoms the radial distribution 
 * functions will be calculated separately from the rest of the atoms 
 * and also that one type of monte carlo step that will be attempted is 
 * the exchange between the first M atoms and the other n-M atoms of the
 * cluster. Notice that the specification of the number of dopant atoms 
 * should be consistent with whatever energy function is specified in 
 * energy.c.If there are no dopant atoms nor it is not necessary to study
 * their RDF separately nor the user does not want to include 
 * Monte Carlo steps that include swapping particles of different kinds 
 * set M 0
*/

/*! \def MAX_LENGTH_CHAR
 * \brief Define the maximum number of characteres that a the file
 * names that the program ptmc.c will read and write.
 */ 
#define CALCULATE_RDF 1
#define NUMBER_OF_DOPANTS 0
#define MAX_LENGTH_CHAR 50

