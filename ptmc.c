/***********************************************************************
 *  Parallel Tempering Monte Carlo.
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

/*! \file ptmc.c
 * \brief Provides the main routine to do Parallel Tempering Monte Carlo
 * See documentation in file ptmc.h
 */

//C Headers
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<malloc.h>

//GSL Headers
#include<gsl/gsl_rng.h>

//MPI Headers
#include<mpi.h>

//Local Headers
#include"energy.h"
#include"lj_params.h"
#include"ptmc.h"

int
main (int argc, char *argv[])
{
#if CALCULATE_RDF
  if (argc != 4)
    {
      fprintf (stdout, "Not enough parameters at input \n");
      return 1;
    }
#else
  if (argc != 3)
    {
      fprintf (stdout, "Not enough parameters at input \n");
      return 1;
    }
#endif
  // Variable Declaration
  char init_conf_name[MAX_LENGTH_CHAR];	//String with the name of the initial conf file

  int i, j, k, l;		//Iterators
  int n;			//Number of particles
  int swap_freq;		//Frequency of swap attemps between replicas
  int mc_steps;			//Number montecarlo steps
  int swi;			//Switch to determine which processes might swap  
  int swip1;			//Value of swi plus one
  int eq_steps;			//Number of equilibration steps
  int rat, ratf;		//Varib
  int gg;			//Variable to receive the output of fscanf
  int zero = 0;

  int *tries;
  int *accepted;

  long random_seed = time (NULL);	//Time when the program starts and random seed.

  double swap_prob;		//Probability of swapping two adjacent replicas
  double step_size;		//Size of the step for MC moves
  double u, us, du, uu;		//Energies
  double xn, yn, zn;		//Temporal vars
  double fb;			//Boltzmann factor
  double ras;			//Ratio of accepted single particle moves
  double box_radius;		//Size of the box;
  double outmoves = 0;		//Number of Moves out of the box
  double tmp;			//Temporal variable
  double t0, tf, t;		//initial, final and temperature of the ith node
  double t_ratio;		//Ratio of the temperatures of the nodes, they are in geometric progression
  double ti;			//Temperature of Each process
  double betai;			//Beta for each process
  double counter = 0;		//Number of MC steps
  double countrej = 0;		//Number of rejected steps
  double myinfo[2];		//temperature and energy at the time of swapping
  double inf[2];		//Array used to communicate between processes
  double *betas;		//Array of inverse temperatures of the nodes kept by master
  double *energies;		//Array of energies of the nodes kept by master
  double *x, *y, *z;		//Coordinates of the particle
  double *xx, *yy, *zz;		//Temporal vectors

  FILE *in;
  FILE *swap_stats_file;
  FILE *parameters;
  FILE *names;
  FILE *nprocs;

  // MPI Variables 
  int my_rank;			// Rank of processes
  int p;			//Number of Processes
  int source;			// Rank of sender
  int dest;			//Rank of receiver
  int tag;			//Tag for messages
  MPI_Status status;		//Return status for receive
  int pm2;			//pm2=p-2

  // GSL Variables for Random Number Generation
  gsl_rng *generator;
  const gsl_rng_type *type;

  // Variables for energy statistics
  char charU[MAX_LENGTH_CHAR];
  int nU;
  int *histU;
  double U = 0, U2 = 0, Uf, U0, dU, normU = 0;
  FILE *fileU;

  // Variables for saving partial configurations and results
  char average_save[MAX_LENGTH_CHAR];
  char conf_save[MAX_LENGTH_CHAR];
  FILE *av_save;
  FILE *cf_save;
  int save_freq;

  // Variables to select atoms to be swapped
#if NUMBER_OF_DOPANTS
  int pe1, pe2;
  double exc_tried = 0, exc_rejected = 0;
#endif

  // Variables for the calculation of radial distribution functions 
#if CALCULATE_RDF
  FILE *rpf1p;
  char rpf1char[MAX_LENGTH_CHAR], ordchar[MAX_LENGTH_CHAR];
  double rpf1a = 0, rpf1a2 = 0;;
  double xxm, yym, zzm, ttmp, xtmp, ytmp, ztmp;
  int inttmp;
  double radialmin, radialmax, dradial;
  int nradial;
  int *rpf1h;
  int normrpf1 = 0;
  int rdf_freq;
  FILE *ord_save;
  FILE *inputs;
#if NUMBER_OF_DOPANTS
  FILE *rpf2p;
  char rpf2char[MAX_LENGTH_CHAR];
  double rpf2a = 0, rpf2a2 = 0;
  int *rpf2h;
  int normrpf2 = 0;
#endif
#endif

  //GSL Vars initialization
  type = gsl_rng_mt19937;
  generator = gsl_rng_alloc (type);

  //The Parallel Program Starts Here:
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &p);
  fprintf (stdout, "%d tt=%d\n", my_rank, (int) random_seed);


  //Reading the input Parameters  
  parameters = fopen (argv[1], "r");
  if (parameters == NULL)
    {
      fprintf (stdout, "The input file %s was not found\n", argv[1]);
      MPI_Abort (MPI_COMM_WORLD, 1);
    }
  else
    {
      gg = fscanf (parameters, "%d %d %d %d %d %lf %lf %lf %lf %lf %lf",
		   &n, &swap_freq, &mc_steps, &save_freq, &eq_steps, &t0, &tf,
		   &U0, &Uf, &dU, &box_radius);
      fclose (parameters);
    }
  nU = (int) ((Uf - U0) / dU);	//Number of Bins in the Energy Histogram

#if CALCULATE_RDF
  inputs = fopen (argv[3], "r");
  if (inputs == NULL)
    {
      fprintf (stdout, "The input file %s was not found\n", argv[3]);
      MPI_Abort (MPI_COMM_WORLD, 1);
    }
  else
    {
      gg =
	fscanf (inputs, "%d %lf %lf %lf", &rdf_freq, &radialmin, &radialmax,
		&dradial);
      fclose (inputs);
    }
  nradial = (int) ((radialmax - radialmin) / dradial);
#endif

  if (my_rank == 0)
    {
      //Random Number Generator Initialization
      //For process zero is time since epoch plus its rank, i.e, zero
      gsl_rng_set (generator, random_seed);

      betas = (double *) malloc (p * sizeof (double));
      energies = (double *) malloc (p * sizeof (double));
      tries = (int *) malloc (p * sizeof (int));
      accepted = (int *) malloc (p * sizeof (int));

      for (i = 0; i < p; i++)
	{
	  tries[i] = 0;
	  accepted[i] = 0;
	}
      tag = 1;
      pm2 = p - 2;
      nprocs = fopen ("nprocs.dat", "w");
      fprintf (nprocs, "%d\n", p);

      //The temperatures will be generated in a GEOMETRIC progression:
      t = t0;
      t_ratio = exp (log (tf / t0) / (p - 2));
      for (i = 1; i < p; i++)
	{
	  dest = i;
	  fprintf (nprocs, "%lf\n", t);
	  MPI_Send (&t, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
	  t = t * t_ratio;
	}
      fclose (nprocs);
      names = fopen (argv[2], "r");
      tag++;
      //Sending to each slave the name of its initial configuration
      for (i = 1; i < p; i++)
	{
	  dest = i;
	  gg = fscanf (names, "%s", init_conf_name);
	  MPI_Send (&init_conf_name, MAX_LENGTH_CHAR, MPI_CHAR, dest, tag,
		    MPI_COMM_WORLD);
	}
      fclose (names);

      for (l = 1; l < mc_steps + 1; l++)
	{
	  //Each swap_freq steps one swap between replicas is tried
	  if (l % swap_freq == 0)
	    {
	      tag = l / swap_freq;
	      //For each cycle the tags are synchronized
	      for (i = 1; i < p; i++)
		{
		  dest = i;
		  MPI_Recv (&inf, 2, MPI_DOUBLE, i, tag, MPI_COMM_WORLD,
			    &status);
		  betas[i] = inf[0];
		  energies[i] = inf[1];
		}
	      //Here node zero selects the nodes that will try to swap
	      swi = ((int) (pm2 * gsl_rng_uniform (generator))) + 1;
	      //swi contains a number 1<=swi<=Number of processes-1
	      //swi and swi+1 will be the nodes that might change configuration
	      tries[swi] = tries[swi] + 1;
	      swip1 = swi + 1;
	      tag++;
	      //In this steps node zero sends zero to all the processes that do not exchange configuration
	      //And does the Metropolis comparison for the 2 subsystems, if possitive sends to each
	      //the rank of the other process
	      for (i = 1; i < swi; i++)
		{
		  MPI_Send (&zero, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
		}
	      swap_prob =
		exp ((betas[swip1] - betas[swi]) * (energies[swip1] -
						    energies[swi]));
	      if (swap_prob > gsl_rng_uniform (generator))
		{
		  MPI_Send (&swi, 1, MPI_INT, swip1, tag, MPI_COMM_WORLD);
		  MPI_Send (&swip1, 1, MPI_INT, swi, tag, MPI_COMM_WORLD);
		  accepted[swi] = accepted[swi] + 1;
		}
	      else
		{
		  MPI_Send (&zero, 1, MPI_INT, swip1, tag, MPI_COMM_WORLD);
		  MPI_Send (&zero, 1, MPI_INT, swi, tag, MPI_COMM_WORLD);
		}
	      for (i = swi + 1; i < p; i++)
		{
		  MPI_Send (&zero, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
		}
	    }
	  //Each save_freq steps partial statistics are saved

	  if ((l % save_freq) == 0)
	    {
	      tag = 0;
	      swap_stats_file = fopen ("exchange.dat", "w");
	      for (i = 1; i < p; i++)
		{
		  fprintf (swap_stats_file, "%d %d\n", tries[i], accepted[i]);
		}
	      fclose (swap_stats_file);
	    }
	}
      free (energies);
      free (betas);
      free (tries);
      free (accepted);
    }


  if (my_rank != 0)
    {
      counter = 0;
      //Variable Allocation
      x = (double *) malloc (n * sizeof (double));
      y = (double *) malloc (n * sizeof (double));
      z = (double *) malloc (n * sizeof (double));
      xx = (double *) malloc (n * sizeof (double));
      yy = (double *) malloc (n * sizeof (double));
      zz = (double *) malloc (n * sizeof (double));

#if CALCULATE_RDF
      rpf1h = (int *) malloc (nradial * sizeof (int));
#if NUMBER_OF_DOPANTS
      rpf2h = (int *) malloc (nradial * sizeof (int));
#endif
      for (i = 0; i < nradial; i++)
	{
	  rpf1h[i] = 0;
#if NUMBER_OF_DOPANTS
	  rpf2h[i] = 0;
#endif
	}
#endif
      histU = (int *) malloc (nU * sizeof (int));
      //Histogram Initialization
      for (i = 0; i < nU; i++)
	{
	  histU[i] = 0;
	}
      //Receiving temperatures from master
      source = 0;
      tag = 1;
      MPI_Recv (&ti, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
      betai = 1 / ti;
      tag++;
      //Receiving name of the initial configuration file from master
      MPI_Recv (&init_conf_name, MAX_LENGTH_CHAR, MPI_CHAR, source, tag,
		MPI_COMM_WORLD, &status);

      //Random number initialization for each node
      gsl_rng_set (generator, random_seed + my_rank);

      //Writing the strings with the names of the output files
      sprintf (charU, "e%d.dat", my_rank);
      sprintf (average_save, "avs%d.dat", my_rank);
      sprintf (conf_save, "cos%d.dat", my_rank);
#if CALCULATE_RDF
      sprintf (rpf1char, "rpfA%d.dat", my_rank);
#if NUMBER_OF_DOPANTS
      sprintf (rpf2char, "rpfB%d.dat", my_rank);
#endif
      sprintf (ordchar, "ord%d.dat", my_rank);
#endif
      step_size = 1;		//The initial step of the random walk is set to one
      box_radius = box_radius * box_radius;	//Takes the square of box_radius to avoid calculating square roots in comparisons.

      //Reading initial configuration
      in = fopen (init_conf_name, "r");
      if (in == NULL)
	{
	  fprintf (stdout, "The input file for process %d was not found \n",
		   my_rank);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
      else
	{
	  for (i = 0; i < n; i++)
	    {
	      gg = fscanf (in, "%lf %lf %lf", &x[i], &y[i], &z[i]);
	      tmp = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
	      if (tmp > box_radius)
		{
		  fprintf (stdout,
			   "The initial configuration for replica %d is out of the box\n",
			   my_rank);
		  return 1;
		}
	    }
	}
      fclose (in);


      //Equilibration steps
      u = cluster_energy (x, y, z, n);
      j = 0;
      rat = 0;
      ratf = 0;
      for (k = 0; k < eq_steps; k++)
	{
	  j++;
	  j = j % n;
	  if (j == 0)
	    {
	      rat = 0;
	      ratf = 0;
	    }
	  //Energy of the cluster before the trial move
	  rat++;
	  u = interaction_energy (x, y, z, j, n);
	  xn = x[j];
	  yn = y[j];
	  zn = z[j];
	  x[j] = xn + step_size * (gsl_rng_uniform (generator) - 0.5);
	  y[j] = yn + step_size * (gsl_rng_uniform (generator) - 0.5);
	  z[j] = zn + step_size * (gsl_rng_uniform (generator) - 0.5);
	  while ((x[j] * x[j] + y[j] * y[j] + z[j] * z[j]) > box_radius)
	    {
	      x[j] = xn + step_size * (gsl_rng_uniform (generator) - 0.5);
	      y[j] = yn + step_size * (gsl_rng_uniform (generator) - 0.5);
	      z[j] = zn + step_size * (gsl_rng_uniform (generator) - 0.5);
	    }

	  //Energy of the cluster after the trial move
	  us = interaction_energy (x, y, z, j, n);
	  du = us - u;
	  fb = exp (-du * betai);
	  if (gsl_rng_uniform (generator) > fb)
	    {
	      ratf++;
	      x[j] = xn;
	      y[j] = yn;
	      z[j] = zn;
	    }
	  //Adapting the steps such that rejection rate is between 40%-60%
	  if (j == 0 && k > 0)
	    {
	      ras = ((double) ratf) / ((double) rat);
	      if (ras < 0.4)
		{
		  step_size = 1.1 * step_size;
		}
	      if (ras > 0.6)
		{
		  step_size = 0.9 * step_size;
		}
	    }
	}
      //MonteCarlo  
      rat = 0;
      ratf = 0;
      uu = cluster_energy (x, y, z, n);

      for (k = 1; k < mc_steps + 1; k++)
	{
	  j = k % n;
	  //Energy of the Cluster before the trial move
	  rat++;
	  u = interaction_energy (x, y, z, j, n);
	  xn = x[j];
	  yn = y[j];
	  zn = z[j];
	  x[j] = xn + step_size * (gsl_rng_uniform (generator) - 0.5);
	  y[j] = yn + step_size * (gsl_rng_uniform (generator) - 0.5);
	  z[j] = zn + step_size * (gsl_rng_uniform (generator) - 0.5);
	  while ((x[j] * x[j] + y[j] * y[j] + z[j] * z[j]) > box_radius)
	    {
	      x[j] = xn + step_size * (gsl_rng_uniform (generator) - 0.5);
	      y[j] = yn + step_size * (gsl_rng_uniform (generator) - 0.5);
	      z[j] = zn + step_size * (gsl_rng_uniform (generator) - 0.5);
	      outmoves++;
	    }

	  //Energy of the Cluster after the trial move
	  us = interaction_energy (x, y, z, j, n);
	  du = us - u;
	  //Temperatures are in energy units
	  fb = exp (-du * betai);
	  counter++;
	  uu = uu + du;

	  //Metropolis comparison
	  if (gsl_rng_uniform (generator) > fb)
	    {
	      x[j] = xn;
	      y[j] = yn;
	      z[j] = zn;
	      uu = uu - du;
	      countrej++;
	      ratf++;
	    }

	  //Adjusting the step of the random walk every n steps (everytime (i%n==0))
	  if (j == 0 && k > 0)
	    {
	      ras = ((double) ratf) / ((double) rat);
	      if (ras < 0.4)
		{
		  step_size = 1.1 * step_size;
		}
	      if (ras > 0.6)
		{
		  step_size = 0.9 * step_size;
		}
	      rat = 0;
	      ratf = 0;
	    }
	  //Energy Histogram and expected values updated
	  l = (int) ((uu - U0) / dU);
	  if (l >= 0 && l < nU)
	    histU[l] = histU[l] + 1;
	  U += uu;
	  U2 += (uu * uu);
	  normU++;

	  //Exchange between different particles
#if NUMBER_OF_DOPANTS
	  pe1 = (int) (n * gsl_rng_uniform (generator));
	  pe2 = (int) (n * gsl_rng_uniform (generator));
	  if (pe1 >= n || pe2 >= n)
	    {
	      fprintf (stdout,
		       "Problems swapping particles tried to get particle %d or %d that doesnt exists\n",
		       pe1, pe2);
	    }
	  if (((pe1 < NUMBER_OF_DOPANTS) && (pe2 >= NUMBER_OF_DOPANTS))
	      || ((pe2 < NUMBER_OF_DOPANTS) && (pe1 >= NUMBER_OF_DOPANTS)))
	    {
	      exc_tried++;
	      xn = x[pe1];
	      yn = y[pe1];
	      zn = z[pe1];
	      x[pe1] = x[pe2];
	      y[pe1] = y[pe2];
	      z[pe1] = z[pe2];
	      x[pe2] = xn;
	      y[pe2] = yn;
	      z[pe2] = zn;
	      du = cluster_energy (x, y, z, n) - uu;
	      uu = uu + du;
	      fb = exp (-du * betai);
	      if (gsl_rng_uniform (generator) > fb)
		{
		  x[pe2] = x[pe1];
		  y[pe2] = y[pe1];
		  z[pe2] = z[pe1];
		  x[pe1] = xn;
		  y[pe1] = yn;
		  z[pe1] = zn;
		  uu = uu - du;
		  exc_rejected++;
		}
	      l = (int) ((uu - U0) / dU);
	      if (l >= 0 && l < nU)
		histU[l] = histU[l] + 1;
	      U += uu;
	      U2 += (uu * uu);
	      normU++;
	    }
	  else
	    {
	      l = (int) ((uu - U0) / dU);
	      if (l >= 0 && l < nU)
		histU[l] = histU[l] + 1;
	      U += uu;
	      U2 += (uu * uu);

	      normU++;
	    }
#endif
	  //Calculating radial distribution functions and mean values of positions
#if CALCULATE_RDF
	  if ((k % rdf_freq) == 0)
	    {
	      //Calculating Geometric Center
	      xxm = yym = zzm = 0.0;
	      for (i = 0; i < n; i++)
		{
		  xxm += x[i];
		  yym += y[i];
		  zzm += z[i];
		}
	      xxm = xxm / n;
	      yym = yym / n;
	      zzm = zzm / n;
	      for (i = NUMBER_OF_DOPANTS; i < n; i++)
		{
		  xtmp = xxm - x[i];
		  ytmp = yym - y[i];
		  ztmp = zzm - z[i];
		  ttmp = sqrt (xtmp * xtmp + ytmp * ytmp + ztmp * ztmp);
		  rpf1a += ttmp;
		  rpf1a2 += ttmp * ttmp;
		  inttmp = (int) ((ttmp - radialmin) / dradial);
		  if (inttmp < nradial)
		    {
		      rpf1h[inttmp]++;
		      normrpf1++;
		    }
		}
	      // For the dopant atoms if any
#if NUMBER_OF_DOPANTS
	      for (i = 0; i < NUMBER_OF_DOPANTS; i++)
		{
		  xtmp = xxm - x[i];
		  ytmp = yym - y[i];
		  ztmp = zzm - z[i];
		  ttmp = sqrt (xtmp * xtmp + ytmp * ytmp + ztmp * ztmp);
		  rpf2a += ttmp;
		  rpf2a2 += ttmp * ttmp;
		  inttmp = (int) ((ttmp - radialmin) / dradial);
		  if (inttmp < nradial)
		    {
		      rpf2h[inttmp]++;
		      normrpf2++;
		    }
		}
#endif
	    }
	  //Saving RDF histogram and mean values
	  if ((k % save_freq) == 0)
	    {
	      ord_save = fopen (ordchar, "a");
	      rpf1p = fopen (rpf1char, "w");
#if NUMBER_OF_DOPANTS
	      rpf2p = fopen (rpf2char, "w");
#endif
	      for (i = 0; i < nradial; i++)
		{
		  fprintf (rpf1p, "%d\n", rpf1h[i]);
#if NUMBER_OF_DOPANTS
		  fprintf (rpf2p, "%d\n", rpf2h[i]);
#endif
		}
	      fprintf (ord_save, "%d %.16e %.16e ", normrpf1, rpf1a, rpf1a2);
#if NUMBER_OF_DOPANTS
	      fprintf (ord_save, "%d %.16e %.16e", normrpf2, rpf2a, rpf2a2);
#endif
	      fprintf (ord_save, "\n");
	      fclose (rpf1p);
#if NUMBER_OF_DOPANTS
	      fclose (rpf2p);
#endif
	      fclose (ord_save);
	    }
#endif
	  //Every swap_freq steps a swap is attempted
	  if ((k % swap_freq) == 0)
	    {
	      tag = k / swap_freq;
	      myinfo[0] = betai;
	      myinfo[1] = uu;
	      dest = 0;
	      MPI_Send (&myinfo, 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	      tag++;
	      MPI_Recv (&swi, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);

	      if (swi != 0)
		{
		  if (swi > my_rank)
		    {
		      tag++;
		      MPI_Send (&(x[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD);
		      tag++;
		      MPI_Send (&(y[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD);
		      tag++;
		      MPI_Send (&(z[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD);
		      tag++;
		      MPI_Recv (&(x[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD, &status);
		      tag++;
		      MPI_Recv (&(y[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD, &status);
		      tag++;
		      MPI_Recv (&(z[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD, &status);
		      uu = cluster_energy (x, y, z, n);
		    }
		  else
		    {
		      for (i = 0; i < n; i++)
			{
			  xx[i] = x[i];
			  yy[i] = y[i];
			  zz[i] = z[i];
			}
		      tag++;
		      MPI_Recv (&(x[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD, &status);
		      tag++;
		      MPI_Recv (&(y[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD, &status);
		      tag++;
		      MPI_Recv (&(z[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD, &status);
		      tag++;
		      MPI_Send (&(xx[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD);
		      tag++;
		      MPI_Send (&(yy[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD);
		      tag++;
		      MPI_Send (&(zz[0]), n, MPI_DOUBLE, swi, tag,
				MPI_COMM_WORLD);
		      uu = cluster_energy (x, y, z, n);
		    }
		}		//if(swi!=0)
	    }			//if(k%swap_freq)

	  //Saving energy histograms and mean values
	  if ((k % save_freq) == 0)
	    {
	      av_save = fopen (average_save, "a");
	      cf_save = fopen (conf_save, "a");
	      fileU = fopen (charU, "w");
	      for (l = 0; l < n; l++)
		{
		  fprintf (cf_save, "%e %e %e\n", x[l], y[l], z[l]);
		}
	      fprintf (av_save, "%.16e %.16e %.16e %.16e %.16e %.16e ",
		       counter, countrej, outmoves, normU, U, U2);
#if NUMBER_OF_DOPANTS
	      fprintf (av_save, "%.16e %.16e", exc_tried, exc_rejected);
#endif
	      fprintf (av_save, "\n");
	      for (i = 0; i < nU; i++)
		{
		  fprintf (fileU, "%d\n", histU[i]);
		}
	      fclose (fileU);
	      fclose (cf_save);
	      fclose (av_save);
	    }
	}			//for(mc_steps)

      free (x);
      free (y);
      free (z);
      free (xx);
      free (yy);
      free (zz);
      free (histU);
    }
  MPI_Finalize ();
  return 0;
}
