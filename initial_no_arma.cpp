//This is the main file that does not use armadillo anymore


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <complex>
#include <mpi.h>
#include <omp.h>

#include "papi_support.h"
//*****************************************
double A[N][M];




//*********************************************
//The LU decomposition fucntions should be defined here






//*********************************************

void init_cond(arma::mat &x,arma::mat &y, arma::mat &u)
{
	//this function return the initial position
	 u=exp(-10.0*(pow(x,2.0) + pow(y,2.0)));
}

 void meshgrid(arma::vec &x_lin, arma::vec &y_lin,arma::mat &xhalf, arma::mat &yhalf)
{
	//this function creates the mesh
	xhalf.each_row() +=x_lin.t();
	yhalf.each_col() +=y_lin;
}

void arange(arma::vec &x,int ini, int fin, int steps)
{
	//figure out the number of elments
	int counter=0;
	for (int i=ini; i<fin; i+=steps)
	{
		counter++;
	}
	int k=0;
	x=arma::zeros(counter);
	//initialize the vector and assign the values
	for(int j=ini; j<fin; j+=steps)
	{
		x(k)=j;
		k++;
	}	

}

void linspace()
{


}




//*********************************************

//This is to test that the LU decomposition functions from serial
//parallel and MPI are the same
void unit_test()
{
  for (int i = 0; i < N; i++) 
  {
    	for (int j = 0; j < N; j++) 
	{
      		if (abs(correct_results[i][j] - C[i][j]) > ERROR_THRESHOLD)
      		{
        		printf("Incorrect result. (%d, %d) %.5f %.5f\n", i, j, correct_results[i][j], C[i][j]);
        		return;
      		}
    	}
  }
  printf("GOOD -- Test passed.\n");
}

int main (int argc, char** argv)
{
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

#ifdef USE_MPI
  if (rank == 0) {
#endif
  papi_init();
  reset_array();

 /**
   * THIS IS WHERE YOU CHOOSE WHICH COUNTERS TO RUN.
   * ONLY ADD AT MOST 2 COUNTERS PER RUN.
   *
   * The current example adds L2 Data Cache Misses and L2 Data Cache Accesses.
   */


#ifndef NOPAPI
  papi_prepare_counter(PAPI_L2_DCM);
  papi_prepare_counter(PAPI_L2_DCA);
#endif

  int nx=4; //This is the number of cell TRY to set it in the Makefile
  int nx=ny;

  //setting the domain of -1 to 1
  double dx=2.0/nx;
  double dy=dx;

  //declaring the mesh grid points

  























printf("Running the basic version.\n");
  papi_start();
  {
    basic_MM();
  }
  papi_stop_and_report();

  memcpy(correct_results, C, sizeof(C));

#ifdef USE_MPI
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  reset_array();

#ifdef USE_OMP
#pragma omp parallel
  {
  #pragma omp master
      printf("Running your openMP version using %d threads\n", omp_get_num_threads());
  }

  papi_start();
  {
    openMP_MM();
  }
  papi_stop_and_report();
  unit_test();
#endif

#ifdef USE_MPI
  reset_array();
  if (rank == 0) {
      printf("Running mpi version with %d processes.\n", size);
      papi_start();
  }
  {
    mpi_MM();
  }

  if (rank == 0) {
      papi_stop_and_report();
  }

  if (rank == 0) unit_test();
  MPI_Finalize();
#endif
}




















