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

// Initialize array to zeros.
void reset_array()
{
  for (int i = 0; i < N; i++)
    for (int j = 0; j < M; j++)
    {
        A[i][j] = 0.0;
    }
}


void mat_create(double **mat, int row, int col)
{
	 mat=new double*[row];
	for (int i=0; i<row; i++)
	{
		mat[i]=new double[col];
	}
}

void array_create(double *array, int len)
{
	array = new double [len];
}


void init_cond(double **mat1,double **mat2, double **mat3, int row, int col)
{
	//this function return the initial position
	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
		{
			mat3[i][j]=exp(-10.0*(pow(mat1[i][j],2) + pow(mat2[i][j],2)));
		}
	}
	 //u=exp(-10.0*(pow(x,2.0) + pow(y,2.0)));
}

//works for square matrices only
void mat_mul(double **A, double **B, double **C, int row, int col)
{ 
    int i, j, k;
    for (i = 0; i < row; ++i) 
    {
        for (j = 0; j < col; ++j) 
	{
            C[i][j] = 0;
            for (k = 0; k < row; ++k) 
	    {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void linspace(double *array, double a, double b, int n)
{
	double step=(a-b)/(n-1);
	int i=0;
	while(a<=b)
	{
		array[i]=a+step;
		i++;
		a+=step;
	}
}

void transpose(double **mat, int row, int col)
{
	for(int i=0; i<row; i++)
	{
		for (int 


 void meshgrid(double *x_lin, double *y_lin, double **xhalf, double **yhalf)
{
	//this function creates the mesh

	xhalf.each_row() +=x_lin.t();
	yhalf.each_col() +=y_lin;
}



//*********************************************
//The LU decomposition fucntions should be defined here








//*********************************************
void linspace()
{


}




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


//************************************
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




















