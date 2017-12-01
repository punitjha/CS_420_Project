#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#include "papi_support.h"

double A[N][M];
double B[M][N];
double C[N][N];

// Used for testing.
double correct_results[N][N];

void basic_MM();
void openMP_MM();
void mpi_MM();

// Initialize array to zeros.
void reset_array()
{
  for (int i = 0; i < N; i++)
    for (int j = 0; j < M; j++) {
        A[i][j] = 1.0;
        B[j][i] = 2.0;
    }
  memset(C, 0.0, sizeof(C));
}

void unit_test()
{
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
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

  printf("Running basic sequential version to save the results.\n");

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

  /** ============================================== */
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
