#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define ln() putchar('\n')
#define GENERIC_TAG (0)

#ifdef USE_PAPI
#include <papi.h>
#endif

//**************************************
//Defining the time function with PAPI
//**************************************

//The  time is in microseconds

#ifdef USE_PAPI
float gettime() 
{
	return((float)PAPI_get_virt_usec());
}

#define NUM_EVENTS 4 
#endif



double *mat_create_init(int size);
void mat_print(double *m, int size);
void mat_print_L(double *m, int size);
void mat_print_U(double *m, int size);

int main(int argc, char**argv) 
{
#ifdef USE_PAPI
	int ret,EventSet= PAPI_NULL;	
	long long values[NUM_EVENTS];
	ret =PAPI_library_init(PAPI_VER_CURRENT);
	if (PAPI_create_eventset(&EventSet) != PAPI_OK){
		   fprintf(stderr, "No hardware counters here, or PAPI not supported.\n");
		      exit(1);
	}

	int events[NUM_EVENTS] = {PAPI_L1_DCM, PAPI_L2_DCM, PAPI_L1_DCA, PAPI_L2_DCA };
	if (PAPI_add_events(EventSet,events,NUM_EVENTS ) != PAPI_OK){
		   fprintf(stderr, "No hardware counters here, or PAPI not supported.\n");
		      exit(1);
	}
#endif


      const int root_p = 0;
         int mx_size = 0, p, id;
	    if (argc < 2) {
	          printf("Matrix size missing in the arguments\n");
		        return EXIT_FAILURE;
			   }
   int size = atoi(argv[1]);
#ifdef USE_PAPI
	float t2 = gettime();
	if ((ret = PAPI_start_counters(events, NUM_EVENTS)) != PAPI_OK) {
		   fprintf(stderr, "PAPI failed to start counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}
#endif

   double *mat = mat_create_init(size);  /* Create and initialize matrix */

   /********************** LU decomposition using MPI ************************/
   int rank, nprocs;

   MPI_Init(NULL, NULL);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   /* Print the origin matrix */
   if (rank == 0) {
      #ifdef PRINT 
      printf("The original matrix:\n");
      mat_print(mat, size);
      ln();
      #endif
   }

   /* Decompose the matrix */
   int i, j, k;
   for (i = 0; i < size - 1; i++) {
       for (j = i + 1; j < size; j++) {
            if (j % nprocs == rank) {
                if(mat[j * size + i] != 0){
                    /* If mat[j][i] is 0, not action required */
                    double factor = mat[j * size + i] / mat[i * size + i];
                    for (k = i + 1; k < size; k++) {
                        mat[j * size + k] = mat[j * size + k] - factor * mat[i * size + k];
                    }
                    mat[j * size + i] = factor;
                }
            }
        }

        for (j = i + 1; j < size; j++) {
            MPI_Bcast(&mat[j * size + i], size - i, MPI_DOUBLE, j % nprocs, MPI_COMM_WORLD);
        }
   }

   /* Print the decomposed L and U */
   if (rank == 0) {
      #ifdef PRINT
      printf("\n Lower matrix: \n");
      mat_print_L(mat, size);
      printf("\n Upper matrix: \n");
      mat_print_U(mat, size);
      #endif
   }
   free(mat);
   if (rank==0){
#ifdef USE_PAPI
	float t4 = gettime();
	if (PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK){
		   fprintf(stderr, "PAPI failed to read counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}

	printf("\n \n This for the Serial  version \n");
	printf("Time in Microsecond %f\n",(t4-t2));
	printf("L1 data cache misses is %lld\n", values[0]);
	printf("L2 data cache misses is %lld\n", values[1]);
	printf("L1 data cache accesses is %lld\n", values[2]);
	printf("L2 data cache accesses is %lld\n", values[3]);
	if (PAPI_stop_counters(values, NUM_EVENTS ) != PAPI_OK){
		fprintf(stderr, "PAPI_stoped_counters - FAILED \n");
		exit(1);
	}
#endif
	}

   MPI_Finalize();
   return 0;
}

/*************************** Supporting funcitons ****************************/
/* Used to create and initialize the matrix */
double *mat_create_init(int size) {
    int total_size = size * size;
    int i, j;
    double *m = (double *)malloc(sizeof(double) * total_size);
    for (i = 0; i < size; i++) {
        j = i;
        for (j = i; j < size; j++) {
            m[i * size + j] = i + 1;
            m[j * size + i] = i + 1;
        }
    }
    return m;
}

/* Used to print out the matrix */
void mat_print(double *m, int size) {
    int total_size = size * size;
    int i, j;
    for (i = 0; i < total_size; i++) {
        printf("%6.3f\t", m[i]);
        if ((i + 1) % size == 0) {
            ln();
        }
    }
}

/* Used to print lower matrix */
void mat_print_L(double *m, int size) {
    int i, j;
    double z = 0;
    double o = 1;
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            if (j > i) {
                printf("%6.3f\t", z);
            } else if (j == i) {
                printf("%6.3f\t", o);
            } else {
                printf("%6.3f\t", m[i * size + j]);
            }
        }
        ln();
    }
}

/* Used to print the upper matrix */
void mat_print_U(double *m, int size) {
    int i, j;
    double z = 0;
    for (i = 0; i < size; i++) {
      for (j = 0; j < size; j++) {
          if (j >= i) {
              printf("%6.3f\t", m[i * size + j]);
          } else {
              printf("%6.3f\t", z);
          }
      }
      ln();
    }
}

