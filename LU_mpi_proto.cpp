#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

void LU_MPI(double **mat, int size)
{
   int i,j,k;
   int rank, nprocs;

   MPI_Init(NULL, NULL);
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   for (i = 0; i < size-1; i++){
        for (j = i+1; j < size; j++){
                if (j % p == id){
                    if(mat[j][i] != 0){
                    // if mat[j][i] is 0, not action required
                        double factor = mat[j][i]/mat[i][i]

                                for (k = i+1; k < size; k++){
                                mat[j][k] = mat[j][k] - factor*mat[i][k]
                                }
                                mat[j][i] = factor;
                    }
                }
        }

        for (j = i+1; j < size; j++){
        MPI_Bcast(&mat[j][i], size-i, MPI_DOUBLE, j % p, MPI_COMM_WORLD);
        }

   }

   MPI_Finalize();
   return;
}
