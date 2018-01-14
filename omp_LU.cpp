#include <omp.h>
#include <stdio.h>
#include <iostream>

void array_create(double *&array, int len);
void array_init(double *array, int len);
void mat_create(double **&mat, int row, int col);
void mat_init(double **mat, int row, int col);

void LU_OMP(double **mat, int size)
{
  int rows,th_min,th_max;
  int pid=0;
  int nprocs;
  // #pragma omp parallel for schedule(static)
     for(int k=0;k<size;k++){
	     for(int j=k+1;j<size;j++){
		     mat[k][j]=mat[k][j]/mat[k][k];
	     }
	     // #pragma omp parallel for schedule(static)
     for(int i=k+1; i<size; i++){
	     for(int j=k+1;j<size;j++){
		     mat[i][j]=mat[i][j]-mat[i][k]*mat[k][j];
		     }
	     }
	}
}

double** omp_LU_decompose(double **mat, int nx, int ny)
{
  int size=nx*ny;
  #pragma omp parallel for schedule(static)
  for(int j=0; j<size; j++){
    for(int i=j+1; i<size; i++){
      mat[i][j] = mat[i][j] / mat[j][j];
	for(int k=j+1; k<size; k++){
	  mat[i][k] = mat[i][k] - mat[j][k] * mat[i][j];
	}
    }
  }
  return mat;
}

double* omp_LU_solve(double **mat, double *vec, int nx, int ny)
{
  int N=nx*ny;
  double *y;
  array_create(y, N); 
  array_init(y, N);	     
  double *x;
  array_create(x, N); 
  array_init(x, N);	     
  double temp=0;

  
  y[0] = vec[0];
  #pragma omp parallel for schedule(static)
  for(int i=1; i<N; i++) {
    temp = 0;
    for(int j=0; j<i; j++) {
      temp += mat[i][j] * y[j];
    }
    y[i] = vec[i] - temp;
  }
  x[N-1] = y[N-1] / mat[N-1][N-1];

  #pragma omp parallel for schedule(static)
  for(int i=N-2; i<N && i>-1; --i) {
    temp = 0;
    for(int j=i+1; j<N; j++) {
      temp += mat[i][j] * x[j];
    }
    x[i] = (y[i] - temp) / mat[i][i];
    }
  return x;
}
