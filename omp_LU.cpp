#include <omp.h>
#include <stdio.h>
#include <iostream>

void array_create(double *&array, int len);
void array_init(double *array, int len);
void mat_create(double **&mat, int row, int col);
void mat_init(double **mat, int row, int col);

double** LU_decompose(double **mat, int nx, int ny)
{
  int size=nx*ny;

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

double* LU_solve(double **mat, double *vec, int nx, int ny)
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
  for(int i=1; i<N; i++) {
    temp = 0;
    for(int j=0; j<i; j++) {
      temp += mat[i][j] * y[j];
    }
    y[i] = vec[i] - temp;
  }
  x[N-1] = y[N-1] / mat[N-1][N-1];

  for(int i=N-2; i>-1; --i) {
    temp = 0;
    for(int j=i+1; j<N; j++) {
      temp += mat[i][j] * x[j];
    }
    x[i] = (y[i] - temp) / mat[i][i];
  }
  return x;
}


double** omp_LU_decompose(double **mat, int nx, int ny)
{
  
#pragma omp parallel shared(mat)
  {
  int size=nx*ny;
  for(int j=0; j<size; j++){
#pragma omp for schedule(static)
    for(int i=j+1; i<size; i++){
      mat[i][j] = mat[i][j] / mat[j][j];
	for(int k=j+1; k<size; k++){
	  mat[i][k] = mat[i][k] - mat[j][k] * mat[i][j];
	}
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


  for(int i=1; i<N; i++) {
  temp = 0.0;
#pragma omp parallel if(i>400) shared(y,mat,temp) 
  {
    #pragma omp for reduction(+:temp) schedule(static)
    for(int j=0; j<i; j++) {
      temp += mat[i][j] * y[j];
    }
  }
    y[i] = vec[i] - temp;
  }

  x[N-1] = y[N-1] / mat[N-1][N-1];

  for(int i=N-2; i>-1; --i) {
    temp = 0;
#pragma omp parallel if(N-i+1>400) shared(x,mat,temp) 
  {
    #pragma omp for reduction(+:temp) schedule(static)
    for(int j=i+1; j<N; j++) {
      temp += mat[i][j] * x[j];
    }
  }
    x[i] = (y[i] - temp) / mat[i][i];
  }
  return x;
}
