#include <stdio.h>


//*********************************************
//Matrix Creation and Initialization functions
//*********************************************

void mat_create(double **&mat, int row, int col)
{
	mat=new double*[row];
	for (int i=0; i<row; i++)
	{
		mat[i]=new double[col];
	}
}

void mat_init(double **mat, int row, int col)
{
	for (int i= 0; i < row; ++i) 
	{
		for (int j= 0; j < col; ++j) 
		{
			mat[i][j] = 0.0;
		}
	}
} 

void array_create(double *&array, int len)
{
	array = new double [len];
}

void array_init(double *array, int len)
{
	for (int i=0; i<len; i++)
	{
		array[i]=0.0;
	}
}
void mat_print(double **mat, int rows, int cols)
{
	for (int i=0; i<rows; i++)
	{
		for (int j=0; j<cols; j++)
		{
			printf( "%7.3f   ",mat[i][j]);
		}
	printf("\n");
	}
}
void array_print(double *array, int len)
{
	for (int i=0; i<len; i++)
	{
		printf("%2.3f   ", array[i]);
	}
	printf("\n\n");
}
void print_gsl_mat(double *m, int rows, int cols)
{
	for (int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
		{
			printf("%7.3f    ",m[j+cols*i]);
		}
	printf("\n");
	}
	printf("\n");
}
void print_mat(double **m, int rows, int cols)
{
	for (int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
		{
			printf("%7.3f    ",m[i][j]);
		}
	printf("\n");
	}
	printf("\n");
}
void copy_gsl_mat(double **cmat, double *gmat, int rows, int cols)
{
	for (int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
		{
			gmat[j+cols*i]=cmat[i][j];
		}
	}
}
void print_vect_in_mat(double *v, int rows, int cols)
{
	for (int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
		{
			printf("%7.3f    ",v[j+cols*i]);
		}
	printf("\n");
	}
	printf("\n");
}
