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
