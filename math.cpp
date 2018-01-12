#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "matrix.cpp"

void array_create(double *&array, int len);
void array_init(double *array, int len);
//*************************************************************************
//Mathematical functions
//*************************************************************************
void init_cond(double **mat1,double **mat2, double **mat3, int row, int col)
{
	//this function return the initial position
	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
		{
			mat3[i][j]=exp(-10.0*(pow(mat1[i][j],2.0) + pow(mat2[i][j],2.0)));
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

void linspace(double *& array, double a, double b, int n)
{
	array_create(array,n);
	array_init(array,n);
	printf("\n \n");
	double epsilon=0.0001;
	double step=(b-a)/(n-1);
	int iti=0;
	if((n==0) || (n==1) || (a==b))
	{
		array[n-1]=b;
	}
	else if (n>1 && step >=0)
	{
		while (a<=b+epsilon)
		{
			array[iti]=a;
			iti++;
			a+=step;
		}
	}
	else
	{
		while(a+epsilon >= b)
		{
		array[iti]=a;
		iti++;
		a+=step;
		}
	}
}



void meshgrid(double *x_lin, double *y_lin, double **xhalf, double **yhalf, int len)
{
	//this function creates the mesh


	for (int i=0; i<len; i++)
	{
		for(int j=0; j<len; j++)
		{
			xhalf[i][j]=x_lin[j];
			yhalf[i][j]=y_lin[i];
		}
	}
}

		
int arange(double *&x,int ini, int fin, int steps)
{
	//figure out the number of elments
	int counter=0;
	for (int i=ini; i<fin; i+=steps)
	{
		counter++;
	}
	array_create(x,counter);
	array_init(x,counter);
	//initialize the vector and assign the values
	int k=0;
	for(int j=ini; j<fin; j+=steps)
	{
		x[k]=j;
		k++;
	}
	return counter;	

}
//Function to add or substract numbers from each element of the matrix
//Pass a negative number if substraction is desired.
void mat_add_subs(double **mat,int rows, int cols, double num)
{
	for(int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
		{
			mat[i][j]+=num;
		}
	}
}

void array_add_subs(double *array, int len, double num)
{
	for (int i=0; i<len; i++)
		array[i]+=num;
}

void diag(double **mat, double *array, int len, double factor, int diag)
{
	if(diag > 0)
	{
		for (int i=0; i<len; i++)
		{
			mat[i][i+diag]=array[i]*factor;
		}
	}
	if(diag ==0)
		{
		for (int i=0; i<len; i++)
		{
			mat[i][i]=array[i]*factor;
		}
	}
	if(diag <0)
	{
		for (int i=0; i<len; i++)
		{
			mat[i+abs(diag)][i]=array[i]*factor;
		}
	}
}
		
void transpose_sq(double **mmat2 ,double **mmat1, int rows, int cols)
{
	for(int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
		{
			mmat2[i][j]=mmat1[j][i];
		}
	}
}
void flatten_gsl(double *mat, double *vect, int rows, int cols, int len)
{
	int k=0;
	for (int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
		{
			vect[k]=mat[i+j*cols];
			k++;
		}
	}
}
