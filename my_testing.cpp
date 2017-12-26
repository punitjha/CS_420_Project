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
#include <omp.h>
#include <armadillo>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

//*****************************************




//******************************************************
//The LU decomposition fucntions should be defined here
//******************************************************




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
			printf( "%7.3f \t ",mat[i][j]);
		}
	printf("\n");
	}
}
void array_print(double *array, int len)
{
	for (int i=0; i<len; i++)
	{
		printf("%2.3f \t", array[i]);
	}
}

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

		
void arange(double *&x,int ini, int fin, int steps)
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


//This is to test that the LU decomposition functions from serial
//parallel and MPI are the same

//void unit_test()
//{
//  for (int i = 0; i < N; i++) 
//  {
//    	for (int j = 0; j < N; j++) 
//	{
//      		if (abs(correct_results[i][j] - C[i][j]) > ERROR_THRESHOLD)
//      		{
//        		printf("Incorrect result. (%d, %d) %.5f %.5f\n", i, j, correct_results[i][j], C[i][j]);
//        		return;
//      		}
//    	}
//  }
//  printf("GOOD -- Test passed.\n");
//}
int main (int argc, char** argv)
{
	int nx;
	printf("\n Enter nx \n");
	scanf("%d",&nx);
	//int nx=4; //This is the number of cells 
	int ny=nx; //No. of grids in both the dimensions are set equal

	//Setting the domain of -1 to 1
	double dx=2.0/nx;
	double dy=dx;

	//declaring mesh grid points
	double **xhalf=NULL,**yhalf=NULL;
	mat_create(xhalf,nx,ny);
	mat_create(yhalf,nx,ny);
	mat_init(xhalf,nx,ny);
	mat_init(yhalf,nx,ny);
	
	double *x_lin, *y_lin;
	linspace(x_lin,-1.0,1.0-dx,nx);	
	linspace(y_lin,-1.0,1.0-dx,nx);	
	//array_print(x_lin,nx);
	meshgrid(x_lin,y_lin,xhalf,yhalf,nx);
//	printf("\n Printing x half \n");
//	mat_print(xhalf,nx,nx);
//	printf("\n Printing y half \n");
//	mat_print(yhalf,nx,nx);	

	double **x, **y;
	mat_create(x,nx,ny);
	mat_create(y,nx,ny);
	mat_init(x,nx,ny);
	mat_init(y,nx,ny);	
	for (int i=0; i<nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
			x[i][j]=xhalf[i][j]+dx/2.0;
			y[i][j]=yhalf[i][j]+dy/2.0;
		}
	}
//	printf("This is x \n ");
//	mat_print(x,nx,ny);
//	printf("This is y\n ");
//	mat_print(y,nx,ny);	
	double T=1.0; // run simulations for T periods

	//Building the initial condition matrix
	double **u;
	mat_create(u,nx,ny);
	mat_init(u,nx,ny);
	init_cond(x,y,u,nx,ny);
	//printf("This is u\n ");
	//mat_print(u,nx,ny);

	double cx=1.0; // flow speed in x
	double cy=1.0; // flow speed in y
	int nt=32; // no. of time steps 
	double dt= T/nt; //the delta t -- time steps
	double lmbda_x = cx*dt/dx; // courant number in x direction
	double lmbda_y = cy*dt/dy; // courant number in y direction
	int svl = nx*ny;// the sol vector length

	double *x_boundary, *y_boundary;
//	array_create(x_boundary,nx);
//	array_create(y_boundary,ny);
//	array_init(x_boundary,nx);
//	array_init(y_boundary,ny);
	arange(x_boundary,nx,svl,nx);
	arange(y_boundary,ny,svl,ny);
	//printf("This is x_boundary \n ");
	//array_print(x_boundary,3);
	array_add_subs(x_boundary,nx,-1);//each row boundary
	array_add_subs(y_boundary,ny,-1);//each col boundary -1 to start from index 0
	
	double *hpb,*vpb;
	arange(hpb,nx,svl+1,nx);
	arange(vpb,ny,svl+1,ny);
	array_add_subs(hpb,nx,-1);//each row boundary
	array_add_subs(vpb,ny,-1);//each col boundary -1 to start from index 0
	//printf("This is x_boundary \n ");
	//array_print(hpb,nx);
	
	//the itirations over time steps start
	//for(int time=0; time <nt; time ++)
	//{
		//creating the matrix m1
		double *l_a1, *l_a2, *I; // last one is identity
		array_create(l_a1,nx*ny-1);
		array_create(l_a2,nx*ny-nx);
		array_create(I,nx*ny);
		array_init(l_a1,nx*ny-1);
		array_init(l_a2,nx*ny-nx);
		array_init(I,nx*ny);
		array_add_subs(l_a1,nx*ny-1,lmbda_x/2.0);//each row boundary
		array_add_subs(l_a2,nx*ny-nx,lmbda_y/2.0);//each row boundary
		array_add_subs(I,nx*ny,1.0);//this is Identity

		double **m1;
		mat_create(m1,nx*ny,nx*ny);
		mat_init(m1,nx*ny,nx*ny);
		diag(m1,l_a1,nx*nx-1,-1.0,-1);
		diag(m1,l_a1,nx*nx-1,1.0,1);		
		diag(m1,l_a2,nx*nx-nx,-1.0,-nx);		
		diag(m1,l_a2,nx*nx-nx,1.0,nx);	
		diag(m1,I,nx*nx,1.0,0);	
		
		//removing the diag lambdas where each row ends and starts
		for(int i=0; i<sizeof(x)/sizeof(x[0]); i++)
		{
			m1[(int)x_boundary[i]][(int)y_boundary[i]+1]=0.0;
			m1[(int)x_boundary[i]+1][(int)y_boundary[i]]=0.0;
		}
		
		//enforce peridic boundary conditions
		
		for(int i=0; i<sizeof(hpb)/sizeof(hpb[0]); i++)
		{
			m1[(int)hpb[i]-(nx-1)][(int)hpb[i]]=-1*lmbda_x/2.0;
			m1[(int)vpb[i]][(int)vpb[i]-(ny-1)]=1*lmbda_x/2.0;
		}


		//mat_print(m1,nx*ny,nx*ny);	
		//assign the vertical periodic boudary conditions
		
		double *lmbda_v;
		array_create(lmbda_v,ny);
		array_init(lmbda_v,ny);
		array_add_subs(lmbda_v,ny,lmbda_y/2.0);	
		//array_print(lmbda_v,ny);

		diag(m1,lmbda_v,ny,-1.0,nx*(ny-1));
		diag(m1,lmbda_v,ny,1.0,-(nx*(ny-1)));
		
		//mat_print(m1,nx*ny,nx*ny);	

		double **m2;
		mat_create(m2,nx*ny,nx*ny);
		mat_init(m2,nx*ny,nx*ny);
		transpose_sq(m2,m1,nx*ny,nx*ny);

		mat_print(m2,nx*ny,nx*ny);	

//we have to code our LU decomposition here

		
	//Using the GNU standard library for solving the linear equation	
		gsl_matrix *mm1=gsl_matrix_alloc(nx*ny,nx*ny);
		gsl_matrix *mm2=gsl_matrix_alloc(nx*ny,nx*ny);

		for(int i=0; i<nx*ny; i++)
		{
			for (int j=0; j<nx*ny; j++)
			{
				gsl_matrix_set(mm1,i,j,m1[i][j]);
			}
		}
		

		//gsl_matrix_view m22= gsl_matrix_view_array( **m2, nx*ny,nx*ny);

	return 0;	


}



