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
#include <mpi.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <hdf5.h>

#ifdef USE_PAPI
#include <papi.h>
#endif

void output_h5(int nx, int ny, gsl_matrix* sol);
void init_h5(double *& lin_x, double *& lin_y, int nx, int ny);
//*****************************************



//******************************************************
//The LU decomposition fucntions should be defined here
//******************************************************







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

//****************************************************************************
//The main function starts here
//****************************************************************************


int main (int argc, char** argv)
{
//****************************
//Declaring and Initialing MPI
//*****************************
//	MPI_Init (&argc, &argv);
//	int rank, size;
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//******************
//PAPI definiations 
//******************

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
	
//****************************
//PAPI is now started initialized
//***************************

//	if (rank==0){
//	printf("This is the MPI Testing");
//	float t0 = gettime();
//	if ((ret = PAPI_start_counters(events, NUM_EVENTS)) != PAPI_OK) {
//		   fprintf(stderr, "PAPI failed to start counters: %s\n", PAPI_strerror(ret));
//		      exit(1);
//	}
//	if (rank ==0){
//	float t1 = gettime();
//	if (PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK){
//		   fprintf(stderr, "PAPI failed to read counters: %s\n", PAPI_strerror(ret));
//		      exit(1);
//	}
//	printf("Time in Microsecond %f\n",(t1-t0));
//	printf("L1 data cache misses is %lld\n", values[0]);
//	printf("L2 data cache misses is %lld\n", values[1]);
//	printf("L1 data cache accesses is %lld\n", values[2]);
//	printf("L2 data cache accesses is %lld\n", values[3]);
//	}
//	}
//	MPI_Finalize();


#ifdef USE_PAPI
	float t2 = gettime();
	if ((ret = PAPI_start_counters(events, NUM_EVENTS)) != PAPI_OK) {
		   fprintf(stderr, "PAPI failed to start counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}
#endif
//************************************
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
	
	meshgrid(x_lin,y_lin,xhalf,yhalf,nx);
//	mat_print(xhalf,nx,nx);
//	mat_print(yhalf,nx,nx);	

	// create h5 file and output coordinates to h5 file	
	double *x_coord, *y_coord;
	array_create(x_coord, nx);
	array_create(y_coord, ny);
	for (int i=0; i<nx; i++) {
	  x_coord[i] = x_lin[i] + dx/2.0;
	  y_coord[i] = y_lin[i] + dy/2.0;
	}
	init_h5(x_coord, y_coord, nx, ny);

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

	double cx=1.0; // flow speed in x
	double cy=1.0; // flow speed in y
	int nt=32; // no. of time steps 
	double dt= T/nt; //the delta t -- time steps
	double lmbda_x = cx*dt/dx; // courant number in x direction
	double lmbda_y = cy*dt/dy; // courant number in y direction
	int svl = nx*ny;// the sol vector length

	double *x_boundary, *y_boundary;
	int x_len=arange(x_boundary,nx,svl,nx);
	int y_len=arange(y_boundary,ny,svl,ny);
	array_add_subs(x_boundary,nx,-1);//each row boundary
	array_add_subs(y_boundary,ny,-1);//each col boundary -1 to start from index 0
	
	double *hpb,*vpb;
	int hpb_len=arange(hpb,nx,svl+1,nx);
	int vpb_len=arange(vpb,ny,svl+1,ny);
	array_add_subs(hpb,nx,-1);//each row boundary
	array_add_subs(vpb,ny,-1);//each col boundary -1 to start from index 0
	
	//the itirations over time steps start
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
	for(int i=0; i<x_len; i++)
	{
		m1[(int)x_boundary[i]][(int)y_boundary[i]+1]=0.0;
		m1[(int)x_boundary[i]+1][(int)y_boundary[i]]=0.0;
	}
	
	//enforce peridic boundary conditions
	for(int i=0; i<hpb_len; i++)
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

	//mat_print(m2,nx*ny,nx*ny);	
	
//********************
//Reading PAPI values
//********************

#ifdef USE_PAPI
	float t2_1 = gettime();
	if (PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK){
		   fprintf(stderr, "PAPI failed to read counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}
	printf("Time in Microsecond %f\n",(t2_1-t2));
	printf("L1 data cache misses is %lld\n", values[0]);
	printf("L2 data cache misses is %lld\n", values[1]);
	printf("L1 data cache accesses is %lld\n", values[2]);
	printf("L2 data cache accesses is %lld\n", values[3]);
	
	float t3 = gettime();
	if ((PAPI_reset(EventSet)) != PAPI_OK){
		   fprintf(stderr, "PAPI failed to reset  counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}
#endif
//******************************************
//we have to code our LU decomposition here*
//******************************************

	//Using the GNU standard library for solving the linear equation	
	gsl_matrix *mm1=gsl_matrix_alloc(nx*ny,nx*ny);
	gsl_matrix *mm2=gsl_matrix_alloc(nx*ny,nx*ny);
	gsl_matrix *uu1=gsl_matrix_alloc(nx,ny);
	
	copy_gsl_mat(m1,gsl_matrix_ptr(mm1,0,0),nx*ny,nx*ny);
	copy_gsl_mat(u,gsl_matrix_ptr(uu1,0,0),nx,ny);
	gsl_matrix_transpose_memcpy(mm2,mm1);

	gsl_vector *u_vec=gsl_vector_alloc(nx*ny);

	//printf("\n Here1 \n");
	for(int i=0; i<10; i++)
	{
		//flatening starts here
	        //printf("\n Here2 \n");
		flatten_gsl(gsl_matrix_ptr(uu1,0,0),gsl_vector_ptr(u_vec,0),nx,ny,nx*ny);	
		//the matrix vector multiplication implemtation.NOTICE NEW VEC USED
		gsl_vector *u_vecc=gsl_vector_alloc (nx*ny);
		flatten_gsl(gsl_matrix_ptr(uu1,0,0),gsl_vector_ptr(u_vecc,0),nx,ny,nx*ny);	
		//printf("\n Here3 \n");
		gsl_blas_dgemv(CblasNoTrans,1.0,mm2,u_vec,0.0,u_vecc);

		//preparing for gsl LU  decompose
		//printf("\n Here4 \n");
		int s; // required signum
		gsl_permutation *p=gsl_permutation_alloc(nx*ny);
		//printf("\n Here5 \n");
		gsl_linalg_LU_decomp (mm1,p,&s);
		gsl_vector *sol=gsl_vector_alloc(nx*ny);
		//printf("\n Here6 \n");
		gsl_linalg_LU_solve (mm1,p,u_vecc,sol);
		//printf("\n Here7 \n");
		gsl_matrix_view sol1=gsl_matrix_view_vector(sol,nx,ny);
		//printf("\n Here8 \n");
		//		print_gsl_mat(gsl_matrix_ptr(&sol1.matrix,0,0),nx,ny);
		//		printf(" \n \n \n");
		//printf("\n Here9 \n");
		gsl_matrix_memcpy(uu1,&sol1.matrix);
		//printf("\n Here10 \n");
	}
	print_gsl_mat(gsl_matrix_ptr(uu1,0,0),nx,ny);
	output_h5(nx,ny,uu1);
//**************
//PAPI
//**************

#ifdef USE_PAPI
	float t4 = gettime();
	if (PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK){
		   fprintf(stderr, "PAPI failed to read counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}
	printf("Time in Microsecond %f\n",(t4-t3));
	printf("L1 data cache misses is %lld\n", values[0]);
	printf("L2 data cache misses is %lld\n", values[1]);
	printf("L1 data cache accesses is %lld\n", values[2]);
	printf("L2 data cache accesses is %lld\n", values[3]);
	if (PAPI_stop_counters(values, NUM_EVENTS ) != PAPI_OK){
		fprintf(stderr, "PAPI_stoped_counters - FAILED \n");
		exit(1);
	}
#endif
	return 0;	


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

