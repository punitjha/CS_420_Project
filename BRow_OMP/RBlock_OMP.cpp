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
#include <ctime>
#include <stdlib.h>
#include <papi.h>
//We new set the num_threads in main or in the run.sub
#ifdef USE_PAPI
float gettime() 
{
	return((float)PAPI_get_virt_usec());
}

#define NUM_EVENTS 4 
#endif



void LU_OMP(double **mat, int size)
{
  int rows,th_min,th_max;
  int pid=0;
  int nprocs;
#pragma omp parallel for schedule(static,100)
     for(int k=0;k<size;k++){
	     for(int j=k+1;j<size;j++){
		     mat[j][k]=mat[j][k]/mat[k][k];
	     }
#pragma omp parallel for schedule(static,100)
     for(int i=k+1; i<size; i++){
	     for(int j=k+1;j<size;j++){
		     mat[i][j]=mat[i][j]-mat[i][k]*mat[k][j];
		     }
	     }
	}
}

//**************************************************************
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
		for (int j= 0; j < col; j++) 
		{
			if (i==0){
			    mat[i][j] =13;
			} 
			else{
			    mat[i][j] =3;
			    mat[j][i] =6;
			}
		}
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

//int check (double **mat, arma::mat A, int size)
//{
//	for(int i=0; i<size; i++)
//	{
//		for(int j=0; j<size; j++)
//		{
//
//			if(fabs(mat[i][j]-A(i,j)) > 0.001)
//			{
//				printf("The mismatch spot  %11.7f   %11.7f  %11.7f  %d   %d ",fabs(mat[i][j]-A(i,j)),mat[i][j],A(i,j),i,j);
//				return 1;
//			}
//		}
//
//	}
//	return 0;
//}



int main(int argc, char **argv)
{
#ifdef USE_PAPI
	int ret,EventSet= PAPI_NULL;	
	long long values[NUM_EVENTS];
	long long start_usec,end_usec;
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
	
	int row=atoi(argv[1]);
	int col=atoi(argv[1]);
	double **mat=NULL;
	mat_create(mat,row,col);
	mat_init(mat,row,col);

#ifdef USE_PAPI
	//float t2 = gettime();
	start_usec=PAPI_get_real_usec();
	if ((ret = PAPI_start_counters(events, NUM_EVENTS)) != PAPI_OK) {
		   fprintf(stderr, "PAPI failed to start counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}
#endif
//************************************


	LU_OMP(mat,row);

//**************
//PAPI
//**************

#ifdef USE_PAPI
	//float t4 = gettime();
	end_usec=PAPI_get_real_usec();
	if (PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK){
		   fprintf(stderr, "PAPI failed to read counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}

	printf("\n \n This for the LAPACK  version \n");
	printf("Time in Microsecond %lld\n",(end_usec-start_usec));
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





