#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <complex>

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



//=============================================================================

extern "C" {

extern void cblas_dswap(int  , double *, int , double *, int);

extern void cblas_daxpy(int , double , double *x, int , double *, int);

extern double cblas_dnrm2(const int, const double *, const int);


extern int  dgemm_(char *, char *, int *, int *, int *, double *, double *,
                      int *, double *, int *, double *, double *, int *);
extern void dgetrf_(int*, int*, double*, int*, int *,int*);

extern int  dtrsm_(char *, char *, char *, char *, int *, int *,
                      double *, double *, int *, double *, int *);
extern int  dsyrk_(char *, char *, int *, int *, double *, double *,
                      int *, double *, double *, int *);
}
//=============================================================================


//Initialize the matrix in column major layout
void mat_init(double *mat, double *A, int  N, int  M){
for(int i=0; i< M; i++){
	for(int j=0; j< N; j++){
	    if(i==j)
       A[j+i*N] = mat[j+i*N] =4.5*j/2;
      else 
       A[j+i*N] = mat[j+i*N] =4.5*i/4;
	  }
	}
}

//Priting the matrix elements
void mat_print(double *mat, int  width, int  M){
printf("\n");
printf("\n");
for(int i=0; i< M; i++){
	for(int j=0; j< width; j++){
	    printf("%6.3f   ",mat[i+j*width]);
	  }
	printf("\n");
	}
}

//This function does the actual copying of the data from source to destination
void convert_copy(int flag, double *source, double *dest, int M, int NBlocs){
	if (flag == 0){
		for (int i=0;i<NBlocs; i++){
			memcpy(dest, source, sizeof(double)*NBlocs);
			source +=M;
			dest+=NBlocs;
		}
	}
	else if(flag ==1){
		for(int i=0; i<NBlocs; i++){
			memcpy (dest, source, sizeof(double)*NBlocs);
			source +=NBlocs;
			dest+=M;
		}
	}
}


//This function converts from columns major layout to column tile layout
void ColToTile(double *mat, double *pmat, int row_chunk, int N, int NBlocs, int M){
	int rt=row_chunk/NBlocs;
	int ct=N/NBlocs;
	double *sptr, *dptr;
	for(int c=0; c<ct; c++){
		for(int r=0; r<rt; r++){
			sptr=pmat+c*NBlocs*M+r*NBlocs;
			dptr=mat+c*NBlocs*M+r*NBlocs*NBlocs;
			convert_copy(0,sptr,dptr,M,NBlocs);
		}
	}
}
//This function converts from column tile layout to column layout
void TileToCol(double *mat, double *pmat, int row_chunk, int N, int NBlocs, int M){
	int rt=row_chunk/NBlocs;
	int ct=N/NBlocs;
	double *sptr, *dptr;
	for(int c=0; c<ct; c++){
		for(int r=0; r<rt; r++){
			sptr=pmat+c*NBlocs*M+r*NBlocs;
			dptr=mat+c*NBlocs*M+r*NBlocs*NBlocs;
			convert_copy(1,dptr,sptr,M,NBlocs);
		}
	}
}


//this function swaps blocks of different sizes with each other based on the pivot element
void exchange_vecs(double *mat, int NBlocs, int k_start, int k_end,int *permut){
       for(int i=k_start; i<k_end; i++){
       		if(permut[i]-1 !=i){
			int val=i;
			int prev=permut[i]-1;
	cblas_dswap(NBlocs, mat+(val/NBlocs)*(NBlocs*NBlocs)+val%NBlocs,NBlocs, 
                            mat+(prev/NBlocs)*(NBlocs*NBlocs)+prev%NBlocs, NBlocs);
		}
       }
}

	 				


void block_LU(double *mat, int  M, int  N, int  NBlocs, int  *permut){
	int info;//required by the lapack interface
	int col_tiles=N/NBlocs;
	int row_tiles=M/NBlocs;
	double alpha=1.0;//required by the lapack functions;
	double neg = -1.0;//required by the lapack functions;
	double *store=(double*)malloc(M*N*sizeof(double));//tmp column tile storage
	//the pointers below store differnt blocks of matrix

	//Convert the matrix from the column major to tile major and store it in store;
	ColToTile(store,mat,M,N,NBlocs,M);


	{
	for (int k=0;k<col_tiles;k++){
	
	int row_chunk=M-k*NBlocs;

	double *source=mat+k*NBlocs*M+k*NBlocs;
	double *store_kk=store+k*NBlocs*M+k*NBlocs*NBlocs;

		//convert back from the tile format
      	TileToCol(store_kk,source,row_chunk,NBlocs,NBlocs,M);

      	dgetrf_(&row_chunk, &NBlocs, source, &M, permut+k*NBlocs, &info);
		
		for (int i=k*NBlocs;i<k*NBlocs+NBlocs; i++){
		  permut[i]+=k*NBlocs;
		}
      		
		ColToTile(store_kk,source,row_chunk,NBlocs,NBlocs,M);		

	for (int j=k+1; j<col_tiles;j++){
		 double *store_kj = store+j*NBlocs*M+k*NBlocs*NBlocs;
		 int rem_row_chunk=M-(k+1)*NBlocs;

		 
		int k_start=k*NBlocs;
	        int k_end=k*NBlocs+NBlocs;
		exchange_vecs(store+j*NBlocs*M, NBlocs,k_start,k_end,permut);

		
		 dtrsm_("l", "l", "n", "u", &NBlocs, &NBlocs,&alpha, store_kk, &NBlocs, store_kj, &NBlocs);
          
          // gemm
          for(int i= k+1; i< row_tiles; i++){
            dgemm_("n", "n", &NBlocs, &NBlocs, &NBlocs, &neg,store+k*NBlocs*M + i*NBlocs*NBlocs, &NBlocs, store_kj, &NBlocs, &alpha, store+j*NBlocs*M + i*NBlocs*NBlocs, &NBlocs);
          }
      }
}
     //pivoting now to left
     for(int ii=1; ii<col_tiles;ii++){
	    int start=ii*NBlocs;
	    int end=N;
		   exchange_vecs(store+(ii-1)*NBlocs*M,NBlocs,start,end,permut);
     }
    }//the parallel region ends here
    TileToCol(store,mat,M,N,NBlocs,M);
    
}
		
int main(int argc, char **argv)
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
	
	int info;
	int N=atoi(argv[1]);
	int M=atoi(argv[2]);
	int RBlocs=atoi(argv[3]);
	int NBlocs=RBlocs/2;
	double *mat=(double*)malloc(M*N*sizeof(double));
	double *A=(double*)malloc(M*N*sizeof(double));
	int *permut=(int*)malloc(N*sizeof(int));//permutation matrix required by lapack libraries
	for (int i=0;i<N;i++){
		permut[i]=i;
	}
	mat_init(mat,A,N,M);
//Starting PAPI counters
#ifdef USE_PAPI
	float t2 = gettime();
	if ((ret = PAPI_start_counters(events, NUM_EVENTS)) != PAPI_OK) {
		   fprintf(stderr, "PAPI failed to start counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}
#endif
//************************************


	block_LU(mat,M,N,NBlocs,permut);
//********************
//Reading PAPI values
//********************

#ifdef USE_PAPI
	float t2_1 = gettime();
	if (PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK){
		   fprintf(stderr, "PAPI failed to read counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}
	printf("This for the block OMP version \n");
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

	//Using PAPI on lapack 
	//mat_print(mat,N,M);
	dgetrf_(&M, &N, A, &M, permut, &info);
//**************
//PAPI
//**************

#ifdef USE_PAPI
	float t4 = gettime();
	if (PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK){
		   fprintf(stderr, "PAPI failed to read counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}

	printf("\n \n This for the LAPACK  version \n");
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





