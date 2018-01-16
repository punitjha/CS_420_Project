#include <iostream>
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
#include <omp.h>
#include <mpi.h>

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


//**************************************************************
void mat_create(double **&mat, int row, int col)
{
	mat=new double*[row];
	for (int i=0; i<row; i++)
	{
		mat[i]=new double[col];
	}
}

void mat_init(double **mat, int row, int col, double **A)
{
	for (int i= 0; i < row; ++i) 
	{
		for (int j= 0; j < col; ++j) 
		{
			    mat[i][j]=A[i][j]=(double)(rand()%10)+1;
		} 

	}
}
void print_matrix_chunk(int dim, int row_start, int col_start, double** matrix) {
	int i,j;
	for(i=row_start;i<(row_start+dim);i++) {
		for (j=col_start;j<(col_start+dim);j++) {
				printf("%7.3f ", matrix[i][j]);
			}
		printf("\n");
	}
}


void mat_initi_zero(double **mat, int row, int col)
{
	for (int i= 0; i < row; ++i) 
	{
		for (int j= 0; j < col; ++j) 
		{
			mat[i][j]=0.0;	
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

//******************************************************************
//LU decomposition with MPI using Cartesian topology
//******************************************************************


void cart_LU(int argc, char **argv, double **matrix, int dim, int nblocks)
{
	int rank;
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;
	MPI_Request request;
	MPI_Status sta[12];
	MPI_Request req[12];

	int cart_rows=sqrt(size);

	int cart_cols=sqrt(size);
	int cart_dims[2];
	int periods[2];
	cart_dims[0]=cart_rows, cart_dims[1]=cart_cols;
	periods[0]=0; 
	periods[1]=0; 

	MPI_Comm comm2D;
	MPI_Cart_create(MPI_COMM_WORLD,2,cart_dims,periods,0, &comm2D);
	//creating the neighbourhood

	int Coords[2];	
	int myrow,mycol;
	MPI_Cart_coords(comm2D,rank,2,Coords);
	myrow=Coords[0];
	mycol=Coords[1];

	
	//Now we determine the neighbours rank
	int right,bottom;
	int left=rank;
	int top=rank;
	MPI_Cart_shift(comm2D, 1,1,&left,&right); //this is right neighbour
	MPI_Cart_shift(comm2D, 0,1,&top,&bottom); //this is right neighbour
	//the size of all the buffers are equal to block_size
	
	//also L is only dim size i.e it can accomodate the whole col/row
	//since its the L matrix -- one of the two matrix needed
	//P stands for Pivot containg buffers
	//L stands for row containing buffers
	double **Lower;
	mat_create(Lower,dim,dim);
	mat_initi_zero(Lower,dim,dim);
	double *LowerBuffSend = (double*) malloc (nblocks* sizeof(double));
	double *LowerBuffRecv = (double*) malloc (nblocks* sizeof(double));
	double *PivotBuffSend = (double*) malloc (nblocks* sizeof(double));
	double *PivotBuffRecv = (double*) malloc (nblocks* sizeof(double));

	int i,j,k;

	//initialize the buffers 
	
	for (i=0;i<nblocks;i++) {
		LowerBuffSend[i] = LowerBuffRecv[i] = PivotBuffSend[i] = PivotBuffRecv[i] = 0;
	} 

	// initialize the Lower matrix  diag with one
	for (i=0;i<dim;i++) {
		Lower[i][i] = 1.0;
	}

	int parts=dim/nblocks;
	int StartCol=(rank*nblocks)%dim;
	int EndCol = StartCol +nblocks -1;
	int StartRow= (rank/parts)*nblocks;
	int EndRow= StartRow+nblocks-1;

	//print_matrix(dim,matrix);
	//This loop runs on all the processors simultaneousely and they pick up there
	//values of rows/cols depending on the partions above
	
	for (k=0; k<dim; k++){
		//***********************************************
		int Row_here= (k>=StartRow && k<=EndRow) ? 1 :0;
		int Row_up= (k <= EndRow-nblocks) ? 1 :0;
		int Row_down= (k>=StartRow +nblocks) ? 1 :0;

		int Col_here= (k>=StartCol && k<=EndCol) ? 1 :0;
		int Col_left= (k<=EndCol-nblocks) ? 1 :0;
		int Col_right= (k>=StartCol +nblocks) ? 1 :0;
                //**********************************************
	
		if (top>=0 && Row_up ==1 && Col_right==0){
			MPI_Recv(PivotBuffRecv,nblocks,MPI_DOUBLE,top,0,MPI_COMM_WORLD,&status);

		for(j=StartCol;j<=EndCol;j++){
			if(j>=k){
				matrix[k][j]=PivotBuffRecv[j-StartCol];
				}
			}
		}
                //**********************************************
		if (bottom>=0 && Col_right ==0){
			if (Row_here ==1 ){
				for (j=StartCol;j<=EndCol;j++){
					if (j>=k){
						PivotBuffSend[j-StartCol]=matrix[k][j];
					}
				}
			}
			else if (Row_up==1){
				for (j=StartCol;j<=EndCol;j++){
					if (j>=k){
						PivotBuffSend[j-StartCol]=PivotBuffRecv[j-StartCol];
					}
				}
			}
			MPI_Isend(PivotBuffSend,nblocks,MPI_DOUBLE,bottom,0,MPI_COMM_WORLD,&request);
		}

                //**********************************************
		//Calculating the lower triangular matrix parts
		if (Col_here ==1){
			for(i=StartRow;i<=EndRow;i++){
				if(i>k){
					Lower[i][k]=matrix[i][k]/matrix[k][k];
				}
			}
		}

                //**********************************************
		//Waiting for the send to complete
		if(bottom >= 0 && Row_here == 1)
			MPI_Wait(&request, &status);
		


                //**********************************************

		if(left >= 0 && Col_left ==1 && Row_down ==0) {
			MPI_Recv(LowerBuffRecv, nblocks, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &status);
			for(i=StartRow;i<=EndRow;i++) {
				if(i>k) {
					Lower[i][k] = LowerBuffRecv[i-StartRow];
				}
			}
		}

                //**********************************************
		if(right >= 0 && Row_down ==0) { 
			if(Col_here == 1) {  
				for(i=StartRow;i<=EndRow;i++) {
					if(i>k) {
						LowerBuffSend[i-StartRow] = Lower[i][k];
					}
				}
			}
			else if(Col_left == 1) { 
				for(i=StartRow;i<=EndRow;i++) {
					if(i>k) {
						LowerBuffSend[i-StartRow] = LowerBuffRecv[i-StartRow];
					}
				}
			}
			MPI_Isend(LowerBuffSend, nblocks, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &request);
		}

                //**********************************************
		//Compute upper triangular matrix
		#pragma omp parallel for private(j,i) firstprivate(k,StartCol,EndCol) 
		for (j=StartCol;j<=EndCol;j++) {
			if (j>=k) {
				for (i=StartRow;i<=EndRow;i++) {
					if (i>k) {
						matrix[i][j] = matrix[i][j]-Lower[i][k]*matrix[k][j];
					}
				}
			}
		}

                //**********************************************
		//Wait for LBuffSend to be usable
		if(right >= 0 && Col_here== 1)
			MPI_Wait(&request, &status);
			}

                
	//**********************************************
		double Lower_chunk[6][6],Upper_chunk[6][6];
		//mat_create(Lower_chunk,nblocks,nblocks);
		//mat_create(Upper_chunk,nblocks,nblocks);
		//mat_initi_zero(Lower_chunk,nblocks,nblocks);
		//mat_initi_zero(Upper_chunk,nblocks,nblocks);
		for (int l=0;l<6;l++){
			for(int h=0;h<6;h++){
				Lower_chunk[l][h]=Upper_chunk[l][h]=0.0;
			}
		}

		if (rank > 0 ){
		int row = 0;
		for(i=StartRow;i<=EndRow;i++) {
			int col = 0;
			for(j=StartCol;j<=EndCol;j++) {
				Lower_chunk[row][col] = Lower[i][j];
				Upper_chunk[row][col] = matrix[i][j];  
				col++;
			}
			row++;
		}
			MPI_Send(&Lower_chunk,nblocks*nblocks,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
			MPI_Send(&Upper_chunk,nblocks*nblocks,MPI_DOUBLE,0,10,MPI_COMM_WORLD);
			MPI_Send(&StartCol,1,MPI_INT,0,1,MPI_COMM_WORLD);
			MPI_Send(&StartRow,1,MPI_INT,0,2,MPI_COMM_WORLD);
			MPI_Send(&EndCol,1,MPI_INT,0,3,MPI_COMM_WORLD);
			MPI_Send(&EndRow,1,MPI_INT,0,4,MPI_COMM_WORLD);
		}
		else if(rank ==0){
			for (int k=1;k<size;k++)
			{
			MPI_Recv(&Lower_chunk,nblocks*nblocks,MPI_DOUBLE,k,0,MPI_COMM_WORLD,&status);
			MPI_Recv(&Upper_chunk,nblocks*nblocks,MPI_DOUBLE,k,10,MPI_COMM_WORLD,&status);
			MPI_Recv(&StartCol,1,MPI_INT,k,1,MPI_COMM_WORLD,&status);
			MPI_Recv(&StartRow,1,MPI_INT,k,2,MPI_COMM_WORLD,&status);
			MPI_Recv(&EndCol,1,MPI_INT,k,3,MPI_COMM_WORLD,&status);
			MPI_Recv(&EndRow,1,MPI_INT,k,4,MPI_COMM_WORLD,&status);
			int row=0;
//			printf("Star Row %d: ,End Row %d:, Start Col %d:, End Col %d: \n", StartRow, EndRow, StartCol, EndCol);
			for(i=StartRow;i<=EndRow;i++) {
				int col = 0;
				for(j=StartCol;j<=EndCol;j++) {
					Lower[i][j]=Lower_chunk[row][col];
					matrix[i][j]=Upper_chunk[row][col] ;  
					col++;
					}
				row++;
				}
			//MPI_Waitall(12,req,sta);
			}
		}

	


}



void serial_LU(double **matrix, int dim)
{
	int i,j,k;
	double **Lower = NULL;
	mat_create(Lower,dim,dim);
	mat_initi_zero(Lower,dim,dim);
	for (i=0; i<dim;i++){
		Lower[i][i]=1;
	}


	for (j=0;j<dim;j++){
		for (i=j+1;i<dim;i++){
			double ratio=matrix[i][j]/matrix[j][j];
			Lower[i][j]=ratio;
			for( k=0;k<dim;k++){
				matrix[i][k]-=ratio*matrix[j][k];
			}
		}
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
//

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int size;
	int rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int dim =atoi(argv[1]);
	//int num_threads=atoi(argv[2]);
	
	int nblocks=dim/(sqrt(size));
	
	int row=dim;
	int col=dim;
	double **mat=NULL;
	double **A=NULL;
	mat_create(mat,row,col);
	mat_create(A,row,col);
	mat_init(mat,row,col,A);
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


#ifdef USE_PAPI
	float t2 = gettime();
	if ((ret = PAPI_start_counters(events, NUM_EVENTS)) != PAPI_OK) {
		   fprintf(stderr, "PAPI failed to start counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}
#endif
	double start=MPI_Wtime();
	cart_LU(argc,argv,mat,dim,nblocks);
	double end=MPI_Wtime();
#ifdef USE_PAPI
	float t2_1 = gettime();
	if (PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK){
		   fprintf(stderr, "PAPI failed to read counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}
	printf("This for the Parallel version \n");
	printf("Time in Microsecond %f\n",(t2_1-t2));
	printf("Time in Microsecond in the MPI time%f\n",(start-end));
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
	if(rank ==0){
		serial_LU(A,dim);
#ifdef USE_PAPI
	float t4 = gettime();
	if (PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK){
		   fprintf(stderr, "PAPI failed to read counters: %s\n", PAPI_strerror(ret));
		      exit(1);
	}

	printf("\n \n This for the Serial  version \n");
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

	}
	MPI_Finalize();
 return 0;
}





