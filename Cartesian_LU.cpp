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
//#include <armadillo>


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
			if (i==j){
			    mat[i][j] =1;
			    A[i][j]=1;
			} 
			else{
			    mat[i][j] =4;
			    A[i][j]=4;
			    mat[j][i] =6;
			    A[j][i]=6;


			}
		}
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

void Cont_create(double**&mat, int size)
{
	mat= new double*[size];
	mat[0] = new double[size*size];
	for (int i = 1; i < size; ++i) {
		mat[i] = mat[0] + i*size;
	}

}

//******************************************************************
//LU decomposition with MPI using Cartesian topology
//******************************************************************


void cart_LU(int argc, char **argv, double **matrix, int dim, int nblocks, int num_threads)
{
	omp_set_num_threads(num_threads);
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
		double *Lower_chunk=(double*)malloc (sizeof(double)*nblocks *nblocks);
		double *Upper_chunk=(double*)malloc(sizeof(double)*nblocks *nblocks);
		//Cont_create(Lower_chunk,nblocks);
		//Cont_create(Upper_chunk,nblocks);
		for (int l=0;l<nblocks;l++){
			for(int h=0;h<nblocks;h++){
				Lower_chunk[l*nblocks+h]=Upper_chunk[l*nblocks+h]=0.0;
			}
		}

		if (rank > 0 ){
		int row = 0;
		for(i=StartRow;i<=EndRow;i++) {
			int col = 0;
			for(j=StartCol;j<=EndCol;j++) {
				Lower_chunk[row*nblocks+col] = Lower[i][j];
				Upper_chunk[row*nblocks+col] = matrix[i][j];  
				col++;
			}
			row++;
		}
			MPI_Send(Lower_chunk,nblocks*nblocks,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
			MPI_Send(Upper_chunk,nblocks*nblocks,MPI_DOUBLE,0,10,MPI_COMM_WORLD);
			MPI_Send(&StartCol,1,MPI_INT,0,1,MPI_COMM_WORLD);
			MPI_Send(&StartRow,1,MPI_INT,0,2,MPI_COMM_WORLD);
			MPI_Send(&EndCol,1,MPI_INT,0,3,MPI_COMM_WORLD);
			MPI_Send(&EndRow,1,MPI_INT,0,4,MPI_COMM_WORLD);
		}
		else if(rank ==0){
			for (int k=1;k<size;k++)
			{
			MPI_Recv(Lower_chunk,nblocks*nblocks,MPI_DOUBLE,k,0,MPI_COMM_WORLD,&status);
			MPI_Recv(Upper_chunk,nblocks*nblocks,MPI_DOUBLE,k,10,MPI_COMM_WORLD,&status);
			MPI_Recv(&StartCol,1,MPI_INT,k,1,MPI_COMM_WORLD,&status);
			MPI_Recv(&StartRow,1,MPI_INT,k,2,MPI_COMM_WORLD,&status);
			MPI_Recv(&EndCol,1,MPI_INT,k,3,MPI_COMM_WORLD,&status);
			MPI_Recv(&EndRow,1,MPI_INT,k,4,MPI_COMM_WORLD,&status);
			int row=0;
//			printf("Star Row %d: ,End Row %d:, Start Col %d:, End Col %d: \n", StartRow, EndRow, StartCol, EndCol);
			for(i=StartRow;i<=EndRow;i++) {
				int col = 0;
				for(j=StartCol;j<=EndCol;j++) {
					Lower[i][j]=Lower_chunk[row*nblocks+col];
					matrix[i][j]=Upper_chunk[row*nblocks+col] ;  
					col++;
					}
				row++;
				}
			//MPI_Waitall(12,req,sta);
			}
		}

	if (rank ==0){
	printf("L from the Cartesian version \n \n ");
	mat_print(Lower,dim,dim);
	printf("U form the Cartesian version\n \n ");
	mat_print(matrix,dim,dim);
	}


}

