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
