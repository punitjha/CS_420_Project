//This is the cpp conversion of the Python code written initially of Heeho
// To compile g++-5 -g -std=c++0x initial.cpp -lblas -llapack -larmadillo
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
#include <armadillo>
#include <complex>

#include <array>   // to be able to use this library compile using this g++ -std=c++0x FileNames


void init_cond(arma::mat &x,arma::mat &y, arma::mat &u)
{
	//this function return the initial position
	 u=exp(-10.0*(pow(x,2.0) + pow(y,2.0)));
}

 void meshgrid(arma::vec &x_lin, arma::vec &y_lin,arma::mat &xhalf, arma::mat &yhalf)
{
	//this function creates the mesh
	xhalf.each_row() +=x_lin.t();
	yhalf.each_col() +=y_lin;
}

void arange(arma::vec &x,int ini, int fin, int steps)
{
	//figure out the number of elments
	int counter=0;
	for (int i=ini; i<fin; i+=steps)
	{
		counter++;
	}
	int k=0;
	x=arma::zeros(counter);
	//initialize the vector and assign the values
	for(int j=ini; j<fin; j+=steps)
	{
		x(k)=j;
		k++;
	}	

}

int main()
{
	int nx=4; //This is the number of cells 
	int ny=nx; //No. of grids in both the dimensions are set equal

	//Setting the domain of -1 to 1
	double dx=2.0/nx;
	double dy=dx;

	//declaring mesh grid points
	
	arma::vec x_linspace=arma::linspace(-1,1-dx,nx);
	arma::vec y_linspace=arma::linspace(-1,1-dx,nx);
	//x_linspace.print();

	//creating the mesh points	
	arma::mat xhalf= arma::zeros(x_linspace.n_rows,x_linspace.n_rows);
	arma::mat yhalf= arma::zeros(y_linspace.n_rows,y_linspace.n_rows);
	meshgrid(x_linspace,y_linspace,xhalf,yhalf);
	//xhalf.print();
	//yhalf.print();

	arma::mat x= xhalf + dx/2.0;	
	arma::mat y= yhalf + dy/2.0;
	//x.print("This is x");
	//y.print("This is y");

	double T=1.0; // run simulations for T periods
	
	arma::mat u= arma::zeros(x.n_rows, y.n_rows);
	init_cond(x,y,u);//building the initial condition matrix
	//u.print("This is the initial condition matrix");
	
	double cx=1.0; // flow speed in x
	double cy=1.0; // flow speed in y
	int nt=32; // no. of time steps 
	double dt= T/nt; //the delta t -- time steps
	double lmbda_x = cx*dt/dx; // courant number in x direction
	double lmbda_y = cy*dt/dy; // courant number in y direction
	int svl = nx*ny;// the sol vector length

	arma::vec x_boundary;	
	arma::vec y_boundary;	
	arange(x_boundary,nx,svl,nx);//equivalent to arange python
	arange(y_boundary,ny,svl,ny);
	x_boundary=x_boundary-1; //each row boundary
	y_boundary=y_boundary-1; //(-1 to start from index 0)
	//x_boundary.print("This is x");	
	//y_boundary.print("This is y");
	arma::vec hpb;
	arma::vec vpb;
	arange(hpb,nx,svl+1,nx);
	arange(vpb,ny,svl+1,ny);
	hpb-=1; //the horizontal periodic boundary
	vpb-=1; //the vertical periodic boundary
	//hpb.print("This is hbp");
	//vpb.print("This is vbp");
	
	//the itirations over time steps start
	for(int time=0; time <nt; time ++)
	{
		//creating the matrix m1
		arma::vec l_a1= arma::ones(nx*ny-1)*(lmbda_x/2.0);
		//l_a1.print();
		arma::vec l_a2= arma::ones(nx*ny-nx)*(lmbda_y/2.0);
		
		//five diagonal matrix
		arma::mat m1= arma::zeros(l_a1.n_rows+1,l_a1.n_rows+1);
	        m1.diag(-1)=-1*l_a1;
		m1.diag(1)=l_a1;
		m1.diag(-nx)=-1*l_a2;
		m1.diag(nx)=l_a2;
		m1.diag()+=1;


		//removing the diag lambdas where each row ends and starts
		for(int i=0; i<x_boundary.n_rows; i++)
		{
			m1(x_boundary(i),y_boundary(i)+1)=0;
			m1(x_boundary(i)+1,y_boundary(i))=0;
		}
		
		//enforce peridic boundary conditions
		
		for(int i=0; i<hpb.n_rows; i++)
		{
			m1(hpb(i)-(nx-1),hpb(i))=-1*lmbda_x/2.0;
			m1(vpb(i),vpb(i)-(ny-1))=1*lmbda_x/2.0;
		}


		//assign the vertical periodic boudary conditions
		
		arma::vec lmbda_v=arma::ones(ny)*(lmbda_y/2.0);

		m1.diag(nx*(ny-1))=-1*lmbda_v;
		m1.diag(-nx*(ny-1))=1*lmbda_v;

		//creating m2
		//m1.print("This is m1 matrix");	
		arma::mat  m2=m1.t();

		//solving the matrix
		//A*x=b		A-matrix 	b-vector
		arma::vec sol=arma::zeros(nx*ny);
		arma::vec u_vec=arma::vectorise(u,0);
		sol=arma::solve(m1,(m2*u_vec));
		u=arma::reshape(sol,nx,ny);
		u.print("This is the solution matrix \n \n");


	}







	return 0;
}
