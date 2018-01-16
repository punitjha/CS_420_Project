#include <hdf5.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <iomanip>
#include <sstream>
#include <iostream>
using namespace std;
#define FILE "out.h5"

void output_h5(int nx, int ny, gsl_matrix* sol, double time) {
  hid_t fid, gid, dsetid, dspaid;
  hsize_t dims[3];
  herr_t status;
  double ans[nx][ny][1];
  const char* time_ptr;
  const char* dset_ptr;
  //printf("hello-1\n");
  stringstream time_stream;
  time_stream << scientific << setprecision(5) << time;
  //printf("hello0\n");
  string time_string = time_stream.str();
  //printf("hello1\n");
  string title = "Time:  ";
  //printf("hello2\n");
  string time_unit = " s";
  //printf("hello3\n");
  string group_name = "";
  //printf("hello4\n");
  string dataset_name = "solution";
  //printf("hello5\n");
  group_name.append(title);
  group_name.append(time_string);
  group_name.append(time_unit);

  time_ptr = group_name.c_str();
  dset_ptr = dataset_name.c_str();
    
  dims[0] = nx;
  dims[1] = ny;
  dims[2] = 1;

  fid = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT);
  gid = H5Gcreate(fid, time_ptr, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  dspaid = H5Screate_simple(3, dims, NULL);  
  dsetid = H5Dcreate(gid, dset_ptr, H5T_NATIVE_DOUBLE, dspaid, \
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      ans[i][j][0] = gsl_matrix_get(sol,i,j);
    }
  }

  status = H5Dwrite(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, \
		    H5P_DEFAULT, ans);
  status = H5Dclose(dsetid);
  status = H5Sclose(dspaid);
  status = H5Gclose(gid);
  status = H5Fclose(fid);
}

void output_mat_h5(int nx, int ny, double* sol, double time) {
  hid_t fid, gid, dsetid, dspaid;
  hsize_t dims[3];
  herr_t status;
  double ans[nx][ny][1];
  const char* time_ptr;
  const char* dset_ptr;

  stringstream time_stream;
  time_stream << scientific << setprecision(5) << time;
  string time_string = time_stream.str();
  string title = "Time:  ";
  string unit = " s";
  string group_name;
  string dataset_name = "u";

  group_name.append(title);
  group_name.append(time_string);
  group_name.append(unit);

  time_ptr = group_name.c_str();
  dset_ptr = dataset_name.c_str();

  dims[0] = nx;
  dims[1] = ny;
  dims[2] = 1;

  fid = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT);
  gid = H5Gcreate(fid, time_ptr, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  dspaid = H5Screate_simple(3, dims, NULL);  
  dsetid = H5Dcreate(gid, dset_ptr, H5T_NATIVE_DOUBLE, dspaid, \
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      ans[i][j][0] = sol[i*ny+j];
    }
  }
  status = H5Dwrite(dsetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, \
		    H5P_DEFAULT, ans);
  status = H5Dclose(dsetid);
  status = H5Sclose(dspaid);
  status = H5Gclose(gid);
  status = H5Fclose(fid);
}


void init_h5(double *& x_lin, double *& y_lin, int nx, int ny) {
  hid_t fid, gid, dsetid_x, dsetid_y, dsetid_z, dspaceid_x, dspaceid_y, dspaceid_z;
  hsize_t x_dims[1], y_dims[1], z_dims[1];
  herr_t status;
  nx = nx + 1;
  ny = ny + 1;
  double x[nx];
  double y[ny];
  double z[2];

  x_dims[0] = nx;
  y_dims[0] = ny;
  z_dims[0] = 2;

  fid = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  gid = H5Gcreate(fid, "Coordinates", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  dspaceid_x = H5Screate_simple(1, x_dims, NULL);
  dspaceid_y = H5Screate_simple(1, y_dims, NULL);
  dspaceid_z = H5Screate_simple(1, z_dims, NULL);
  dsetid_x = H5Dcreate(gid, "X [m]", H5T_NATIVE_DOUBLE, dspaceid_x, \
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dsetid_y = H5Dcreate(gid, "Y [m]", H5T_NATIVE_DOUBLE, dspaceid_y, \
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dsetid_z = H5Dcreate(gid, "Z [m]", H5T_NATIVE_DOUBLE, dspaceid_z, \
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  for (int i = 0; i < nx; i++) {
    x[i] = x_lin[i];
  }
  for (int i = 0; i < ny; i++) {
    y[i] = y_lin[i];
  }
  for (int i = 0; i < 2; i++) {
    z[i] = double(i);
  }

  status = H5Dwrite(dsetid_x, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, \
		    H5P_DEFAULT, x);
  status = H5Dwrite(dsetid_y, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, \
		    H5P_DEFAULT, y);
  status = H5Dwrite(dsetid_z, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, \
		    H5P_DEFAULT, z);
  status = H5Dclose(dsetid_x);
  status = H5Sclose(dspaceid_x);
  status = H5Dclose(dsetid_y);
  status = H5Sclose(dspaceid_y);
  status = H5Dclose(dsetid_z);
  status = H5Sclose(dspaceid_z);
  status = H5Gclose(gid);
  status = H5Fclose(fid);

}
