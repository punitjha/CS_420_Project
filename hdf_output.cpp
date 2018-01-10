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
  hsize_t dims[2];
  herr_t status;
  double ans[nx][ny];
  const char* time_ptr;
  const char* dset_ptr;

  stringstream time_stream;
  time_stream << fixed << setprecision(5) << time;
  string time_string = time_stream.str();
  string title = "Time: ";
  string unit = " [sec]";
  string group_name;
  string dataset_name = "u";

  group_name.append(title);
  group_name.append(time_string);
  group_name.append(unit);

  time_ptr = group_name.c_str();
  dset_ptr = dataset_name.c_str();

  dims[0] = nx;
  dims[1] = ny;

  fid = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT);
  gid = H5Gcreate(fid, time_ptr, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  dspaid = H5Screate_simple(2, dims, NULL);  
  dsetid = H5Dcreate(gid, dset_ptr, H5T_NATIVE_DOUBLE, dspaid, \
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      ans[i][j] = gsl_matrix_get(sol,i,j);
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
  hid_t fid, gid, dsetid_x, dsetid_y, dspaceid_x, dspaceid_y;
  hsize_t x_dims[1], y_dims[1];
  herr_t status;
  double x[nx];
  double y[ny];

  x_dims[0] = nx;
  y_dims[0] = ny;

  fid = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  gid = H5Gcreate(fid, "Coordinates", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  dspaceid_x = H5Screate_simple(1, x_dims, NULL);
  dspaceid_y = H5Screate_simple(1, y_dims, NULL);
  dsetid_x = H5Dcreate(fid, "Coordinates/X", H5T_NATIVE_DOUBLE, dspaceid_x, \
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dsetid_y = H5Dcreate(fid, "Coordinates/Y", H5T_NATIVE_DOUBLE, dspaceid_y, \
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  for (int i = 0; i < nx; i++) {
    x[i] = x_lin[i];
  }
  for (int i = 0; i < ny; i++) {
    y[i] = y_lin[i];
  }

  status = H5Dwrite(dsetid_x, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, \
		    H5P_DEFAULT, x);
  status = H5Dwrite(dsetid_y, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, \
		    H5P_DEFAULT, y);
  status = H5Dclose(dsetid_x);
  status = H5Sclose(dspaceid_x);
  status = H5Dclose(dsetid_y);
  status = H5Sclose(dspaceid_y);
  status = H5Gclose(gid);
  status = H5Fclose(fid);

}
