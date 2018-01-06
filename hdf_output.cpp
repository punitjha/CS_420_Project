#include <hdf5.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#define FILE "out.h5"

void output_h5(int nx, int ny, gsl_matrix* sol) {
  hid_t fid, gid;
  herr_t status;

  fid = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT);
  //  gid = H5Gcreate(fid, "Coordinates", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  //  status = H5Gclose(gid);
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
