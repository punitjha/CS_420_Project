#include <hdf5.h>
#define FILE "out.h5"

void output_h5() {
  hid_t fid;
  herr_t status;

  fid = H5Fcreate(FILE, H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);

  status = H5Fclose(fid);
}
