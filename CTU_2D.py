import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import os
import sys
from pyvisfile.vtk import write_structured_grid

nx = 4  # number of cells
ny = nx  # assume equal (however the rest of the code is written generally,
         # but not guaranteed to run correctly. not tested.

# domain on -1 to 1
dx = 2.0 / nx
dy = dx

# needs for plotting 
# set x i-1/2
xhalf, yhalf = np.meshgrid(np.linspace(-1, 1-dx, nx), np.linspace(-1, 1-dy, ny))
# set x i
x = xhalf + dx/2
y = yhalf + dy/2
mesh = np.rollaxis(np.dstack((x, y)), 2)


def u0c(x, y):
    # Gaussian initial condition
    return np.exp(-10.0 * (x**2.0 + y**2.0))

# remove all old vts output files
os.system("rm *.vts")

T = 1.0  # run simulations for T periods

u = u0c(x, y)  # set initial condition
cx = 1.0  # flow speed in x
cy = 1.0  # flow speed in y

write_structured_grid('test0.vts', mesh,
                      point_data=[('u', u[np.newaxis, :, :])])

nt = 32  # number of timesteps
dt = T/nt  # delta t
lmbda_x = cx*dt/dx  # courant number in x direction
lmbda_y = cy*dt/dy  # courant number in y direction

v = 0  # vis file index
svl = nx*ny  # solution vector length
x_boundary = np.arange(nx, svl, nx)-1  # each row boundary
y_boundary = np.arange(ny, svl, ny)-1  # (-1 to start index from 0)

hpb = np.arange(nx, svl+1, nx)-1  # horizontal periodic boundary
vpb = np.arange(ny, svl+1, ny)-1  # vertical periodic boundary

# loop over timesteps
for time in range(0, nt):
    # CREATE THE MATRIX m1
    lmbda_array_1 = np.ones(nx*ny-1)*[lmbda_x/2.0]
    lmbda_array_2 = np.ones(nx*ny-nx)*[lmbda_y/2.0]

    # five diagonal matrix
    m1 = np.diag(-1*lmbda_array_1, -1) + np.diag(lmbda_array_1, 1) + \
         np.diag(-1*lmbda_array_2, -nx) + np.diag(lmbda_array_2, nx) + \
         np.identity(nx*ny)

    # remove diag lambdas where each row ends and starts
    for i in range(len(x_boundary)):
        m1[x_boundary[i], y_boundary[i]+1] = 0
        m1[x_boundary[i]+1, y_boundary[i]] = 0
            
    # enforce periodic boundary conditions
    for i in range(len(hpb)):
        # assign horizontal periodic boundary conditions
        m1[hpb[i]-(nx-1), hpb[i]] = -1*lmbda_x/2.0
        m1[vpb[i], vpb[i]-(ny-1)] = 1*lmbda_x/2.0
    # assign vertical periodic boundary conditons
    lmbda_v = np.ones(ny)*[lmbda_y/2.0]

    m1 = m1 + np.diag(-1*lmbda_v,(nx*(ny-1))) + \
         np.diag(1*lmbda_v,-(nx*(ny-1)))

    # create m2:
    m2 = m1.transpose()
    # solve the matrix:
    # A*x=b          A-matrix       b-vector
    u = np.linalg.solve(m1, np.dot(m2, u.flatten()))
    # vector solution back into matrix
    u = u.reshape([nx,ny])
    # output the result which can be viewed by paraview
    write_structured_grid('test%d.vts' % (v+1), mesh,
                          point_data=[('u', u[np.newaxis, :, :])])
    v += 1
