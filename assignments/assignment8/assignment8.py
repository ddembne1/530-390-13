#!/usr/bin/env python
import numpy as np
import pde
import matplotlib.pyplot as plt

Lx = 10
Ly = 4
dx = 0.1
dy = 0.1
Nx = int(1 + (Lx/dx))
Ny = int(1 + (Ly/dy))
nx = Nx-2
ny = Ny-2
d2dx = 1/(dx*dx)

u = np.zeros(Nx*Ny)
f = np.zeros(Nx*Ny)

i = 0
for j in range(Ny):
  u[i+j*Nx] = 1

u = pde.jacobi(u,f,dx,Nx,Ny,1000)

u_int = np.zeros(nx*ny)
for j in range(ny):
  for i in range(nx):
    u_int[i+j*nx] = u[(i+1) + (j+1)*Nx]

u_plot = np.zeros([Nx,Ny])
for j in range(Ny):
  for i in range(Nx):
    u_plot[i,j] = u[i+j*Nx]

x = np.zeros([Nx+1,Ny+1])
y = np.zeros([Nx+1,Ny+1])
for j in range(Ny+1):
  for i in range(Nx+1):
    x[i,j] = -0.5*dx + i*dx
    y[i,j] = -0.5*dy + j*dy

plt.pcolormesh(x,y,u_plot)
plt.colorbar()
plt.show()
