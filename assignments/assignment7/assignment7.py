#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pde

L = 1
U = 1
mu = 1
T = 0.1
Ny = 50
F = 0
dy = L/(Ny-1)
dt = 0.0001
Nt = round(T / dt + 1)

t = np.zeros(Nt)
y = np.zeros(Ny)
u = np.zeros([Nt,Ny])
u_sol = np.zeros(Ny)

for i in range(Ny):
  y[i] = -L/2 + i*dy
  u[0,i] = 0
  u_sol[i] = 2*U*y[i]/L

i = 1
while t[i-1]+dt <= T:
  t[i] = t[i-1] + dt
  u[i,:] = pde.euler_diffusion(u[i-1,:],mu,dt,dy,F)
  u[i,0] = -U
  u[i,Ny-1] = U
  i = i+1

plt.plot(y,u_sol,'o',y,u[0,:],y,u[np.floor(Nt/8),:],y,u[np.floor(Nt/4),:],y,u[np.floor(Nt/2),:],y,u[Nt-2,:])
plt.legend(["Exact solution","0","T/8","T/4","T/2","T"],loc = "upper left")
plt.show()

