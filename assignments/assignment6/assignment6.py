#!/usr/bin/env python
import numpy as np
import integrator as intr
import rng
import matplotlib.pyplot as plt

def hist(A,N):
  L = len(A)
  d = (max(A)-min(A))/N
  H = np.zeros(2*N)
  X = np.zeros(2*N)
  for i in range(N):
    mi = min(A) + i*d
    mx = min(A) + (i+1) * d
    X[2*i] = i*d/N
    X[2*i+1] = (i+1)*d/N
    if i == 0:
      for k in range(L):
        if A[k] >= mi and A[k] <= mx:
          H[2*i] = H[2*i] + 1
          H[2*i+1] = H[2*i+1] + 1
    else:
      for j in range(L):
        if A[j] > mi and A[j] <= mx:
          H[2*i] = H[2*i] + 1
          H[2*i+1] = H[2*i+1] + 1
    H[2*i] = H[2*i] / L
    H[2*i+1] = H[2*i+1] / L
  plt.plot(X,H)
  plt.ylim(0,1)
  plt.show()

A = np.zeros(100000)
for m in range(100000):
  A[m] = intr.exponential(2)

hist(A,100)

#The histogram looks very similar to the expected probability density function.

[w,Sw] = intr.monte_carlo_3d(intr.density,[-0.5,-0.5,-1],[1,1,1],intr.sphere,100000)
[x,Sx] = intr.monte_carlo_3d(intr.xmoment,[-0.5,-0.5,-1],[1,1,1],intr.sphere,100000)
[y,Sy] = intr.monte_carlo_3d(intr.ymoment,[-0.5,-0.5,-1],[1,1,1],intr.sphere,100000)
[z,Sz] = intr.monte_carlo_3d(intr.zmoment,[-0.5,-0.5,-1],[1,1,1],intr.sphere,100000)

print(x/w,y/w,z/w)

