#!/usr/bin/env python
import numpy as np
import fft

def correlation(A,B):
  fft.fft(A,1.)
  fft.fft(B,1.)
  N = len(A) >> 1
  corr = np.zeros(2*N)
  for i in range(2*N):
    B[i] = -B[i]
  for j in range(N):
    corr[2*j] = A[2*j]*B[2*j] - A[2*j+1]*B[2*j+1]
    corr[2*j+1] = A[2*j+1]*B[2*j] + A[2*j]*B[2*j+1]
  fft.fft(corr,-1.)
  return corr

N1 = 32
N2 = 128
L = 1
dx1 = L/N1
dx2 = L/N2
x1 = np.zeros(N1)
x2 = np.zeros(N2)
g1 = np.zeros(2*N1)
h1 = np.zeros(2*N1)
g2 = np.zeros(2*N2)
h2 = np.zeros(2*N2)
k = 0
l = 0
m = 0
n = 0

for i in range(N1):
  x1[i] = i * dx1
  if x1[i] <= 0.1:
    g1[2*i] = 1
    k = k + 1
  elif x1[i] >= 0.4 and x1[i] <= 0.6:
    h1[2*i] = 1
    l = l + 1
for i in range(N1):
  g1[i] = g1[i] / k
  h1[i] = h1[i] / l

for j in range(N2):
  x2[j] = j * dx2
  if x2[j] <= 0.1:
    g2[2*j] = 1
    m = m+1
  elif x2[j] >= 0.4 and x2[j] <= 0.6:
    h2[2*j] = 1
    n = n+1
for j in range(N2):
  g2[j] = g2[j] / m
  h2[j] = h2[j] / n

corr1 = correlation(g1,h1)
corr2 = correlation(g2,h2)

fft.plot_c(x1,g1)
fft.plot_c(x1,h1)
fft.plot_c(x1,corr1)
fft.plot_c(x2,g2)
fft.plot_c(x2,h2)
fft.plot_c(x2,corr2)
