#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import improc
import time
import fft

N = [100,1000,10000]
N1 = np.random.rand(100)
N2 = np.random.rand(1000)
N3 = np.random.rand(10000)

t0 = time.time()
N1_sel = improc.selectionsort(N1,100)
t1 = time.time()
N1_mer = improc.mergesort(N1,100)
t2 = time.time()
N2_sel = improc.selectionsort(N2,1000)
t3 = time.time()
N2_mer = improc.mergesort(N2,1000)
t4 = time.time()
N3_sel = improc.selectionsort(N3,10000)
t5 = time.time()
N3_mer = improc.mergesort(N3,10000)
t6 = time.time()

t_sel = [t1-t0,t3-t2,t5-t4]
t_mer = [t2-t1,t4-t3,t6-t5]
plt.figure()
plt.plot(N,t_sel,N,t_mer)
plt.legend(["selection sort","merge sort"])
plt.show()

tr100 = t_sel[0]/t_mer[0]
tr1000 = t_sel[1]/t_mer[1]
tr10000 = t_sel[2]/t_mer[2]
print(tr100,tr1000,tr10000)

#The data does appear to match the performance expectations for the two difference 
#sorting methods.  Also, looking at the ratios we can see that the selection sort
#is actually faster for the size 100 array however it gets progressibely worse as
#the size of the array increases.

def Fibonacci(n):
  if n == 0:
    return 0
  elif n == 1:
    return 1
  else:
    return Fibonacci(n-1) + Fibonacci(n-2)

print(Fibonacci(24))

b = [4, 13, 28, 39, 50]
b = np.array(b)
N = 64
PI = 2*np.arcsin(1)
dx = 2*PI / (N-1)

x = np.zeros(N)
f = np.zeros(N)
y_1 = np.zeros(N)
y_2 = np.zeros(N)
y_3 = np.zeros(N)
y_4 = np.zeros(N)
y_5 = np.zeros(N)

for i in range(N):
  x[i] = i*dx
  f[i] = i
  y_1[i] = np.sin(x[i]) + np.sin(b[0]*x[i])
  y_2[i] = np.sin(x[i]) + np.sin(b[1]*x[i])
  y_3[i] = np.sin(x[i]) + np.sin(b[2]*x[i])
  y_4[i] = np.sin(x[i]) + np.sin(b[3]*x[i])
  y_5[i] = np.sin(x[i]) + np.sin(b[4]*x[i])

y_1 = fft.fft_slow(y_1,1.)
fft.plot_c(f[:0.5*N],y_1[:N])
y_2 = fft.fft_slow(y_2,1.)
fft.plot_c(f[:0.5*N],y_2[:N])
y_3 = fft.fft_slow(y_3,1.)
fft.plot_c(f[:0.5*N],y_3[:N])
y_4 = fft.fft_slow(y_4,1.)
fft.plot_c(f[:0.5*N],y_4[:N])
y_5 = fft.fft_slow(y_5,1.)
fft.plot_c(f[:0.5*N],y_5[:N])

#At the different spectra, the peaks moved to the higher frequencies until
#they reached b=39.  
