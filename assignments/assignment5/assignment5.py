#!/usr/bin/env python
import numpy as np
import integrator as intr
import matplotlib.pyplot as plt

PI = 2. * np.arcsin(1)
N = 1000
a = 0
b = 2. * PI
h = (b - a) / (N - 1)
x = np.zeros(N)
y1 = np.zeros(N)
y2 = np.zeros(N)
y3 = np.zeros(N)
for i in range(N):
  x[i] = a + h * i
  y1[i] = np.sin(x[i])
  y2[i] = np.sin(2. * x[i])
  y3[i] = np.sin(8. * x[i])

print(intr.trap_d(x,intr.multiply(y1,y1)))
print(intr.trap_d(x,intr.multiply(y1,y2)))
print(intr.trap_d(x,intr.multiply(y2,y3)))
print(intr.trap_d(x,intr.multiply(y3,y3)))

#We find that when m=n the scalar product is pi and when m=/n the scalar product
#is 0. This is very similar to the scalar product of the legendre functions.

def powr(x,n):
  m = 1
  if n == 0:
    return 1
  else:
    for l in range(n):
      m = m * x
    return m

def x11(x):
  return x*x*x*x*x*x*x*x*x*x*x

def x12(x):
  return x*x*x*x*x*x*x*x*x*x*x*x

N1 = 10
N2 = 100
N3 = 1000
n1 = 1
n2 = 2
n3 = 6
H1 = 1 / (N1 -1)
H2 = 1 / (N2 -1)
H3 = 1 / (N3 -1)

# rather than having a whole lot of functions I just swatched the N values for each
[P1,W1] = intr.gauss_leg(0,1,n3)
X1 = np.zeros(N3)
Y1 = np.zeros(N3)
Y2 = np.zeros(N3)
Y3 = np.zeros(N3)
Y4 = np.zeros(N3)
Y5 = np.zeros(N3)
for i in range(N3):
  X1[i] = H3*i
  Y1[i] = powr(X1[i],1)
  Y2[i] = powr(X1[i],2)
  Y3[i] = powr(X1[i],11)
  Y4[i] = powr(X1[i],12)
  Y5[i] = np.exp(X1[i])

print(intr.trap_d(X1,Y1))
print(intr.trap_d(X1,Y2))
print(intr.trap_d(X1,Y3))
print(intr.trap_d(X1,Y4))
print(intr.trap_d(X1,Y5))
print(intr.gauss_quad(intr.lin,P1,W1))
print(intr.gauss_quad(intr.quad,P1,W1))
print(intr.gauss_quad(x11,P1,W1))
print(intr.gauss_quad(x12,P1,W1))
print(intr.gauss_quad(np.exp,P1,W1))
