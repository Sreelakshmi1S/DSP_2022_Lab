import numpy as np
import matplotlib.pyplot as plt

def dft(x):
   x1=np.array([x]).T
   N=len(x)
   D=np.empty((N,N),dtype=complex)
   W=np.exp(-2j*np.pi/N)
   for k in np.arange(N):
       for n in np.arange(N):
           D[k,n]=W**(k*n)
   X = D @ x1
   return X

def idft(X):
    X1=np.array([X]).T
    N=len(X)
    D=np.empty((N,N),dtype=complex)
    W=np.exp(2j*np.pi/N)
    for k in np.arange(N):
       for n in np.arange(N):
           D[k,n]=W**(k*n)
    a = (1/N)*(D @ X1)
    return a

sr=100
ts=1/sr
t=np.arange(0,1,ts)
x = 2*np.sin(20*np.pi*t) + np.cos(40*np.pi*t)

plt.show()
X=dft(x)
h=[1,-1,1]
y=np.convolve(x,h)
Y=dft(y)
N=len(Y)
n=np.arange(N)
T=N/sr
f=n/T
plt.stem(f,abs(Y),markerfmt=" ",basefmt="-b")
plt.show()
