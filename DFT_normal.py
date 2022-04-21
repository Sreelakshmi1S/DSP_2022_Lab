import numpy as np 
import matplotlib.pyplot as plt  
N=len(x) 
D=np.empty((N,N),dtype=complex) 
W=np.exp(-1j*2*np.pi/N) 
for k in np.arange(N): 
     for n in np.arange(N): 
         D[k,n]=W**(k*n) 
np.round(D) 
x=np.array([[1.+0.j, 3.+0.j, 5.+0.j, 7.+0.j, 9.+0.j]]).T 
X=D@x 
np.round(X) 
print(X) 
