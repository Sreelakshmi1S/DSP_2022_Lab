import numpy as np
import matplotlib.pyplot as plt 
X=np.array([[25.+0.j,-5.+6.8819096j,-5.+1.62459848j,-5.-1.62459848j,-5.-6.8819096j ]]).T
N=len(X)
D=np.empty((N,N),dtype=np.cdouble)
W=np.exp(1j*2*np.pi/N)
k=np.arange(N)
for n in np.arange(N):
    for n in np.arange(N):
        D[:,n]=W**(k*n)
       
np.round(D)

x=(D/N)@X
np.round(x)
print(x)
