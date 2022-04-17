"""1) Implement convolution of two sequences of length K and M (K ≥ M) using python. (Do not use
 available library function for convolution.). Use the program to compute [1 2 3] ∗ [1 1 1 1]."""
   
  
import numpy as np
import matplotlib.pyplot as plt
def dirconv(x,h):
   k=len(x)
   m=len(h)
   print(m)
   print(k)
   l=k+m-1
   y=np.zeros(l)
   for i in range(0,k-1):
       s=0
       for j in range(0,m):
           s=s+np.dot(x[i+j],h[m-j-1])
       y[i+1]=s
    y[l-1]=np.dot(h[1],x[k-1])
    y[0]=np.dot(x[0],h[0])
    return y
x=np.random.randint(1,4,4)
h=np.random.randint(1,4,2)
y=dirconv(x,h)
print("Convolution sum=",y)
n=np.arange(len(y))
plt.title("linear convolution")
plt.xlabel('n')
plt.ylabel('y[n]')
plt.stem(n,y)
plt.xticks(n)
plt.show()
