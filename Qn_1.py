"""1) Implement convolution of two sequences of length K and M (K ≥ M) using python. (Do not use
 available library function for convolution.). Use the program to compute [1 2 3] ∗ [1 1 1 1]."""
   
import numpy as np
import matplotlib.pyplot as plt
#impulse response
h = [1,1,1,1];
#input response
x = [1,2,3];
def dirconv(x,h):
    N1 = len(x)
    N2 = len(h)
    N = N1+N2-1
    y = np.zeros(N)
    m = N-N1
    n = N-N2
#Padding zeros to x and h to make their length to N
    x =np.pad(x,(0,m),'constant')
    h =np.pad(h,(0,n),'constant')

#Linear convolution using convolution sum formula
    for n in range (N):
        for k in range (N):
            if n >= k:
               y[n] = y[n]+x[n-k]*h[k]
    return y
y=dirconv(x,h)
print('Linear convolution using convolution sum formula output response y =\n',y)
plt.title("linear convolution")
n=np.arange(len(y))
plt.xlabel('n')
plt.ylabel('y[n]')
plt.stem(n,y)
plt.xticks(n)
plt.show()         
