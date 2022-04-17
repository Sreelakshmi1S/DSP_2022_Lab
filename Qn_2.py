"""A system with input x[n] produces output

y[n] = 2x[n] + x[n − 1] + 0.5x[n − 2]

Three such systems are cascaded as a chain to produce output z[n]. Write a python program that
takes x[n] as input and produces z[n] as output."""

import numpy as np
import matplotlib.pyplot as plt
#impulse response
h= [2,1,0.5]
#input response
x = []
 
# number of elements as input
n = int(input("Enter number of elements : "))
 
# iterating till the range
for i in range(0, n):
    ele = int(input())
 
    x.append(ele) # adding the element
     
print(x)
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
z=dirconv(h,dirconv(h,dirconv(x,h)))
print('Linear convolution of cascaded inputs, z =\n',z)
plt.title("linear convolution")
n=np.arange(len(z))
plt.xlabel('n')
plt.ylabel('z[n]')
plt.stem(n,z)
plt.xticks(n)
plt.show()     
