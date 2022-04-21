""""6) Write a python program to compute discrete Fourier transform (DFT) of a signal, without using
library functions for DFT. Let x(t) = 3 sin(4πt) + cos(10πt) + 0.5 cos(16πt) containing frequencies
2Hz, 5Hz and 8Hz. Use the program to plot the one-sided magnitude spectrum of y(t) = x(t) +
x(t) cos(32πt) showing analog frequency on the X-axis. Assume sampling frequency of Fs =
128Hz. Identify the frequency components in y(t) by plotting its magnitude spectrum."""



mport numpy as np
import matplotlib.pyplot as plt 
fs = 128
# sampling interval
ts = 1.0/fs
t = np.arange(0,1,ts)

freq = 1.
y = 3*np.sin(4*np.pi*t)+np.cos(10*np.pi*t)+0.5*np.cos(16*np.pi*t)

x=y+y*(np.cos(32*np.pi*t))

x=x.T
N=len(x)
D=np.empty((N,N),dtype=np.cdouble)
W=np.exp(-1j*2*np.pi/N)
for k in np.arange(N):
    for n in np.arange(N):
        D[k,n]=W**(k*n)
np.round(D)

X=D@x
np.round(X)
print(X)
N = len(X)
k = np.arange(N)
T = N/fs
freq = k/T  


n_oneside = N//2
# get the one side frequency
f_oneside = freq[:n_oneside]

# normalize the amplitude
X_oneside =X[:n_oneside]/n_oneside


plt.stem(f_oneside, abs(X_oneside), 'b',  markerfmt=" ")#, basefmt="-b")
plt.xlabel('Freq (Hz)')
plt.ylabel('Normalized FFT Amplitude |X(freq)|')
plt.tight_layout()
plt.show()

