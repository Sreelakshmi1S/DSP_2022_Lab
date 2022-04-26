import numpy as np
import matplotlib.pyplot as plt

"""
A recursive implementation of the 1D Cooley-Tukey FFT, the 
input should have a length of power of 2. 
"""
def FFT(x):
	N = len(x)    
	if N == 1:
		return x
	else:
		X_even = FFT(x[::2])
		X_odd = FFT(x[1::2])
		factor = np.exp(-2j*np.pi*np.arange(N)/ N)
        
	X = np.concatenate(\
		[X_even+factor[:int(N/2)]*X_odd, X_even+factor[int(N/2):]*X_odd])
	return X

# sampling rate
fs = 128
# sampling interval
ts = 1.0/fs
t = np.arange(0,1,ts)

freq = 1.
x = 3*np.sin(2*np.pi*freq*t)

freq = 4
x += np.sin(2*np.pi*freq*t)

freq = 7
x += 0.5* np.sin(2*np.pi*freq*t)

plt.figure(figsize = (8, 6))
plt.plot(t, x, 'r')
plt.ylabel('Amplitude')
plt.show()

# Compute FFT of x
X=FFT(x)

# calculate the frequency
N = len(X)
k = np.arange(N)
T = N/fs
freq = k/T  

plt.figure(figsize = (12, 6))
plt.subplot(121)
plt.stem(freq, abs(X), 'b',markerfmt=" ")#, basefmt="-b")
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |X(freq)|')

# Get the one-sided specturm
n_oneside = N//2
# get the one side frequency
f_oneside = freq[:n_oneside]

# normalize the amplitude
X_oneside =X[:n_oneside]/n_oneside

plt.subplot(122)
plt.stem(f_oneside, abs(X_oneside), 'b',  markerfmt=" ")#, basefmt="-b")
plt.xlabel('Freq (Hz)')
plt.ylabel('Normalized FFT Amplitude |X(freq)|')
plt.tight_layout()
plt.show()








import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

def sinc_filter(M, fc):
	"""Return an M + 1 point symmetric point sinc kernel with normalised cut-off 
    		frequency fc 0->0.5."""
	if M%2:
		raise Exception('M must be odd')
	return np.sinc(2*fc*(np.arange(M + 1) - M/2))

def hamming(M):
	"""Return an M + 1 point symmetric hamming window."""
	if M%2:
		raise Exception('M must be odd')
	return 0.54 - 0.46*np.cos(2*np.pi*np.arange(M + 1)/M)

def build_filter(M, fc, window=None):
	"""Construct filter using the windowing method for filter parameters M
	number of taps, cut-off frequency fc and window. Window defaults to None 
	i.e. a rectangular window."""
	if window is None:
		h = sinc_filter(M, fc)
	else:
		h = sinc_filter(M, fc)*window(M)
	return h/h.sum()

f0 = 20 #20Hz
ts = 0.01 # i.e. sampling frequency is 1/ts = 100Hz
sr = 1/ts
x = np.arange(-10, 10, ts)
signal = (np.cos(2*np.pi*f0*x) + np.sin(2*np.pi*2*f0*x) + 
                np.cos(2*np.pi*0.5*f0*x) + np.cos(2*np.pi*1.5*f0*x))

#build up some filters
#Low pass
M = 100 #number of taps in filter
fc = 0.25 #i.e. normalised cutoff frequency 1/4 of sampling rate i.e. 25Hz
ham_lp = build_filter(M, fc, window=hamming)
y_ham = np.convolve(signal, ham_lp)

X = fft(signal)
Nx = len(X)
n = np.arange(Nx)
T = Nx/sr
freq = n/T 

plt.figure(figsize = (12, 6))
plt.subplot(121)
plt.stem(freq, abs(X), 'b', \
         markerfmt=" ", basefmt="-b")
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude of Input |X(freq)|')

Y=fft(y_ham)
Ny = len(Y)
n = np.arange(Ny)
T = Ny/sr
freq = n/T 

plt.subplot(122)
plt.stem(freq, abs(Y), 'b', \
         markerfmt=" ", basefmt="-b")
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude of Output |Y(freq)|')
plt.tight_layout()
plt.show()
