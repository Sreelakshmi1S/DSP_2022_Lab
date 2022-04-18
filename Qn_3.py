"""Write a python program to simulate a low-pass FIR filter using window method with parameters:

Order of filter, M = 100
Normalized cut-off frequency, ωc = 0.2

The sampling frequency Fs = 200Hz. You may use Hamming window with parameter a whose
equation is given below:

wH[n] = a − (1 − a) cos(2 pi n/M), 0 ≤ n ≤ M.

Use the program to simulate filtering of x(t) = sin(20πt) + cos(30πt) + sin(80πt). Plot magnitude
spectrum of both x(t) and filtered output y(t) showing analog frequency on the X-axis. You may
use scipy.fft.fft to compute DFT. Use values a = 0.54 and a = 0.35. Verify which one performs
better."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

def sinc_filter(M, fc):
    if M%2:
       raise Exception('M must be odd')
    return np.sinc(2*fc*(np.arange(M + 1) - M/2))

def hamming(M):
    a=0.35
    if M%2:
       raise Exception('M must be odd')
    return (a - (1-a)*np.cos(2*np.pi*np.arange(M + 1)/M))

def build_filter(M, fc, window):
    if window is None:
       h = sinc_filter(M, fc)
    else:
       h = sinc_filter(M, fc)*window(M)
    return h/h.sum()

ts = 1/200
sr = 1/ts
x = np.arange(-10, 10, ts)
signal = (np.cos(np.pi*30*x) + np.sin(2*np.pi*10*x) + 
                np.sin(np.pi*80*x))

#build up some filters
#Low pass
M = 100 #number of taps in filter
fc = 0.2/(np.pi*2) #i.e. normalised cutoff frequency 1/4 of sampling rate i.e. 25Hz
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
