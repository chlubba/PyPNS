"""
Code from https://blancosilva.wordpress.com/teaching/mathematical-imaging/denoising-wavelet-thresholding/
"""

import pywt
import numpy as np
import matplotlib.pylab as plt


def wden(signal, threshold=0.5, mode='soft', wavelet='Sym7', level=6):
    wavelet = pywt.Wavelet(wavelet)
    waveletCoeffs = pywt.wavedec(signal, wavelet, level=level)

    newWaveletCoeffs = map(lambda x: pywt.threshold(x, threshold, mode=mode),
                           waveletCoeffs)

    newSignal = pywt.waverec(newWaveletCoeffs, wavelet)

    return newSignal

# x = np.arange(0,10,0.01)
# signal = np.sin(x)
#
# signal += np.random.normal(0, 0.3, size=np.shape(signal))
#
# newSignal = wden(signal)
#
# plt.plot(x, signal)
# plt.plot(x, newSignal, color=[1,0,0], linewidth=2.0)
# plt.show()

