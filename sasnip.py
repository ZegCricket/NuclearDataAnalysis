import numpy as np
import matplotlib.pyplot as plt

def smoothSignal(signal):
    """
    This function smooths the signal with the one which was considered as the best smoothing method for low-statistics spectra. The signal is padded in both sides with the values of the edges in order to return an array with the same size of the signal.

    Parameters
    ----------
    signal : array-like
        The values of the measured data, i.e. the spectrum itself.

    Returns
    -------
    smoothedSignal : numpy.ndarray
        The smoothed signal.
        
    References
    ----------
    .. [1] M. Morháč, “Multidimensional peak searching algorithm for low-statistics nuclear spectra,” Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, vol. 581, no. 3, pp. 821-830, 2007.
    """
    paddedSignal = np.pad(signal, 2, mode = 'edge')
    smoothedSignal = np.zeros(paddedSignal.size)

    for channel in range(2, smoothedSignal.size - 2):
        smoothedSignal[channel] = 1/9 * (paddedSignal[channel - 2] + 2 * paddedSignal[channel - 1] + 3 * paddedSignal[channel] + 2 * paddedSignal[channel + 1] + paddedSignal[channel + 2])

    return smoothedSignal[2: smoothedSignal.size - 2]

def firstDerivative(signal):
    """
    This function calculates de first derivative of the signal. The signal is padded in both sides with the values of the edges in order to return an array with the same size of the signal.

    Parameters
    ----------
    signal : array-like
        The values of the measured data, i.e. the spectrum itself.

    Returns
    -------
    firstDerivative : numpy.ndarray
        The first derivative of the signal.
    """
    paddedSignal = np.pad(signal, 1, mode = 'edge')
    firstDerivative = np.zeros(paddedSignal.size)

    for channel in range(1, firstDerivative.size - 1):
        firstDerivative[channel] = 1/2 * (- paddedSignal[channel - 1] + paddedSignal[channel + 1])

    return firstDerivative[1: firstDerivative.size - 1]

def findPeaks(signal, peakMinimum = 4, peakMaximum = 100, derivativeThreshold = 0):
    """
    This function finds the regions of the peaks and fills an array with the values of the FWHM (full width at half maximum) in those regions of the array. To do so, this function calculates the first derivative of the input signal, evaluates where the derivative crosses y = 0 and determines the regions of the peaks based on whether the crossing happens from a positive value to a negative value or vice-versa. Since the noise of the signal can have a big impact in the peak identification, it is possible to apply some constraints to peak identification, such as minimum and maximum peak sizes (base width, not FWHM) and a threshold value the derivative must reach to be considered a peak (Warning: this parameter may impact greatly the results).

    Parameters
    ----------
    signal : array-like
        The values of the measured data, i.e. the spectrum itself.
    peakMinimum : int, optional
        The minimum base width (not FWHM) that a peak must have.
    peakMaximum : int, optional
        The maximum base width (not FWHM) that a peak must have.
    derivativeThreshold : float, optional
        The value that the derivative must cross (both in the positive and negative sides of the derivative of the peak) in order to consider a certain region as a valid peak. (Warning: this parameter may impact greatly the results)
    
    Returns
    -------
    fwhmArray : numpy.ndarray
        Array containing the values of the FWHM in the base width regions of the peak.
    """
    smoothedSignal = smoothSignal(signal)
    #smoothedSignal = smoothSignal(smoothedSignal)
    #smoothedSignal = smoothSignal(smoothedSignal)
    signalPrime = firstDerivative(smoothedSignal)
    zeroCrossing = np.zeros(signal.size)
    fwhmArray = np.zeros(signal.size)

    for i in range(signal.size - 1):
        if (signalPrime[i] > 0 and signalPrime[i + 1] <= 0):
            zeroCrossing[i] = -1
            zeroCrossing[i + 1] = -1
        elif (signalPrime[i] < 0 and signalPrime[i + 1] >= 0):
            zeroCrossing[i] = 1
            zeroCrossing[i + 1] = 1

    for i in range(signal.size - 1):
        if (zeroCrossing[i] == -1 and zeroCrossing[i + 1] == -1):
            left = 0
            right = 0
            for j in range(1, int(peakMaximum/2)):
                if (i - j == 0 or i + 1 + j == signal.size):
                    break
                if (zeroCrossing[i - j] == 1 and left == 0):
                    left = j
                if (zeroCrossing[i + 1 + j] == 1 and right == 0):
                    right = j
                if (left != 0 and right != 0):
                    #VER O EFEITO DESTE IF
                    if (left + right <= peakMinimum):
                        for k in range(i - left, i + 1 + right):
                            zeroCrossing[k] = 0
                    elif (max(signalPrime[i - left: i]) > derivativeThreshold and min(signalPrime[i + 1: i + 1 + right]) < - derivativeThreshold):
                        fwhmLeft = int((signalPrime[i - left: i].size - np.argmax(signalPrime[i - left: i])) * np.sqrt(2 * np.log(2)))
                        fwhmRight = int((np.argmin(signalPrime[i + 1: i + 1 + right])) * np.sqrt(2 * np.log(2)))
                        for k in range(left + 1):
                            fwhmArray[i - k] = fwhmLeft + fwhmRight
                        for k in range(right + 1):
                            fwhmArray[i + 1 + k] = fwhmLeft + fwhmRight
                    break

    return fwhmArray

def stopCondition(background, tolerance, fwhmArray, signalSum, previousB):
    """
    This function evaluates if the desired background has been estimated according with the tolerance given.

    Parameters
    ----------
    background : array-like
        The background of the spectrum estimated so far.
    tolerance : float
        The value whose evaluation should be inferior to.
    fwhmArray : numpy.ndarray
        Array containing the values of the FWHM in the base width regions of the peak. Used to identify the regions of peaks and regions of no-peaks.
    signalSum : int
        The total sum of counts of the original spectrum.
    previousB : float
        The previous value of parameter B whose new calculated B will be compared to.

    Returns
    -------
    (percentage > tolerance) : bool
        If true, the SASNIP algorithm will continue, if false, it will stop.
    parameterB : float
        The value of the new parameter B calculated in the function.
    """    
    yin = 0
    yout = 0
    for i in range(len(background)):
        if (fwhmArray[i] == 1 or fwhmArray[i] == 0):
            yout += background[i]
        else:
            yin += background[i]

    omega = (yin + yout) / signalSum
    parameterB =  yin * omega + yout * (1 - omega)
    percentage = abs(parameterB - previousB)/previousB
    
    return (percentage > tolerance), parameterB

def sasnip(signal, t = 1, tolerance = 0.005, decrease = True, peakMinimum = 4, peakMaximum = 100, derivativeThreshold = 0, smooth = True):
    """
    This function executes the SASNIP algorithm and calls all the functions needed to do so.

    Parameters
    ----------
    signal : array-like
        The values of the measured data, i.e. the spectrum itself.
    t : int, optional
        A scalar that can be multiplied to the array of the FWHM of the peaks. The default value is recommended.
    tolerance : float, optional
        The value whose evaluation made in the stopCondition function should be inferior to. The default value is recommended.
    decrease : bool, optional
        This variable determines whether the background estimation in a point starts in the closest or in the furthest neighbours. To a better understanding, read [2].
    peakMinimum : int, optional
        The minimum base width (not FWHM) that a peak must have. Called in the findPeaks function.
    peakMaximum : int, optional
        The maximum base width (not FWHM) that a peak must have. Called in the findPeaks function.
    derivativeThreshold : float, optional
        The value that the derivative must cross (both in the positive and negative sides of the derivative of the peak) in order to consider a certain region as a valid peak. Called in the findPeaks function. (Warning: this parameter may impact greatly the results)
    smooth : bool, optional
        If true, an initial smooth will be applied to the raw signal. If false, the smooth will only be applied in the functions that need it to obtain better results (example: findPeaks function).

    Returns
    -------
    background : np.ndarray
        Array with the final background estimated for the input signal.
    
    References
    ----------
    .. [1] R. Shi, X. Tuo, H. Zheng et al., “Step-approximation SNIP background-elimination algorithm for HPGe gamma spectra,” Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, vol. 885, pp. 60-66, 2018.
    .. [2] M. Morháč and V. Matoušek, “Peak clipping algorithms for background estimation in spectroscopic data,” Applied spectroscopy, vol. 62, no. 1, pp. 91-106, 2008.
    """   
    parameterB = 1
    fwhm = findPeaks(signal, peakMinimum, peakMaximum, derivativeThreshold)
    r = t * fwhm
    m = int(max(r))
    continueCondition = True
    signalSum = signal.sum()
    if smooth:
        signal = smoothSignal(signal)
    llsSignal = np.log(np.log(np.sqrt(signal + 1) + 1) + 1)
    baseline = np.zeros(llsSignal.size)

    if decrease:
        range_args = (m, 0, -1)
    else:
        range_args = (1, m + 1, 1)

    while (continueCondition):
        for i in range(*range_args):
            for channel in range(signal.size):
                if (i <= r[channel]):
                    try:
                        baseline[channel] = min(llsSignal[channel], (llsSignal[channel - i] + llsSignal[channel + i])/2)
                    except:
                        baseline[channel] = llsSignal[channel]
                else:
                    baseline[channel] = llsSignal[channel]
            llsSignal = baseline.copy()

        background = -1 + (np.exp(np.exp(baseline) - 1) - 1) ** 2
        continueCondition, parameterB = stopCondition(background, tolerance, fwhm, signalSum, parameterB)

    return background