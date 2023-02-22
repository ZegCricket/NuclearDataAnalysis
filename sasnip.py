import numpy as np
import matplotlib.pyplot as plt

def smoothSignal(signal):
    """
    DESCRIPTION

    Parameters
    ----------
    signal : 
        DESCRIPTION

    Returns
    -------
    smoothedSignal : 
        DESCRIPTION
        
    References
    ----------
    .. [1] 
    """
    paddedSignal = np.pad(signal, 2, mode = 'edge')
    smoothedSignal = np.zeros(paddedSignal.size)

    for channel in range(2, smoothedSignal.size - 2):
        smoothedSignal[channel] = 1/9 * (paddedSignal[channel - 2] + 2 * paddedSignal[channel - 1] + 3 * paddedSignal[channel] + 2 * paddedSignal[channel + 1] + paddedSignal[channel + 2])

    return smoothedSignal[2: smoothedSignal.size - 2]

def firstDerivative(signal):
    """
    DESCRIPTION

    Parameters
    ----------
    signal : 
        DESCRIPTION

    Returns
    -------
    firstDerivative : 
        DESCRIPTION
    """
    paddedSignal = np.pad(signal, 1, mode = 'edge')
    firstDerivative = np.zeros(paddedSignal.size)

    for channel in range(1, firstDerivative.size - 1):
        firstDerivative[channel] = 1/2 * (- paddedSignal[channel - 1] + paddedSignal[channel + 1])

    return firstDerivative[1: firstDerivative.size - 1]

def findPeaks(signal, peakMinimum = 4, peakMaximum = 100, derivativeThreshold = 0):
    """
    DESCRIPTION

    Parameters
    ----------
    signal : 
        DESCRIPTION

    Returns
    -------
    firstDerivative : 
        DESCRIPTION
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
    DESCRIPTION

    Parameters
    ----------
    signal : 
        DESCRIPTION

    Returns
    -------
    firstDerivative : 
        DESCRIPTION
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

def sasnip(signal, t = 1, tolerance = 0.005, decrease = True, peakMinimum = 4, peakMaximum = 100, derivativeThreshold = 0):
    """
    DESCRIPTION

    Parameters
    ----------
    signal : 
        DESCRIPTION

    Returns
    -------
    firstDerivative : 
        DESCRIPTION
    """   
    parameterB = 1
    fwhm = findPeaks(signal, peakMinimum, peakMaximum, derivativeThreshold)
    r = t * fwhm
    m = int(max(r))
    continueCondition = True
    signalSum = signal.sum()
    #signal = smoothSignal(signal)
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



# def drawTurningPoints(crossingPoints):
#     for channel in range(crossingPoints.size):
#         if (crossingPoints[channel] == 1):
#             plt.axvspan(channel, channel, color = 'red', alpha = 0.3)
#         elif (crossingPoints[channel] == -1):
#             plt.axvspan(channel, channel, color = 'blue', alpha = 0.3)

# if __name__ == '__main__':
#     data = np.genfromtxt(r"RBS1LiF_Ag2_660keV_8uC_182.dat", comments = '$')
#     data = np.reshape(data, -1)
#     data = data[10:]
#     zeros = np.zeros(data.size)
    
#     plt.plot(data)
#     plt.xlim(0, data.size)

#     background = sasnip(data)
#     data_no_bg = data - background
#     plt.plot(background)
#     plt.plot(data_no_bg)

#     plt.show()