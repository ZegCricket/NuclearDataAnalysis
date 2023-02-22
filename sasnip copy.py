import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage import uniform_filter1d

def smoothedDerivative(data, smooth_half_window = 3, pad = True, pad_width = 2):
    if (pad or pad_width < 2):
        data = np.pad(data, pad_width, mode = 'linear_ramp')
    dataDerivative = np.zeros(len(data))

    for i in range(pad_width, len(data) - pad_width):
        dataDerivative[i] = 1/12 * (data[i - 2] - 8 * data[i - 1] + 8 * data[i + 1] - data[i + 2])
    
    dataDerivative = uniform_filter1d(dataDerivative, 2 * smooth_half_window + 1)
    return dataDerivative

def findPeaks(data, max_half_fwhm = 10, derivativeThreshold = 10):
    dataDerivative = smoothedDerivative(data, pad = False, pad_width = max_half_fwhm)
    fwhmArray = np.zeros(len(data))

    for i in range(max_half_fwhm, len(data) - max_half_fwhm):
        if (dataDerivative[i] > 0 and dataDerivative[i + 1] < 0):
            if (max(dataDerivative[i - max_half_fwhm: i]) > derivativeThreshold and min(dataDerivative[i: i + max_half_fwhm]) < -derivativeThreshold):
                leftLimit = np.argmax(dataDerivative[i - max_half_fwhm: i]) + i - max_half_fwhm
                rightLimit = np.argmin(dataDerivative[i: i + max_half_fwhm]) + i
                twosigma = rightLimit - leftLimit + 1
                fwhm = twosigma * np.sqrt(2 * np.log(2))
                leftLimit = i - int(fwhm/2)
                rightLimit = i + int(fwhm/2)
                for j in range(leftLimit, rightLimit + 1):
                    if (fwhm < 2 * max_half_fwhm):
                        fwhmArray[j] = int(fwhm)
    return fwhmArray

def sumy(background, r):
    yin = 0
    yout = 0
    for i in range(len(background)):
        if (r[i] == 0):
            yout += background[i]
        else:
            yin += background[i]
    return yin, yout

def sasnip(data, t = 1, tolerance = 0.005, max_half_fwhm = 10, derivativeThreshold = 10, smooth_half_window = 1, decrease = True):
    newB = 1
    percentage = 1
    spectrumSum = data.sum()
    paddedData = np.pad(data, max_half_fwhm, mode = 'edge')
    #paddedData = uniform_filter1d(paddedData, 2 * smooth_half_window + 1)
    fwhm = findPeaks(paddedData, max_half_fwhm = max_half_fwhm, derivativeThreshold = derivativeThreshold)
    transformed_data = np.log(np.log(np.sqrt(paddedData + 1) + 1) + 1)
    background = np.zeros(transformed_data.size)
    flag = 0

    while (percentage > tolerance):
        if flag != 0:
            background = np.log(np.log(np.sqrt(background + 1) + 1) + 1)
        previousB = newB
        r = t * fwhm
        m = int(max(r))

        if decrease:
            range_args = (m, 0, -1)
        else:
            range_args = (1, m + 1, 1)

        for i in range(*range_args):
            for channel in range(max_half_fwhm, len(data) - max_half_fwhm):
                if (i <= r[channel] or r[channel] == 0): #or r[channel] == 0
                    background[channel] = min(transformed_data[channel], (transformed_data[channel - i] + transformed_data[channel + i])/2)
                # elif(r[channel] == 0):
                #     background[channel] = min(transformed_data[channel], (transformed_data[channel - 1] + transformed_data[channel + 1])/2)
                else:
                    background[channel] = transformed_data[channel]
            transformed_data = background.copy()

        background = -1 + (np.exp(np.exp(background) - 1) - 1) ** 2

        yin, yout = sumy(background, r)
        omega = (yin + yout) / spectrumSum

        newB =  yin * omega + yout * (1 - omega)
        percentage = abs(newB - previousB)/previousB

        fwhm = findPeaks(background, max_half_fwhm = max_half_fwhm, derivativeThreshold = derivativeThreshold)
        flag += 1
    if smooth_half_window != 0:
        background = uniform_filter1d(background, 2 * smooth_half_window + 1)

    return background[max_half_fwhm:-max_half_fwhm]