import numpy as np
import matplotlib.pyplot as plt

def findPeaks(data, peakMinimum = 4, peakMaximum = 100, derivativeThreshold = 20):
    firstDerivative = np.zeros(data.size)
    zeroCrossing = np.zeros(data.size)
    teste = np.zeros(data.size)
    r = np.ones(data.size)

    for i in range(4, data.size - 4):
        teste[i] = 1/9 * (data[i - 2] + 2 * data[i - 1] + 3 * data[i] + 2 * data[i + 1] + data[i + 2])

    # for i in range(4, data.size - 4):
    #     teste[i] = data[i]

    for i in range(4, data.size - 4):
        firstDerivative[i] = 1/2 * (- teste[i - 1] + teste[i + 1])

    for i in range(data.size - 1):
        if (firstDerivative[i] > 0 and firstDerivative[i + 1] < 0):
            zeroCrossing[i] = -1
            zeroCrossing[i + 1] = -1
        elif (firstDerivative[i] < 0 and firstDerivative[i + 1] > 0):
            zeroCrossing[i] = 1
            zeroCrossing[i + 1] = 1

    for i in range(data.size):
        if (zeroCrossing[i] == -1 and zeroCrossing[i + 1] == -1):
            left = 0
            right = 0
            j = 1
            while (j < int(peakMaximum/2)):
                if (i - j == 0 or i + 1 + j == data.size):
                    break
                if (zeroCrossing[i - j] == 1 and left == 0):
                    left = j
                if (zeroCrossing[i + 1 + j] == 1 and right == 0):
                    right = j
                if (left != 0 and right != 0):
                    if (left + right <= peakMinimum):
                        for k in range(i - left, i + 1 + right):
                            zeroCrossing[k] = 0
                    elif (max(firstDerivative[i - left: i]) > derivativeThreshold and min(firstDerivative[i + 1: i + 1 + right]) < - derivativeThreshold):
                        left = firstDerivative[i - left: i].size - np.argmax(firstDerivative[i - left: i])
                        right = np.argmin(firstDerivative[i + 1: i + 1 + right])
                        left = int(left * np.sqrt(2 * np.log(2)))
                        right = int(right * np.sqrt(2 * np.log(2)))
                        for k in range(left + 1):
                            r[i - k] = left + right
                        for k in range(right + 1):
                            r[i + 1 + k] = left + right
                    break
                j += 1

    return r

def sumy(background, r):
    yin = 0
    yout = 0
    for i in range(len(background)):
        if (r[i] == 1 or r[i] == 0):
            yout += background[i]
        else:
            yin += background[i]
    return yin, yout

def sasnip(data, t = 1, tolerance = 0.005, max_half_fwhm = 10, derivativeThreshold = 10, smooth_half_window = 1, decrease = True):
    newB = 1
    percentage = 1
    spectrumSum = data.sum()
    paddedData = np.zeros(data.size)

    # for i in range(4, data.size - 4):
    #     paddedData[i] = 1/9 * (data[i - 2] + 2 * data[i - 1] + 3 * data[i] + 2 * data[i + 1] + data[i + 2])
    for i in range(4, data.size - 4):
        paddedData[i] = data[i]

    paddedData = np.pad(paddedData, max_half_fwhm, mode = 'edge')
    #paddedData = uniform_filter1d(paddedData, 2 * smooth_half_window + 1)
    fwhm = findPeaks(paddedData)
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
                if (i <= r[channel]): #or r[channel] == 0
                    try:
                        background[channel] = min(transformed_data[channel], (transformed_data[channel - i] + transformed_data[channel + i])/2)
                    except:
                        background[channel] = transformed_data[channel]
                else:
                    background[channel] = transformed_data[channel]
            transformed_data = background.copy()

        background = -1 + (np.exp(np.exp(background) - 1) - 1) ** 2
        #plt.plot(background[max_half_fwhm:-max_half_fwhm])

        yin, yout = sumy(background, r)
        omega = (yin + yout) / spectrumSum

        newB =  yin * omega + yout * (1 - omega)
        percentage = abs(newB - previousB)/previousB
        print(percentage, newB, previousB, omega)

        #fwhm = findPeaks(background)
        flag += 1

    return background[max_half_fwhm:-max_half_fwhm]

if __name__ == '__main__':
    data = np.genfromtxt(r"ERD_7Li_imp_Al_1800_500uC_1.dat", comments = '$')
    data = np.reshape(data[10:], -1)
    firstDerivative = np.zeros(data.size)
    zeroCrossing = np.zeros(data.size)
    teste = np.zeros(data.size)
    zeros = np.zeros(data.size)
    r = np.zeros(data.size)
    background = sasnip(data, decrease = True, derivativeThreshold=1)

    plt.plot(data)
    plt.plot(background)

    plt.show()