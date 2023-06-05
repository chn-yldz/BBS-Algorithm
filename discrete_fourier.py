import sympy
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.special
from scipy.fft import fft


def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array


def DFT_test(bits):
    dft_result = fft(bits)
    magnitude = np.abs(dft_result)
    return magnitude