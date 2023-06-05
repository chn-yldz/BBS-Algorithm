import sympy
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.special 

def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array

def runs_test(bit_array):
    n = len(bit_array)
    s = sum(bit_array)
    pi = s / n # 1'lerin oranı
    tau = 2 / np.sqrt(n)

    if abs(pi - 0.5) >= tau:
        return False
    else:
        runs = 1 + sum(bit_array[i] != bit_array[i+1] for i in range(n-1))
        R_obs = abs(runs - 2*n*pi*(1-pi)) / np.sqrt(2*n*pi*(1-pi))  # pi*(1-pi) olarak düzeltildi
        p_value = 2 * stats.norm.cdf(-R_obs)
        if p_value < 0.01:
            return False
        else:
            return True
