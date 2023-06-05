import sympy
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.special 

def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array


def frequency(bit_array):
    M = 200
    n = len(bit_array)
    N = n // M
    blocks = [bit_array[i*M:(i+1)*M] for i in range(N)]
    X = [sum(block) - M/2 for block in blocks]
    chi_squared = 4/M * sum(xi*xi for xi in X)
    p_value = 1 - stats.chi2.cdf(chi_squared,N)
    if p_value < 0.01:
        return False
    else:
        return True        
