import sympy
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.special 
from scipy.special import gammaincc

def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array

def overlapping_template_matching_test(binary_data, m=9, M=1032, N=968, a=[0,0,0,0,0,0,0,0,1]):
    n = len(binary_data)
    
    # Compute the number of occurrences of the template in each block
    blocks = [binary_data[i * M : (i + 1) * M] for i in range(N)]
    counts = [block.count(''.join(map(str, a))) for block in blocks]
    
    # Compute the test statistic
    mu = (M - m + 1) / 2 ** m
    var = M * (1 / 2 ** m - (2 * m - 1) / 2 ** (2 * m))
    chi_squared = sum([(c - mu) ** 2 / var for c in counts])
    
    # Compute the p-value
    p_value = gammaincc(N / 2, chi_squared / 2)
    
    return p_value
