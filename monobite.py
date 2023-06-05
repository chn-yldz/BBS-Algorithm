import sympy
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.special 

def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array



def monobit_test(bit_sequence):
    if not bit_sequence:
        raise ValueError("Bit sequence is empty")
    
    count = 0
    for bit in bit_sequence:
        if bit == 1:
            count += 1
    n = len(bit_sequence)
    S_obs = abs(count - n/2)
    S = S_obs / (n**0.5)
    return S < 2.16