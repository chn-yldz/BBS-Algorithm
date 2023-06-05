import sympy
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.special 
import itertools

def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array

def longest(bit_array):
    M = 10 # Blok Boyutu
    blocks = [bit_array[i:i+M] for i in range(0,len(bit_array),M)]

    longest_series = [max((len(list(g)) for k,g in itertools.groupby(block) if k == 1), default=0) for block in blocks ]
    K = 6 # Katsayı İndex'î
    pi = [0.2148, 0.3672, 0.2305, 0.1250, 0.0527, 0.0102] # Katsayılar
    nu = [sum(1 for serie in longest_series if serie == i ) for i in range(1, K+2)]
    N = len(blocks)
    X_obs = sum((nu[i] - N*pi[i])**2 / (N*pi[i]) for i in range(K))
    p_value = scipy.special.gammaincc(K/2,X_obs/2)
    if p_value < 0.01:
        return False
    else:
        return True


