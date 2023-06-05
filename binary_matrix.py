import sympy
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.special 
from scipy.stats import chi2
def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array

def binaryrank(bit_array):
    M = 3
    Q = 3
    if len(bit_array) % (M*Q) != 0:
        extra_bits = M*Q - len(bit_array) % (M*Q)
        bit_array = np.concatenate((bit_array, np.zeros(extra_bits)))

    # Matris Boyutu
    matrixes = [np.reshape(bit_array[i:i+M*Q], (M,Q)) for i in range(0,len(bit_array),M*Q)] # Matris Oluşturma
    # Rank  Hesaplama
    ranklar = [np.linalg.matrix_rank(matrix) for matrix in matrixes]
    # Rank Dağılımlarını Hesaplama
    rank1_sayisi = ranklar.count(1)
    rank2_sayisi = ranklar.count(2)
    rank3_sayisi = ranklar.count(3)
    # Test İstatistiği
    N = len(matrixes)
    pi = [0.288, 0.5776, 0.1336] # Beklenen rank değerleri
    X_obs = ((rank1_sayisi - N*pi[0])**2/(N*pi[0])) + ((rank2_sayisi - N*pi[1])**2/(N*pi[1])) + ((rank3_sayisi - N*pi[2])**2/(N*pi[2]))
    p_value = chi2.sf(X_obs,2) # 2 Derece Serbestlik

    if p_value < 0.01:
        return False
    else:
        return True

