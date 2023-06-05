import sympy
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.special 
from scipy.spatial.distance import cdist
def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array

import numpy as np

def ApEn_fast(U, m, r):
    """Approximate Entropy hesaplama fonksiyonu (hızlı versiyon)"""
    def _phi(m):
        # Verinin m-uzunluklu histogramını hesapla
        histogram = np.bincount([int(''.join(map(str, U[i:i+m])), 2) for i in range(len(U)-m+1)])
        # Sadece frekansı 0'dan büyük olan değerleri al
        histogram = histogram[histogram > 0]
        # Frekansları olasılıklara dönüştür
        p = histogram / float(sum(histogram))
        return -sum(p * np.log(p))

    return abs(_phi(m+1) - _phi(m))

