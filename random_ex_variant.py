import sympy
import numpy 
import pandas as pd
import scipy.stats as stats
import scipy.special 
import matplotlib.pyplot as plt
from scipy import special as spc


def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array

import numpy
from scipy import special as spc

def get_frequency(list_data, trigger):
    frequency = 0
    for (x, y) in list_data:
        if x == trigger:
            frequency = y
    return frequency

def random_excursions_variant(bin_data):
    int_data = numpy.zeros(len(bin_data))
    for i in range(len(bin_data)):
        int_data[i] = int(bin_data[i])
    sum_int = (2 * int_data) - numpy.ones(len(int_data))
    cumulative_sum = numpy.cumsum(sum_int)

    li_data = []
    for xs in sorted(set(cumulative_sum)):
        if numpy.abs(xs) <= 9:
            li_data.append([xs, len(numpy.where(cumulative_sum == xs)[0])])

    j = get_frequency(li_data, 0) + 1
    p_values = []
    for xs in range(-9, 9 + 1):
        if not xs == 0:
            den = numpy.sqrt(2 * j * (4 * numpy.abs(xs) - 2))
            p_values.append(spc.erfc(numpy.abs(get_frequency(li_data, xs) - j) / den))
    return p_values
