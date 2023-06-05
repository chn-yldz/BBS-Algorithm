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


def get_pik_value(k, x):
    if k == 0:
        out = 1 - 1.0 / (2 * numpy.abs(x))
    elif k >= 5:
        out = (1.0 / (2 * numpy.abs(x))) * (1 - 1.0 / (2 * numpy.abs(x))) ** 4
    else:
        out = (1.0 / (4 * x * x)) * (1 - 1.0 / (2 * numpy.abs(x))) ** (k - 1)
    return out

def random_excursions(bin_data):
    int_data = numpy.zeros(len(bin_data))
    for i in range(len(bin_data)):
        if bin_data[i] == '0':
            int_data[i] = -1.0
        else:
            int_data[i] = 1.0

    cumulative_sum = numpy.cumsum(int_data)
    cumulative_sum = numpy.append(cumulative_sum, [0])
    cumulative_sum = numpy.append([0], cumulative_sum)

    x_values = numpy.array([-4, -3, -2, -1, 1, 2, 3, 4])

    position = numpy.where(cumulative_sum == 0)[0]
    cycles = []
    for pos in range(len(position) - 1):
        cycles.append(cumulative_sum[position[pos]:position[pos + 1] + 1])
    num_cycles = len(cycles)

    state_count = []
    for cycle in cycles:
        state_count.append(([len(numpy.where(cycle == state)[0]) for state in x_values]))
    state_count = numpy.transpose(numpy.clip(state_count, 0, 5))

    su = []
    for cycle in range(6):
        su.append([(sct == cycle).sum() for sct in state_count])
    su = numpy.transpose(su)

    piks = ([([get_pik_value(uu, state) for uu in range(6)]) for state in x_values])
    inner_term = num_cycles * numpy.array(piks)
    chi = numpy.sum(1.0 * (numpy.array(su) - inner_term) ** 2 / inner_term, axis=1)
    p_values = ([spc.gammaincc(2.5, cs / 2.0) for cs in chi])

    return p_values
