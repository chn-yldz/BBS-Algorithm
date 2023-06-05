import sympy
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.special 

def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array


def serial_test(bits, m):
    # Calculate all possible permutations of m-bits
    bit_permutations = [bin(i)[2:].zfill(m) for i in range(2**m)]

    # Count occurrences of each permutation in the bit string
    counts = {bp: bits.count(bp) for bp in bit_permutations}

    # Number of blocks
    n = len(bits) // m

    # Expected count under null hypothesis
    expected_count = n / 2**m

    # Chi-square test
    chi_square = sum((count - expected_count)**2 / expected_count for count in counts.values())

    # Calculate p-value
    p_value = scipy.special.gammaincc((2**(m-1)) / 2, chi_square / 2)

    return p_value

