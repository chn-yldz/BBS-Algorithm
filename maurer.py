import sympy
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.special 

def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array

def maurers_universal_test(bitstring, L=7, Q=1280):
    """
    Maurer's Universal Statistical Test

    Parameters
    ----------
    bitstring : str
        The bitstring to be tested for randomness.
    L : int, optional
        The length of each block. The default is 7.
    Q : int, optional
        The number of blocks. The default is 1280.

    Returns
    -------
    float
        The p-value from the test.
    """

    # Convert the bitstring to an array of integers
    bitstring = np.array([int(bit) for bit in bitstring])
    
    # Calculate the number of blocks
    n = len(bitstring)
    K = n // L - Q
    if K <= 0 or Q <= 0:
        raise ValueError("The length of the bitstring is not appropriate for this test.")
    
    # Initialize the array of blocks
    blocks = np.zeros((2**L,), dtype=int)
    
    # Populate the array with the indices of the blocks
    for i in range(Q):
        block = bitstring[i*L:(i+1)*L]
        j = int("".join(str(b) for b in block), 2)
        blocks[j] = i + 1
    
    # Calculate the sum of the log2 values
    sum_log2 = 0
    for i in range(Q, Q+K):
        block = bitstring[i*L:(i+1)*L]
        j = int("".join(str(b) for b in block), 2)
        sum_log2 += np.log2(i + 1 - blocks[j])
        blocks[j] = i + 1
    
    # Calculate the test statistic
    fn = sum_log2 / K
    c = 0.7 - 0.8 / L + (4 + 32 / L) * (K**(-3 / L)) / 15
    sigma = c * np.sqrt((1 / L) * (32 / L - 2))
    test_statistic = (fn - 0.95) / sigma
    
    # Calculate the p-value
    p_value = scipy.special.ndtr(test_statistic)
    
    return p_value