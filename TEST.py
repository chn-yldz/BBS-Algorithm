import numpy as np
import itertools
import matplotlib.pyplot as plt
import scipy.special as sps
from io import StringIO
import scipy.stats as stats
from sympy import isprime
import random
from scipy.stats import chi2
import pandas as pd
from scipy.stats import chisquare
import scipy.special
from scipy.fft import fft
from scipy.special import gammaincc

class MonobitTest:
    def __init__(self, bit_sequence):
        self.bit_sequence = bit_sequence

    def perform_test(self):
        if not self.bit_sequence:
            raise ValueError("Bit sequence is empty")
        
        count = 0
        for bit in self.bit_sequence:
            if bit == 1:
                count += 1
        n = len(self.bit_sequence)
        S_obs = abs(count - n/2)
        S = S_obs / (n**0.5)
        return f"P Değeri : {S}"


class FrequencyTest:
    def __init__(self, bit_array):
        self.bit_array = bit_array

    def perform_test(self):
        M = 200
        n = len(self.bit_array)
        N = n // M
        blocks = [self.bit_array[i*M:(i+1)*M] for i in range(N)]
        X = [sum(block) - M/2 for block in blocks]
        chi_squared = 4/M * sum(xi*xi for xi in X)
        p_value = 1 - stats.chi2.cdf(chi_squared, N)
        if p_value < 0.01:
            return f"P Değeri : {p_value}"
        else:
            return f"P Değeri : {p_value}"
     

class RunsTest:
    def __init__(self, bit_array):
        self.bit_array = bit_array

    def perform_test2(self):
        n = len(self.bit_array)
        s = sum(self.bit_array)
        pi = s / n  # 1'lerin oranı
        tau = 2 / np.sqrt(n)

        if abs(pi - 0.5) >= tau:
            return False
        else:
            runs = 1 + sum(self.bit_array[i] != self.bit_array[i + 1] for i in range(n - 1))
            R_obs = abs(runs - 2 * n * pi * (1 - pi)) / np.sqrt(2 * n * pi * (1 - pi))  # pi*(1-pi) olarak düzeltildi
            p_value = 2 * stats.norm.cdf(-R_obs)
            if p_value < 0.01:
                return f"P Değeri : {p_value}"
            else:
                return f"P Değeri : {p_value}"


class LongestSeriesTest:
    def __init__(self, bit_array):
        self.bit_array = bit_array

    def perform_test3(self):
        M = 10  # Blok Boyutu
        blocks = [self.bit_array[i:i + M] for i in range(0, len(self.bit_array), M)]

        longest_series = [max((len(list(g)) for k, g in itertools.groupby(block) if k == 1), default=0) for block in blocks]
        K = 6  # Katsayı İndex'i
        pi = [0.2148, 0.3672, 0.2305, 0.1250, 0.0527, 0.0102]  # Katsayılar
        nu = [sum(1 for serie in longest_series if serie == i) for i in range(1, K + 2)]
        N = len(blocks)
        X_obs = sum((nu[i] - N * pi[i]) ** 2 / (N * pi[i]) for i in range(K))
        p_value = scipy.special.gammaincc(K / 2, X_obs / 2)
        if p_value < 0.01:
            return f"P Değeri : {p_value}"
        else:
            return f"P Değeri : {p_value}"


class BinaryRankTest:
    def __init__(self, bit_array):
        self.bit_array = bit_array

    def perform_test5(self):
        M = 3
        Q = 3
        if len(self.bit_array) % (M*Q) != 0:
            extra_bits = M*Q - len(self.bit_array) % (M*Q)
            self.bit_array = np.concatenate((self.bit_array, np.zeros(extra_bits)))

        # Matris Boyutu
        matrixes = [np.reshape(self.bit_array[i:i+M*Q], (M,Q)) for i in range(0, len(self.bit_array), M*Q)]
        # Rank Hesaplama
        ranks = [np.linalg.matrix_rank(matrix) for matrix in matrixes]
        # Rank Dağılımlarını Hesaplama
        rank1_count = ranks.count(1)
        rank2_count = ranks.count(2)
        rank3_count = ranks.count(3)
        # Test İstatistiği
        N = len(matrixes)
        pi = [0.288, 0.5776, 0.1336] # Beklenen rank değerleri
        X_obs = ((rank1_count - N*pi[0])**2/(N*pi[0])) + ((rank2_count - N*pi[1])**2/(N*pi[1])) + ((rank3_count - N*pi[2])**2/(N*pi[2]))
        p_value = 1 - chi2.cdf(X_obs, 2) # 2 Derece Serbestlik

        if p_value < 0.01:
            return f"P Değeri: {p_value}  -->  False"
        else:
            return f"P Değeri: {p_value}  -->  True"



class DFTTest:
    def __init__(self, bits):
        self.bits = bits

    def perform_test6(self):
        dft_result = fft(self.bits)
        magnitude = np.abs(dft_result)
        return f"P Değeri : {magnitude}"


class NonOverlappingTemplateMatching:
    def __init__(self, bits, m, b):
        self.bits = bits
        self.m = m
        self.b = b

    def perform_test7(self):
        n = len(self.bits)
        M = self.m
        N = n // M
        blocks = [self.bits[i * M:(i + 1) * M] for i in range(N)]
        counts = [block.count(self.b) for block in blocks]

        # Statistic
        mu = (M - len(self.b) + 1) / 2**len(self.b)
        var = M * (1 / 2**len(self.b) - (2 * len(self.b) - 1) / 2**(2 * len(self.b)))
        chi2 = sum([(c - mu)**2 / var for c in counts])

        # Compute the p-value
        p_value = sps.gammaincc(N / 2, chi2 / 2)

        return f"P Değeri : {p_value}"




class OverlappingTemplateMatchingTest:
    def __init__(self, binary_data, m=9, M=1032, N=968, a=[0, 0, 0, 0, 0, 0, 0, 0, 1]):
        self.binary_data = binary_data
        self.m = m
        self.M = M
        self.N = N
        self.a = a

    def perform_test(self):
        n = len(self.binary_data)

        # Compute the number of occurrences of the template in each block
        blocks = [self.binary_data[i * self.M: (i + 1) * self.M] for i in range(self.N)]
        counts = [block.count(''.join(map(str, self.a))) for block in blocks]

        # Compute the test statistic
        mu = (self.M - self.m + 1) / 2 ** self.m
        var = self.M * (1 / 2 ** self.m - (2 * self.m - 1) / 2 ** (2 * self.m))
        chi_squared = sum([(c - mu) ** 2 / var for c in counts])

        # Compute the p-value
        p_value = gammaincc(self.N / 2, chi_squared / 2)

        return f"P Değeri : {p_value}"


class MaurersUniversalTest:
    def __init__(self, bitstring, L=7, Q=1280):
        self.bitstring = bitstring
        self.L = L
        self.Q = Q

    def perform_test(self):
        bitstring = np.array([int(bit) for bit in self.bitstring])

        n = len(bitstring)
        K = n // self.L - self.Q
        if K <= 0 or self.Q <= 0:
            raise ValueError("The length of the bitstring is not appropriate for this test.")

        blocks = np.zeros((2**self.L,), dtype=int)

        for i in range(self.Q):
            block = bitstring[i*self.L:(i+1)*self.L]
            j = int("".join(str(b) for b in block), 2)
            blocks[j] = i + 1

        sum_log2 = 0
        for i in range(self.Q, self.Q+K):
            block = bitstring[i*self.L:(i+1)*self.L]
            j = int("".join(str(b) for b in block), 2)
            sum_log2 += np.log2(i + 1 - blocks[j])
            blocks[j] = i + 1

        fn = sum_log2 / K
        c = 0.7 - 0.8 / self.L + (4 + 32 / self.L) * (K**(-3 / self.L)) / 15
        sigma = c * np.sqrt((1 / self.L) * (32 / self.L - 2))
        test_statistic = (fn - 0.95) / sigma

        p_value = scipy.special.ndtr(test_statistic)

        return f"P Değeri : {p_value}"

class Result:
    def __init__(self, name, pass_test, p_value):
        self.name = name
        self.pass_test = pass_test
        self.p_value = p_value

class Test:
    def __init__(self, name, significance_value):
        self.name = name
        self.significance_value = significance_value

class LinearComplexityTest(Test):  
    def __init__(self):  
        self._sequence_size_min: int = 1000000  
        self._pattern_length: int = 512  
        self._freedom_degrees: int = 6  
        self._probabilities: np.ndarray = np.array([0.010417, 0.03125, 0.125, 0.5, 0.25, 0.0625, 0.020833])  
        self._mu: float = (self._pattern_length / 2.0) + (((-1) ** (self._pattern_length + 1)) + 9.0) / 36.0 - ((self._pattern_length / 3.0) + (2.0 / 9.0)) / (2 ** self._pattern_length)  
        self._last_bits_size: int = -1  
        self._blocks_number: int = -1  
        super(LinearComplexityTest, self).__init__("Linear Complexity", 0.01)  

    def _execute(self, bits: np.ndarray) -> Result:  
        if self._last_bits_size == -1 or self._last_bits_size != bits.size:  
            blocks_number: int = int(bits.size // self._pattern_length)  
            self._last_bits_size = bits.size  
            self._blocks_number = blocks_number  
        else:  
            blocks_number: int = self._blocks_number  
        blocks_linear_complexity: np.ndarray = np.zeros(blocks_number, dtype=int)  
        for i in range(blocks_number):  
            blocks_linear_complexity[i] = self._berlekamp_massey(bits[(i * self._pattern_length):((i + 1) * self._pattern_length)])  
        tickets: np.ndarray = ((-1.0) ** self._pattern_length) * (blocks_linear_complexity[:] - self._mu) + (2.0 / 9.0)  
        frequencies: np.ndarray = np.zeros(self._freedom_degrees + 1, dtype=int)  
        for ticket in tickets:  
            frequencies[min(self._freedom_degrees, int(max(-2.5, ticket) + 2.5))] += 1  
        chi_square: float = float(np.sum(((frequencies[:] - (blocks_number * self._probabilities[:])) ** 2.0) / (blocks_number * self._probabilities[:])))  
        score: float = scipy.special.gammaincc((self._freedom_degrees / 2.0), (chi_square / 2.0))  
        if score >= self.significance_value:  
            return Result(self.name, True, np.array(score))  
        return Result(self.name, False, np.array(score))  

    def is_eligible(self, bits: np.ndarray) -> bool:  
        if bits.size < self._sequence_size_min:  
            return False  
        return True  

    @staticmethod  
    def _berlekamp_massey(sequence: np.ndarray) -> int:  
        b: np.ndarray = np.zeros(sequence.size, dtype=int)  
        c: np.ndarray = np.zeros(sequence.size, dtype=int)  
        b[0] = 1  
        c[0] = 1  
        generator_length: int = 0  
        m: int = -1  
        n: int = 0  
        while n < sequence.size:  
            discrepancy = sequence[n]  
            for j in range(1, generator_length + 1):  
                discrepancy: int = discrepancy ^ (c[j] & sequence[n - j])  
            if discrepancy != 0:  
                t = c[:]  
                for j in range(0, sequence.size - n + m):  
                    c[n - m + j] = c[n - m + j] ^ b[j]  
                if generator_length <= n / 2:  
                    generator_length = n + 1 - generator_length  
                    m = n  
                    b = t  
                n = n + 1  
        return generator_length



class SerialTest:
    def __init__(self, m):
        self.m = m

    def perform_test(self, bits):
        bit_permutations = [bin(i)[2:].zfill(self.m) for i in range(2**self.m)]
        counts = {bp: bits.count(bp) for bp in bit_permutations}
        n = len(bits) // self.m
        expected_count = n / 2**self.m
        chi_square = sum((count - expected_count)**2 / expected_count for count in counts.values())
        p_value = scipy.special.gammaincc((2**(self.m-1)) / 2, chi_square / 2)
        return p_value




class ApEnFast:
    def __init__(self, m, r):
        self.m = m
        self.r = r

    def perform_calculation(self, U):
        def _phi(m):
            histogram = np.bincount([int(''.join(map(str, U[i:i+m])), 2) for i in range(len(U)-m+1)])
            histogram = histogram[histogram > 0]
            p = histogram / float(sum(histogram))
            return -sum(p * np.log(p))

        return abs(_phi(self.m+1) - _phi(self.m))


class CusumTest:
    def __init__(self):
        self.mu = 0.5

    def perform_test(self, array):
        array = np.array(array)
        cusum = np.cumsum(array - self.mu)
        return f"P Değerleri : {cusum}"
    
import numpy as np
import scipy.special as spc

class RandomExcursionsTest:
    def __init__(self):
        self.x_values = np.array([-4, -3, -2, -1, 1, 2, 3, 4])

    def perform_test(self, bin_data):
        int_data = np.zeros(len(bin_data))
        for i in range(len(bin_data)):
            if bin_data[i] == '0':
                int_data[i] = -1.0
            else:
                int_data[i] = 1.0

        cumulative_sum = np.cumsum(int_data)
        cumulative_sum = np.append(cumulative_sum, [0])
        cumulative_sum = np.append([0], cumulative_sum)

        position = np.where(cumulative_sum == 0)[0]
        cycles = []
        for pos in range(len(position) - 1):
            cycles.append(cumulative_sum[position[pos]:position[pos + 1] + 1])
        num_cycles = len(cycles)

        state_count = []
        for cycle in cycles:
            state_count.append(([len(np.where(cycle == state)[0]) for state in self.x_values]))
        state_count = np.transpose(np.clip(state_count, 0, 5))

        su = []
        for cycle in range(6):
            su.append([(sct == cycle).sum() for sct in state_count])
        su = np.transpose(su)

        piks = ([([self.get_pik_value(uu, state) for uu in range(6)]) for state in self.x_values])
        inner_term = num_cycles * np.array(piks)
        chi = np.sum(1.0 * (np.array(su) - inner_term) ** 2 / inner_term, axis=1)
        p_values = ([spc.gammaincc(2.5, cs / 2.0) for cs in chi])

        return p_values

    @staticmethod
    def get_pik_value(k, x):
        if k == 0:
            out = 1 - 1.0 / (2 * np.abs(x))
        elif k >= 5:
            out = (1.0 / (2 * np.abs(x))) * (1 - 1.0 / (2 * np.abs(x))) ** 4
        else:
            out = (1.0 / (4 * x * x)) * (1 - 1.0 / (2 * np.abs(x))) ** (k - 1)
        return out

class RandomExcursionsVariantTest:
    def __init__(self):
        self.x_values = list(range(-9, 10))
        self.x_values.remove(0)

    def perform_test(self, bin_data):
        int_data = np.array([int(bit) for bit in bin_data])
        sum_int = (2 * int_data) - np.ones(len(int_data))
        cumulative_sum = np.cumsum(sum_int)

        li_data = []
        for xs in sorted(set(cumulative_sum)):
            if np.abs(xs) <= 9:
                li_data.append([xs, len(np.where(cumulative_sum == xs)[0])])

        j = self.get_frequency(li_data, 0) + 1
        p_values = []
        for xs in self.x_values:
            den = np.sqrt(2 * j * (4 * np.abs(xs) - 2))
            p_values.append(spc.erfc(np.abs(self.get_frequency(li_data, xs) - j) / den))
        return p_values

    @staticmethod
    def get_frequency(list_data, trigger):
        frequency = 0
        for (x, y) in list_data:
            if x == trigger:
                frequency = y
        return frequency