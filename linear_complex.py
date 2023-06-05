import numpy  
import scipy.special  

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
        self._probabilities: numpy.ndarray = numpy.array([0.010417, 0.03125, 0.125, 0.5, 0.25, 0.0625, 0.020833])  
        self._mu: float = (self._pattern_length / 2.0) + (((-1) ** (self._pattern_length + 1)) + 9.0) / 36.0 - ((self._pattern_length / 3.0) + (2.0 / 9.0)) / (2 ** self._pattern_length)  
        self._last_bits_size: int = -1  
        self._blocks_number: int = -1  
        super(LinearComplexityTest, self).__init__("Linear Complexity", 0.01)  

    def _execute(self, bits: numpy.ndarray) -> Result:  
        if self._last_bits_size == -1 or self._last_bits_size != bits.size:  
            blocks_number: int = int(bits.size // self._pattern_length)  
            self._last_bits_size = bits.size  
            self._blocks_number = blocks_number  
        else:  
            blocks_number: int = self._blocks_number  
        blocks_linear_complexity: numpy.ndarray = numpy.zeros(blocks_number, dtype=int)  
        for i in range(blocks_number):  
            blocks_linear_complexity[i] = self._berlekamp_massey(bits[(i * self._pattern_length):((i + 1) * self._pattern_length)])  
        tickets: numpy.ndarray = ((-1.0) ** self._pattern_length) * (blocks_linear_complexity[:] - self._mu) + (2.0 / 9.0)  
        frequencies: numpy.ndarray = numpy.zeros(self._freedom_degrees + 1, dtype=int)  
        for ticket in tickets:  
            frequencies[min(self._freedom_degrees, int(max(-2.5, ticket) + 2.5))] += 1  
        chi_square: float = float(numpy.sum(((frequencies[:] - (blocks_number * self._probabilities[:])) ** 2.0) / (blocks_number * self._probabilities[:])))  
        score: float = scipy.special.gammaincc((self._freedom_degrees / 2.0), (chi_square / 2.0))  
        if score >= self.significance_value:  
            return Result(self.name, True, numpy.array(score))  
        return Result(self.name, False, numpy.array(score))  

    def is_eligible(self, bits: numpy.ndarray) -> bool:  
        if bits.size < self._sequence_size_min:  
            return False  
        return True  

    @staticmethod  
    def _berlekamp_massey(sequence: numpy.ndarray) -> int:  
        b: numpy.ndarray = numpy.zeros(sequence.size, dtype=int)  
        c: numpy.ndarray = numpy.zeros(sequence.size, dtype=int)  
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
def read_file_to_array(file):
        data = file.read().decode()  # Dosyanın içeriğini oku ve bir stringe dönüştür
        array = [int(char) for char in data if char in '01']  # Dosyanın içeriğini bir bit dizisine dönüştür
        return array


