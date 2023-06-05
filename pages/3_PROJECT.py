import streamlit as st
import TEST
import numpy
import itertools
import matplotlib.pyplot as plt
import scipy.special as sps
from io import StringIO
import scipy.stats as stats
from sympy import isprime
import random
import pandas as pd
from scipy.special import gammaincc
from scipy.stats import chisquare
import time
st.set_page_config(page_title="PROJECT",
                   page_icon=":computer:",
                   layout="wide")

st.title("BBS ALGORİTMASI NIST TESTLERİ")
st.header("BBS Ekranında Oluşturulan Algoritma Aynı Anda 15 Adet NIST Testinden Geçecektir.")
st.markdown(""":green[Algoritma Aşağıdaki Gibidir.]
""")
code = '''import random
from sympy import isprime

def asal_bul():
    while True:
        p = random.getrandbits(512)
        if isprime(p) and p % 4 == 3:
            return p

p = asal_bul()
q = asal_bul()
M = p * q
x = random.randint(1,M)
bit_array = list()

for i in range(1000000):
    x = (x*x) % M 
    b = x % 2
    bit_array.append(b)

dosya = open("bit_array.txt","w")
dosya.write(''.join(map(str, bit_array)))
dosya.close()
'''
st.code(code, language='python')
st.subheader(":red[TEST SONUÇLARINI GÖRÜNTÜLEMEK İÇİN TIKLAYINIZ]")
button = st.button("TEST",key="primary")
st.write(button)
if button:
    def asal_bul():
        while True:
            p = random.getrandbits(512)
            if isprime(p) and p % 4 == 3:
                return p

    p = asal_bul()
    q = asal_bul()
    M = p * q
    x = random.randint(1,M)
    bit_array = list()

    for i in range(1000000):
        x = (x*x) % M 
        b = x % 2
        bit_array.append(b)
    st.subheader(":blue[Monobit Test Sonucu]")
    bbs = TEST.MonobitTest(bit_array)
    result = bbs.perform_test()
    st.markdown(result)
    st.subheader(":blue[Frequency Test within a Block Test Sonucu]")
    bbs2 = TEST.FrequencyTest(bit_array)
    result2 = bbs2.perform_test()
    st.markdown(result2)
    st.subheader(":blue[Runs Test Sonucu]")
    bbs3 = TEST.RunsTest(bit_array)
    result3 = bbs3.perform_test2()
    st.markdown(result3)
    st.subheader(":blue[Test for the Longest Run of Ones in a Block Test Sonucu]")
    bbs4 = TEST.LongestSeriesTest(bit_array)
    result4 = bbs4.perform_test3()
    st.markdown(result4)
    st.subheader(":blue[Binary Matrix Rank Test Sonucu]")
    bbs5 = TEST.BinaryRankTest(bit_array)
    result5 = bbs5.perform_test5()
    st.markdown(result5)
    st.subheader(":blue[Discrete Fourier Transform (Spectral) Test Sonucu]")
    bbs6 = TEST.DFTTest(bit_array)
    result6 = bbs6.perform_test6()
    st.markdown(result6)
    st.subheader(":blue[Non-overlapping Template Matching Test Sonucu]")
    bits = ''.join(map(str, bit_array))
    m = 10000  # you can change this to the desired block length
    b = '101'
    bbs7 = TEST.NonOverlappingTemplateMatching(bits,m,b)
    result7 = bbs7.perform_test7()
    st.markdown(result7)
    st.subheader(":blue[Overlapping Template Matching Test Sonucu]")
    template_tester = TEST.OverlappingTemplateMatchingTest(bit_array, m=4, M=16, N=8, a=[1, 0, 1])
    p_value = template_tester.perform_test()
    st.markdown(p_value)
    st.subheader(":blue[Maurer’s “Universal Statistical” Test Sonucu]")
    test = TEST.MaurersUniversalTest(bit_array, L=7, Q=1280)
    p_value = test.perform_test()
    st.markdown(p_value)
    st.subheader(":green[⚠️ Linear Complexity Test HATALI ⚠️]")
    st.subheader(":blue[Serial Test Sonucu]")
    m = 2
    serial_test = TEST.SerialTest(m=3)
    p_value = serial_test.perform_test(bit_array)
    st.markdown(p_value)
    st.subheader(":blue[Approximate Entropy Test Sonucu]")
    apen_fast = TEST.ApEnFast(m=2, r=0.2)
    result = apen_fast.perform_calculation(bit_array)
    st.markdown(result)
    st.subheader(":blue[Cumulative Sums (Cusum) Test Sonucu]")
    cusum_test = TEST.CusumTest()
    result = cusum_test.perform_test(bit_array)
    st.markdown(result)
    st.subheader(":blue[Random Excursions Test Sonucu]")
    random_excursions_test = TEST.RandomExcursionsTest()
    result = random_excursions_test.perform_test(bit_array)
    st.markdown(f"P Değerleri : {result}")
    st.subheader(":blue[Random Excursions Variant Test Sonucu]")
    random_excursions_variant_test = TEST.RandomExcursionsVariantTest()
    result = random_excursions_variant_test.perform_test(bit_array)
    st.markdown(f"P Değerleri : {result}")
    col1,col2,col3 = st.columns([1,2,1])
    progress_bar = col2.progress(0)
    for perc_completed in range(100):
        time.sleep(.1)
        progress_bar.progress(perc_completed + 1)
    col2.success("Bütün Testler Tamamlandı")
st.subheader("Oluşturulacak Array Sayısını Seçiniz")
array_custom = st.slider(':green[Array]', 0, 100, 1)
st.write(f"Seçilen Array Sayısı : {array_custom}")
if array_custom == 100:
    code="""import random

def asal_bul():
    while True:
        p = random.getrandbits(512)
        if isprime(p) and p % 4 == 3:
            return p

bit_arrays = []
for _ in range(100):
    p = asal_bul()
    q = asal_bul()
    M = p * q
    x = random.randint(1, M)
    bit_array = []
    
    for i in range(1000000):
        x = (x * x) % M
        b = x % 2
        bit_array.append(b)
    
    bit_arrays.append(bit_array)

# 100 farklı liste oluşturuldu ve bit_arrays listesine eklendi.

# Test amaçlı, oluşturulan listelerin uzunluğunu kontrol edelim
for i in range(100):
    print(f"Liste {i+1} uzunluğu: {len(bit_arrays[i])}")

"""
    st.code(code, language='python')

st.sidebar.success("Gitmek İstediğiniz Sayfayı Seçiniz")