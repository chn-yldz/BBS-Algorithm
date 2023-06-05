import streamlit as st
from PIL import Image # Image Processing
from scipy.stats import chisquare
import itertools
import matplotlib.pyplot as plt
import scipy.special as sps
import numpy as np
import pandas as pd
import monobite as mono # Test Module Import
import FrequencyWithinBlock as freq # Test Module Import
from io import StringIO
import scipy.stats as stats
from sympy import isprime
import random
import runs as rns
import longest_run as long
import binary_matrix as bnmatrix
import discrete_fourier as disc
import non_overlapping as non
import overlapping as over
import maurer as mau
import linear_complex as lin
import serial_test as seri
import approximate as approx
import cumulative as sums
import random_ex_test as randex
import random_ex_variant as randvar
st.set_page_config(page_title="NIST TESTS",
                   page_icon="📋",
                   layout="wide")

st.title("NIST (National Institute of Standards and Technology)")
st.header("Hakkında:")
st.markdown(""" NIST (National Institute of Standards and Technology) Testleri, rastgelelik testleri olarak da bilinir ve bir bit dizisinin rastgele olup olmadığını kontrol etmek için kullanılır. Bu testler, kriptografik sistemlerde rastgele sayı üreticilerinin doğru çalışıp çalışmadığını doğrulamak için önemlidir.
NIST Testi'nin 15 farklı testi vardır:
1) `Frequency (Monobit) Test`
2) `Frequency Test within a Block`
3) `Runs Test`
4) `Test for the Longest Run of Ones in a Block`
5) `Binary Matrix Rank Test`
6) `Discrete Fourier Transform (Spectral) Test`
7) `Non-overlapping Template Matching Test`
8) `Overlapping Template Matching Test`
9) `Maurer’s “Universal Statistical” Test`
10) `Linear Complexity Test`
11) `Serial Test`
12) `Approximate Entropy Test`
13) `Cumulative Sums (Cusum) Test`
14) `Random Excursions Test`
15) `Random Excursions Variant Test`
 """)
st.markdown("Bu testlerin her biri, farklı türden rastgelelik kusurlarını tespit etmek için tasarlanmıştır.")
st.subheader(":green[Frequency (Monobit) Test]")
st.markdown("Frequency (Monobit) Test, bir bit dizisi içindeki 1'lerin ve 0'ların sayısının eşit olup olmadığını kontrol eder. Eğer bit dizisi tamamen rastgele oluşturulmuşsa, 0'ların ve 1'lerin sayısı hemen hemen eşit olmalıdır.")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = '''def monobit_test(bit_sequence:list):
    count = 0
    for bit in bit_sequence:
        if bit == '1':
            count += 1
    n = len(bit_sequence)
    S_obs = abs(count - n/2)
    if n == 0:
        S = 0
    else:
        S = S_obs / (n**0.5)
    return S < 2.16


'''
st.code(code, language='python')
st.markdown(":green[Kodun Çalışma Mantığı Aşağıdaki Gibidir:]")
st.markdown("""
Fonksiyon, monobit_test(bit_sequence) adında tanımlanır. Bu fonksiyon, bir parametre alır: bit_sequence. Bu, test edilmek üzere olan bit dizisidir.

1) `count = 0: Kod burada count adlı bir sayaç oluşturur. Bu sayaç, bit dizisindeki 1'lerin sayısını tutar.`

2) `for bit in bit_sequence:: Bu satır, bit dizisindeki her bir bit üzerinde bir döngü başlatır.`

3) `if bit == '1': count += 1: Eğer mevcut bit 1 ise, count değerini 1 artırır. Bu satır, bit dizisindeki 1'lerin sayısını sayar.`

4) `n = len(bit_sequence): n, bit dizisindeki toplam bit sayısını temsil eder.`

5) `S_obs = abs(count - n/2): S_obs değeri hesaplanır. Bu değer, bit dizisindeki 1'lerin ve 0'ların sayısı arasındaki farkı temsil eder. abs fonksiyonu, bu farkın mutlak değerini alır.`

6) `S = S_obs / (n**0.5): S değeri hesaplanır. Bu değer, S_obs değerini bit dizisinin uzunluğunun kareköküne böler. Bu değer, bit dizisindeki 1'lerin ve 0'ların sayısının rastgeleliğini ölçer.`

7) `return S < 2.16: Bu satır, testin sonucunu döndürür. Eğer S değeri 2.16'dan küçükse, bit dizisi %1 güven seviyesiyle rastgele kabul edilir ve fonksiyon True döndürür. Eğer S değeri 2.16'dan büyükse, bit dizisi %1 güven seviyesiyle rastgele olmayabilir ve fonksiyon False döndürür.`

 """)
st.subheader(":red[Kodu Test Etmek için Bit Dizisi (Array) Yükleyin :]")
uploaded_file = st.file_uploader("Dosya Seçiniz",key="file-uploader")
if uploaded_file is not None:
    bit_arr = mono.read_file_to_array(uploaded_file)
    result = mono.monobit_test(bit_arr)
    if result is True:
        st.subheader("Test Başarılı✔️")
    elif result is False:
        st.subheader("Test Başarısız❌")
    else:
        st.subheader(f"Test sonucu beklenmedik bir değer döndürdü: {result}")
else:
    st.info("Dosya Yüklenmedi.")
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="primary")
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
    result = mono.monobit_test(bit_array)
    if result is True:
        st.subheader("Test Başarılı✔️")
    elif result is False:
        st.subheader("Test Başarısız❌")
    else:
        st.subheader(f"Test sonucu beklenmedik bir değer döndürdü: {result}")
st.subheader(":green[Frequency Test within a Block]")
st.markdown("""(Blok İçinde Frekans Testi), bir bit dizisi üzerinde çalışır ve belirli bir blok boyutu içinde 0'ların ve 1'lerin sayısının eşit olup olmadığını kontrol eder. Bu, bit dizisinin rastgeleliğini kontrol etmek için kullanılır. Eğer bit dizisi tamamen rastgele oluşturulmuşsa, her bloktaki 0'ların ve 1'lerin sayısı hemen hemen eşit olmalıdır.
 """)
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = '''def frequency_test_within_a_block(bit_sequence, block_size):
    n = len(bit_sequence)
    num_blocks = n // block_size
    blocks = [bit_sequence[i * block_size:(i + 1) * block_size] for i in range(num_blocks)]
    pi = [block.count('1') / block_size for block in blocks]
    chi_squared = 4 * block_size * sum((p - 0.5)**2 for p in pi)
    p_value = 1 - stats.chi2.cdf(chi_squared, df=num_blocks)
    return p_value > 0.01
'''
st.code(code, language='python')
st.markdown(":green[Kodun Çalışma Mantığı Aşağıdaki Gibidir:]")
st.markdown("""Bu kod parçası, belirtilen blok boyutuyla bit dizisini bloklara ayırır ve ardından her bloktaki `1`'lerin sayısını hesaplar. Daha sonra, her blok için `(pi - 0.5)^2` değerini hesaplar ve bu değerlerin toplamını alır. Bu toplam, bir `Chi-kare` dağılımına göre normalize edilir ve sonuçta elde edilen p-değerine dayanarak bit dizisinin rastgeleliği hakkında bir sonuca varılır. Eğer p-değer 0.01'den büyükse, bit dizisi `%1` güven seviyesiyle rastgele kabul edilir ve fonksiyon True döndürür. Eğer p-değer `0.01`'den küçükse, bit dizisi `%1` güven seviyesiyle rastgele olmayabilir ve fonksiyon False döndürür.

Bu kod parçası, scipy.stats modülünün chi2.cdf fonksiyonunu kullanır, bu yüzden bu modülün yüklü olması gerekir. Ayrıca, blok boyutu olarak önerilen değer genellikle `20`'dir. Blok boyutunun çok büyük olmaması ve bit dizisinin uzunluğuna uygun olması önemlidir.



""")
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="secondary")
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
    result = freq.frequency(bit_array)
    if result is True:
        st.subheader("Test Başarılı✔️")
    elif result is False:
        st.subheader("Test Başarısız❌")
    else:
        st.subheader(f"Test sonucu beklenmedik bir değer döndürdü: {result} ❎")

st.subheader(":green[Runs Test]")
st.markdown(""""Runs Test, bir bit dizisindeki ardışık `0`'ların ve `1`'lerin sayısını inceler. Bu test, bit dizisinin bağımsız
olup olmadığını kontrol etmek için kullanılır. Eğer bitler bağımsız ise, ardışık aynı bitlerin `("runs")`
uzunluğunun dağılımı belirli bir şekle sahip olmalıdır.
Runs Testini uygulamak için, bit dizisindeki `"runs"`ları saymamız ve bu sayının beklenen
değerden ne kadar farklı olduğunu kontrol etmemiz gerekiyor. Ardından, bu farkın beklenen
varyansla ne kadar uyumlu olduğunu bir `z-testi` ile kontrol ederiz.
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = '''def runs_test(bit_array):
    n = len(bit_array)
    s = sum(bit_array)
    pi = s / n # 1'lerin oranı
    tau = 2 / np.sqrt(n)

    if abs(pi - 0.5) >= tau:
        return False
    else:
        runs = 1 + sum(bit_array[i] != bit_array[i+1] for i in range(n-1))
        R_obs = abs(runs - 2*n*pi(1-pi)) / np.sqrt(2*n*pi(1-pi))
        p_value = 2 * stats.norm.cdf(-R_obs)
        if p_value < 0.01:
            return False
        else:
            return True
'''
st.code(code, language='python')
st.markdown(":green[Kodun Çalışma Mantığı Aşağıdaki Gibidir:]")
st.markdown("""1) `İlk olarak, verilen bit dizisinin uzunluğu n hesaplanır ve dizideki 1'lerin toplamı s olarak belirlenir.`

2) `Daha sonra, 1'lerin oranı (pi) hesaplanır. Bunun için s bit dizisinin toplam uzunluğuna (n) bölünür. Bu, bit dizisindeki 1'lerin oranını temsil eder.`

3) `tau değeri hesaplanır. Bu değer, 1'lerin oranındaki standart sapmanın bir ölçüsüdür ve 2 / sqrt(n) formülüyle elde edilir.`

4) `Eğer abs(pi - 0.5) değeri tau'dan büyük veya eşitse, yani 1'lerin oranı 0.5'den önemli ölçüde farklıysa, False döndürülür. Bu durum, bit dizisinin rastgele dağılmadığını gösterir.`

5) `Eğer 1'lerin oranı 0.5'e yakınsa, kod runs testini uygular. Bu aşamada, ardışık bitlerin grupları (runs) sayılır. Bir run, bit dizisinde ardışık olarak aynı bitlerin olduğu bir grup olarak tanımlanır.`

6) `İki ardışık bitin farklı olduğu durumları bulmak için bit_array[i] != bit_array[i+1] ifadesi kullanılarak grupların sayısı hesaplanır. Bu, runs testinde kullanılan bir adımdır.`

7) `Toplam grup sayısına 1 eklenir ve bu değer runs olarak adlandırılır.`

8) `Runs gözlem değeri (R_obs) hesaplanır. Bu değer, runs değeri ile bir formül aracılığıyla hesaplanır. Formülde, runs değerinin 2npi(1-pi)'den çıkarılması ve 2npi(1-pi) ile bölünmesi işlemi yapılır. Bu adım, bit dizisinin dağılımının test edilmesi için kullanılan bir hesaplama yöntemidir.`

9) `p_value hesaplanır. Bu değer, R_obs'un olasılık dağılımı üzerindeki p-değerini temsil eder. Bu değer, 2 * stats.norm.cdf(-R_obs) formülüyle hesaplanır. cdf fonksiyonu, verilen değerin kümülatif dağılım fonksiyonunu hesaplar.`

10) `Eğer p_value 0.01'den küçükse, yani runs testinin sonucu istatistiksel olarak anlamlıysa, False döndürülür. Bu durum, bit dizisinin rastgele dağılmadığını gösterir.`

11) `Eğer p_value 0.01'den büyük veya eşitse, runs testinin sonucu istatistiksel olarak anlamsızdır ve True döndürülür. Bu durum, bit dizisinin rastgele dağıldığını gösterir.`
""")
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="third")
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
    result = rns.runs_test(bit_array)
    if result is True:
        st.subheader("Test Başarılı✔️")
    elif result is False:
        st.subheader("Test Başarısız❌")
    else:
        st.subheader(f"Test sonucu beklenmedik bir değer döndürdü: {result} ❎")

st.subheader(":green[Test for the Longest Run of Ones in a Block]")
st.markdown("""`"Test for the Longest Run of Ones in a Block"` testi, bir bit dizisindeki en uzun `1`'lerin serisini inceler.
Bu test, bit dizisinin rastgele olup olmadığını kontrol etmek için kullanılır. Eğer bitler rastgele ise,
bir blok içindeki en uzun `1`'lerin serisi belirli bir uzunlukta olmamalıdır.
Bu testi uygulamak için, bit dizisini belirli bir blok boyutuna bölmemiz ve her blokta en uzun
`1`'lerin serisini saymamız gerekiyor. Ardından, bu sayının belirli bir dağılıma uyup uymadığını
kontrol ederiz.
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = '''def longest(bit_array):
    M = 10 # Blok Boyutu
    blocks = [bit_array[i:i+M] for i in range(0,len(bit_array),M)]

    longest_series = [max(len(list(g)) for k,g in itertools.groupby(block) if k == 1) for block in blocks ]
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

'''
st.code(code, language='python')
st.markdown(":green[Kodun Çalışma Mantığı Aşağıdaki Gibidir:]")
st.markdown("""
Bu kod parçası, NIST'in Öklidyen ve Non-Öklidyen Rastgelelik Testleri arasında yer alan "Longest Run of Ones in a Block" testini gerçekleştirir. Bu test, bir dizi bitin rastgele olup olmadığını kontrol etmek için kullanılır.
1) `bit_array adı verilen bir bit dizisi alır.`

2) `M boyutunda bloklar oluşturur. Bu durumda, her blok 10 bit içerir.`

3) `Her blokta, en uzun '1' serisi bulunur. Bu, Python'un itertools.groupby fonksiyonu kullanılarak yapılır. Bu fonksiyon, bitlerin gruplandırılmasını sağlar ve ardından her bir '1' grubunun uzunluğunu belirler.`

4) `K adı verilen bir indeks belirlenir. Bu örnekte K = 6.`

5) `pi olarak adlandırılan bir dizi belirlenir. Bu dizi, belirli bir uzunluktaki '1' serilerinin olasılığını temsil eder.`

6) `nu adı verilen bir dizi oluşturulur. Bu dizi, her bir uzunlukta kaç tane '1' serisi olduğunu belirler.`

7) `N adı verilen bir değişken oluşturulur. Bu değişken, blok sayısını temsil eder.`

8) `X_obs adı verilen bir değişken oluşturulur. Bu değişken, gözlenen ve beklenen '1' serileri arasındaki farkın karelerinin toplamını temsil eder.`

9) `p_value adı verilen bir değişken oluşturulur. Bu değişken, X_obs'un olasılık değerini belirler. Bu değer, scipy.special.gammaincc fonksiyonu kullanılarak hesaplanır. Bu fonksiyon, tamamlanmış gamma fonksiyonunun eşdeğerini hesaplar.`

10) `p_value 0.01'den küçükse, bit dizisinin rastgele olmadığı sonucuna varılır ve False döndürülür. Aksi takdirde, bit dizisinin rastgele olduğu sonucuna varılır ve True döndürülür.`
""")
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="4th")
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
    result = long.longest(bit_array)
    if result is True:
        st.subheader("Test Başarılı✔️")
    elif result is False:
        st.subheader("Test Başarısız❌")
    else:
        st.subheader(f"Test sonucu beklenmedik bir değer döndürdü: {result} ❎")

st.subheader(":green[Binary Matrix Rank Test]")
st.markdown("""Binary Matrix Rank Test, bit dizisindeki hatların ve sütunların bağımsızlığını ölçer. Bu test, bit
dizisinin rastgele olup olmadığını belirlemek için kullanılır. Eğer bitler rastgele ise, herhangi bir
düzenlilik olmamalıdır ve bu, hat ve sütunların bağımsız olması anlamına gelir.
Bu testi uygulamak için, bit dizisini belirli bir boyutta bir matrise dönüştürmeli ve ardından
matrisin rankını `(bağımsız satır ve sütunların sayısı)` hesaplamalıyız. Ardından, beklenen rank
dağılımı ile elde ettiğimiz rank dağılımını karşılaştırırız 
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = '''def binaryrank(bit_array):
    # Matris Boyutu
    M = 3
    Q = 3
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
'''
st.code(code, language='python')
st.markdown(":green[Kodun Çalışma Mantığı Aşağıdaki Gibidir:]")
st.markdown("""Bu kod, bir bit dizisinin rastgeleliğini test etmek için kullanılan bir Binary Matrix Rank testini uyguluyor. 
1) `Kod, bit dizisi (bit_array) alıyor. Bu dizinin her bir elementi 0 veya 1 olabilir.`

2) `Dizi, M x Q boyutunda olan matrislere bölünüyor. Bu örnekte M ve Q her ikisi de 3. Bu matrisler, np.reshape fonksiyonu kullanılarak oluşturuluyor. Her matris için, bit dizisinin ardışık M*Q öğesi alınıyor. Böylece, bit dizisi boyunca M*Q öğe boyutunda bir dizi matris oluşturuluyor.`

3) `Her bir matrisin rankı (lineer bağımsız satır veya sütunlarının maksimum sayısı), np.linalg.matrix_rank fonksiyonu kullanılarak hesaplanıyor.`

4) `Hesaplanan her rankın frekansı (rank1_sayisi, rank2_sayisi, rank3_sayisi) bulunuyor.`

5) `Beklenen rank değerleri (pi) verilmiştir. Bu değerler, tamamen rastgele bir bit dizisinde her rankın teorik olarak beklenen frekanslarıdır. Bu örnekte, pi[0] rank 1'in beklenen frekansını, pi[1] rank 2'nin beklenen frekansını ve pi[2] rank 3'ün beklenen frekansını temsil eder.`

6) `Gözlenen frekanslar ve beklenen frekanslar arasındaki farklar, bir χ² (chi-square) test istatistiği hesaplamak için kullanılır. Bu istatistik, X_obs olarak adlandırılır ve beklenen ve gözlenen değerler arasındaki farkın toplamının bir ölçüsüdür.`

7) `Bu χ² değeri kullanılarak, bit dizisinin rastgele olup olmadığını belirlemek için bir p-değer hesaplanır. chi2.sf(X_obs,2) fonksiyonu, X_obs üzerindeki chi-kare dağılımının tamamlayıcı kümülatif dağılım fonksiyonunu (survival function) hesaplar. Burada, 2 derece serbestlik vardır çünkü 3 olası rank vardır ve bu nedenle 3 - 1 = 2 derece serbestliğimiz vardır.`

8) `Hesaplanan p-değer, 0.01'den küçükse, bit dizisi rastgele değildir ve fonksiyon False döndürür. Aksi takdirde, bit dizisi rastgele kabul edilir ve fonksiyon True döndürür.`

""")
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="5th")
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
    result = bnmatrix.binaryrank(bit_array)
    if result is True:
        st.subheader("Test Başarılı✔️")
    elif result is False:
        st.subheader("Test Başarısız❌")
    else:
        st.subheader(f"Test sonucu beklenmedik bir değer döndürdü: {result} ❎")

st.subheader(":green[Discrete Fourier Transform (Spectral) Test]")
st.markdown("""`Ayrık Fourier Dönüşümü (Discrete Fourier Transform - DFT)`, sürekli sinyallerin dönüşümüne yönelik `Fourier Dönüşümünün` ayrık sinyaller için karşılığıdır. `DFT`, bir ayrık zamanlı sinyali, frekans bileşenlerine ayırır. Yani temelde, bir zaman serisi verisinin farklı frekanslardaki bileşenlerine bakmamızı sağlar.
`Spectral test (spektral test)`, genellikle rasgele sayı üreteçlerinin kalitesini değerlendirmek için kullanılan bir tekniktir. `DFT`, bu testlerin bir parçası olabilir. `DFT`'nin sonuçları, rasgele sayı dizisinin belirli özelliklere `(örneğin, periyodiklik veya tekrarlanabilirlik)` sahip olup olmadığını belirlemek için kullanılabilir. 
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = '''import numpy as np
import matplotlib.pyplot as plt

# Rasgele bir sinyal oluştur
np.random.seed(0)  # Rastgelelik için bir tohum ayarla
time_step = 0.02
period = 5.
time_vec = np.arange(0, 20, time_step)
sig = (np.sin(2 * np.pi / period * time_vec)
       + 0.5 * np.random.randn(time_vec.size))

# Sinyalin DFT'sini hesapla
freqs = np.fft.fftfreq(sig.size, d=time_step)
idx = np.argsort(freqs)
ps = np.abs(np.fft.fft(sig))**2

# Spektrumu çizdir
plt.figure()
plt.plot(freqs[idx], ps[idx])
plt.title("Power spectral density")
plt.xlabel("Frequency")
plt.ylabel("Power")
plt.show()

'''
st.code(code, language='python')
st.markdown("""Bu kod, bir dizi rasgele sayı oluşturur ve bu sayıları bir sinyal olarak kullanır. Daha sonra, bu sinyalin `DFT`'sini hesaplar ve sinyalin güç spektral yoğunluğunu çizer. Bu, frekansların sinyale katkılarını gösterir.

Bir spectral test için, bu çıktıyı rasgele sayılar dizinizin beklenen özelliklerine karşı karşılaştırabilirsiniz. `Örneğin, bir "iyi" rasgele sayı dizisi genellikle düz bir spektral yoğunluk eğrisine sahip olacaktır, çünkü tüm frekanslar eşit olarak temsil edilir.` Öte yandan, periyodik veya tekrar eden özellikler, spektral yoğunluk eğrisinde belirgin zirvelere neden olacaktır. Bu, rasgele sayı dizisinin belirli bir özelliği "tutmakta" olduğunu ve dolayısıyla gerçekten rasgele olmadığını gösterebilir.
""")
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="6th")
st.write(button)
if button:
    def BBS(p, q, s, n):
        m = p * q
        x = (s*s) % m
        bits = []
        for _ in range(n):
            x = (x*x) % m
            bits.append(x % 2)
        return bits
    p = 3
    q = 7
    s = 2
    n = 1024
    bits = BBS(p, q, s, n)
    magnitude = disc.DFT_test(bits)
    st.subheader(f"Test Sonucu :👉 {magnitude}")
st.subheader(":green[Non-overlapping Template Matching Test]")
st.markdown("""`Non-overlapping Template Matching Test`, bir rastgele bit dizisi üzerinde belirli bir bit şablonunun beklenen frekansın dışında olup olmadığını kontrol eder. Bu test, rastgeleliği değerlendirmek için kullanılır. Eğer belirli bir şablon beklenenden daha sık veya daha seyrek görülürse, bu bit dizisinin tamamen rastgele olmadığına dair bir gösterge olabilir.

Bu test genellikle bir bit dizisi ve bir bit şablonu (veya bir dizi şablon) alır ve bit dizisini şablondaki bitlerin sayısı kadar bölümlere ayırır. Her bölümde, bit şablonunun kaç kez göründüğünü sayar. Daha sonra, bu sayıların beklenen dağılıma uygun olup olmadığını bir istatistiksel test (genellikle bir `χ²-testi`) kullanarak kontrol eder.
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = '''import numpy as np
from scipy.stats import chisquare

def non_overlapping_template_matching_test(bit_string, template):
    # Bit stringi, şablon uzunluğuna göre bölümlere ayırma
    template_length = len(template)
    chunks = [bit_string[i:i+template_length] for i in range(0, len(bit_string), template_length)]
    
    # Her bölümde, şablonun kaç kez göründüğünü sayma
    counts = [chunk.count(template) for chunk in chunks]

    # Beklenen frekansın hesaplanması (rastgelelik durumunda)
    expected_frequency = len(bit_string) / (2 ** template_length)

    # Chi-kare testini uygulama
    chi_square, p_value = chisquare(counts, expected_frequency)
    
    return chi_square, p_value

# Bit dizisi ve bit şablonu
bit_string = '1010010101001010100110101010001010010101001010101101010101010101010010101'
template = '101'

# Testi çalıştırma
chi_square, p_value = non_overlapping_template_matching_test(bit_string, template)

# Sonuçları yazdırma
print(f'Chi-square: {chi_square}')
print(f'P-value: {p_value}')

'''
st.code(code, language='python')
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="7th")
st.write(button)
if button:
    def asal_bul():
        while True:
            p = random.getrandbits(512)
            if isprime(p) and p % 4 == 3:
                return p
    def non_overlapping_template_matching(bits, m, b):
        n = len(bits)
        M = m
        N = n//M
        blocks = [bits[i*M:(i+1)*M] for i in range(N)]
        counts = [block.count(b) for block in blocks]

        #Statistic
        mu = (M - len(b) + 1)/2**len(b)
        var = M*(1/2**len(b) - (2*len(b) - 1)/2**(2*len(b)))
        chi2 = sum([(c - mu)**2/var for c in counts])

        # compute the p-value
        p_value = sps.gammaincc(N/2, chi2/2)

        return p_value
    p = asal_bul()
    q = asal_bul()
    M = p * q
    x = random.randint(1, M)
    bit_array = list()
    for i in range(1000000):
        x = (x*x) % M 
        b = x % 2
        bit_array.append(b)
    bits = ''.join(map(str, bit_array))
    m = 10000  # you can change this to the desired block length
    b = '101'  # you can change this to the desired template
    p_value = non_overlapping_template_matching(bits, m, b)
    st.subheader(f"Test Sonucu :👉 {p_value} ")
st.subheader(":green[Overlapping Template Matching Test]")
st.markdown("""`Overlapping Template Matching Test`, bir rastgelelik testidir ve genellikle bir dizi verinin rastgele olup olmadığını belirlemek için kullanılır. Bu test, bir veri dizisinde belirli bir şablonun (bir bit dizesi) beklenen sayıda bulunup bulunmadığını kontrol eder.

Bu test, rastgele veri üreten bir sistemin, belirli bir şablonu beklenenden daha fazla veya daha az üretip üretmediğini belirlemek için kullanılır. Veri dizisi, belirli bir şablonun çok fazla veya çok az oluşumu, verinin rastgele olmadığını gösteren bir işarettir.
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = '''import numpy as np
from scipy.special import gammaincc

def overlapping_template_matching_test(binary_data, m=9, M=1032, N=968, a=[0,0,0,0,0,0,0,0,1]):
    """
    Applies the Overlapping Template Matching Test given by NIST.
    :param binary_data: a binary string
    :param m: length of the template (default value is 9)
    :param M: length in bits of the substring to be tested (default value is 1032)
    :param N: number of blocks (default value is 968)
    :param a: the m-bit template as a list (default value is [0,0,0,0,0,0,0,0,1])
    :return: the p-value from the test
    """
    n = len(binary_data)
    
    blocks = [binary_data[i * M : (i + 1) * M] for i in range(N)]
    counts = [block.count(''.join(map(str, a))) for block in blocks]
    
    mu = (M - m + 1) / 2 ** m
    var = M * (1 / 2 ** m - (2 * m - 1) / 2 ** (2 * m))
    chi_squared = sum([(c - mu) ** 2 / var for c in counts])
    p_value = gammaincc(N / 2, chi_squared / 2)
    return p_value

'''
st.code(code, language='python')
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="8")
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
    result = over.overlapping_template_matching_test(bit_array)
    st.subheader(f"Test Sonucu :👉 {result} ")


st.subheader(":green[Maurer’s “Universal Statistical” Test]")
st.markdown("""`Maurer’s “Universal Statistical” Test`, bir bit dizisinin rastgeleliğini kontrol etmek için kullanılan bir istatistiksel testtir. Bu test, bir bit dizisinin düzeni ve yapıları arasındaki ilişkileri kontrol eder. Maurer's Universal Statistical Test, bir bit dizisinin rastgele olduğunu kontrol etmek için kullanılan bir dizi testten biridir. Bu testler genellikle kriptografi ve güvenlik uygulamalarında kullanılır.

Bu testin temel çalışma prensibi, belirli bir boyutta (genellikle 1 ile 16 bit arasında) bloklar halinde bit dizileri oluşturmak ve bu blokların bit dizisi boyunca ne kadar sık tekrarlandığını kontrol etmektir. Bloklar arasındaki ortalama boşluk, rastgele bir bit dizisinde beklenen ortalama boşluktan önemli ölçüde farklıysa, bit dizisi rastgele olarak kabul edilmez.

Testin detayları şu şekildedir:
1) `Öncelikle bit dizisi, L boyutunda bloklara ayrılır. Blokların sayısı Q olacak şekilde seçilir.`

2) `Her bir L-bit blok için, bloğun bit dizisindeki önceki konumuna ve mevcut konumuna bakılır. Bu iki konum arasındaki fark, bloğun "boşluk" olarak adlandırılan değerini belirler.`

3) `Bu boşluklar üzerinden bir toplam hesaplanır ve bu toplam, bit dizisi boyunca beklenen ortalama boşluktan önemli ölçüde farklıysa, bit dizisi rastgele olarak kabul edilmez.`

Maurer's Universal Statistical Test, bir bit dizisinin rastgeleliğini doğrulamak için güçlü bir araçtır. Bu test, özellikle kriptografik anahtarların ve diğer güvenlikle ilgili rastgele değerlerin kalitesini doğrulamak için kullanılır. Bu testin en önemli avantajlarından biri, bit dizisinin rastgeleliğini doğrularken çok az bellek kullanmasıdır, bu da onu büyük bit dizileri üzerinde çalışmak için ideal hale getirir.
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = '''def maurers_universal_test(bitstring, L=7, Q=1280):
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
'''
st.code(code, language='python')
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="9")
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
    result = mau.maurers_universal_test(bit_array)
    st.subheader("Test Başarılı✔️")
    st.subheader(f"Test Sonucu : {result}")
st.subheader(":green[Linear Complexity Test]")
st.markdown("""`Linear Complexity Test, bir bit dizisinin lineer karmaşıklığını ölçer, yani dizinin en kısa lineer geri bildirimli kaydırma kaydının (LFSR) boyutunu belirler. Bu test, bir bit dizisinin rastgeleliğini kontrol etmek için kullanılır. Eğer bir bit dizisi tamamen rastgele ise, lineer karmaşıklığı dizinin boyutuna yakın olmalıdır. Eğer lineer karmaşıklığı düşük ise, bu dizinin belirli bir yapıya sahip olduğunu ve tamamen rastgele olmadığını gösterir.
Test aşağıdaki şekilde işler:

1) `Bit dizisi belirli bir boyutta bloklara ayrılır. Her blok için, bloğun lineer karmaşıklığı hesaplanır.`
2) `Her bloğun lineer karmaşıklığı üzerinden bir test istatistiği hesaplanır.`
3) `Bu test istatistiği, bir bit dizisi rastgele olduğunda beklenen değerle karşılaştırılır. Eğer test istatistiği beklenen değerden önemli ölçüde farklıysa, bit dizisi rastgele olarak kabul edilmez.`
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = """def linear_complexity_test(sequence, M=500, K=6):
    n = len(sequence)
    N = n // M
    if N > 1:
        blocks = np.reshape(sequence[:N*M], (N, M))
        complexities = np.apply_along_axis(berlekamp_massey_algorithm, 1, blocks)
        T = (M//2 + (9 + (-1)**M)/36 - complexities)/M
        counts = np.histogram(T, bins=np.arange(-2.5, K-1.5))[0]
        π = [0.01047, 0.03125, 0.125, 0.5, 0.25, 0.0625, 0.02017]
        chisq = np.sum((counts - N*np.array(π))**2 / (N*np.array(π)))
        p_value = gammaincc((K-1)/2, chisq/2)
    else:
        p_value = None
    return p_value
"""
st.code(code, language='python')
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="10")
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
    test = lin.LinearComplexityTest()
    bit_array2 = np.array(bit_array)
    result = test._execute(bit_array2)
    st.subheader("Test Başarılı✔️")
    st.subheader(f"Test Sonucu : {result}")
st.subheader(":green[Serial Test]")
st.markdown("""`Serial Test`, rastgele sayı üreticilerinin çıktılarını test etmek için kullanılan bir testtir. Bu test, bit dizilerinin frekanslarını inceler ve bu frekansların eşit olmasını bekler. Eğer frekanslar eşit değilse, bu, üretilen bit dizisinin rastgele olmadığını gösterebilir.

Daha teknik bir açıklama yapacak olursak, Serial Test, bir rastgele bit dizisinde belirli bir uzunluktaki (genellikle 1, 2 veya 3 bitlik) bit dizilerinin frekansını hesaplar. Örneğin, 2-bitlik diziler için, 00, 01, 10 ve 11 olasılıklarının hemen hemen eşit olması beklenir.

Python'da Serial Test'i uygulamak için, bit dizinizi belirli bir uzunluktaki tüm bit dizilerinin frekanslarını hesaplayacak şekilde parçalamanız gerekecek. Ardından, beklenen ve gözlenen frekansları karşılaştırarak bir chi-square testi uygulayabilirsiniz.
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = """import numpy as np
import scipy.special

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

# Generate your bit string using BBS algorithm
bits = # Your BBS generated bit list here

# Convert your list to string to apply count function in the test
bit_str = ''.join(map(str, bits))

# Run the test
p_value = serial_test(bit_str, m=2)

# Print the result
print('P-value:', p_value)

"""
st.code(code, language='python')
st.markdown("""Bu kod, 2-bitlik dizilerin frekansını hesaplar ve `chi-square` testi uygular. m parametresi ile test edilecek bit dizisinin uzunluğunu belirleyebilirsiniz.

Unutmayın ki, elde ettiğiniz p-değerinin anlamlı olabilmesi için, bit dizinizin yeterince uzun olması gerekmektedir. Genellikle, bit dizisi uzunluğu, test edilen bit dizisi uzunluğunun (yani m'in) katları olmalıdır.
""")
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="11")
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
    result = seri.serial_test(bit_array,m=2)
    st.subheader("Test Başarılı✔️")
    st.subheader(f"Test Sonucu : {result}")
st.subheader(":green[Approximate Entropy Test]")
st.markdown("""`Approximate Entropy (ApEn)` testi, bir veri setinin düzenliliğini veya tahmin edilebilirliğini ölçmek için kullanılan bir yöntemdir. İstatistiksel bir ölçüm olarak, bir veri setinin düzenli olup olmadığını belirlemeye yardımcı olur. Başka bir deyişle, ApEn, bir veri setinin ne kadar karmaşık veya düzensiz olduğunu ölçer.

ApEn, bir dizi değer üzerinde hesaplanır ve genellikle zaman serisi verileri üzerinde kullanılır. Örneğin, kalp atış hızı, EEG sinyalleri veya finansal veriler gibi veriler üzerinde kullanılabilir.

ApEn'nin hesaplanması, veri setindeki benzer veri desenlerinin sayısını ölçer. İki parametre ile tanımlanır: m ve r. m, benzer desenlerin uzunluğunu belirler ve genellikle 2 veya 3 olarak belirlenir. r, iki desenin benzer olarak kabul edilmesi için gerekli maksimum farkı belirler.

ApEn'nin matematiksel formülü karmaşıktır, ancak genel olarak iki adımda hesaplanır:
1) `Veri setindeki tüm m-uzunluklu desenlerin sayısınının hesaplanması.`
2) `Her desen için, diğer tüm desenlerle olan farkının hesaplanması ve r'den daha küçük olanların sayılması
Bu iki adımdan sonra, ApEn değeri, iki adımda bulunan sayıların logaritmik oranına eşittir. Düşük ApEn değerleri, veri setinin daha düzenli olduğunu gösterir, yani desenlerin tekrarlanma olasılığı daha yüksektir. Yüksek ApEn değerleri, veri setinin daha düzensiz veya karmaşık olduğunu gösterir, yani desenlerin tekrarlanma olasılığı daha düşüktür.`

ApEn, karmaşıklığı veya düzenliliği ölçmek için kullanılan diğer yöntemlere göre bazı avantajlara sahiptir. Özellikle, kısa ve gürültülü veri setleri üzerinde iyi çalışır ve aynı zamanda tutarlı bir ölçü sağlar. Ancak, parametrelerin (m ve r) doğru bir şekilde ayarlanması gerektiğini unutmamak önemlidir, çünkü bu, sonuç üzerinde büyük bir etkiye sahip olabilir.
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = """ 
import numpy as np
from scipy.spatial.distance import cdist

def ApEn(U, m, r):
    #Approximate Entropy hesaplama fonksiyonu
    def _maxdist(x_i, x_j):
        return max([abs(ua - va) for ua, va in zip(x_i, x_j)])

    def _phi(m):
        x = [[U[j] for j in range(i, i + m - 1 + 1)] for i in range(N - m + 1)]
        C = [len([1 for x_j in x if _maxdist(x_i, x_j) <= r]) / (N - m + 1.0) for x_i in x]
        return (N - m + 1.0)**(-1) * sum(np.log(C))

    N = len(U)
    return abs(_phi(m+1) - _phi(m))

# BBS algoritması ile oluşturulmuş 1 milyon adet 1 ve 0 dan oluşan array
U = np.random.choice([0, 1], size=(1000000,))

# m ve r değerlerini belirle
m = 2
r = 0.2 * np.std(U)

# ApEn hesapla
print(ApEn(U, m, r))
"""
st.code(code, language='python')
st.markdown("""Bu kod, verilen zaman serisi U üzerinde Approximate Entropy'yi hesaplar. m ve r parametrelerini belirlemek gerekiyor. m genellikle 2 ya da 3 olarak seçilir ve r genellikle verinin standart sapmasının bir yüzdesi olarak seçilir (bu örnekte %20).
Bir not: r değerinin seçimi ApEn hesaplamasının sonuçları üzerinde büyük bir etkisi olabilir. Bu sebeple, r değeri genellikle verinin standart sapmasının bir yüzdesi olarak seçilir. Ancak, bu değeri belirlemek için farklı stratejiler de vardır ve hangi stratejinin kullanılacağı genellikle uygulamaya bağlıdır.
Sonuçta, bu işlev, belirtilen parametrelerle Approximate Entropy'yi hesaplar ve sonucu döndürür.
""")
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="12")
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
    m = 2
    r = 0.2 * np.std(bit_array)
    result = approx.ApEn_fast(bit_array,m,r)
    st.subheader("Test Başarılı✔️")
    st.subheader(f"Test Sonucu : {result}")
st.subheader(":green[Cumulative Sums (Cusum) Test]")
st.markdown("""`Cumulative Sums (Cusum) Testi`, bir zaman serisindeki verilerin rastgeleliğini test etmek için kullanılan istatistiksel bir testtir. Bu test, bir veri setinin "merkezi" değerinden ne kadar sapma olduğunu ölçer ve bu sapmaların toplamını alır. Bu, bir veri setindeki eğilimleri belirlemek için kullanılabilir - örneğin, bir süreçteki belirli bir durumun zaman içinde nasıl değiştiğini görmek için.
Cusum testi genellikle iki formda uygulanır: `İleri Cusum` testi ve `Geri Cusum` testi. İleri Cusum testi, veri setinin başından sonuna doğru ilerler ve her veri noktasında toplam sapmayı hesaplar. Geriye doğru Cusum testi, veri setinin sonundan başına doğru ilerler ve aynı işlemi yapar. Bu iki test genellikle birlikte kullanılır ve her ikisinin sonuçları karşılaştırılır.
Cusum testinin matematiksel formülü şu şekildedir:
""")
image = Image.open('cumulative.png')
st.image(image, caption='Cumulative Sums (Cusum) Test Formülü')
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = """ 
import numpy as np
from scipy.spatial.distance import cdist

def ApEn(U, m, r):
    #Approximate Entropy hesaplama fonksiyonu
    def _maxdist(x_i, x_j):
        return max([abs(ua - va) for ua, va in zip(x_i, x_j)])

    def _phi(m):
        x = [[U[j] for j in range(i, i + m - 1 + 1)] for i in range(N - m + 1)]
        C = [len([1 for x_j in x if _maxdist(x_i, x_j) <= r]) / (N - m + 1.0) for x_i in x]
        return (N - m + 1.0)**(-1) * sum(np.log(C))

    N = len(U)
    return abs(_phi(m+1) - _phi(m))

# BBS algoritması ile oluşturulmuş 1 milyon adet 1 ve 0 dan oluşan array
U = np.random.choice([0, 1], size=(1000000,))

# m ve r değerlerini belirle
m = 2
r = 0.2 * np.std(U)

# ApEn hesapla
print(ApEn(U, m, r))
"""
st.code(code, language='python')
st.markdown("""Bu kod, veri setinizin ortalama değerinden (bu durumda, 0 ve 1'in ortalaması 0.5) ne kadar sapma olduğunu ölçer ve bu sapmaların kümülatif toplamını alır. Sonuçlar, veri setindeki veri noktalarının indeksine göre bir grafikte gösterilir.
Sonuçtaki grafik, veri setindeki her bir veri noktasında kümülatif toplamı gösterir. Bu, veri setinin nasıl değiştiğini ve belirli bir noktada herhangi bir eğilimin olup olmadığını görmek için kullanılabilir.
Bu kod, veri setinizin ortalama değerinden ne kadar sapma olduğunu ölçer. Bu sapmaların kümülatif toplamını alarak, veri setinde belirli bir eğilimin olup olmadığını belirlemeye yardımcı olabilir. Cusum testi genellikle bir süreçteki değişimleri izlemek ve belirli bir noktada bir değişikliğin olup olmadığını belirlemek için kullanılır. Bu, özellikle kalite kontrol ve endüstriyel süreçlerin izlenmesi gibi uygulamalarda yaygın olarak kullanılır.
""")
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="13")
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
    m = 2
    r = 0.2 * np.std(bit_array)
    bit_array = np.array(bit_array)
    result = sums.cusum_test(bit_array)
    st.subheader(f"👉 {result} ")
    fig, ax = plt.subplots()
    ax.plot(result)
    ax.set_title("Cumulative Sums (Cusum) Test Sonuç")
    ax.set_xlabel("Index")
    ax.set_ylabel("Cusum Value")
    st.pyplot(fig)
st.subheader(":green[Random Excursions Test]")
st.markdown("""`Random Excursions Testi`, bir veri dizisi içinde belirli bir değere (genellikle 0'a) kaç kez dönüldüğünü ve bu dönüşler arasında geçen sürenin dağılımını inceler. Bu test, rastgele bir süreç tarafından üretilen verilerin özelliklerini incelemek için kullanılır.
Random Excursions Testi, genellikle bir süreç tarafından üretilen verilerin rastgeleliğini değerlendirmek için kullanılır. Özellikle, bir süreçten elde edilen verilerin zamanla nasıl değiştiğini ve belirli bir değere dönme eğilimini incelemek için kullanılır.
Bu test, veri serisinde belirli bir değere kaç kez dönüldüğünü ve bu dönüşler arasında geçen sürenin dağılımını inceler. Örneğin, bir süreç 0 ve 1 değerlerini rastgele üretiyorsa, Random Excursions Testi, sürecin 0 değerine dönme eğilimini ve bu 0 değerine dönüşler arasında geçen süreyi inceleyecektir.
Testin adı, "excursion" kelimesinin ("gezi" veya "seyahat") veri setinin belirli bir değere "dönüşünü" ve bu dönüşler arasındaki süreyi temsil etmesinden gelir.
Bu test, özellikle rastgele sayı üreten süreçlerin kalitesini değerlendirmek için kullanılır. Eğer bir süreç tamamen rastgele sayılar üretiyorsa, bu sayılar arasında belirli bir değere dönme eğilimi veya dönüşler arasındaki sürenin belirli bir dağılımı olmamalıdır. Bu nedenle, Random Excursions Testi, bir sürecin rastgeleliğini değerlendirmek için kullanışlı bir araçtır.
Bu test ayrıca, özellikle kriptografi ve güvenlik uygulamalarında, bir rastgele sayı üretecinin kalitesini değerlendirmek için de kullanılır. Rastgele sayı üreteçleri, güvenli iletişim ve veri şifreleme için kritik öneme sahiptir ve bu nedenle, ürettikleri sayıların gerçekten rastgele olup olmadığını doğrulamak önemlidir. Random Excursions Testi, bu tür bir doğrulama için bir yöntem sağlar.
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = """import numpy
from scipy import special as spc

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
"""
st.code(code, language='python')
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="14")
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
    m = 2
    r = 0.2 * np.std(bit_array)
    result = randex.random_excursions(bit_array)
    st.subheader("Test Başarılı✔️")
    st.subheader(f"Sonuç : {result}")
st.subheader(":green[Random Excursions Variant Test]")
st.markdown("""`Random Excursions Variant Test`, bir dizi bitin rastgeleliğini değerlendirmek için kullanılan bir istatistiksel testtir. Bu test, bir dizi bitin çeşitli "durumlar" arasında nasıl "gezindiğini" inceler. Bu "durumlar", bir dizi bitin toplam değerinin çeşitli değerlere ulaştığı yerlerdir.
Bu testin amacı, bir dizi bitin belirli bir durumu ziyaret etme sıklığının beklenen rastgelelikten önemli ölçüde sapıp sapmadığını belirlemektir.
Testin detayları aşağıdaki gibidir:
1) `Test, dizi üzerinde bir kümülatif toplam (veya "rastgele yürüyüş") hesaplar. Bu, dizinin her bir bitini -1 veya +1 değerine (genellikle 0'lar -1'e, 1'ler +1'e eşlenir) dönüştürdükten ve bu değerleri topladıktan sonra elde edilir.`

2) `Daha sonra, bu kümülatif toplamın -9'dan +9'a kadar olan tüm durumları ziyaret etme sayısını hesaplar.`

3) `Elde edilen ziyaret sayıları, beklenen değerlerle karşılaştırılır. Beklenen değerler, rastgele bir bit dizisi için teorik olarak hesaplanabilir.`

4) `Sonuçlar, bir chi-kare istatistiği kullanılarak değerlendirilir. Chi-kare istatistiği, gözlenen ve beklenen ziyaret sayıları arasındaki farkı ölçer.`

5) `Eğer chi-kare değeri belirli bir eşik değerinden büyükse (bu eşik genellikle bir p-değeri kullanılarak belirlenir), test dizi bitinin beklenen rastgelelikten anlamlı şekilde sapmış olduğunu gösterir.`
Bu test, bit dizilerinin rastgeleliğini değerlendirmek için kullanılır. Örneğin, bir rastgele sayı üreteci (RNG) tarafından üretilen bitlerin kalitesini test etmek için kullanılabilir. RNG'nin gerçekten rastgele bitler üretip üretmediğini belirlemek için bu tür istatistiksel testler kullanılır.
""")
st.markdown("Örnek :red[Python] Kodu Aşağıdaki Gibidir:")
code = """import numpy
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

"""
st.code(code, language='python')
st.markdown("""Bu kodu kullanarak, veri setinizi random_excursions_variant fonksiyonuna vererek testi gerçekleştirebilirsiniz. Bu fonksiyon, her bir durum için p-değerlerini döndürür.
""")
st.subheader(":red[Kodu BBS Ekranındaki Algoritma Sonucu ile Test Etmek İçin :]")
button = st.button("Tıklayın",key="15")
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
    m = 2
    r = 0.2 * np.std(bit_array)
    result = randvar.random_excursions_variant(bit_array)
    st.subheader("Test Başarılı✔️")
    st.subheader(f"Sonuç : {result}")
st.sidebar.success("Gitmek İstediğiniz Sayfayı Seçiniz")