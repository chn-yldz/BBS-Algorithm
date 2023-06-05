import streamlit as st
from PIL import Image
import random
from sympy import isprime
from io import StringIO
import scipy.stats as stats

st.set_page_config(page_title="Blum Blum Shum (BBS) Algorithm",
                   page_icon=":computer:",
                   layout="wide")
st.title("Blum Blum Shum (BBS) Algoritması")


st.header("Algoritma Hakkında:")
st.markdown("""Blum Blum Shub (BBS) algoritması, 1986 yılında Manuel Blum,
Lenore Blum ve Michael Shub tarafından geliştirilen bir
rastgele sayı üreteci algoritmasıdır. Bu algoritma, özellikle
kriptografik uygulamalar için güvenli rastgele sayılar üretmek
amacıyla tasarlanmıştır. BBS, genellikle güvenlik açısından
yüksek gereksinimlere sahip uygulamalarda kullanılır. BBS
algoritması, iki büyük asal sayının çarpımı olan bir sayı
üzerinde çalışır. Bu sayılar genellikle p ve q olarak adlandırılır
ve çarpımları M (M = p*q) olarak adlandırılır. Formül
aşağıdaki gibidir :
""")
image = Image.open('formul.png')
st.image(image, caption='BBS Formülü')
st.markdown("""BBS algoritması, genel olarak bu adımları takip eder:
İlk olarak, iki büyük asal sayı seçilir. Bu asal sayıların 4 ile
bölündüğünde 3 kalanını vermesi gerekmektedir. Bu asal
sayıları p ve q olarak adlandıralım.
Bu iki asal sayının çarpımı olan M = p*q hesaplanır. M
sayısı, BBS algoritmasının modülü olarak kullanılır.
Ardından, M'den küçük ve M ile aralarında asal olan bir
başlangıç değeri (seed) seçilir. Bu başlangıç değerini x
olarak adlandıralım.
Her yeni rastgele biti üretmek için, x
'in karesi alınır ve M'
ye
göre modülü alınır. Yani, yeni x değeri, x = x^2 mod M
formülü ile hesaplanır.
Üretilen rastgele bit, x
'in en düşük anlamlı biti (LSB) olarak
alınır. Bu adımlar, yeni rastgele bitler üretmek için tekrarlanır.
""")
st.header("Örnek BBS Algoritması Aşağıdaki Gibidir:")
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
st.subheader("Algoritmanın Çalışma Mantığı :")
st.markdown("""Bu Python kodu, büyük bir boyutta rastgele asal sayılar üretir ve ardından bir Blum Blum Shub rastgele bit üreteci kullanır. Aşağıda her bir kod parçasının ne yaptığı açıklanmıştır:

1. `import random` ve `from sympy import isprime`: Bu satırlar, Python'da rastgele sayılar oluşturmayı ve bir sayının asal olup olmadığını kontrol etmeyi sağlayan kütüphaneleri dahil eder.

2. `asal_bul` fonksiyonu: Bu fonksiyon, 512 bitlik bir rastgele sayı üretir ve bu sayının asal olup olmadığını kontrol eder. Ayrıca üretilen asal sayının 4 ile bölümünden kalanın 3 olup olmadığını kontrol eder. Eğer bu koşullar sağlanıyorsa, asal sayıyı döndürür. Bu, Blum asal sayısını bulmak için kullanılan bir metottur. Bir asal sayının Blum asalı olabilmesi için 4 ile bölündüğünde kalanın 3 olması gerekir.

3. `p = asal_bul()` ve `q = asal_bul()`: Bu satırlar, `asal_bul` fonksiyonunu çağırarak iki farklı Blum asal sayısı üretir.

4. `M = p * q`: Bu satır, `p` ve `q` asal sayılarının çarpımını hesaplar. Bu değer, Blum Blum Shub üretecinin modülo değeri olarak kullanılacaktır.

5. `x = random.randint(1,M)`: Bu satır, 1 ve `M` arasında rastgele bir tam sayı seçer.

6. Bit dizisi oluşturma döngüsü: Bu döngü, Blum Blum Shub rastgele bit üreteci algoritmasını uygular. Algoritma, `x` değerinin karesini `M` ile modulo alır ve sonucu yeni `x` değeri olarak atar. Ardından, `x` değerini 2 ile modulo alır ve sonucu bit dizisine ekler. Bu işlem, bir milyon kez tekrarlanır.

7. `dosya = open("bit_array.txt","w")`: Bu satır, yazma modunda "bit_array.txt" adlı bir dosya açar.

8. `dosya.write(''.join(map(str, bit_array)))`: Bu satır, bit dizisini bir metin dizesine dönüştürür ve bu diziyi dosyaya yazar.

9. `dosya.close()`: Bu satır, dosyayı kapatır. Bu, dosyanın düzgün bir şekilde kaydedilmesini ve bellek kaynaklarının serbest bırakılmasını sağlar.

Sonuç olarak, bu kod, kriptografide sıklıkla kullanılan güçlü rastgelelik özelliklerine sahip bir bit dizisi oluşturur.""")
st.subheader(":red[Algoritmayı Denemek İçin :]")
button = st.button("Tıklayın")
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
    dosya = open("bit_array.txt","w")
    dosya.write(''.join(map(str, bit_array)))
    dosya.close()
    with open("bit_array.txt", "r") as dosya:
        data = dosya.read()

    st.download_button('Kodun Çıktısını İndirmek İçin Tıklayın', data)

#st.sidebar.header("NIST Test Menüsü")
st.subheader("Proje Hakkında Görüş Ve Önerileriniz İçin :mailbox:")
contact_form = """
<form action="https://formsubmit.co/chnx35x@gmail.com" method="POST">
     <input type="hidden" name="_captcha" value="false">
     <input type="text" name="name" placeholder="İsminiz" required>
     <input type="email" name="email" placeholder = "E-Posta Adresiniz" required>
     <textarea name="message" placeholder="Mesajınız"></textarea>
     <button type="submit">Gönder</button>
</form>
"""
st.markdown(contact_form,unsafe_allow_html=True)
def local_css(file_name):
    with open(file_name) as f:
        st.markdown(f"<style>{f.read()}</style>",unsafe_allow_html=True)
local_css("style/style.css")
st.sidebar.success("Gitmek İstediğiniz Sayfayı Seçiniz")