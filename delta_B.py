import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import *
from scipy.special import *
def sb(l, n, prob, output) :
    tol = 0.00001
    sig = l ** 0.5
    k = l + n * sig
    p0 = l + n * sig
    p1 = l + (n + 1) * sig
    q0 = gammaincc(k+1, p0) - (1 - prob)
    q1 = gammaincc(k+1, p1) - (1 - prob)
    while True :
        p = p1 - q1 * ((p1 - p0) / (q1 - q0))
        if abs(p - p1) < tol :
            break
        q = gammaincc(k+1, p) - (1 - prob)
        if q * q1 < 0 :
            p0 = p1
            q0 = q1

        p1 = p
        q1 = q


    if output == 's' or output == 'S':
        return p - l
    elif output == 'total':
        return p        
    elif output == 'ceil':
        return math.ceil(p - l)

def delta_B(B0, dB, acc, n, prob, mod):
    r = np.arange(0.0001,100,acc)
    rl = r.tolist()
    sbv = np.vectorize(sb)
    B_ac = sbv(rl, n, prob,'s')
    
    lower_limit = (B0 - dB)
    upper_limit = (B0 + dB)

    index_l = np.argmin(abs(r - lower_limit))
    index_u = np.argmin(abs(r - upper_limit))

    interval = r[index_l : index_u + 1]
    #print(interval)
    G = norm.pdf(interval,loc = B0, scale = (dB/3))
    #print(G)
    Gn = np.asarray(G)
    B_ac_n = np.asarray(B_ac)
    B_interval = B_ac_n[index_l : index_u + 1]
    result = Gn * B_interval
    #print(result)
    value = sum(result)/sum(Gn)

    if mod == 'minus':
        return (value - B0)
    if mod == 'origin':
        return value
    if mod == 'ceil' :
        return math.ceil(value - B0)

print(delta_B(50, 20, 0.01, 3, 0.5, 'minus'))

delta_Bv = np.vectorize(delta_B)

a = np.arange(0.01,5,0.01)
aa = a.tolist()

plt.plot(aa, delta_Bv(aa, 0.5, 0.01, 3, 0.9, 'ceil'))
'''
plt.subplot(511)
plt.plot(aa, delta_Bv(5, aa, 0.01, 3, 0.5, 'minus'), label = 'B0 = 5')
plt.legend(loc = 'best')

plt.subplot(512)
plt.plot(aa, delta_Bv(10, aa, 0.01, 3, 0.5, 'minus'), label = 'B0 = 10')
plt.legend(loc = 'best')

plt.subplot(513)
plt.plot(aa, delta_Bv(30, aa, 0.01, 3, 0.5, 'minus'), label = 'B0 = 30')
plt.legend(loc = 'best')

plt.subplot(514)
plt.plot(aa, delta_Bv(50, aa, 0.01, 3, 0.5, 'minus'), label = 'B0 = 50')
plt.legend(loc = 'best')

plt.subplot(515)
plt.plot(aa, delta_Bv(70, aa, 0.01, 3, 0.5, 'minus'), label = 'B0 = 70')
plt.legend(loc = 'best')
'''
plt.show()

