from scipy.optimize import curve_fit
import numpy as np
import math



# ys = np.array([0.0322428941726685,0.00215978622436523,0.00177663564682007,0.00187919139862061])
# ys = np.array([float('6.9301894732884e-05'), float('8.52587677183605e-06'),
#                float('5.09910640262422e-06'), float('4.42791552770706e-06')])
# xs = np.array([1,25,64,100])
# #ys = np.array([0.0334906578063965,0.00186474323272705,0.00111260414123535,0.00100641250610352])

# def f(xs, a, b, c):
#     return a + b / (xs ** c)
#     fanin, fanout = xs[:,0], xs[:,1]
#     return a + b*fanin# + c*fanout
# popt, pcov = curve_fit(f, xs, ys)
# print map(float,popt)
# print pcov


