#!/usr/bin/env python3
import matplotlib.pyplot as plt
import math
import numpy as np

#<functions>
def F (x_data):
    assert (isinstance (x_data, np.ndarray))
    return np.sin (x_data)

def dF (x_data):
    assert (isinstance (x_data, np.ndarray))
    return np.cos (x_data)

def d2F (x_data):
    assert (isinstance (x_data, np.ndarray))
    return -np.sin (x_data)
#</functions>

#<dimensions>
x_low = -5
x_high = 5
y_low = -1.5
y_high = 1.5
x_num = 100
x_num_max = 100 #min is 30
x_data = np.linspace (x_low, x_high, x_num)
x_delta = (x_high - x_low) / (x_num - 1)
#</dimensions>

#<style>
graph_title = "Graph"
x_label = ""
y_label = ""
plt.title (graph_title)
plt.xlabel (x_label)
plt.ylabel (y_label)
#</style>

#<dF/dx>
plt.subplot (221) # 2 line 2 graphs in line 
plt.grid (True)

#<presize dF/dx>
plt.plot (x_data, dF (x_data), 'b-')
#</presize dF/dx>
#</dF/dx>

#<O(h) right (Fi+1 - Fi) / h>
def Oright_dy (x_data, y_data): # differentiating net-function (сеточная ф-ция)
    assert (isinstance (x_data, np.ndarray))
    assert (isinstance (y_data, np.ndarray))
    assert (x_data.size == y_data.size)
    len = x_data.size
    dy = np.empty (len)
    for i in range (0, len - 1):
        dy[i] = (y_data[i + 1] - y_data[i]) / (x_data[i + 1] - x_data[i])
    dy[len - 1] = (y_data[len - 1] - y_data [len - 2]) / (x_data[len - 1] - x_data[len - 2])
    return dy

plt.plot (x_data, Oright_dy (x_data, F (x_data)), 'g-')
#</O(h) right (Fi+1 - Fi) / h>

#<O(h) left (Fi - Fi-1) / h>
def Oleft_dy (x_data, y_data): # differentiating net-function (сеточная ф-ция)
    assert (isinstance (x_data, np.ndarray))
    assert (isinstance (y_data, np.ndarray))
    assert (x_data.size == y_data.size)
    len = x_data.size
    dy = np.empty (len)
    for i in range (1, len):
        dy[i] = (y_data[i] - y_data[i - 1]) / (x_data[i] - x_data[i - 1])
    dy[0] = (y_data[1] - y_data [0]) / (x_data[1] - x_data[0])
    return dy

plt.plot (x_data, Oleft_dy (x_data, F (x_data)), 'r-')
#</O(h) left (Fi - Fi-1) / h>

#<O(h2) (Fi+1 - Fi-1) / 2h; F'0 = (-3Fi + 4Fi+1 - Fi+2) / 2h; F'n = -(-3Fi + 4Fi-1 - Fi-2) / 2h>
def O2_dy (x_data, y_data): # differentiating net-function (сеточная ф-ция)
    assert (isinstance (x_data, np.ndarray))
    assert (isinstance (y_data, np.ndarray))
    assert (x_data.size == y_data.size)
    len = x_data.size
    dy = np.empty (len)
    for i in range (1, len - 1):
        dy[i] = (y_data[i + 1] - y_data[i - 1]) / (x_data[i + 1] - x_data[i - 1])
    dy[0] = (-3 * y_data[0] + 4 * y_data[1] - y_data[2]) / (x_data[2] - x_data[0])
    dy[len-1] = -(-3 * y_data[len - 1] + 4 * y_data[len - 2] - y_data[len - 3]) / (x_data[len - 1] - x_data[len - 3])
    return dy

plt.plot (x_data, O2_dy (x_data, F (x_data)), 'm-')
#</O(h2) (Fi+1 - Fi-1) / 2h; F'0 = (-3Fi + 4Fi+1 - Fi+2) / 2h; F'n = -(-3Fi + 4Fi-1 - Fi-2) / 2h>
#</dF/dx>

#<d2F/dx>
#<presize d2F/dx>
plt.subplot (222)
plt.grid (True)
plt.plot (x_data, d2F (x_data), 'b-')
#</presize d2F/dx>

#<O2(h) (dFi+1 - 2dFi + 2dFi-1) / h^2>
def O2_d2y_fixed_h (x_data, y_data):
    assert (isinstance (x_data, np.ndarray))
    assert (isinstance (y_data, np.ndarray))
    assert (x_data.size == y_data.size)
    len = x_data.size
    h2 = (x_data[1] - x_data[0]) * (x_data[1] - x_data[0])
    d2y = np.empty (len)
    for i in range (1, len - 1):
        d2y[i] = (y_data[i + 1] + y_data[i - 1] - 2 * y_data[i]) / h2
    d2y[len - 1] = (2 * y_data[len - 1] - 5 * y_data[len - 2] + 4 * y_data[len - 3] - y_data[len - 4]) / h2
    d2y[0] = (2 * y_data[0] - 5 * y_data[1] + 4 * y_data[2] - y_data[3]) / h2
    return d2y

plt.subplot (222)
plt.grid (True)
plt.plot (x_data, O2_d2y_fixed_h (x_data, F (x_data)), 'm-')
#</O2(h) (dFi+1 - 2dFi + 2dFi-1) / h^2>
#</d2F/dx>

#<F mistake D = max |F(i) - Fcalculated(i)|>
def D_y (x_data, y_data, y0_data):
    assert (isinstance (x_data, np.ndarray))
    assert (isinstance (y0_data, np.ndarray))
    assert (isinstance (y_data, np.ndarray))
    assert (x_data.size == y0_data.size == y_data.size)
    len = x_data.size
    D_y = -1
    for i in range (0, len):
        if D_y < abs (y_data[i] - y0_data[i]):
            D_y = abs (y_data[i] - y0_data[i])
    return D_y
#</F mistake D = max |F(i) - Fcalculated(i)|>

#<mistake graphs>
n_data = np.arange (4, x_num_max) #4 is minimal due to o2
dy_data = np.empty (x_num_max - 4)
h_data = (x_high - x_low) / n_data

#<dF mistake graph>
plt.subplot (223)
plt.grid (True)
plt.yscale ('log')
plt.xscale ('log')
for i in range (4, x_num_max):
    temp_x_data = np.linspace (x_low, x_high, i)
    dy_data[i - 4] = D_y (temp_x_data, Oleft_dy (temp_x_data, F (temp_x_data)), dF (temp_x_data))
plt.plot (n_data, dy_data, 'ro')
for i in range (4, x_num_max):
    temp_x_data = np.linspace (x_low, x_high, i)
    dy_data[i - 4] = D_y (temp_x_data, Oright_dy (temp_x_data, F (temp_x_data)), dF (temp_x_data))
plt.plot (n_data, dy_data, 'g.')
kdy = (math.log (dy_data[x_num_max - 5]) - math.log (dy_data[math.floor (0.7 * (x_num_max - 5))])) / (math.log (n_data[x_num_max - 5]) - math.log (n_data[math.floor (0.7 * (x_num_max - 5))]))
plt.text (n_data[math.floor(0.7 * (x_num_max - 5))], dy_data[math.floor(0.7 * (x_num_max - 5))],  '{:.6f}'.format(kdy))
for i in range (4, x_num_max):
    temp_x_data = np.linspace (x_low, x_high, i)
    dy_data[i - 4] = D_y (temp_x_data, O2_dy (temp_x_data, F (temp_x_data)), dF (temp_x_data))
plt.plot (n_data, dy_data, 'm-')
plt.plot (n_data, h_data, 'c-')
kh = (math.log (h_data[x_num_max - 5]) - math.log (h_data[math.floor (0.7 * (x_num_max - 5))])) / (math.log (n_data[x_num_max - 5]) - math.log (n_data[math.floor (0.7 * (x_num_max - 5))]))
plt.text (n_data[0], dy_data[10] / 10, 'kdy_O1 / kh = ' + '{:.6f}'.format(kdy / kh))
kdy = (math.log (dy_data[x_num_max - 5]) - math.log (dy_data[math.floor (0.7 * (x_num_max - 5))])) / (math.log (n_data[x_num_max - 5]) - math.log (n_data[math.floor (0.7 * (x_num_max - 5))]))
plt.text (n_data[math.floor(0.7 * (x_num_max - 5))], dy_data[math.floor(0.7 * (x_num_max - 5))],  '{:.6f}'.format(kdy))
plt.text (n_data[math.floor(0.7 * (x_num_max - 5))], h_data[math.floor(0.7 * (x_num_max - 5))], '{:.6f}'.format(kh))
plt.text (n_data[0], dy_data[10] / 20, 'kdy_O2 / kh = ' + '{:.6f}'.format(kdy / kh))
#</dF mistake graph>

#<d2F mistake graph>
plt.subplot (224)
plt.grid (True)
plt.xscale ('log')
plt.yscale ('log')
for i in range (4, x_num_max):
    temp_x_data = np.linspace (x_low, x_high, i)
    dy_data[i - 4] = D_y (temp_x_data, O2_d2y_fixed_h (temp_x_data, F (temp_x_data)), d2F (temp_x_data))
plt.plot (n_data, dy_data, 'm-')
plt.plot (n_data, h_data, 'c-')
kh = (math.log (h_data[x_num_max - 5]) - math.log (h_data[math.floor (0.7 * (x_num_max - 5))])) / (math.log (n_data[x_num_max - 5]) - math.log (n_data[math.floor (0.7 * (x_num_max - 5))]))
plt.text (n_data[math.floor(0.7 * (x_num_max - 5))], h_data[math.floor(0.7 * (x_num_max - 5))],  '{:.6f}'.format(kh))
kd2y = (math.log (dy_data[x_num_max - 5]) - math.log (dy_data[math.floor (0.7 * (x_num_max - 5))])) / (math.log (n_data[x_num_max - 5]) - math.log (n_data[math.floor (0.7 * (x_num_max - 5))]))
plt.text (n_data[math.floor(0.7 * (x_num_max - 5))], dy_data[math.floor(0.7 * (x_num_max - 5))],  '{:.6f}'.format(kd2y))
plt.text (n_data[0], dy_data[10] / 20, 'k2dy / kh = ' + '{:.6f}'.format(kd2y / kh))
#</d2F mistake graph>

plt.show ()