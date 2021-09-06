#!/usr/local/bin/python3
import matplotlib.pyplot as plt
import numpy as np

#<functions>
def F (x_data):
    assert (isinstance (x_data, np.ndarray))
    return np.exp (x_data)

def dF (x_data):
    assert (isinstance (x_data, np.ndarray))
    return np.exp (x_data)
#</functions>

#<dimensions>
x_low = -5
x_high = 5
y_low = -0.5
y_high = 100
x_num = 20
x_data = np.linspace (x_low, x_high, x_num)
x_delta = (x_high - x_low) / (x_num - 1)
#</dimensions>

#<style>
graph_title = "Graph"
x_label = ""
y_label = ""
plt.axis ([x_low, x_high, y_low, y_high])
plt.title (graph_title)
plt.xlabel (x_label)
plt.ylabel (y_label)
plt.grid (True)
#</style>

#<presize d(sin(x))/dx>
plt.plot (x_data, dF (x_data), 'b-')
#</presize d(sin(x))/dx>

#<O(n) right (Fi+1 - Fi) / h>
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
#</O(n) right (Fi+1 - Fi) / h>

#<O(n) left (Fi - Fi-1) / h>
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
#</O(n) left (Fi - Fi-1) / h>

#<O(n2) (Fi+1 - Fi-1) / 2h; F'0 = (-3Fi + 4Fi+1 - Fi+2) / 2h; F'n = -(-3Fi + 4Fi-1 - Fi-2) / 2h>
def Oleft_dy (x_data, y_data): # differentiating net-function (сеточная ф-ция)
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

plt.plot (x_data, Oleft_dy (x_data, F (x_data)), 'm-')
#</O(n2) (Fi+1 - Fi-1) / 2h; F'0 = (-3Fi + 4Fi+1 - Fi+2) / 2h; F'n = -(-3Fi + 4Fi-1 - Fi-2) / 2h>

plt.show ()