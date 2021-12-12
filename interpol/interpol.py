#!/usr/bin/env python3
import matplotlib.pyplot as plt
import math
import numpy as np
import sys

#<functions>
function_name = "exp(x)"
def F (x_data):
    assert (isinstance (x_data, np.ndarray))
    return np.exp (x_data)
#</functions>

#<dimensions and data>
x_low = -5
x_high = 5
y_low = -1.5
y_high = 1.5
x_num = 5
x_num_big = 1000

x_data = np.linspace (x_low, x_high, x_num)
x_data_big = np.linspace (x_low, x_high, x_num_big)
x_delta = (x_high - x_low) / (x_num - 1)
x_delta_big = (x_high - x_low) / (x_num_big - 1)

f_data = F(x_data)
interpol_data = np.empty (x_num_big)
#</dimensions and data>

#<lagrange>
def calculate_lagrange ():
    lagrange_data = np.empty (x_num_big)
    for i in range (0, x_num_big):
        for j in range (0, x_num):
            temp = f_data[j]
            for k in range (0, j):
                temp = temp * (x_data_big[i] - x_data[k]) / (x_data[j] - x_data[k])
            for k in range (j + 1, x_num):
                temp = temp * (x_data_big[i] - x_data[k]) / (x_data[j] - x_data[k])
            lagrange_data[i] += temp
    return lagrange_data

def change_knot ():
#requires interpol_data to be precalculated by lagrange
    knot_num = int (input ("enter knot num [0-%d]: " % (x_num - 1)))
    while knot_num < 0 or knot_num > x_num - 1:
        knot_num = int (input ("please, enter knot num [0-%d]: " % (x_num - 1)))
    knot_new_val = np.float64 (input ("enter new knot value: "))
    for i in range (0, x_num_big):
        temp = 1
        for k in range (0, knot_num):
            temp = temp * (x_data_big[i] - x_data[k]) / (x_data[knot_num] - x_data[k])
        for k in range (knot_num + 1, x_num):
            temp = temp * (x_data_big[i] - x_data[k]) / (x_data[knot_num] - x_data[k])
        interpol_data[i] += temp * knot_new_val - temp * f_data[knot_num]
    f_data[knot_num] = knot_new_val
    return
#</lagrange>

#<newton>
divsubs_all = []
def calculate_newton ():
    #<calculate divsubs>
    for i in range (0, x_num - 1):
        divsubs_i = []
        for j in range (0, x_num - 1 - i):
            if i == 0:
                divsubs_i.append ((f_data[j + 1] - f_data[j]) / (x_data[j + 1] - x_data[j]))
            else:
                divsubs_i.append ((divsubs_all[i - 1][j + 1] - divsubs_all[i - 1][j]) / (x_data[j + i + 1] - x_data[j]))
        divsubs_all.append (divsubs_i)
    #</calculate divsubs>
    newton_data = np.empty (x_num_big)
    for i in range (0, x_num_big):
        newton_data[i] = f_data[0]
        for j in range (0, x_num - 1):
            temp = divsubs_all[j][0]
            for k in range (0, j + 1):
                temp *= (x_data_big[i] - x_data[k])
            newton_data[i] += temp
    return newton_data

def add_knot ():
    global x_num
    global x_data
    global f_data
    global x_data_big
    global x_low
    global x_high
    knot_x = np.float64 (input ("enter knot coord on the edge [<%f or >%f]: " % (x_low, x_high)))
    while (knot_x >= x_low and knot_x <= x_high):
        knot_x = np.float64 (input ("please, enter knot coord on the edge [<%f or >%f]: " % (x_low, x_high)))
    knot_val = np.float64 (input ("enter knot value: "))
    ind = 0
    if knot_x > x_high:
        ind = x_num
        x_data = np.append (x_data, [knot_x])
        f_data = np.append (f_data, [knot_val])
        x_data_big = np.linspace (x_low, knot_x, x_num_big)
        x_high = knot_x
    else:
        x_data = np.insert (x_data, ind, [knot_x])
        f_data = np.insert (f_data, ind, [knot_val])
        x_data_big = np.linspace (knot_x, x_high, x_num_big)
        x_low = knot_x
    x_num += 1
    #<recalculate divsubs>, max N + N + 1 + N - 1, no matter where insertion takes place
    if ind == 0:
        for i in range (0, x_num - 1):
            if i == x_num - 2:
                divsubs_last = []
                divsubs_last.append ((divsubs_all[i - 1][1] - divsubs_all[i - 1][0]) / (x_data[i + 1] - x_data[0]))
                divsubs_all.append (divsubs_last)
            else:
                if i == 0:
                    divsubs_all[i].insert (0, (f_data[1] - f_data[0]) / (x_data[1] - x_data[0]))
                else:
                    divsubs_all[i].insert (0, (divsubs_all[i - 1][1] - divsubs_all[i - 1][0]) / (x_data[i + 1] - x_data[0]))
    if ind == x_num - 1:
        for i in range (0, x_num - 1):
            if i == x_num - 2:
                divsubs_last = []
                divsubs_last.append ((divsubs_all[i - 1][1] - divsubs_all[i - 1][0]) / (x_data[i + 1] - x_data[0]))
                divsubs_all.append (divsubs_last)
            else:
                if i == 0:
                    divsubs_all[i].append ((f_data[ind] - f_data[ind - 1]) / (x_data[ind] - x_data[ind - 1]))
                else:
                    divsubs_all[i].append ((divsubs_all[i - 1][ind - i] - divsubs_all[i - 1][ind - i - 1]) / (x_data[ind] - x_data[ind - i - 1]))
    #</recalculate divsubs>
    for i in range (0, x_num_big):
        interpol_data[i] = f_data[0]
        for j in range (0, x_num - 1):
            temp = divsubs_all[j][0]
            for k in range (0, j + 1):
                temp *= (x_data_big[i] - x_data[k])
            interpol_data[i] += temp
    return
#</newton>

#<spline>
def calculate_spline ():
    # natural cubic splines: f'' (ends) == 0
    # si = fi + bi * (x - xi) + ci * (x - xi)^2 + di * (x - xi)^3
    #<calculate coefficients>
    spline_data = np.empty (x_num_big)
    a = np.empty (x_num - 1)
    a[0] = 0 #not used
    c = np.empty (x_num)
    l = np.empty (x_num)
    l[0] = 1
    m = np.empty (x_num)
    m[0] = 0
    z = np.empty (x_num)
    z[0] = 0
    for i in range (1, x_num - 1):
        a[i] = (3 * (f_data[i + 1] - f_data[i]) / (x_data[i + 1] - x_data[i]) - 3 * (f_data[i] - f_data[i - 1]) / (x_data[i] - x_data[i - 1]))
        l[i] = (2 * (x_data[i + 1] - x_data[i - 1]) - (x_data[i] - x_data[i - 1]) * m[i - 1])
        m[i] = ((x_data[i + 1] - x_data[i]) / l[i])
        z[i] = ((a[i] - (x_data[i] - x_data[i - 1]) * z[i]) / l[i])
    l[x_num - 1] = 1
    z[x_num - 1] = 0
    c[x_num - 1] = 0
    b = np.empty (x_num - 1)
    d = np.empty (x_num - 1)
    for i in range (x_num - 2, -1, -1):
        c[i] = z[i] - m[i] * c[i + 1]
        b[i] = (f_data[i + 1] - f_data[i]) / (x_data[i + 1] - x_data[i]) - (x_data[i + 1] - x_data[i]) * (c[i + 1] + 2 * c[i]) / 3
        d[i] = (c[i + 1] - c[i]) / 3 / (x_data[i + 1] - x_data[i])
    #</calculate coefficients>
    #<calculate net spline_data>
    k = 0
    for i in range (0, x_num_big):
        if k <= x_num - 1 and x_data_big[i] > x_data[k + 1]:
            k += 1
        spline_data[i] = f_data[k] + b[k] * (x_data_big[i] - x_data[k]) + c[k] * (x_data_big[i] - x_data[k])**2 + d[k] * (x_data_big[i] - x_data[k])**3
    #</calculate net spline_data>
    return spline_data
#</spline>

#<plotting>
def plt_init ():
    plt.clf ()
    x_label = ""
    y_label = ""
    plt.xlabel (x_label)
    plt.ylabel (y_label)
    plt.grid (True)
    plt.ion ()
    plt.show (block = False)
#</plotting>

#<main>
if __name__ == "__main__":
    if len (sys.argv) != 2:
        print ("Usage: " + sys.argv[0] + "<lagrange, newton, spline, combine>")
        sys.exit (1)
    if not (sys.argv[1] == "lagrange" or sys.argv[1] == "newton" or sys.argv[1] == "spline" or sys.argv[1] == "combine"):
        print ("Usage: " + sys.argv[0] + "<lagrange, newton, spline, combine>")
        sys.exit (1)
    if sys.argv[1] == "lagrange":
        interpol_data = calculate_lagrange ()
        cmd = '0'
        while cmd != 'q':
            if (cmd == 'c'):
                change_knot ()
            plt_init ()
            plt.title ("Lagrange interpolation " + function_name)
            plt.plot (x_data_big, F (x_data_big), "b-", label = 'reference f')
            plt.plot (x_data, f_data, "ro", label = 'net f')
            plt.plot (x_data_big, interpol_data, "g-", label = 'lagrange')
            plt.legend ()
            cmd = input ("enter c for knot change, q for quit: ")
        sys.exit (0)
    if sys.argv[1] == "newton":
        interpol_data = calculate_newton ()
        cmd = '0'
        while cmd != 'q':
            if (cmd == 'a'):
                add_knot ()
            plt_init ()
            plt.title ("Newton interpolation " + function_name)
            plt.plot (x_data_big, F (x_data_big), "b-", label = 'reference f')
            plt.plot (x_data, f_data, "ro", label = 'net f')
            plt.plot (x_data_big, interpol_data, "g-", label = 'newton')
            plt.legend ()
            cmd = input ("enter a for knot add, q for quit: ")
        sys.exit (0)
    if sys.argv[1] == "spline":
        interpol_data = calculate_spline ()
        cmd = '0'
        while cmd != 'q':
            plt_init ()
            plt.title ("Cubic natural spline interpolation " + function_name)
            plt.plot (x_data_big, F (x_data_big), "b-", label = 'reference f')
            plt.plot (x_data, f_data, "ro", label = 'net f')
            plt.plot (x_data_big, interpol_data, "g-", label = 'spline')
            plt.legend ()
            cmd = input ("enter q for quit: ")
        sys.exit (0)
    if sys.argv[1] == "combine":
        interpol_data_newton = calculate_newton ()
        interpol_data_lagrange = calculate_lagrange ()
        interpol_data_spline = calculate_spline ()
        cmd = '0'
        while cmd != 'q':
            plt_init ()
            plt.title ("Lagrange Newton and spline interpolation " + function_name)
            plt.plot (x_data_big, F (x_data_big), "b-", label = 'reference f')
            plt.plot (x_data, f_data, "ro", label = 'net f')
            plt.plot (x_data_big, interpol_data_lagrange, "g-", label = 'lagrange')
            plt.plot (x_data_big, interpol_data_newton, "c--", label = 'newton')
            plt.plot (x_data_big, interpol_data_spline, "m--", label = 'spline')
            plt.legend ()
            cmd = input ("enter q for quit: ")
        sys.exit (0)
#</main>