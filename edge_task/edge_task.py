#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import subprocess

#<dimensions and data>
x_low = 0
x_high = math.pi / 8
x_num = 20

eq_x0 = x_low
eq_x1 = x_high

#<1 kind>
eq_y0 = 1
eq_y1 = 0
#</1 kind>

#<2 kind>
eq_dy0 = -math.sqrt (2) * math.sin (math.pi / 8) * math.exp (-math.pi / 8) + 0.5
eq_dy1 = (-(1 + (1 / math.sqrt (2)) * math.sin (math.pi / 8) * math.exp (-math.pi / 8)) * 2 * math.sqrt (2) * math.exp (math.pi / 4)) + math.exp (math.pi / 8) * (math.sin (math.pi / 8) + math.cos (math.pi / 8)) / 2
#</2 kind>

mistake_num_max = 100
mistake_num_min = 20

x_data = np.linspace (x_low, x_high, x_num)
#</dimensions and data>

#<equation>
def eq_g (x, y, z):
    return z

def eq_f (x, y, z):
    return np.exp (x) * (2 * np.sin (x) - np.cos (x)) + 4 * z - 8 * y

def eq_ref_y (x):
    return (np.cos (2 * x) - (1 + 1 / math.sqrt (2) * math.exp (-math.pi / 8) * math.sin (math.pi / 8)) * np.sin (2 * x)) * np.exp (2 * x) + 1 / 2 * np.exp (x) * np.sin (x)

# z' = f (x, y, z)
# y' = z = g (x, y, z), z0, y0
class difequation:
    def __init__ (self, f, g, x0, y0, z0, name = None, ref_y = None, dif_F = None):
        self.f = f
        self.ref_y = ref_y
        self.g = g
        self.y0 = y0
        self.x0 = x0
        self.z0 = z0
        self.name = name
        self.dif_F = dif_F

    def solve_runge_cutt_4 (self, x_data):
        assert (isinstance (x_data, np.ndarray))
        y_data = np.empty (x_data.size)
        z_data = np.empty (x_data.size)
        i = 0
        #while x_data[i] < self.x0:
        #    i += 1
        y_data[0] = self.y0
        z_data[0] = self.z0
        for i in range (1, x_data.size):
            h = (x_data[i] - x_data[i - 1])

            k1 = self.f (x_data[i - 1], y_data[i - 1], z_data[i - 1])
            q1 = self.g (x_data[i - 1], y_data[i - 1], z_data[i - 1])

            k2 = self.f (x_data[i - 1] + h / 2, y_data[i - 1] + q1 * h / 2, z_data[i - 1] + k1 * h / 2)
            q2 = self.g (x_data[i - 1] + h / 2, y_data[i - 1] + q1 * h / 2, z_data[i - 1] + k1 * h / 2)

            k3 = self.f (x_data[i - 1] + h / 2, y_data[i - 1] + q2 * h / 2, z_data[i - 1] + k2 * h / 2)
            q3 = self.g (x_data[i - 1] + h / 2, y_data[i - 1] + q2 * h / 2, z_data[i - 1] + k2 * h / 2)

            k4 = self.f (x_data[i - 1] + h, y_data[i - 1] + q3 * h, z_data[i - 1] + k3 * h)
            q4 = self.g (x_data[i - 1] + h, y_data[i - 1] + q3 * h, z_data[i - 1] + k3 * h)
            z_data[i] = z_data[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6
            y_data[i] = y_data[i - 1] + (q1 + 2 * q2 + 2 * q3 + q4) * h / 6
        return y_data
#</equation>

#<shoot 1st kind>
z1 = 3
z2 = -3
def solve_1_kind_shoot_linear (x_data):
    eq1 = difequation (eq_f, eq_g, eq_x0, eq_y0, z1, ref_y = eq_ref_y)
    y1 = eq1.solve_runge_cutt_4 (x_data)
    eq2 = difequation (eq_f, eq_g, eq_x0, eq_y0, z2)
    y2 = eq2.solve_runge_cutt_4 (x_data)
    C1 = (y2[x_data.size - 1] - eq_y1) / (y2[x_data.size - 1] - y1[x_data.size - 1])
    C2 = (y1[x_data.size - 1] - eq_y1) / (y1[x_data.size - 1] - y2[x_data.size - 1])
    return C1 * y1 + C2 * y2
#</shoot 1st king>

#<fin dif 1st kind>
def solve_1_kind_fin_dif (x_data):
    result = subprocess.run (["./edge_task", "1", str (eq_x0), str (eq_y0), str (eq_x1), str (eq_y1), str (x_data.size)], capture_output=True, text=True)
    res_y = np.empty (x_data.size)
    res_split = result.stdout.split ()
    for i in range (0, x_data.size):
        res_y[i] = res_split[i]
    return res_y
#</fin dif 1st kind>

#<fin dif 2st kind>
def solve_2_kind_fin_dif (x_data):
    result = subprocess.run (["./edge_task", "23", str (eq_x0), str (eq_dy0), "0", str (eq_x1),  str (eq_dy1), "0", str (x_data.size)], capture_output=True, text=True)
    res_y = np.empty (x_data.size)
    res_split = result.stdout.split ()
    for i in range (0, x_data.size):
        res_y[i] = res_split[i]
    return res_y
#</fin dif 2st kind>

#<mistakes graph calculation>

def mistakes_graph (method, reference, num_start, num_end):
    mistake_graph = np.empty (num_end - num_start)
    for num in range (num_start, num_end):
        temp_x_data = np.linspace (x_low, x_high, num)
        temp_f_data = method (temp_x_data)
        temp_ref_data = reference (temp_x_data)
        max = 0
        for i in range (0, temp_x_data.size):
            if abs (temp_f_data[i] - temp_ref_data[i]) > max:
                max = abs (temp_f_data[i] - temp_ref_data[i])
        mistake_graph[num - num_start] = max
    return mistake_graph

#</mistakes graph calculation>

if __name__ == "__main__":
    #if len (sys.argv) != 2:
    #    print ("Usage: " + sys.argv[0] + "<solution, mistake>")
    #    sys.exit (1)
    #if not (sys.argv[1] == "solution" or sys.argv[1] == "mistake"):
    #    print ("Usage: " + sys.argv[0] + "<solution, mistake>")
    #    sys.exit (1)
    y_shoot = solve_1_kind_shoot_linear (x_data)
    y_1_kind = solve_1_kind_fin_dif (x_data)
    y_2_kind = solve_2_kind_fin_dif (x_data)
    plt.subplot (121)
    plt.grid (True)
    plt.title ("Solution")
    plt.plot (x_data, eq_ref_y (x_data), "r-")
    plt.plot (x_data, y_shoot, "b--")
    plt.plot (x_data, y_1_kind, "g--")
    plt.plot (x_data, y_2_kind, "c.")

    #mistake
    mistake_x_data = np.linspace (mistake_num_min, mistake_num_max, mistake_num_max - mistake_num_min)
    mistake_graph_1_kind_shoot = mistakes_graph (solve_1_kind_shoot_linear, eq_ref_y, mistake_num_min, mistake_num_max)
    mistake_graph_1_kind_fin_dif = mistakes_graph (solve_1_kind_fin_dif, eq_ref_y, mistake_num_min, mistake_num_max)
    mistake_graph_2_kind_fin_dif = mistakes_graph (solve_2_kind_fin_dif, eq_ref_y, mistake_num_min, mistake_num_max)
    plt.subplot (122)
    plt.grid (True)
    plt.title ("Mistake")
    plt.xscale ('log')
    plt.yscale ('log')

    plt.plot (mistake_x_data, mistake_graph_1_kind_shoot, "r-")
    k_1_kind_shoot = (math.log (mistake_graph_1_kind_shoot[mistake_num_max - mistake_num_min - 1]) - math.log (mistake_graph_1_kind_shoot[0])) / (math.log (mistake_x_data[mistake_num_max - mistake_num_min - 1]) - math.log (mistake_x_data[0]))
    plt.text (mistake_x_data[5], mistake_graph_1_kind_shoot[0], 'k_1_kind_shoot = {:.6f}'.format(k_1_kind_shoot))

    plt.plot (mistake_x_data, mistake_graph_1_kind_fin_dif, "b-")
    k_1_kind_fin_dif = (math.log (mistake_graph_1_kind_fin_dif[mistake_num_max - mistake_num_min - 1]) - math.log (mistake_graph_1_kind_fin_dif[0])) / (math.log (mistake_x_data[mistake_num_max - mistake_num_min - 1]) - math.log (mistake_x_data[0]))
    plt.text (mistake_x_data[5], mistake_graph_1_kind_fin_dif[0], 'k_1_kind_fin_dif = {:.6f}'.format(k_1_kind_fin_dif))

    plt.plot (mistake_x_data, mistake_graph_2_kind_fin_dif, "g-")
    k_2_kind_fin_dif = (math.log (mistake_graph_2_kind_fin_dif[mistake_num_max - mistake_num_min - 1]) - math.log (mistake_graph_2_kind_fin_dif[0])) / (math.log (mistake_x_data[mistake_num_max - mistake_num_min - 1]) - math.log (mistake_x_data[0]))
    plt.text (mistake_x_data[5], mistake_graph_2_kind_fin_dif[0], 'k_1_kind_fin_dif = {:.6f}'.format(k_2_kind_fin_dif))
    plt.show ()